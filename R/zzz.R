
# Package environment for cached lookup tables
.xrf_cache <- new.env(parent = emptyenv())

# Keyed FIFO store for xrf_energies() results (multi-condition batches alternate beam energies, so a
# last-call memo missed on every switch; see .xrf_kv_get/.xrf_kv_put below).
.xrf_energies_store <- new.env(parent = emptyenv())

# Keyed cache of deconvolution design matrices, reused across a batch of spectra that share the
# same energy grid, peak list and line-shape (only the response vector changes). Holds SEVERAL
# entries (FIFO), because real batches interleave beam conditions (e.g. 50/12/10 kV cycles per
# sample): a single-slot cache missed on every condition switch and rebuilt every template. See
# xrf_deconvolute_gaussian_least_squares(cache_templates = TRUE).
.xrf_template_cache <- new.env(parent = emptyenv())

# tiny keyed FIFO store on an environment: entries = list of list(key, val); identical() key match.
# max_entries bounds memory (a cached design matrix is ~1-2 MB); linear scan over <= 8 keys is ~us.
.xrf_kv_get <- function(env, key) {
  for (e in env$entries) if (identical(e$key, key)) return(e$val)
  NULL
}
.xrf_kv_put <- function(env, key, val, max_entries = 8L) {
  env$entries <- c(env$entries, list(list(key = key, val = val)))
  if (length(env$entries) > max_entries) env$entries <- env$entries[-1]
  invisible()
}

# Collapse duplicated grid energies (the EPDL tables double each absorption edge: one row just below,
# one just at the edge) ONCE here, so the interpolators can use the fast ties = "ordered" path (the
# default ties handling routes through tapply()/mean() per call, which dominated the physics-lookup
# profile). The collapse MUST average the LOGS (geometric mean): the interpolators have always called
# approx() on log-transformed grids, whose ties = "mean" averaged log values -- averaging the linear
# values instead changes edge-point results enormously (the O K-edge pair spans 0.29 -> 2.3e5 cm^2/g:
# arithmetic mean 1.1e5 vs the legacy geometric 255, a 440x difference that leaked into
# efficiency-weighted areas for lines near absorber edges). Geometric keeps every lookup byte-identical
# to the pre-cache behaviour. Also precomputes the log grids so nothing re-logs the table per call.
.xrf_collapse_loglog <- function(x, y) {
  ux <- unique(x)                       # x arrives sorted ascending, so ux is sorted too
  ly <- log(y)
  if (length(ux) != length(x)) {
    idx <- match(x, ux)
    ly <- as.vector(rowsum(ly, idx) / rowsum(rep(1, length(ly)), idx))
    x <- ux
  }
  data.frame(energy_kev = x, value = exp(ly), loge = log(x), logv = ly)
}

# Lazy initialization of lookup tables using named vectors
# Named vector lookup with match() is faster than left_join for repeated lookups
.init_physics_cache <- function() {
  if (is.null(.xrf_cache$initialized)) {
    # Build named vector lookup for emission probabilities
    ep <- x_ray_emission_probabilities
    ep_keys <- paste(ep$element, ep$trans, sep = ":")
    .xrf_cache$emission_prob_vec <- setNames(ep$emission_probability, ep_keys)

    # Build named vector lookup for fluorescence yields
    fy <- x_ray_fluorescence_yields
    fy_keys <- paste(fy$element, fy$shell, sep = ":")
    .xrf_cache$fluorescence_yield_vec <- setNames(fy$fluorescence_yield, fy_keys)

    # Build named vector lookup for Coster-Kronig probabilities
    ck <- x_ray_coster_kronig_probabilities
    ck_keys <- paste(ck$element, ck$shell, ck$coster_kronig_trans, sep = ":")
    .xrf_cache$coster_kronig_vec <- setNames(ck$coster_kronig_prob, ck_keys)

    # Optional MODERN Coster-Kronig source (Elam, Ravel & Sieber 2002; opt-in via xrf_set_ck_source("Elam")
    # or option xrftools.ck_source). Built only if the data-raw-generated table is present, else NULL.
    .xrf_cache$coster_kronig_vec_elam <- tryCatch({
      cke <- x_ray_coster_kronig_probabilities_elam
      setNames(cke$coster_kronig_prob, paste(cke$element, cke$shell, cke$coster_kronig_trans, sep = ":"))
    }, error = function(e) NULL)

    # Build per (element:shell) photoionization cross-section grids for log-log interpolation.
    # Below an edge the tabulated cross section is NA; drop those so each grid starts at the edge.
    # Duplicate (edge-doubled) energies are collapsed and the log grids precomputed (.xrf_collapse_loglog).
    xs <- x_ray_cross_sections
    xs <- xs[!is.na(xs$cross_section) & xs$cross_section > 0 & xs$energy_kev > 0, ]
    xs <- xs[order(xs$energy_kev), ]
    .xrf_cache$cross_section_split <- lapply(
      split(data.frame(energy_kev = xs$energy_kev, cross_section = xs$cross_section),
            paste(xs$element, xs$shell, sep = ":")),
      function(g) .xrf_collapse_loglog(g$energy_kev, g$cross_section)
    )

    # Per-element mass attenuation grids (cm^2/g) from the full EPDL table: photoelectric (for the
    # active-volume absorption / escape mu-ratios) and total (for window/dead-layer transmission),
    # plus coherent (Rayleigh) / incoherent (Compton) scatter for the scatter templates. All collapsed
    # + pre-logged like the cross sections.
    ma <- x_ray_mass_attenuation
    ma <- ma[is.finite(ma$energy_kev) & ma$energy_kev > 0, ]
    ma <- ma[order(ma$element, ma$energy_kev), ]
    mu_split <- function(col) {
      ok <- is.finite(ma[[col]]) & ma[[col]] > 0
      lapply(split(data.frame(energy_kev = ma$energy_kev[ok], mu = ma[[col]][ok]), ma$element[ok]),
             function(g) .xrf_collapse_loglog(g$energy_kev, g$mu))
    }
    .xrf_cache$pe_mu_split <- mu_split("photoelectric")
    .xrf_cache$total_mu_split <- mu_split("total")
    .xrf_cache$rayleigh_mu_split <- mu_split("rayleigh")
    .xrf_cache$compton_mu_split <- mu_split("compton")

    # Atomic form factors F(q, Z) and incoherent scattering functions S(q, Z) for the fixed-angle
    # scatter shapes, split per element with pre-logged grids (reusing .xrf_collapse_loglog; the
    # "energy_kev" slot holds q in 1/angstrom). The q = 0 / value <= 0 rows are dropped for the
    # log-log grids: F below the first positive-q point is flat (= Z), S extrapolates down its
    # low-q power law (~q^2), S above the grid saturates to Z -- handled in .xrf_scatter_factor_at.
    ffd <- x_ray_form_factors
    ok_f <- ffd$q_inv_angstrom > 0 & is.finite(ffd$form_factor) & ffd$form_factor > 0
    .xrf_cache$ff_split <- lapply(
      split(data.frame(q = ffd$q_inv_angstrom[ok_f], v = ffd$form_factor[ok_f]), ffd$element[ok_f]),
      function(g) .xrf_collapse_loglog(g$q, g$v)
    )
    sfd <- x_ray_incoherent_functions
    ok_s <- sfd$q_inv_angstrom > 0 & is.finite(sfd$incoherent_function) & sfd$incoherent_function > 0
    .xrf_cache$sf_split <- lapply(
      split(data.frame(q = sfd$q_inv_angstrom[ok_s], v = sfd$incoherent_function[ok_s]), sfd$element[ok_s]),
      function(g) .xrf_collapse_loglog(g$q, g$v)
    )

    # Biggs Compton profiles J(p_z) per element, for the impulse-approximation Doppler shape of the
    # Compton scatter peaks: linear p_z grid with log(J) for interpolation (J falls ~exponentially).
    cpd <- x_ray_compton_profiles
    ok_c <- is.finite(cpd$j_total) & cpd$j_total > 0
    .xrf_cache$cp_split <- lapply(
      split(data.frame(pz = cpd$pz_au[ok_c], logj = log(cpd$j_total[ok_c])), cpd$element[ok_c]),
      function(g) g[order(g$pz), , drop = FALSE]
    )

    # Per-subshell Compton profiles with occupancy + binding energy (bound-electron thresholds for
    # the Doppler cluster): per element a compact list(pz, logj[pz x shell], occ, bind_kev).
    cpsd <- x_ray_compton_profiles_shells
    .xrf_cache$cp_shell_split <- lapply(split(cpsd, cpsd$element), function(g) {
      pzv <- sort(unique(g$pz_au))
      shv <- sort(unique(g$shell))
      J <- matrix(0, length(pzv), length(shv))
      occ <- numeric(length(shv)); bind <- numeric(length(shv))
      for (k in seq_along(shv)) {
        gk <- g[g$shell == shv[k], ]
        J[, k] <- gk$j_shell[match(pzv, gk$pz_au)]
        occ[k] <- gk$occupancy[1]; bind[k] <- gk$binding_kev[1]
      }
      list(pz = pzv, logj = log(pmax(J, 1e-300)), occ = occ, bind_kev = bind)
    })

    # Subshell absorption-edge energies keyed "element:shell" (from the line table), used by the
    # electron-impact ionization model for the overvoltage U = E0 / E_edge.
    xe <- x_ray_xrf_energies
    ek <- unique(xe[, c("element", "edge", "edge_kev")])
    .xrf_cache$edge_kev_vec <- setNames(ek$edge_kev, paste(ek$element, ek$edge, sep = ":"))
    # per-element split of the line table: xrf_fp_sensitivity and the secondary-fluorescence
    # exciter-line lists subset by element once per element per call, so a pre-split list lookup
    # beats scanning the full table each time
    .xrf_cache$xe_by_element <- split(xe, xe$element)

    # Bote-Salvat (2008) electron-impact ionization coefficients, keyed "z:subshell".
    bs <- x_ray_bote_salvat
    .xrf_cache$bote_salvat <- bs
    .xrf_cache$bote_salvat_key <- paste(bs$z, bs$subshell, sep = ":")

    .xrf_cache$initialized <- TRUE
  }
  invisible()
}
