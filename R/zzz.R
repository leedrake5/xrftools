
# Package environment for cached lookup tables
.xrf_cache <- new.env(parent = emptyenv())

# Last-call cache of the deconvolution design matrix, reused across a batch of spectra that share
# the same energy grid, peak list and line-shape (only the response vector changes). See
# xrf_deconvolute_gaussian_least_squares(cache_templates = TRUE).
.xrf_template_cache <- new.env(parent = emptyenv())

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

    # Build per (element:shell) photoionization cross-section grids for log-log interpolation.
    # Below an edge the tabulated cross section is NA; drop those so each grid starts at the edge.
    xs <- x_ray_cross_sections
    xs <- xs[!is.na(xs$cross_section) & xs$cross_section > 0 & xs$energy_kev > 0, ]
    xs <- xs[order(xs$energy_kev), ]
    .xrf_cache$cross_section_split <- split(
      data.frame(energy_kev = xs$energy_kev, cross_section = xs$cross_section),
      paste(xs$element, xs$shell, sep = ":")
    )

    # Per-element mass attenuation grids (cm^2/g) from the full EPDL table: photoelectric (for the
    # active-volume absorption / escape mu-ratios) and total (for window/dead-layer transmission).
    ma <- x_ray_mass_attenuation
    ma <- ma[is.finite(ma$energy_kev) & ma$energy_kev > 0, ]
    ma <- ma[order(ma$element, ma$energy_kev), ]
    pe_ok <- ma[is.finite(ma$photoelectric) & ma$photoelectric > 0, ]
    .xrf_cache$pe_mu_split <- split(
      data.frame(energy_kev = pe_ok$energy_kev, mu = pe_ok$photoelectric), pe_ok$element
    )
    tot_ok <- ma[is.finite(ma$total) & ma$total > 0, ]
    .xrf_cache$total_mu_split <- split(
      data.frame(energy_kev = tot_ok$energy_kev, mu = tot_ok$total), tot_ok$element
    )

    # Subshell absorption-edge energies keyed "element:shell" (from the line table), used by the
    # electron-impact ionization model for the overvoltage U = E0 / E_edge.
    xe <- x_ray_xrf_energies
    ek <- unique(xe[, c("element", "edge", "edge_kev")])
    .xrf_cache$edge_kev_vec <- setNames(ek$edge_kev, paste(ek$element, ek$edge, sep = ":"))

    # Bote-Salvat (2008) electron-impact ionization coefficients, keyed "z:subshell".
    bs <- x_ray_bote_salvat
    .xrf_cache$bote_salvat <- bs
    .xrf_cache$bote_salvat_key <- paste(bs$z, bs$subshell, sep = ":")

    .xrf_cache$initialized <- TRUE
  }
  invisible()
}
