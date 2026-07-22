
# Momentum-transfer variable q = sin(theta/2)/lambda = E[keV] sin(theta/2) / 12.39842 (1/angstrom),
# the argument of the atomic form factor F(q, Z) and incoherent scattering function S(q, Z).
.xrf_q_momentum_transfer <- function(energy_kev, theta_rad) {
  energy_kev * sin(theta_rad / 2) / 12.39842
}

# Interpolate a cached F/S grid (log-log; grids from .xrf_collapse_loglog with q in the energy_kev
# slot). below = "clamp" (F: flat at Z toward q -> 0) or "extrapolate" (S: continue the low-q power
# law ~ q^2); above the grid F extrapolates its power-law tail, S clamps at its last value (-> Z).
.xrf_scatter_factor_at <- function(g, q, below = c("clamp", "extrapolate"),
                                   above = c("extrapolate", "clamp")) {
  below <- match.arg(below); above <- match.arg(above)
  if (is.null(g) || nrow(g) < 2) return(rep(NA_real_, length(q)))
  m <- nrow(g)
  lo <- g$energy_kev[1]; hi <- g$energy_kev[m]
  out <- rep(NA_real_, length(q))
  within <- is.finite(q) & q >= lo & q <= hi
  if (any(within)) out[within] <- exp(stats::approx(g$loge, g$logv, log(q[within]), ties = "ordered")$y)
  bel <- is.finite(q) & q < lo
  if (any(bel)) {
    out[bel] <- if (below == "clamp") g$value[1] else {
      sl <- (g$logv[2] - g$logv[1]) / (g$loge[2] - g$loge[1])
      exp(g$logv[1] + sl * (log(pmax(q[bel], 1e-12)) - g$loge[1]))
    }
  }
  abv <- is.finite(q) & q > hi
  if (any(abv)) {
    out[abv] <- if (above == "clamp") g$value[m] else {
      sl <- (g$logv[m] - g$logv[m - 1]) / (g$loge[m] - g$loge[m - 1])
      exp(g$logv[m] + sl * (log(q[abv]) - g$loge[m]))
    }
  }
  out
}

# Composition-weighted fixed-angle scattering strengths of a material, per unit mass (relative --
# the constant N_A drops into the free template amplitudes): sum_i (w_i/A_i) F_i(q)^2 for Rayleigh
# and sum_i (w_i/A_i) S_i(q) for Compton.
.xrf_material_F2 <- function(material, q) {
  .init_physics_cache()
  info <- .xrf_material_info(material)
  acc <- rep(0, length(q))
  for (el in names(info$composition)) {
    A <- unname(.xrf_atomic_mass[el]); if (!is.finite(A)) next
    f <- .xrf_scatter_factor_at(.xrf_cache$ff_split[[el]], q, below = "clamp", above = "extrapolate")
    f[!is.finite(f)] <- 0
    acc <- acc + (info$composition[[el]] / A) * f^2
  }
  acc
}
.xrf_material_S <- function(material, q) {
  .init_physics_cache()
  info <- .xrf_material_info(material)
  acc <- rep(0, length(q))
  for (el in names(info$composition)) {
    A <- unname(.xrf_atomic_mass[el]); if (!is.finite(A)) next
    sv <- .xrf_scatter_factor_at(.xrf_cache$sf_split[[el]], q, below = "extrapolate", above = "clamp")
    sv[!is.finite(sv)] <- 0
    acc <- acc + (info$composition[[el]] / A) * sv
  }
  acc
}

# Ribberfors projection p_z (atomic units): the component of the electron's initial momentum along
# the scattering vector that shifts a photon scattered at theta from E1 to E2. p_z = 0 reproduces
# the standard Compton line; the Biggs profile J(p_z) weights the rest.
.xrf_pz_au <- function(E1, E2, theta) {
  mc2 <- 510.999
  R <- sqrt(pmax(E1^2 + E2^2 - 2 * E1 * E2 * cos(theta), 1e-20))
  137.036 * (E1 - E2 - E1 * E2 * (1 - cos(theta)) / mc2) / R
}

# Composition-weighted total Compton profile of a material, per unit mass (relative):
# sum_i (w_i/A_i) J_i(|p_z|). Beyond the tabulated p_z grid (100 a.u.) J is ~0.
.xrf_material_J <- function(material, pz) {
  .init_physics_cache()
  info <- .xrf_material_info(material)
  acc <- rep(0, length(pz))
  apz <- abs(pz)
  for (el in names(info$composition)) {
    A <- unname(.xrf_atomic_mass[el]); if (!is.finite(A)) next
    g <- .xrf_cache$cp_split[[el]]; if (is.null(g) || nrow(g) < 2) next
    j <- exp(stats::approx(g$pz, g$logj, xout = pmin(apz, max(g$pz)), ties = "ordered")$y)
    j[apz > max(g$pz) | !is.finite(j)] <- 0
    acc <- acc + (info$composition[[el]] / A) * j
  }
  acc
}

# Bound-electron momentum density of a material: like .xrf_material_J, but summed per SUBSHELL with
# the kinematic binding threshold -- subshell k contributes only where the energy TRANSFER exceeds
# its binding energy B_k (an electron cannot be ejected with less). This is what displaces the
# core-shell contributions to the low-energy side of the Compton line (at Rh Kalpha / 135 deg the
# transfer at the peak is ~1.3 keV, so Si K electrons, B = 1.84 keV, are excluded from the hump's
# core entirely) and produces the small "Compton defect" shift. transfer_kev recycles with pz.
.xrf_material_J_bound <- function(material, pz, transfer_kev) {
  .init_physics_cache()
  info <- .xrf_material_info(material)
  apz <- abs(pz)
  acc <- rep(0, length(apz))
  for (el in names(info$composition)) {
    A <- unname(.xrf_atomic_mass[el]); if (!is.finite(A)) next
    g <- .xrf_cache$cp_shell_split[[el]]
    if (is.null(g)) {                     # element beyond the per-shell table: unthresholded total
      acc <- acc + (info$composition[[el]] / A) * .xrf_material_J(el, pz)
      next
    }
    pin <- pmin(apz, max(g$pz))
    for (k in seq_along(g$occ)) {
      open_k <- transfer_kev >= g$bind_kev[k]
      if (!any(open_k)) next
      j <- exp(stats::approx(g$pz, g$logj[, k], xout = pin, ties = "ordered")$y)
      j[apz > max(g$pz) | !is.finite(j)] <- 0
      acc <- acc + (info$composition[[el]] / A) * g$occ[k] * j * open_k
    }
  }
  acc
}

# Impulse-approximation Doppler cluster for one Compton-scattered line: the scattered-energy
# distribution rho(E2) ~ J(|p_z(E2)|) |dp_z/dE2|, sampled as n_nodes pseudo-lines whose weights sum
# to 1 (the caller multiplies by the line's total Compton weight KN * S * kernel, so upgrading the
# SHAPE leaves the template's total intensity untouched). The region kept holds ~99.9% of the
# profile mass; the high side is kinematically truncated at E0. This replaces the "Gaussian at the
# Compton shift x compton_broadening" heuristic with the physical, asymmetric line shape.
.xrf_compton_doppler_cluster <- function(E0, theta, scatterer, n_nodes = 31L, bound = TRUE) {
  mc2 <- 510.999
  E_C <- E0 / (1 + (E0 / mc2) * (1 - cos(theta)))
  Jat <- function(e2, pz) {
    if (isTRUE(bound)) .xrf_material_J_bound(scatterer, pz, transfer_kev = E0 - e2)
    else .xrf_material_J(scatterer, pz)
  }
  probe <- seq(max(E_C * 0.55, 0.05), E0 * (1 - 1e-9), length.out = 600L)
  pz <- .xrf_pz_au(E0, probe, theta)
  rho <- Jat(probe, pz)
  grad <- abs(c(diff(pz) / diff(probe), NA)); grad[length(grad)] <- grad[length(grad) - 1L]
  rho <- rho * grad
  ok <- is.finite(rho) & rho > 0
  if (!any(ok) && isTRUE(bound)) {        # fully bound-suppressed (no open shell): fall back to total
    return(.xrf_compton_doppler_cluster(E0, theta, scatterer, n_nodes, bound = FALSE))
  }
  if (!any(ok)) return(tibble::tibble(energy_kev = E_C, weight = 1))
  keep <- ok & rho >= 1e-3 * max(rho[ok])
  nodes <- seq(min(probe[keep]), max(probe[keep]), length.out = n_nodes)
  pzn <- .xrf_pz_au(E0, nodes, theta)
  w <- Jat(nodes, pzn)
  gn <- abs(c(diff(pzn) / diff(nodes), NA)); gn[length(gn)] <- gn[length(gn) - 1L]
  w <- w * gn
  w[!is.finite(w) | w < 0] <- 0
  if (sum(w) <= 0) return(tibble::tibble(energy_kev = E_C, weight = 1))
  tibble::tibble(energy_kev = nodes, weight = w / sum(w))
}

# Thick-sample scatter self-absorption kernel: detected scatter from a thick scatterer is
# N(E) * sigma_scatter(E) / (mu(E)/sin psi_in + mu(E_out)/sin psi_out) -- the same Sherman-type
# denominator as fluorescence, with the incident leg at the tube energy and the exit leg at the
# detected (Rayleigh: same; Compton: shifted) energy. Without it the soft scatter components are
# wildly over-weighted (a 3 keV photon scatters more per gram but almost none of it escapes the
# sample). Uses the representative scatterer's TOTAL attenuation and the geometry's incidence /
# takeoff angles. Relative shape only (templates carry a free amplitude).
.xrf_scatter_sample_kernel <- function(scatterer, e_in, e_out, geometry) {
  sin_in <- sin(geometry$incidence_deg * pi / 180)
  sin_out <- sin(geometry$takeoff_deg * pi / 180)
  a <- xrf_mass_attenuation(scatterer, e_in, type = "total") / max(sin_in, 1e-3) +
    xrf_mass_attenuation(scatterer, e_out, type = "total") / max(sin_out, 1e-3)
  k <- 1 / a
  k[!is.finite(k) | k < 0] <- 0
  k
}

#' Build Rayleigh and Compton tube-scatter peak templates
#'
#' The two largest features in most real XRF spectra are the tube's own characteristic lines
#' scattered off the sample: coherently (Rayleigh, at the tube-line energy) and incoherently
#' (Compton, at a lower, Doppler-broadened energy). They are not element lines, so unless they are
#' modelled explicitly the deconvolution assigns the residual bump to whatever element sits under
#' it (false positives), and they dominate at high tube voltages. This constructs peak templates for
#' both, suitable for binding onto the \code{peaks} passed to
#' \link{xrf_deconvolute_gaussian_least_squares} (which is done automatically when a \code{tube} is
#' supplied to the deconvolution).
#'
#' The Compton (incoherent) peak is shifted to
#' \deqn{E' = E / (1 + (E/m_e c^2)(1 - \cos\theta))}
#' with \eqn{m_e c^2 = 511} keV and \eqn{\theta} the scattering angle (\code{geometry$scatter_angle_deg}),
#' and is broadened relative to the detector resolution (Doppler broadening of the Compton profile);
#' \code{compton_broadening} scales the Compton width relative to the detector width at \eqn{E'}.
#'
#' The tube characteristic lines and their fluxes are the \emph{Ebel} anode-line intensities of
#' \link{xrf_tube_spectrum} (\code{discrete_lines = TRUE}; tube \code{filter} applied), and each
#' line is weighted by the \strong{exact fixed-angle differential} scattering strength of the
#' \code{scatterer} at the instrument's scattering angle: \eqn{(1+\cos^2\theta)\,F^2(q, Z)} for
#' Rayleigh and \eqn{KN(E,\theta)\,S(q, Z)} for Compton, using the Hubbell et al. (1975) atomic form
#' factors and incoherent scattering functions (\link{x_ray_form_factors}) with
#' \eqn{q = E\sin(\theta/2)/12.39842} -- so the template's internal line balance (e.g. anode L vs K
#' components) follows what the detector actually sees at its angle, not angle-integrated averages. Rayleigh and Compton are returned as
#' two separate "elements" (each a fixed-shape multi-line template with one free amplitude), because
#' their ratio varies with the sample (light matrices scatter more incoherently). Note the detector
#' efficiency is \emph{not} applied here; \link{xrf_deconvolute_gaussian_least_squares} applies its
#' efficiency model to these rows together with the element templates when \code{efficiency = TRUE}.
#'
#' @param tube An \link{xrf_tube} object (needs \code{anode} and \code{kv}).
#' @param geometry An \link{xrf_geometry} object (uses \code{scatter_angle_deg}).
#' @param detector_type,fano,epsilon_ev,noise_fwhm_ev Detector resolution parameters (see
#'   \link{xrf_detector_sigma_kev}) used for the scatter-peak widths. If none are supplied the
#'   \code{default_sigma} constant width is used.
#' @param default_sigma Fallback constant peak width (keV) when no detector model is given.
#' @param compton_broadening Multiplier on the detector width approximating the Doppler-broadened
#'   Compton peak -- used by \code{compton_shape = "gaussian"} only (the default "profile" mode
#'   computes the physical Doppler shape instead).
#' @param compton_shape "profile" (default): each anode line's Compton feature is the exact
#'   impulse-approximation Doppler shape -- \eqn{\rho(E_2) \propto J(|p_z|)\,|dp_z/dE_2|} from the
#'   Biggs Compton profiles (\link{x_ray_compton_profiles}) -- rendered as a pseudo-line cluster at
#'   detector width: the asymmetric hump with its long low-energy tail and kinematic high-side
#'   cutoff. "gaussian": the legacy single Gaussian at the Compton shift, widened by
#'   \code{compton_broadening} (optionally with \code{compton_tail}). The cluster weights sum to 1
#'   per line, so both modes carry the same total Compton intensity; note the reported
#'   \code{peak_area} of the fitted template is its strongest single component, so Compton
#'   normalizations calibrated in one mode must not be mixed with the other.
#' @param min_relative_intensity Smallest scatter line to include (relative to the strongest of
#'   either component; the Rayleigh and Compton templates keep a common line set).
#' @param scatterer Representative sample scatterer (element symbol or supported material) whose
#'   fixed-angle form factor \eqn{F^2(q)} / incoherent function \eqn{S(q)} and self-absorption weight
#'   the line intensities; default \code{"Si"}, matching \link{xrf_scatter_continuum}.
#' @param compton_tail,compton_beta Optional low-energy exponential tail (relative amplitude and decay
#'   length, keV) on the Compton rows, approximating the asymmetric Compton profile of bound
#'   electrons; emitted as per-line \code{tail}/\code{beta} columns consumed by the deconvolution's
#'   per-line lineshape support. Default 0 (symmetric Gaussian). \code{compton_beta = NULL} uses the
#'   broadened Compton width.
#'
#' @return A tibble of peak templates with columns \code{element}
#'   ("scatter_rayleigh"/"scatter_compton"), \code{trans}, \code{energy_kev}, \code{sigma},
#'   \code{relative_peak_intensity}.
#' @export
#'
#' @examples
#' xrf_scatter_peaks(xrf_tube("Rh", kv = 40), xrf_geometry(scatter_angle_deg = 135),
#'                   detector_type = "SDD")
#'
xrf_scatter_peaks <- function(tube, geometry = xrf_geometry(),
                              detector_type = NULL, fano = NULL, epsilon_ev = NULL,
                              noise_fwhm_ev = NULL, default_sigma = 0.07,
                              compton_broadening = 2, min_relative_intensity = 0.05,
                              scatterer = "Si", compton_tail = 0, compton_beta = NULL,
                              compton_shape = c("profile", "gaussian")) {
  compton_shape <- match.arg(compton_shape)
  stopifnot(inherits(tube, "xrf_tube"))
  if (is.null(geometry)) geometry <- xrf_geometry()
  stopifnot(inherits(geometry, "xrf_geometry"))

  empty <- tibble::tibble(
    element = character(0), trans = character(0), energy_kev = numeric(0),
    sigma = numeric(0), relative_peak_intensity = numeric(0)
  )
  if (is.na(tube$anode) || is.na(tube$kv)) {
    stop("xrf_scatter_peaks() needs the tube's anode and kv.")
  }

  # tolerate messy anode strings from vendor metadata (e.g. "Ag 50"): take the leading token
  anode <- strsplit(trimws(as.character(tube$anode)), "\\s+")[[1]][1]
  if (!(anode %in% all_elements)) {
    stop("tube anode '", tube$anode, "' is not a recognized element symbol.")
  }

  # Ebel anode-line fluxes (tube filter applied) x the scatterer's coherent/incoherent cross section
  # at each line energy: what actually reaches the detector from the sample, up to the free template
  # amplitude. (Previously the line weights were the anode's own fluorescence-production ratios from
  # xrf_energies -- no emitted-flux physics, no filter, no scatter cross section.)
  tube_lines <- xrf_tube_spectrum(tube, 1, discrete_lines = TRUE)
  if (is.null(tube_lines) || nrow(tube_lines) == 0) return(empty)

  m_e_c2 <- 510.999                                  # electron rest energy (keV)
  theta <- geometry$scatter_angle_deg * pi / 180
  E <- tube_lines$energy_kev
  E_compton <- E / (1 + (E / m_e_c2) * (1 - cos(theta)))

  # EXACT fixed-angle differentials (Hubbell 1975 F/S grids, x_ray_form_factors /
  # x_ray_incoherent_functions): Rayleigh dsigma/dOmega ~ (1 + cos^2 theta) F^2(q, Z) and Compton
  # dsigma/dOmega ~ KN(E, theta) S(q, Z), with q = E sin(theta/2)/12.39842. This replaces the earlier
  # angle-INTEGRATED mu_coh/mu_inc heuristics: at backscatter angles the form factor kills coherent
  # scattering of the hard components far faster than mu_coh(E) implies, and S(q) suppresses soft
  # Compton -- both directly shape what the detector actually sees at ITS angle.
  qm <- .xrf_q_momentum_transfer(E, theta)
  F2 <- .xrf_material_F2(scatterer, qm)
  Sq <- .xrf_material_S(scatterer, qm)
  P <- E_compton / E
  knF <- P^2 * (P + 1 / P - sin(theta)^2)          # fixed-angle Klein-Nishina (constants dropped)
  knF[!is.finite(knF) | knF < 0] <- 0

  rpi_ray <- tube_lines$flux * (1 + cos(theta)^2) * F2 *
    .xrf_scatter_sample_kernel(scatterer, E, E, geometry)
  rpi_com <- tube_lines$flux * knF * Sq *
    .xrf_scatter_sample_kernel(scatterer, E, E_compton, geometry)
  # one COMMON keep set (a line survives if either component clears the threshold), so the Rayleigh
  # and Compton templates stay row-aligned
  mr <- max(rpi_ray, 0); mc <- max(rpi_com, 0)
  keep <- (mr > 0 & rpi_ray >= min_relative_intensity * mr) |
    (mc > 0 & rpi_com >= min_relative_intensity * mc)
  if (!any(keep)) return(empty)
  E <- E[keep]; E_compton <- E_compton[keep]
  rpi_ray <- rpi_ray[keep]; rpi_com <- rpi_com[keep]

  detector_given <- !is.null(detector_type) || !is.null(fano) ||
    !is.null(epsilon_ev) || !is.null(noise_fwhm_ev)
  sigma_at <- function(e) {
    if (detector_given) {
      xrf_detector_sigma_kev(e, detector_type = detector_type, fano = fano,
                             epsilon_ev = epsilon_ev, noise_fwhm_ev = noise_fwhm_ev)
    } else {
      rep(default_sigma, length(e))
    }
  }

  ray_rows <- tibble::tibble(
    element = "scatter_rayleigh",
    trans = paste0("rayleigh_", round(E, 3)),
    energy_kev = E,
    sigma = sigma_at(E),
    relative_peak_intensity = rpi_ray
  )
  if (compton_shape == "profile") {
    # EXACT Doppler shape (impulse approximation, Biggs profiles): each anode line's Compton feature
    # is a cluster of pseudo-lines sampling rho(E2) ~ J(|p_z|)|dp_z/dE2| at DETECTOR width -- the
    # physical asymmetric hump (long low-energy tail, kinematic cutoff at E0) replaces the broadened
    # Gaussian. Cluster weights sum to 1 per line, so the line's total Compton intensity
    # (KN * S * kernel) is unchanged; compton_broadening / compton_tail are gaussian-mode-only.
    com_rows <- purrr::map_dfr(seq_along(E), function(i) {
      cl <- .xrf_compton_doppler_cluster(E[i], theta, scatterer)
      spc <- if (nrow(cl) > 1) cl$energy_kev[2] - cl$energy_kev[1] else 0
      tibble::tibble(
        element = "scatter_compton",
        trans = paste0("compton_", round(E[i], 3), "_", seq_len(nrow(cl))),
        energy_kev = cl$energy_kev,
        sigma = pmax(sigma_at(cl$energy_kev), 0.8 * spc),
        relative_peak_intensity = rpi_com[i] * cl$weight
      )
    })
    return(dplyr::bind_rows(ray_rows, com_rows))
  }
  com_sigma <- sigma_at(E_compton) * compton_broadening
  out <- dplyr::bind_rows(
    ray_rows,
    tibble::tibble(
      element = "scatter_compton",
      trans = paste0("compton_", round(E, 3)),
      energy_kev = E_compton,
      sigma = com_sigma,
      relative_peak_intensity = rpi_com
    )
  )
  # Optional asymmetric Compton shape for GAUSSIAN mode: a per-line Hypermet tail on the Compton rows
  # (consumed by the deconvolution's per-line lineshape columns) as a first-order stand-in for the
  # Compton-profile asymmetry; `compton_beta` defaults to the broadened Compton width. Rayleigh pure.
  if (is.finite(compton_tail) && compton_tail > 0) {
    out$tail <- ifelse(out$element == "scatter_compton", compton_tail, 0)
    out$step <- 0
    out$beta <- ifelse(out$element == "scatter_compton",
                       if (is.null(compton_beta)) com_sigma else compton_beta, NA_real_)
  }
  out
}

#' Build scattered-tube-continuum peak templates (Rayleigh + Compton of the bremsstrahlung)
#'
#' \link{xrf_scatter_peaks} scatters only the tube's discrete anode \emph{lines}. The tube
#' \strong{bremsstrahlung continuum} also scatters off the sample -- coherently (Rayleigh, at the same
#' energy) and incoherently (Compton, Doppler-broadened and shifted) -- producing the broad scattered-
#' continuum background (with its Compton edge near the shifted endpoint) that dominates real spectra at
#' high tube voltage and for light matrices. Left unmodelled it is absorbed into the SNIP baseline (which
#' cannot reproduce the Compton edge) and biases every net area. This samples the Ebel continuum
#' (\link{xrf_tube_spectrum}, characteristic lines excluded so they are not double-counted with
#' \link{xrf_scatter_peaks}) and returns it as two fixed-shape multi-line templates (one free amplitude
#' each), suitable to bind onto the \code{peaks} of \link{xrf_deconvolute_gaussian_least_squares}.
#'
#' \strong{Use against an un-baselined spectrum, with a fitted background.} These templates ARE the smooth
#' scatter background, so they must replace -- not stack on top of -- a subtracted baseline. The validated
#' recipe is to fit \code{xrf_deconvolute_gaussian_least_squares(..., scatter_continuum = TRUE, background =
#' 6..12)} against \code{.values = cps} (no SNIP baseline): the continuum + the smooth \code{background}
#' basis together model the background, and real trace areas are recovered cleanly (r^2 = 1 on synthetics;
#' SNIP-only over-estimates traces 100-300\% under a Compton hump). Fitting them alongside a
#' SNIP-baseline-subtracted response instead double-counts the background and biases real peaks.
#'
#' @inheritParams xrf_scatter_peaks
#' @param scatterer Representative sample scatterer (element symbol or a supported material) whose
#'   fixed-angle \eqn{F^2(q, Z)} / \eqn{S(q, Z)} scattering strengths and self-absorption set the
#'   continuum \emph{shape} (the overall amplitude is free in the fit). Defaults to \code{"Si"} (a
#'   light-matrix stand-in).
#' @param n_grid Number of energies sampled across the continuum (pseudo-lines per scatter type).
#' @param energy_min_kev Low-energy start of the sampled continuum (keV).
#' @param second_order If TRUE (default), add a third template ("scatter_compton_continuum2"): the
#'   collapsed \emph{second-order} Compton continuum -- tube photons (continuum + anode lines) pushed
#'   through the Klein-Nishina Compton operator twice with independently-drawn angles. Doubly-scattered
#'   photons fill the valley below the single-scatter Compton features and the broad low-energy
#'   pedestal that single-scatter models systematically miss; the direction correlation between the
#'   two scatters and the two internal path legs are approximated (see the source), so the template is
#'   a physically-shaped basis function with a free amplitude, not an absolute prediction.
#' @param active_thickness_um,be_window_um,dead_layer_um,air_path_cm,atmosphere,window Detector /
#'   measurement-path overrides for the built-in efficiency weighting of the continuum shape
#'   (\link{xrf_detector_efficiency}); \code{NULL}/\code{"Air"} use the \code{detector_type} preset.
#'   \link{xrf_deconvolute_gaussian_least_squares} threads its own efficiency settings through, so
#'   the continuum shape sees the same window/air path as the element templates.
#'
#' @return A tibble of templates with \code{element} ("scatter_rayleigh_continuum"/
#'   "scatter_compton_continuum"), \code{trans}, \code{energy_kev}, \code{sigma},
#'   \code{relative_peak_intensity}, or zero rows.
#' @export
#'
#' @examples
#' xrf_scatter_continuum(xrf_tube("Rh", kv = 45), xrf_geometry(scatter_angle_deg = 135),
#'                       detector_type = "SDD")
#'
xrf_scatter_continuum <- function(tube, geometry = xrf_geometry(),
                                  detector_type = NULL, fano = NULL, epsilon_ev = NULL,
                                  noise_fwhm_ev = NULL, default_sigma = 0.07,
                                  compton_broadening = 2, scatterer = "Si", n_grid = 250L,
                                  energy_min_kev = 1, second_order = TRUE,
                                  active_thickness_um = NULL, be_window_um = NULL,
                                  dead_layer_um = NULL, air_path_cm = NULL, atmosphere = "Air",
                                  window = NULL) {
  stopifnot(inherits(tube, "xrf_tube"))
  if (is.null(geometry)) geometry <- xrf_geometry()
  stopifnot(inherits(geometry, "xrf_geometry"))
  empty <- tibble::tibble(element = character(0), trans = character(0), energy_kev = numeric(0),
                          sigma = numeric(0), relative_peak_intensity = numeric(0))
  if (is.na(tube$anode) || is.na(tube$kv)) stop("xrf_scatter_continuum() needs the tube's anode and kv.")
  E0 <- tube$kv
  if (!is.finite(E0) || E0 <= energy_min_kev) return(empty)

  grid <- seq(energy_min_kev, E0 * 0.999, length.out = n_grid)
  spacing <- grid[2] - grid[1]
  N <- xrf_tube_spectrum(tube, grid, char_peak_ratio = 0)      # CONTINUUM only (lines handled elsewhere)

  theta <- geometry$scatter_angle_deg * pi / 180
  m_e_c2 <- 510.999
  E_compton <- grid / (1 + (grid / m_e_c2) * (1 - cos(theta)))
  # exact fixed-angle scattering strengths at THIS angle (see xrf_scatter_peaks): F^2(q) for
  # Rayleigh, S(q) for Compton, with q = E sin(theta/2)/12.39842
  qg <- .xrf_q_momentum_transfer(grid, theta)
  F2g <- .xrf_material_F2(scatterer, qg)
  Sg <- .xrf_material_S(scatterer, qg)

  detector_given <- !is.null(detector_type) || !is.null(fano) ||
    !is.null(epsilon_ev) || !is.null(noise_fwhm_ev)
  # render width: at least ~1.5 grid steps so adjacent pseudo-lines overlap into a SMOOTH continuum
  # (the background is intrinsically smooth, so this sampling broadening is harmless).
  sigma_at <- function(e, extra = 1) {
    s <- if (detector_given) {
      xrf_detector_sigma_kev(e, detector_type = detector_type, fano = fano,
                             epsilon_ev = epsilon_ev, noise_fwhm_ev = noise_fwhm_ev)
    } else rep(default_sigma, length(e))
    pmax(s * extra, 1.5 * spacing)
  }

  # The deconvolution renders each template line at HEIGHT proportional to intensity/sigma (it treats
  # relative_peak_intensity as a line AREA). A densely-sampled sum of such Gaussians has envelope
  # height(E) ~ (intensity/sigma) * sigma * sqrt(2*pi) / spacing = intensity * sqrt(2*pi) / spacing, so to
  # make the rendered HEIGHT follow the physical value N(E)*sigma_scatter(E) we pre-multiply by the constant
  # grid SPACING (NOT the per-line sigma). Using sigma cancels only while sigma is clamped to ~1.5*spacing;
  # once sigma un-clamps (finer n_grid or a coarse detector) the sigma(E) growth would distort the shape.
  # (For the Compton branch the E->E' Jacobian cancels between the physical density and the sampling density,
  # so the same constant spacing is correct there too.)
  s_ray <- sigma_at(grid)
  s_com <- sigma_at(E_compton, compton_broadening)
  # Fixed-angle shape corrections on top of the angle-integrated scatter cross sections:
  #  - Compton: the exact fixed-angle differential KN(E, theta) * S(q, Z); Rayleigh: (1+cos^2) F^2(q, Z)
  #    (Hubbell 1975 grids). The (1+cos^2) factor is constant at fixed angle but kept for clarity.
  #  - detector efficiency at the DETECTED energy (grid for Rayleigh, E' for Compton) shapes both continua --
  #    notably the low-energy window cutoff and the high-energy active-layer falloff -- when a preset is known.
  P <- E_compton / grid
  kn <- P^2 * (P + 1 / P - sin(theta)^2)
  kn[!is.finite(kn) | kn < 0] <- 0
  eff_at <- function(e) {
    xrf_detector_efficiency(e, detector_type = detector_type,
                            active_thickness_um = active_thickness_um,
                            be_window_um = be_window_um, dead_layer_um = dead_layer_um,
                            air_path_cm = air_path_cm, atmosphere = atmosphere, window = window)
  }
  eff_ray <- if (!is.null(detector_type)) eff_at(grid) else rep(1, length(grid))
  eff_com <- if (!is.null(detector_type)) eff_at(E_compton) else rep(1, length(E_compton))
  eff_ray[!is.finite(eff_ray)] <- 0; eff_com[!is.finite(eff_com)] <- 0
  # thick-sample self-absorption of the scattered photon (incident leg at E, exit leg at the detected
  # energy): the same kernel as xrf_scatter_peaks, so line and continuum scatter share one physics
  ker_ray <- .xrf_scatter_sample_kernel(scatterer, grid, grid, geometry)
  ker_com <- .xrf_scatter_sample_kernel(scatterer, grid, E_compton, geometry)
  out <- dplyr::bind_rows(
    tibble::tibble(element = "scatter_rayleigh_continuum", trans = paste0("raycont_", seq_along(grid)),
                   energy_kev = grid, sigma = s_ray,
                   relative_peak_intensity = N * (1 + cos(theta)^2) * F2g * ker_ray * eff_ray * spacing),
    tibble::tibble(element = "scatter_compton_continuum", trans = paste0("comcont_", seq_along(grid)),
                   energy_kev = E_compton, sigma = s_com,
                   relative_peak_intensity = N * kn * Sg * ker_com * eff_com * spacing)
  )

  # ---- collapsed SECOND-ORDER Compton continuum ------------------------------------------------
  # Photons that Compton-scatter TWICE before escaping fill the "valley" below the single-scatter
  # Compton features and the broad pedestal toward low energy -- a single-scatter-only model
  # systematically under-predicts the measured background there (the effect Monte Carlo codes get
  # for free). Collapsed 2nd-order model: tube photons (continuum + anode lines) are pushed through
  # the Compton operator twice -- energy shift E -> E'(theta) with the Klein-Nishina angular pdf
  # (normalized per energy; the incoherent mass attenuation carries the total probability at each
  # step) -- with the direction CORRELATION between the two scatters ignored (angles drawn
  # independently), the two internal legs collapsed into the standard thick-sample kernel evaluated
  # at (flux-weighted mean source energy) -> E'', and the detector efficiency at E''. One more
  # fixed-shape, free-amplitude template ("scatter_compton_continuum2"); its ^scatter_compton prefix
  # keeps it inside the Compton normalizer denominator, where 2nd-order incoherent scatter belongs.
  if (isTRUE(second_order)) {
    th <- seq(0.15, pi - 0.15, length.out = 12L)
    # per-step angular weight: the exact differential KN(e, theta) * S(q(e, theta), Z) * sin(theta)
    # (unnormalized -- the template amplitude is free)
    kn_pdf <- function(e, t) {
      Pp <- 1 / (1 + (e / m_e_c2) * (1 - cos(t)))
      Pp^2 * (Pp + 1 / Pp - sin(t)^2) * sin(t) *
        .xrf_material_S(scatterer, .xrf_q_momentum_transfer(e, t))
    }
    cshift <- function(e, t) e / (1 + (e / m_e_c2) * (1 - cos(t)))
    src_E <- grid; src_w <- N * spacing
    tl <- tryCatch(xrf_tube_spectrum(tube, 1, discrete_lines = TRUE), error = function(e) NULL)
    if (!is.null(tl) && nrow(tl) > 0) { src_E <- c(src_E, tl$energy_kev); src_w <- c(src_w, tl$flux) }
    ok0 <- is.finite(src_w) & src_w > 0 & is.finite(src_E) & src_E > 0
    src_E <- src_E[ok0]; src_w <- src_w[ok0]
    if (length(src_E)) {
      two_pass <- function(eV, wV) {                     # one Compton step: returns c(E', w') pairs
        W <- outer(eV, th, kn_pdf)                        # exact differential weight per (E, theta)
        W[!is.finite(W) | W < 0] <- 0
        Ei <- outer(eV, th, cshift)
        Fi <- wV * W                                      # row-recycled source weights
        list(e = as.vector(Ei), w = as.vector(Fi))
      }
      s1 <- two_pass(src_E, src_w)
      k1 <- s1$w > 0 & is.finite(s1$e)
      s2 <- two_pass(s1$e[k1], s1$w[k1])
      bin <- as.integer(round((s2$e - grid[1]) / spacing)) + 1L
      okb <- is.finite(s2$w) & s2$w > 0 & bin >= 1L & bin <= length(grid)
      acc <- numeric(length(grid))
      if (any(okb)) {
        agg <- rowsum(s2$w[okb], bin[okb])
        acc[as.integer(rownames(agg))] <- as.vector(agg)
      }
      if (any(acc > 0)) {
        e_bar <- sum(src_E * src_w) / sum(src_w)
        rpi_o2 <- acc * .xrf_scatter_sample_kernel(scatterer, rep(e_bar, length(grid)), grid, geometry) *
          (if (!is.null(detector_type)) eff_at(grid) else rep(1, length(grid)))
        out <- dplyr::bind_rows(out, tibble::tibble(
          element = "scatter_compton_continuum2", trans = paste0("comcont2_", seq_along(grid)),
          energy_kev = grid, sigma = sigma_at(grid, compton_broadening),
          relative_peak_intensity = rpi_o2
        ))
      }
    }
  }

  out[is.finite(out$relative_peak_intensity) & out$relative_peak_intensity > 0 &
        is.finite(out$energy_kev) & is.finite(out$sigma), , drop = FALSE]
}

#' Infer the effective scatter angle (and Compton width) from a measured spectrum
#'
#' The tube's characteristic anode K-line scatters off the sample both coherently (Rayleigh, pinned to the
#' anode line energy \eqn{E_0}) and incoherently (Compton, shifted DOWN to
#' \eqn{E_0 / (1 + (E_0/510.999)(1 - \cos\theta))} and Doppler-broadened to ~2-3x the detector width). Inverting
#' the \emph{measured} Rayleigh-to-Compton separation \eqn{r = (E_{Ray} - E_{Comp})/E_{Ray}} gives the
#' EFFECTIVE scatter angle \eqn{\theta} of the real instrument geometry -- a free self-calibration from any
#' spectrum with visible anode scatter. Because both energies are read off the SAME (possibly miscalibrated)
#' axis, a common energy gain/offset error cancels; only \eqn{k = E_0/m_ec^2} uses the theoretical anode
#' energy (the incident photon energy is a calibration-independent physical constant). \eqn{\theta}
#' then pins the Compton peak position/width for the scatter templates (\link{xrf_scatter_peaks}) and the
#' downstream phantom gate. (It is not automatically wired into \link{xrf_quantify}, which takes explicit
#' \code{incidence_deg} / \code{takeoff_deg}; a single scatter angle only constrains their sum via
#' \eqn{\theta \approx 180 - \mathrm{incidence} - \mathrm{takeoff}}, not each separately.)
#'
#' Robust by construction: returns \code{confident = FALSE} (the caller should keep its prior angle) when
#' neither the anode K nor L lines are excited, the scatter band carries too few counts, or the candidate
#' Compton feature is too NARROW to be Compton -- i.e. it is a fluorescence /
#' contaminant line sitting in the window (the width test: its FWHM must be
#' \code{min_width_ratio}-\code{max_width_ratio}x the detector FWHM; e.g. an Hf K\eqn{\alpha} contaminant
#' landing on a high-kV W-Compton). The Rayleigh energy is a baseline-subtracted centroid in a tight window pinned to
#' \eqn{E_0} rather than a global maximum (the Compton peak is frequently the stronger of the two in low-Z
#' matrices), and the Compton centroid/width windows are kept a safe margin below the Rayleigh peak so its
#' tail cannot bias the angle.
#'
#' @param energy_kev,counts Measured spectrum (channel energies in keV and counts/CPS). Non-finite dropped.
#' @param anode Tube anode element symbol (e.g. "Rh", "Ag", "W").
#' @param tube_kv Tube voltage (kV). Above the anode K-edge the K\eqn{\alpha} scatter anchors the inference;
#'   below it (low-kV tubes) the strongest excited anode L line is used instead (less precise -- small shift).
#' @param detector_type Detector type for the width test (see \link{xrf_detector_sigma_kev}); default "SDD".
#' @param angle_prior_deg Fallback scatter angle returned when inference is not confident (default 135).
#' @param min_band_counts Minimum total signal in the scatter band required to attempt inference (default
#'   200). In the units of \code{counts} (raw counts or CPS) -- raise it for raw-count spectra, lower it for
#'   short-livetime / CPS data.
#' @param min_width_ratio Minimum Compton-FWHM / detector-FWHM to accept the feature as Compton (default 1.6).
#' @param max_width_ratio Maximum Compton-FWHM / detector-FWHM accepted before the feature is rejected as
#'   noise or merged peaks (default 4).
#' @return A list: \code{scatter_angle_deg}, \code{compton_broadening} (data-measured Compton width / detector
#'   width, clamped to [1,4]), \code{e_rayleigh}, \code{e_compton}, \code{e_anode} (the anode anchor line
#'   energy) and \code{anode_line} ("K" or "L"), \code{width_ratio}, \code{band_counts},
#'   \code{confident} (logical), and \code{reason}.
#' @examples
#' \dontrun{
#' g <- xrf_infer_scatter_geometry(spec$energy_kev, spec$cps, anode = "Rh", tube_kv = 40)
#' if (g$confident) geometry <- xrf_geometry(scatter_angle_deg = g$scatter_angle_deg)
#' }
#' @export
xrf_infer_scatter_geometry <- function(energy_kev, counts, anode, tube_kv,
                                       detector_type = "SDD", angle_prior_deg = 135,
                                       min_band_counts = 200, min_width_ratio = 1.6, max_width_ratio = 4) {
  out <- list(scatter_angle_deg = angle_prior_deg, compton_broadening = 2,
              e_rayleigh = NA_real_, e_compton = NA_real_, e_anode = NA_real_, anode_line = NA_character_,
              width_ratio = NA_real_, band_counts = NA_real_, confident = FALSE, reason = "")
  E <- suppressWarnings(as.numeric(energy_kev)); C <- suppressWarnings(as.numeric(counts))
  keep <- is.finite(E) & is.finite(C); E <- E[keep]; C <- C[keep]
  if (length(E) < 16) { out$reason <- "too few channels"; return(out) }
  en <- tryCatch(xrf_energies(anode, beam_energy_kev = max(tube_kv, 120)), error = function(e) NULL)
  if (is.null(en) || !nrow(en)) { out$reason <- "no anode lines"; return(out) }
  # Choose the anode scatter anchor line: the K-alpha centroid if the anode K shell is excited (Ka ~ 0.87 x
  # K-edge -> require a comfortable overvoltage), otherwise the strongest anode L line. The L fallback lets a
  # low-kV tube below its anode K-edge (e.g. a W tube at 40-50 kV, scattering off the W L-series) still self-
  # calibrate -- but L-line inference is inherently less precise (the Compton shift is small at ~8-12 keV and
  # sample fluorescence crowds the band), so the width/count gates reject it more readily; treat a returned
  # L-based theta as approximate.
  kl <- en[grepl("^K", en$trans_siegbahn) & is.finite(en$energy_kev) &
             is.finite(en$relative_peak_intensity), , drop = FALSE]
  E0k <- if (nrow(kl)) {
    d <- kl[kl$relative_peak_intensity >= 0.5 * max(kl$relative_peak_intensity), , drop = FALSE]  # Ka1,Ka2
    sum(d$energy_kev * d$relative_peak_intensity) / sum(d$relative_peak_intensity)
  } else NA_real_
  if (is.finite(E0k) && is.finite(tube_kv) && tube_kv >= E0k * 1.2) {
    E0 <- E0k; out$anode_line <- "K"
  } else {
    ll <- en[grepl("^L", en$trans_siegbahn) & is.finite(en$energy_kev) &
               is.finite(en$relative_peak_intensity), , drop = FALSE]
    if (!nrow(ll)) { out$reason <- "no excited anode K or L lines"; return(out) }
    sel <- which.max(ll$relative_peak_intensity)
    E0 <- ll$energy_kev[sel]
    E0edge <- if ("edge_kev" %in% names(ll)) suppressWarnings(as.numeric(ll$edge_kev[sel])) else NA_real_
    excited <- if (is.finite(E0edge)) isTRUE(tube_kv >= E0edge) else isTRUE(is.finite(tube_kv) && tube_kv >= E0 * 1.15)
    if (!excited) { out$reason <- "anode K not excited and L not excited (kv too low)"; return(out) }
    out$anode_line <- "L"
  }
  out$e_anode <- E0
  k <- E0 / 510.999; E180 <- E0 / (1 + 2 * k)
  bandw <- E >= E180 - 1 & E <= E0 + 1; out$band_counts <- sum(C[bandw], na.rm = TRUE)
  if (!is.finite(out$band_counts) || out$band_counts < min_band_counts) { out$reason <- "insufficient scatter counts"; return(out) }
  # detector width near the scatter band -- used for baseline-subtracted centroids and shift-aware windows
  detsig0 <- tryCatch(xrf_detector_sigma_kev(E0, detector_type), error = function(e) NA_real_)
  if (!is.finite(detsig0) || detsig0 <= 0) detsig0 <- 0.05 + 0.0025 * E0
  # baseline-subtracted intensity-weighted centroid over channel indices `idx`: subtract the window minimum
  # so the smooth scatter continuum under a peak does not pull the centroid; NA if there is no net signal.
  bcen <- function(idx) {
    if (length(idx) < 1L) return(NA_real_)
    b <- min(C[idx]); net <- pmax(C[idx] - b, 0)
    if (sum(net) <= 0) return(NA_real_)
    sum(E[idx] * net) / sum(net)
  }
  # Rayleigh peak, pinned to E0 (NOT the global max -- Compton is often stronger in low-Z matrices): its
  # baseline-subtracted centroid in a tight +/-0.5 keV window. This MEASURED energy -- not the theoretical E0
  # -- anchors the angle inversion below, so a common energy gain/offset miscalibration cancels (Rayleigh and
  # Compton are read off the same axis). Falls back to the local argmax if the window has no net signal.
  rw <- which(E >= E0 - 0.5 & E <= E0 + 0.5); if (!length(rw)) { out$reason <- "no Rayleigh window"; return(out) }
  E_ray <- bcen(rw); if (!is.finite(E_ray)) E_ray <- E[rw][which.max(C[rw])]
  out$e_rayleigh <- E_ray
  # Keep every Compton-side measurement a safe margin (3 detector sigma) BELOW the Rayleigh peak, so the
  # Rayleigh tail cannot masquerade as, or bias, the Compton feature (the old raw +/-1.2 keV centroid ingested
  # the Rayleigh tail and biased the angle systematically low, ~8-14 deg at 120-135 deg geometries).
  c_hi <- E_ray - 3 * detsig0
  # Compton search window: from ~180 deg (E180) up to the ~90 deg (least-shifted) Compton energy E0/(1+k),
  # clamped below Rayleigh. The shift-aware upper bound fixes the old fixed E0-0.6, which excluded the ~90 deg
  # Compton peak for light anodes (e.g. Mo at 90 deg: E' = 16.86 > E0-0.6 = 16.84).
  lo <- E180 - 0.8; hi <- min(E0 / (1 + k) + 2 * detsig0, c_hi); cw <- which(E >= lo & E <= hi)
  if (length(cw) < 3) { out$reason <- "no Compton window"; return(out) }
  Ecp <- E[cw][which.max(C[cw])]
  # Compton centroid: baseline-subtracted, window clamped below Rayleigh (c_hi).
  cc <- which(E >= Ecp - 1.2 & E <= min(Ecp + 1.2, c_hi)); Ecen <- bcen(cc)
  if (!is.finite(Ecen)) { out$reason <- "no Compton centroid"; return(out) }
  out$e_compton <- Ecen
  # width test on a Compton-scale window (also clamped below Rayleigh): reject a narrow fluorescence /
  # contaminant line sitting in the window.
  ww <- which(E >= Ecp - 0.9 & E <= min(Ecp + 0.9, c_hi))
  if (length(ww) < 3) { out$reason <- "no Compton width window"; return(out) }
  base <- min(C[ww]); net <- pmax(C[ww] - base, 0)
  cFWHM <- if (sum(net) > 0) { mu <- sum(E[ww] * net) / sum(net); 2.355 * sqrt(sum(net * (E[ww] - mu)^2) / sum(net)) } else NA_real_
  detsig <- tryCatch(xrf_detector_sigma_kev(Ecen, detector_type), error = function(e) NA_real_)
  if (!is.finite(detsig) || detsig <= 0) detsig <- 0.05 + 0.0025 * Ecen
  out$width_ratio <- cFWHM / (2.355 * detsig)
  if (!is.finite(out$width_ratio) || out$width_ratio < min_width_ratio) { out$reason <- "Compton candidate too narrow (likely an element line)"; return(out) }
  if (out$width_ratio > max_width_ratio) { out$reason <- "Compton region too broad (noise / merged features)"; return(out) }
  # Angle from the MEASURED Rayleigh->Compton shift r = (E_ray - Ecen)/E_ray (gain-cancelling); k = E0/mc^2
  # stays theoretical (the incident photon energy is a physical constant, independent of the spectrum's cal).
  r <- (E_ray - Ecen) / E_ray; cth <- max(-1, min(1, 1 - r / (k * (1 - r))))
  out$scatter_angle_deg <- acos(cth) * 180 / pi
  out$compton_broadening <- max(1, min(out$width_ratio, 4))
  out$confident <- TRUE; out$reason <- "ok"
  out
}
