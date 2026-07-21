
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
#' line is weighted by the \code{scatterer}'s coherent (Rayleigh) or incoherent (Compton, with the
#' fixed-angle Klein-Nishina factor) scattering cross section at its energy -- so the template's
#' internal line balance (e.g. anode L vs K components) follows the physical emitted-and-scattered
#' intensity, not the anode's fluorescence-production ratios. Rayleigh and Compton are returned as
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
#' @param compton_broadening Multiplier on the detector width to approximate the Doppler-broadened
#'   Compton peak (default 2).
#' @param min_relative_intensity Smallest scatter line to include (relative to the strongest of
#'   either component; the Rayleigh and Compton templates keep a common line set).
#' @param scatterer Representative sample scatterer (element symbol or supported material) whose
#'   coherent/incoherent cross sections weight the line intensities; default \code{"Si"}, matching
#'   \link{xrf_scatter_continuum}.
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
                              scatterer = "Si") {
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

  sig_coh <- .xrf_material_scatter_mu(scatterer, E, "rayleigh"); sig_coh[!is.finite(sig_coh)] <- 0
  sig_inc <- .xrf_material_scatter_mu(scatterer, E, "compton");  sig_inc[!is.finite(sig_inc)] <- 0
  # fixed-angle Klein-Nishina reshaping, normalized by the angle-integrated KN total already inside
  # mu_inc -- identical composition to xrf_scatter_continuum, so line and continuum scatter agree.
  P <- E_compton / E
  eps <- pmax(E / m_e_c2, 1e-6)
  kn_tot <- 0.75 * (((1 + eps) / eps^2) * (2 * (1 + eps) / (1 + 2 * eps) - log(1 + 2 * eps) / eps) +
                      log(1 + 2 * eps) / (2 * eps) - (1 + 3 * eps) / (1 + 2 * eps)^2)
  kn <- (P^2 * (P + 1 / P - sin(theta)^2) / (1 + cos(theta)^2)) / kn_tot
  kn[!is.finite(kn) | kn < 0] <- 0

  rpi_ray <- tube_lines$flux * sig_coh * .xrf_scatter_sample_kernel(scatterer, E, E, geometry)
  rpi_com <- tube_lines$flux * sig_inc * kn * .xrf_scatter_sample_kernel(scatterer, E, E_compton, geometry)
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

  dplyr::bind_rows(
    tibble::tibble(
      element = "scatter_rayleigh",
      trans = paste0("rayleigh_", round(E, 3)),
      energy_kev = E,
      sigma = sigma_at(E),
      relative_peak_intensity = rpi_ray
    ),
    tibble::tibble(
      element = "scatter_compton",
      trans = paste0("compton_", round(E, 3)),
      energy_kev = E_compton,
      sigma = sigma_at(E_compton) * compton_broadening,
      relative_peak_intensity = rpi_com
    )
  )
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
#'   coherent/incoherent mass-attenuation energy dependence sets the continuum \emph{shape} (the overall
#'   amplitude is free in the fit). Defaults to \code{"Si"} (a light-matrix stand-in).
#' @param n_grid Number of energies sampled across the continuum (pseudo-lines per scatter type).
#' @param energy_min_kev Low-energy start of the sampled continuum (keV).
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
                                  energy_min_kev = 1,
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
  sig_coh <- .xrf_material_scatter_mu(scatterer, grid, "rayleigh")
  sig_inc <- .xrf_material_scatter_mu(scatterer, grid, "compton")
  sig_coh[!is.finite(sig_coh)] <- 0; sig_inc[!is.finite(sig_inc)] <- 0

  theta <- geometry$scatter_angle_deg * pi / 180
  m_e_c2 <- 510.999
  E_compton <- grid / (1 + (grid / m_e_c2) * (1 - cos(theta)))

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
  #  - Compton: mu_inc(E) already carries the ANGLE-INTEGRATED Klein-Nishina energy falloff (times the
  #    incoherent scattering function), so multiplying by the raw fixed-angle KN differential would count
  #    the KN energy dependence twice. The correct fixed-angle reshaping replaces the angle-integrated KN
  #    with the fixed-angle one: multiply by [dsigma_KN/dOmega(E, theta) / sigma_KN_tot(E)], normalized by
  #    its Thomson limit [ (1+cos^2)/2 / sigma_T ] so the factor -> 1 at low energy. (The residual is the
  #    fixed-angle vs angle-integrated difference of the incoherent function S(q, Z), whose table is not
  #    carried; the co-fit background absorbs that smooth remainder.)
  #  - Rayleigh's angular factor (1+cos^2) is energy-independent at fixed angle, so it does not reshape the
  #    coherent continuum; its q-dependent atomic form factor F(q,Z) would, but that table is not carried
  #    here (the co-fit background absorbs the smooth residual).
  #  - detector efficiency at the DETECTED energy (grid for Rayleigh, E' for Compton) shapes both continua --
  #    notably the low-energy window cutoff and the high-energy active-layer falloff -- when a preset is known.
  P <- E_compton / grid
  eps <- pmax(grid / m_e_c2, 1e-6)
  # total Klein-Nishina cross section relative to Thomson: sigma_KN(eps)/sigma_T = (3/4) * [ ... ] -> 1 - 2 eps
  kn_tot <- 0.75 * (((1 + eps) / eps^2) * (2 * (1 + eps) / (1 + 2 * eps) - log(1 + 2 * eps) / eps) +
                      log(1 + 2 * eps) / (2 * eps) - (1 + 3 * eps) / (1 + 2 * eps)^2)
  kn <- (P^2 * (P + 1 / P - sin(theta)^2) / (1 + cos(theta)^2)) / kn_tot
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
                   relative_peak_intensity = N * sig_coh * ker_ray * eff_ray * spacing),
    tibble::tibble(element = "scatter_compton_continuum", trans = paste0("comcont_", seq_along(grid)),
                   energy_kev = E_compton, sigma = s_com,
                   relative_peak_intensity = N * sig_inc * kn * ker_com * eff_com * spacing)
  )
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
