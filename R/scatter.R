
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
#' The tube characteristic lines and their relative intensities are taken from \link{xrf_energies}
#' for the anode element at the tube voltage. Rayleigh and Compton are returned as two separate
#' "elements" (each a fixed-shape multi-line template with one free amplitude), because their ratio
#' varies with the sample (light matrices scatter more incoherently).
#'
#' @param tube An \link{xrf_tube} object (needs \code{anode} and \code{kv}).
#' @param geometry An \link{xrf_geometry} object (uses \code{scatter_angle_deg}).
#' @param detector_type,fano,epsilon_ev,noise_fwhm_ev Detector resolution parameters (see
#'   \link{xrf_detector_sigma_kev}) used for the scatter-peak widths. If none are supplied the
#'   \code{default_sigma} constant width is used.
#' @param default_sigma Fallback constant peak width (keV) when no detector model is given.
#' @param compton_broadening Multiplier on the detector width to approximate the Doppler-broadened
#'   Compton peak (default 2).
#' @param min_relative_intensity Smallest anode line to include (relative to the strongest).
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
                              compton_broadening = 2, min_relative_intensity = 0.05) {
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

  tube_lines <- xrf_energies(anode, beam_energy_kev = tube$kv,
                             min_relative_intensity = min_relative_intensity)
  if (nrow(tube_lines) == 0) return(empty)

  m_e_c2 <- 510.999                                  # electron rest energy (keV)
  theta <- geometry$scatter_angle_deg * pi / 180
  E <- tube_lines$energy_kev
  E_compton <- E / (1 + (E / m_e_c2) * (1 - cos(theta)))

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
      trans = paste0("rayleigh_", tube_lines$trans_siegbahn),
      energy_kev = E,
      sigma = sigma_at(E),
      relative_peak_intensity = tube_lines$relative_peak_intensity
    ),
    tibble::tibble(
      element = "scatter_compton",
      trans = paste0("compton_", tube_lines$trans_siegbahn),
      energy_kev = E_compton,
      sigma = sigma_at(E_compton) * compton_broadening,
      relative_peak_intensity = tube_lines$relative_peak_intensity
    )
  )
}
