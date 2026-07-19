
#' Energy-dependent detector resolution
#'
#' Semiconductor (EDS) detector resolution is not constant: the peak width grows with photon
#' energy because the Fano-limited statistical broadening adds in quadrature to the (energy
#' independent) electronic noise. The full width at half maximum is
#'
#' \deqn{FWHM(E) = \sqrt{ N^2 + (2.3548)^2 \, F \, \epsilon \, E }}
#'
#' where \eqn{N} is the electronic-noise FWHM (eV), \eqn{F} the Fano factor, \eqn{\epsilon} the
#' mean electron-hole pair creation energy (eV), and \eqn{E} the photon energy (eV). The Gaussian
#' standard deviation is \eqn{\sigma = FWHM / 2.3548}. A single constant \code{sigma} (the historic
#' default of 0.07 keV, i.e. ~165 eV FWHM) is only correct near ~8-10 keV; it is ~3x too wide at
#' the C K\eqn{\alpha} line (0.28 keV) and ~3x too narrow at the U K\eqn{\alpha} line (98 keV),
#' which is why it must vary with energy for wide-range (SEM-EDS to ~100 keV) work.
#'
#' \code{detector_type} selects a set of material defaults (see \link{xrf_detector_presets}); any of
#' \code{fano}, \code{epsilon_ev}, \code{noise_fwhm_ev} passed explicitly override the preset. If
#' \code{detector_type} is \code{NULL} but some scalars are supplied, an \code{"SDD"} (silicon) base
#' is assumed and the supplied scalars override it.
#'
#' @param energy_kev Photon energy (keV). May be a vector.
#' @param detector_type One of the names in \link{xrf_detector_presets} (e.g. "SDD", "SiLi",
#'   "SiPIN", "HPGe", "CdTe"). Case-insensitive.
#' @param fano Fano factor (dimensionless). \code{NULL} uses the preset.
#' @param epsilon_ev Mean electron-hole pair creation energy (eV). \code{NULL} uses the preset
#'   (Si 3.64, Ge 2.96, CdTe 4.43).
#' @param noise_fwhm_ev Electronic-noise contribution as a FWHM in eV (the width at zero energy).
#'   \code{NULL} uses the preset.
#'
#' @return \code{xrf_detector_fwhm_ev}: FWHM in eV. \code{xrf_detector_sigma_kev}: Gaussian sigma
#'   in keV (suitable for the \code{sigma} column consumed by
#'   \link{xrf_deconvolute_gaussian_least_squares}).
#' @export
#'
#' @examples
#' # resolution grows with energy
#' xrf_detector_fwhm_ev(c(0.28, 5.9, 59.3, 98.4), "SDD")
#' xrf_detector_sigma_kev(c(0.28, 5.9, 59.3, 98.4), "SDD")
#'
xrf_detector_sigma_kev <- function(energy_kev, detector_type = NULL, fano = NULL,
                                   epsilon_ev = NULL, noise_fwhm_ev = NULL) {
  xrf_detector_fwhm_ev(energy_kev, detector_type, fano, epsilon_ev, noise_fwhm_ev) / 2.3548 / 1000
}

#' @rdname xrf_detector_sigma_kev
#' @export
xrf_detector_fwhm_ev <- function(energy_kev, detector_type = NULL, fano = NULL,
                                 epsilon_ev = NULL, noise_fwhm_ev = NULL) {
  p <- resolve_detector_params(detector_type, fano, epsilon_ev, noise_fwhm_ev)
  e_ev <- energy_kev * 1000
  e_ev[!is.finite(e_ev) | e_ev < 0] <- 0
  sqrt(p$noise_fwhm_ev^2 + (2.3548^2) * p$fano * p$epsilon_ev * e_ev)
}

#' Detector material presets
#'
#' Physical constants used by the energy-dependent resolution model
#' (\link{xrf_detector_sigma_kev}) and, in future, the detector efficiency / escape-peak models.
#' These are sensible defaults only -- override any value with the explicit \code{fano},
#' \code{epsilon_ev}, or \code{noise_fwhm_ev} arguments, or seed \code{noise_fwhm_ev}/\code{fano}
#' from an instrument's own calibration (e.g. the \code{NoiseReference}/\code{FanoReference} fields
#' parsed by \link{read_xrf_panalytical}).
#'
#' The \code{"SDD-windowless"} preset (no Be window, thin dead layer) is for light-element (B/C/N/O)
#' SEM-EDS work: an 8 um Be window blocks sub-0.5 keV lines entirely, so pass a thin polymer or
#' silicon-nitride window instead (e.g. \code{window = "Si3N4 0.15"} or \code{"polypropylene 0.4"} to
#' \link{xrf_detector_efficiency}). \strong{Note on CdTe:} the FWHM here is the Fano-limited electronic
#' value only; real CdTe/CZT resolution is dominated by incomplete charge collection (hole tailing),
#' so the tabulated width is a best case -- enable the Hypermet tail in the deconvolution
#' (\link{xrf_lineshape}) and expect broader, low-energy-tailed peaks in practice.
#'
#' @return A tibble of presets: \code{type}, \code{material}, \code{fano}, \code{epsilon_ev},
#'   \code{noise_fwhm_ev}, \code{active_thickness_um}, \code{be_window_um}, \code{dead_layer_um}.
#' @export
#'
#' @examples
#' xrf_detector_presets()
#'
xrf_detector_presets <- function() {
  tibble::tribble(
    ~type,            ~material, ~fano,  ~epsilon_ev, ~noise_fwhm_ev, ~active_thickness_um, ~be_window_um, ~dead_layer_um,
    "SDD",            "Si",      0.115,  3.64,        50,             0.45e3,               8.0,           0.10,
    "SDD-windowless", "Si",      0.115,  3.64,        50,             0.45e3,               0.0,           0.05,
    "SiLi",           "Si",      0.115,  3.64,        80,             3.0e3,                12.0,          0.10,
    "SiPIN",          "Si",      0.115,  3.64,        90,             0.5e3,                12.5,          0.10,
    "Si",             "Si",      0.115,  3.64,        50,             0.45e3,               8.0,           0.10,
    "HPGe",           "Ge",      0.110,  2.96,        90,             5.0e3,                25.0,          0.30,
    "Ge",             "Ge",      0.110,  2.96,        90,             5.0e3,                25.0,          0.30,
    "CdTe",           "CdTe",    0.150,  4.43,        90,             1.0e3,                4.0,           0.10
  )
}

#' Energy-dependent detector line-shape tailing (incomplete charge collection)
#'
#' A real semiconductor-detector peak has a low-energy exponential \strong{tail} and a flat \strong{shelf}
#' (step) from incomplete charge collection (ICC), whose amplitudes vary with photon energy and detector
#' material -- unlike the single global \code{tail}/\code{step} scalars of
#' \link{xrf_deconvolute_gaussian_least_squares}. This returns a physically-motivated, \strong{tunable}
#' estimate of \code{(tail, step, beta)} versus energy that you can feed to the deconvolution (or
#' \link{xrf_lineshape}) instead of guessing a constant: for silicon the tail grows toward \emph{low}
#' energy (front-surface charge loss), while for germanium / CdTe the hole-tailing grows toward \emph{high}
#' energy. These are approximate defaults (real ICC is instrument-specific) -- scale them with
#' \code{tail_ref}/\code{step_ref} to match a measured line shape.
#'
#' @param energy_kev Photon energy (keV), vectorized.
#' @param detector_type A preset name (\link{xrf_detector_presets}); sets the active material.
#' @param tail_ref,step_ref Reference tail / shelf amplitudes (relative to the Gaussian peak) at the
#'   pivot energy. \code{NULL} uses a material default (Si ~0.01/0.005, Ge ~0.02/0.008, CdTe ~0.05/0.02).
#' @param fano,epsilon_ev,noise_fwhm_ev Detector resolution overrides for the returned \code{beta}
#'   (defaults to the Gaussian sigma at each energy).
#'
#' @return A list with numeric vectors \code{tail}, \code{step}, \code{beta} (keV), one per
#'   \code{energy_kev}.
#' @export
#'
#' @examples
#' xrf_detector_tailing(c(2, 6, 20), "SDD")
#' xrf_detector_tailing(c(20, 60, 100), "CdTe")   # hole-tailing grows with energy
#'
xrf_detector_tailing <- function(energy_kev, detector_type = NULL, tail_ref = NULL, step_ref = NULL,
                                 fano = NULL, epsilon_ev = NULL, noise_fwhm_ev = NULL) {
  material <- resolve_detector_params(detector_type)$material
  # material-family defaults (relative amplitude at the pivot energy) and the energy trend of the tail
  d <- switch(material,
    "Si"   = list(tail0 = 0.010, step0 = 0.005, pivot = 6,  trend = "low"),   # surface loss: grows at low E
    "Ge"   = list(tail0 = 0.020, step0 = 0.008, pivot = 60, trend = "high"),  # hole tailing: grows at high E
    "CdTe" = list(tail0 = 0.050, step0 = 0.020, pivot = 60, trend = "high"),
    list(tail0 = 0.010, step0 = 0.005, pivot = 6, trend = "low")
  )
  if (!is.null(tail_ref)) d$tail0 <- tail_ref
  if (!is.null(step_ref)) d$step0 <- step_ref
  e <- pmax(energy_kev, 0.1)
  f <- if (d$trend == "low") pmin(d$pivot / e, 5) else pmin(e / d$pivot, 5)  # bounded energy trend
  beta <- xrf_detector_sigma_kev(energy_kev, detector_type = detector_type, fano = fano,
                                 epsilon_ev = epsilon_ev, noise_fwhm_ev = noise_fwhm_ev)
  # both the tail and the shelf (step) arise from incomplete charge collection, so give the shelf the same
  # energy trend as the tail (rather than holding it constant) for internal consistency.
  list(tail = d$tail0 * f, step = d$step0 * f, beta = beta)
}

# Resolve (detector_type + explicit overrides) -> concrete (fano, epsilon_ev, noise_fwhm_ev).
# Not exported.
resolve_detector_params <- function(detector_type = NULL, fano = NULL, epsilon_ev = NULL,
                                     noise_fwhm_ev = NULL) {
  presets <- xrf_detector_presets()
  key <- if (is.null(detector_type)) "SDD" else detector_type
  idx <- match(tolower(key), tolower(presets$type))
  if (is.na(idx)) {
    stop("Unknown detector_type '", key, "'. Known types: ",
         paste(presets$type, collapse = ", "))
  }
  base <- presets[idx, ]
  out <- list(
    fano          = if (is.null(fano)) base$fano else fano,
    epsilon_ev    = if (is.null(epsilon_ev)) base$epsilon_ev else epsilon_ev,
    noise_fwhm_ev = if (is.null(noise_fwhm_ev)) base$noise_fwhm_ev else noise_fwhm_ev,
    material      = base$material
  )
  stopifnot(
    is.numeric(out$fano), out$fano > 0,
    is.numeric(out$epsilon_ev), out$epsilon_ev > 0,
    is.numeric(out$noise_fwhm_ev), out$noise_fwhm_ev >= 0
  )
  out
}

#' X-ray tube description (reserved)
#'
#' A lightweight container describing the excitation source (X-ray tube). It is carried through the
#' analysis so that a single object can hold all instrument conditions, but it is \strong{not yet
#' consumed by the deconvolution fit} -- it is reserved for the forthcoming tube-spectrum,
#' Rayleigh/Compton scatter, and fundamental-parameters work. Supplying it today is harmless and
#' future-proofs your scripts.
#'
#' @param anode Anode element symbol (e.g. "Rh", "W", "Ag", "Mo", "Au").
#' @param kv Tube accelerating voltage (kV); the bremsstrahlung endpoint energy in keV.
#' @param ma Tube current (mA), optional.
#' @param filter Primary beam filter description, optional.
#' @param ... Additional named fields to carry.
#'
#' @return An \code{xrf_tube} object (a classed list).
#' @export
#'
#' @examples
#' xrf_tube("Rh", kv = 40)
#'
xrf_tube <- function(anode = NA_character_, kv = NA_real_, ma = NA_real_,
                     filter = NA_character_, ...) {
  stopifnot(length(anode) == 1, length(kv) == 1)
  structure(
    list(anode = as.character(anode), kv = as.numeric(kv), ma = as.numeric(ma),
         filter = filter, ...),
    class = "xrf_tube"
  )
}

#' Measurement geometry (reserved)
#'
#' A lightweight container describing the source-sample-detector geometry. Like \link{xrf_tube} it is
#' carried through the analysis but \strong{not yet consumed by the deconvolution fit}; it is
#' reserved for the scatter (the Compton shift depends on \code{scatter_angle_deg}) and
#' fundamental-parameters (path lengths depend on the incidence/takeoff angles) work.
#'
#' @param incidence_deg Angle of the incident beam to the sample surface (degrees).
#' @param takeoff_deg Angle of the emergent (detected) beam to the sample surface (degrees).
#' @param scatter_angle_deg Source-sample-detector scattering angle (degrees); sets the Compton
#'   energy shift \eqn{E' = E / (1 + (E/511)(1 - \cos\theta))}.
#' @param ... Additional named fields to carry.
#'
#' @return An \code{xrf_geometry} object (a classed list).
#' @export
#'
#' @examples
#' xrf_geometry(incidence_deg = 45, takeoff_deg = 45, scatter_angle_deg = 135)
#'
xrf_geometry <- function(incidence_deg = 45, takeoff_deg = 45, scatter_angle_deg = 135, ...) {
  stopifnot(
    is.numeric(incidence_deg), is.numeric(takeoff_deg), is.numeric(scatter_angle_deg)
  )
  structure(
    list(incidence_deg = incidence_deg, takeoff_deg = takeoff_deg,
         scatter_angle_deg = scatter_angle_deg, ...),
    class = "xrf_geometry"
  )
}
