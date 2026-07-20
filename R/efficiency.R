
# Elemental densities (g/cm^3) for materials used as detectors, windows, anodes and beam filters.
.xrf_element_density <- c(
  Be = 1.848, C = 2.100, Na = 0.971, Mg = 1.738, Al = 2.700, Si = 2.330, P = 1.823, S = 2.067,
  Ti = 4.506, V = 6.11, Cr = 7.19, Mn = 7.47, Fe = 7.874, Co = 8.90, Ni = 8.908, Cu = 8.96,
  Zn = 7.14, Ga = 5.91, Ge = 5.323, Se = 4.81, Zr = 6.52, Nb = 8.57, Mo = 10.28, Rh = 12.41,
  Pd = 12.02, Ag = 10.49, Cd = 8.65, In = 7.31, Sn = 7.31, Sb = 6.68, Te = 6.24, W = 19.25,
  Ta = 16.65, Hf = 13.31, Re = 21.02, Pt = 21.45, Au = 19.30, Hg = 13.53, Pb = 11.34, Bi = 9.78,
  U = 19.05
)

# Material properties used by the efficiency, escape-peak and tube-spectrum models. Densities in
# g/cm^3; composition is a named vector of mass fractions (elements sum to 1). Any element symbol is
# accepted; a few compounds are special-cased.
.xrf_material_info <- function(material) {
  m <- as.character(material)
  cmp <- switch(
    m,
    "CdTe" = list(density = 5.850, composition = c(Cd = 0.4684, Te = 0.5316)),
    NULL
  )
  if (!is.null(cmp)) return(cmp)
  # Common measurement-path materials: the atmosphere (air / helium) and thin polymer snout/detector
  # windows. Compositions are element mass fractions; densities in g/cm^3. Matched case-insensitively with
  # a few aliases. (Air drops negligible Kr; polymers include H.)
  cmp <- switch(
    tolower(m),
    "air"                        = list(density = 0.0012048, composition = c(N = 0.7553, O = 0.2318, Ar = 0.0128, C = 0.000124)),
    "he" = , "helium"           = list(density = 0.0001786, composition = c(He = 1.0)),
    "vacuum" = , "none"         = list(density = 0,         composition = c(He = 1.0)),   # no attenuation
    "polypropylene" = , "pp" = , "prolene" = list(density = 0.905, composition = c(C = 0.8563, H = 0.1437)),
    "kapton" = , "polyimide"    = list(density = 1.42,  composition = c(C = 0.6911, H = 0.0264, N = 0.0733, O = 0.2092)),
    "mylar" = , "pet"           = list(density = 1.40,  composition = c(C = 0.625017, H = 0.041959, O = 0.333025)),
    # silicon nitride: common thin, light-element-transparent detector window (Moxtek AP / SiN membranes)
    "si3n4" = , "sin" = , "silicon nitride" = list(density = 3.17, composition = c(Si = 0.6006, N = 0.3994)),
    NULL
  )
  if (!is.null(cmp)) return(cmp)
  if (m %in% all_elements) {
    return(list(density = unname(.xrf_element_density[m]), composition = stats::setNames(1, m)))
  }
  stop("Unknown material '", m, "'. Use an element symbol, CdTe, or a known path material ",
       "(Air, He, polypropylene, Kapton, Mylar, Si3N4).")
}

# log-log interpolate a single element's mu/rho grid (photoelectric or total), with flat-clamp below
# the grid and power-law extrapolation above it. The EPDL grid runs from ~5 eV to 200 keV, so the
# clamp/extrapolation are essentially never reached for real analytical lines (sub-keV light-element
# lines are interpolated directly from tabulated values, not extrapolated).
.xrf_interp_mu <- function(g, energy_kev) {
  if (is.null(g) || nrow(g) < 2) return(rep(NA_real_, length(energy_kev)))
  lo <- min(g$energy_kev); hi <- max(g$energy_kev); m <- nrow(g)
  e <- energy_kev
  out <- rep(NA_real_, length(e))
  within <- is.finite(e) & e >= lo & e <= hi
  if (any(within)) {
    out[within] <- exp(suppressWarnings(   # edge-doubled grid points are collapsed by approx()
      stats::approx(log(g$energy_kev), log(g$mu), log(e[within]))$y))
  }
  below <- is.finite(e) & e < lo
  if (any(below)) out[below] <- g$mu[1]
  above <- is.finite(e) & e > hi
  if (any(above)) {
    slope <- (log(g$mu[m]) - log(g$mu[m - 1])) / (log(g$energy_kev[m]) - log(g$energy_kev[m - 1]))
    out[above] <- exp(log(g$mu[m]) + slope * (log(e[above]) - log(g$energy_kev[m])))
  }
  out
}

# photoelectric / total mass attenuation of one element vs energy (cm^2/g)
.xrf_element_pe_mu <- function(element, energy_kev) {
  if (is.numeric(element)) element <- all_elements[element]
  .init_physics_cache()
  .xrf_interp_mu(.xrf_cache$pe_mu_split[[element]], energy_kev)
}
.xrf_element_total_mu <- function(element, energy_kev) {
  if (is.numeric(element)) element <- all_elements[element]
  .init_physics_cache()
  .xrf_interp_mu(.xrf_cache$total_mu_split[[element]], energy_kev)
}
# coherent (Rayleigh) and incoherent (Compton) scatter mass attenuation (cm^2/g) of one element vs energy
.xrf_element_rayleigh_mu <- function(element, energy_kev) {
  if (is.numeric(element)) element <- all_elements[element]
  .init_physics_cache()
  .xrf_interp_mu(.xrf_cache$rayleigh_mu_split[[element]], energy_kev)
}
.xrf_element_compton_mu <- function(element, energy_kev) {
  if (is.numeric(element)) element <- all_elements[element]
  .init_physics_cache()
  .xrf_interp_mu(.xrf_cache$compton_mu_split[[element]], energy_kev)
}
# mass-weighted scatter mass attenuation of a material (element / compound / named composition)
.xrf_material_scatter_mu <- function(material, energy_kev, type = c("rayleigh", "compton")) {
  type <- match.arg(type)
  fn <- if (type == "rayleigh") .xrf_element_rayleigh_mu else .xrf_element_compton_mu
  info <- .xrf_material_info(material)
  mu <- rep(0, length(energy_kev))
  for (el in names(info$composition)) {
    m <- fn(el, energy_kev); m[!is.finite(m)] <- 0
    mu <- mu + info$composition[[el]] * m
  }
  mu
}

#' Mass attenuation of a material
#'
#' The mass attenuation coefficient \eqn{\mu/\rho} (cm^2/g) of an element or simple compound, from
#' the EPDL97 \link{x_ray_mass_attenuation} table (mass-weighted over compound constituents). Used by
#' \link{xrf_detector_efficiency} and \link{xrf_escape_peaks}. \code{type = "photoelectric"} (default)
#' returns only the photoelectric component -- the one that feeds the full-energy peak (detector
#' absorption) -- while \code{type = "total"} includes coherent + incoherent scatter, appropriate for
#' window / dead-layer transmission.
#'
#' @param material An element symbol ("Si", "Ge", "Be", "Au", "C") or a supported compound ("CdTe").
#' @param energy_kev Photon energy (keV), vectorized.
#' @param type "photoelectric" (default) or "total".
#'
#' @return \eqn{\mu/\rho} in cm^2/g.
#' @export
#'
#' @examples
#' xrf_mass_attenuation("Si", c(2, 10, 60))
#' xrf_mass_attenuation("CdTe", c(30, 60, 100), type = "total")
#'
xrf_mass_attenuation <- function(material, energy_kev, type = c("photoelectric", "total")) {
  type <- match.arg(type)
  info <- .xrf_material_info(material)
  fn <- if (type == "total") .xrf_element_total_mu else .xrf_element_pe_mu
  mu <- rep(0, length(energy_kev))
  for (el in names(info$composition)) {
    mu <- mu + info$composition[[el]] * fn(el, energy_kev)
  }
  mu
}

# Resolve detector geometry (active material/thickness, window, dead layer) from a preset plus
# explicit overrides.
.resolve_detector_geometry <- function(detector_type = NULL, active_material = NULL,
                                        active_thickness_um = NULL, be_window_um = NULL,
                                        dead_layer_um = NULL) {
  presets <- xrf_detector_presets()
  key <- if (is.null(detector_type)) "SDD" else detector_type
  idx <- match(tolower(key), tolower(presets$type))
  if (is.na(idx)) stop("Unknown detector_type '", key, "'.")
  base <- presets[idx, ]
  list(
    active_material     = if (is.null(active_material)) base$material else active_material,
    active_thickness_um = if (is.null(active_thickness_um)) base$active_thickness_um else active_thickness_um,
    be_window_um        = if (is.null(be_window_um)) base$be_window_um else be_window_um,
    dead_layer_um       = if (is.null(dead_layer_um)) base$dead_layer_um else dead_layer_um
  )
}

#' Energy-dependent detector efficiency
#'
#' The fraction of photons of energy \code{energy_kev} that are recorded in the full-energy
#' (photo)peak, combining transmission through the (beryllium) window and the dead layer with
#' absorption in the active volume:
#' \deqn{eff(E) = e^{-(\mu/\rho)_{win}\rho_{win} t_{win}} \; e^{-(\mu/\rho)_{dead}\rho_{dead} t_{dead}}
#'   \; (1 - e^{-(\mu/\rho)_{active}\rho_{active} t_{active}}).}
#' The active-volume term uses the photoelectric \eqn{\mu/\rho} (only photoelectric absorption feeds
#' the full-energy peak) and captures the transparency of thin Si detectors at high energy (e.g. a
#' few-hundred-micron Si detector detects only a small fraction of 60-100 keV photons, so
#' heavy-element K-lines are strongly under-counted); the window and dead-layer terms use the total
#' \eqn{\mu/\rho} (any interaction removes the photon). The EPDL mass-attenuation grid runs from a
#' few eV to 200 keV, so the low-energy window absorption (which matters most for light elements) is
#' now computed directly rather than extrapolated. Geometry defaults come from the
#' \code{detector_type} preset (\link{xrf_detector_presets}); override any of them explicitly.
#'
#' By default the active-volume term is the total photoabsorption probability. Set
#' \code{full_energy_peak = TRUE} to instead return the \strong{full-energy-peak} fraction, i.e.
#' multiplied by \eqn{(1 - f_{escape}(E))} where \eqn{f_{escape}} is the detector-fluorescence escape
#' probability (\link{xrf_escape_peaks}); this is the quantity that predicts the measured photopeak area
#' and is what \link{xrf_fp_sensitivity} uses. The difference is negligible for Si (escape ~1\%) but
#' reaches ~10-14\% for Ge/CdTe just above their K edges.
#'
#' @param energy_kev Photon energy (keV), vectorized.
#' @param detector_type A preset name (\link{xrf_detector_presets}); sets the material and default
#'   geometry.
#' @param active_material,active_thickness_um Active-volume material and thickness (microns).
#' @param be_window_um Beryllium window thickness (microns).
#' @param dead_layer_um Dead-layer thickness (microns), assumed to be the active material.
#' @param full_energy_peak If TRUE, return the full-energy-(photo)peak fraction, i.e. multiply the active
#'   term by \eqn{(1 - f_{escape})}. Default FALSE returns the total photoabsorption. See Details.
#' @param active_absorption Active-volume absorption model: "photoelectric" (default; single-interaction
#'   full-peak, a lower bound for thick detectors at high energy) or "total" (any interaction, an upper
#'   bound assuming full containment). See Details. Note: combining \code{active_absorption = "total"} with
#'   \code{full_energy_peak = TRUE} is approximate -- the escape correction \eqn{(1-f_{escape})} is derived
#'   for photoelectric absorption, so applying it to the total (photoelectric + scatter) interaction
#'   probability over-applies escape to the scatter events; use "photoelectric" for the self-consistent
#'   full-energy-peak fraction.
#'
#' @return Efficiency in [0, 1], same length as \code{energy_kev}.
#' @export
#'
#' @examples
#' # thin Si becomes transparent at high energy; a Be window cuts low energy
#' xrf_detector_efficiency(c(0.5, 6, 30, 60, 100), "SDD")
#' xrf_detector_efficiency(c(30, 60, 100), "CdTe")   # CdTe stays efficient at high energy
#'
#' @param air_path_cm,atmosphere,window Measurement path for the efficiency model: air/He path length (cm), atmosphere ("Air"/"He"/"Vacuum") and an optional snout window; NULL / "Air" reproduce the historical no-path behaviour.
xrf_detector_efficiency <- function(energy_kev, detector_type = NULL, active_material = NULL,
                                    active_thickness_um = NULL, be_window_um = NULL,
                                    dead_layer_um = NULL,
                                    air_path_cm = NULL, atmosphere = "Air", window = NULL,
                                    full_energy_peak = FALSE,
                                    active_absorption = c("photoelectric", "total")) {
  active_absorption <- match.arg(active_absorption)
  geom <- .resolve_detector_geometry(detector_type, active_material, active_thickness_um,
                                     be_window_um, dead_layer_um)
  cm_per_um <- 1e-4

  # Measurement-path attenuation between the sample and the detector's Be window: the atmosphere over the
  # air path (air, or ~transparent helium), and any polymer snout/detector window(s) such as a 4 um
  # polypropylene foil ("window" is a filter-style spec in microns, and may be a stack "Kapton 8;
  # polypropylene 4"). Both mainly cut low-energy lines, on top of the Be window / dead layer / active
  # absorption below.
  t_path <- rep(1, length(energy_kev))
  if (!is.null(air_path_cm) && is.finite(air_path_cm) && air_path_cm > 0) {
    atm <- if (is.null(atmosphere) || is.na(atmosphere) || !nzchar(trimws(as.character(atmosphere)))) "Air" else atmosphere
    ainfo <- tryCatch(.xrf_material_info(atm), error = function(e) NULL)
    if (!is.null(ainfo) && is.finite(ainfo$density) && ainfo$density > 0) {
      t_path <- t_path * exp(-xrf_mass_attenuation(atm, energy_kev, type = "total") * ainfo$density * air_path_cm)
    }
  }
  t_window <- .xrf_filter_transmission(window, energy_kev)   # 1 when no window given

  # window and dead layer attenuate by ANY interaction -> total mu; the active volume records only
  # photoelectric events in the full-energy peak -> photoelectric mu.
  win_rho <- .xrf_material_info("Be")$density
  t_win <- exp(-xrf_mass_attenuation("Be", energy_kev, type = "total") * win_rho *
                 geom$be_window_um * cm_per_um)

  dead_rho <- .xrf_material_info(geom$active_material)$density
  t_dead <- exp(-xrf_mass_attenuation(geom$active_material, energy_kev, type = "total") * dead_rho *
                  geom$dead_layer_um * cm_per_um)

  act_rho <- .xrf_material_info(geom$active_material)$density
  # Active-volume absorption. "photoelectric" (default) counts only photons that photo-absorb on their
  # FIRST interaction -- a LOWER bound on the full-energy-peak efficiency for a thick detector at high
  # energy, where a growing share of full-peak events are Compton-scatter(s)-then-photoabsorb (E8).
  # "total" counts every ENERGY-DEPOSITING interaction -- an UPPER bound, assuming multiple-interaction
  # events are fully contained in the peak (they are not: some scattered photons escape). We use
  # photoelectric + incoherent (Compton), i.e. total MINUS coherent (Rayleigh) scatter, since Rayleigh
  # deposits no energy and cannot contribute to the peak. The true full-peak efficiency lies between the two;
  # for thin Si they nearly coincide, for thick HPGe/CdTe above ~50 keV they bracket a real uncertainty.
  act_mu <- if (active_absorption == "total") {
    pmax(xrf_mass_attenuation(geom$active_material, energy_kev, type = "total") -
           .xrf_material_scatter_mu(geom$active_material, energy_kev, "rayleigh"), 0)
  } else {
    xrf_mass_attenuation(geom$active_material, energy_kev, type = "photoelectric")
  }
  a_active <- 1 - exp(-act_mu * act_rho * geom$active_thickness_um * cm_per_um)

  # Full-energy-peak fraction: remove the events whose detector characteristic X-ray escapes the crystal
  # (recorded in an escape peak, not the photopeak). Negligible for Si, ~10-14% for Ge/CdTe near their K edges.
  if (isTRUE(full_energy_peak)) {
    a_active <- a_active * (1 - .xrf_detector_escape_fraction(energy_kev, geom$active_material))
  }

  pmin(pmax(t_path * t_window * t_win * t_dead * a_active, 0), 1)
}
