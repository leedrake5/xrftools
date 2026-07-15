
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
  if (m %in% all_elements) {
    return(list(density = unname(.xrf_element_density[m]), composition = stats::setNames(1, m)))
  }
  stop("Unknown material '", m, "'. Use an element symbol or the compound CdTe.")
}

# Total photoelectric mass attenuation (cm^2/g) of a single element vs energy, with log-log
# interpolation of the packaged cross-section grid and power-law extrapolation beyond it. NOTE: the
# grid starts at 5 keV, so values below ~5 keV are power-law extrapolated and only approximate.
# log-log interpolate a single element's mu/rho grid (photoelectric or total), with flat-clamp below
# the grid and power-law extrapolation above it. The EPDL grid runs from ~5 eV to 200 keV, so the
# clamp/extrapolation are essentially never reached for real analytical lines.
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
#' @param energy_kev Photon energy (keV), vectorized.
#' @param detector_type A preset name (\link{xrf_detector_presets}); sets the material and default
#'   geometry.
#' @param active_material,active_thickness_um Active-volume material and thickness (microns).
#' @param be_window_um Beryllium window thickness (microns).
#' @param dead_layer_um Dead-layer thickness (microns), assumed to be the active material.
#'
#' @return Efficiency in [0, 1], same length as \code{energy_kev}.
#' @export
#'
#' @examples
#' # thin Si becomes transparent at high energy; a Be window cuts low energy
#' xrf_detector_efficiency(c(0.5, 6, 30, 60, 100), "SDD")
#' xrf_detector_efficiency(c(30, 60, 100), "CdTe")   # CdTe stays efficient at high energy
#'
xrf_detector_efficiency <- function(energy_kev, detector_type = NULL, active_material = NULL,
                                    active_thickness_um = NULL, be_window_um = NULL,
                                    dead_layer_um = NULL) {
  geom <- .resolve_detector_geometry(detector_type, active_material, active_thickness_um,
                                     be_window_um, dead_layer_um)
  cm_per_um <- 1e-4

  # window and dead layer attenuate by ANY interaction -> total mu; the active volume records only
  # photoelectric events in the full-energy peak -> photoelectric mu.
  win_rho <- .xrf_material_info("Be")$density
  t_win <- exp(-xrf_mass_attenuation("Be", energy_kev, type = "total") * win_rho *
                 geom$be_window_um * cm_per_um)

  dead_rho <- .xrf_material_info(geom$active_material)$density
  t_dead <- exp(-xrf_mass_attenuation(geom$active_material, energy_kev, type = "total") * dead_rho *
                  geom$dead_layer_um * cm_per_um)

  act_rho <- .xrf_material_info(geom$active_material)$density
  a_active <- 1 - exp(-xrf_mass_attenuation(geom$active_material, energy_kev, type = "photoelectric") *
                        act_rho * geom$active_thickness_um * cm_per_um)

  pmin(pmax(t_win * t_dead * a_active, 0), 1)
}
