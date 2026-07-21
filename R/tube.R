
# Ebel characteristic-line fluxes of the anode (photons/sr/mA/s), K + L series. K: Ebel (1999) with
# const_K = 5.0e13 (the value unified with the L extension; PyMCA's K-only implementation uses 6.0e13,
# a ~20% calibration spread), z_K = 2, b_K = 0.35. L: Ebel (2003) per-subshell constants 0.71e13 (L1),
# 2.70e13 (L2), 4.94e13 (L3), with the empirical fcorr(Z) on L1/L2 below Z = 80 (as in XMI-MSIM's
# xmi_ebel). Per line: flux = const_sub * (z b / Z)(U0 lnU0 + 1 - U0)[1 + 16.05 sqrt(J/E_edge)(...)] *
# R_backscatter(U0) * (omega p)_line * f_abs(E_line), with U0 = E0/E_edge of the line's SUBSHELL and
# the same Ebel depth model as the continuum evaluated at that overvoltage. We use the per-transition
# emission probabilities (omega_subshell * branching) where the reference multiplies RadRate by a mean
# omega_L -- more granular, and within the L model's ~10-20% calibration accuracy. Anode M lines are
# not parametrized by Ebel and are omitted (sub-2-keV; they vanish through any real window/filter).
.xrf_ebel_char_lines <- function(anode, Z, E0, eta, rhozmax, sinfac) {
  xe <- xrftools::x_ray_xrf_energies
  ln <- xe[xe$element == anode & xe$edge %in% c("K", "L1", "L2", "L3") &
             xe$edge_kev < E0 & xe$energy_kev < E0, , drop = FALSE]
  if (nrow(ln) == 0) return(NULL)
  wp <- xrf_transition_probability(anode, ln$trans)
  ok <- is.finite(wp) & wp > 0
  ln <- ln[ok, , drop = FALSE]; wp <- wp[ok]
  if (nrow(ln) == 0) return(NULL)
  u0 <- E0 / ln$edge_kev
  lu <- log(u0)
  oneovers <- (sqrt(u0) * lu + 2 * (1 - sqrt(u0))) / (u0 * lu + 1 - u0)   # stopping factor 1/S
  oneovers <- 1 + 16.05 * sqrt(0.0135 * Z / ln$edge_kev) * oneovers
  is_k <- ln$edge == "K"
  oneovers <- (ifelse(is_k, 2.0 * 0.35, 8.0 * 0.25) / Z) * (u0 * lu + 1 - u0) * oneovers
  Rbk <- 1 - 0.0081517 * Z + 3.613e-5 * Z^2 + 0.009583 * Z * exp(-u0) + 0.001141 * E0
  fc <- if (Z >= 80) 1 else max(-0.4814 + 0.03781 * Z - 2.413e-4 * Z^2, 0)
  const_s <- ifelse(is_k, 5.0e13,
                    ifelse(ln$edge == "L1", 0.71e13 * fc,
                           ifelse(ln$edge == "L2", 2.70e13 * fc, 4.94e13)))
  p1 <- lu * (0.49269 - 1.0987 * eta + 0.78557 * eta^2)                   # Ebel depth at the SUBSHELL U0
  p2 <- 0.70256 - 1.09865 * eta + 1.0046 * eta^2 + lu
  tau <- 2 * xrf_mass_attenuation(anode, ln$energy_kev, type = "photoelectric") * rhozmax * (p1 / p2) * sinfac
  f_abs <- ifelse(is.finite(tau) & tau > 1e-9, (1 - exp(-tau)) / tau, 1)
  out <- data.frame(energy_kev = ln$energy_kev, flux = const_s * oneovers * Rbk * wp * f_abs)
  out <- out[is.finite(out$flux) & out$flux > 0, , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  out[out$flux >= 1e-4 * max(out$flux), , drop = FALSE]                    # drop negligible satellites
}

# Transmission through a primary-beam filter, parsed from a tube$filter string. Accepts a single filter
# ("Al" or "Al 100") OR a stacked tri-filter like "Cu 100; Ti 25; Al 300" (as handhelds report): filters
# are separated by ';' or ',', each "Element thickness_um" (thickness in microns, default 100), and their
# transmissions multiply. Returns 1 (no attenuation) for any unparseable/absent spec.
.xrf_filter_transmission <- function(filter, energy_kev) {
  if (is.null(filter) || length(filter) != 1 || is.na(filter) || !nzchar(trimws(filter))) {
    return(rep(1, length(energy_kev)))
  }
  parts <- strsplit(trimws(as.character(filter)), "\\s*[;,]\\s*")[[1]]
  parts <- parts[nzchar(parts)]
  trans <- rep(1, length(energy_kev))
  for (p in parts) {
    toks <- strsplit(trimws(p), "\\s+")[[1]]
    mat <- toks[1]
    # Any material .xrf_material_info knows: an element, CdTe, or a polymer window (polypropylene/Kapton/...).
    info <- tryCatch(.xrf_material_info(mat), error = function(e) NULL)
    if (is.null(info) || !is.finite(info$density) || info$density <= 0) next   # unknown/vacuum -> no attenuation
    thick_um <- suppressWarnings(as.numeric(toks[2])); if (is.na(thick_um)) thick_um <- 100
    trans <- trans * exp(-xrf_mass_attenuation(mat, energy_kev, type = "total") * info$density * thick_um * 1e-4)
  }
  trans
}

#' X-ray tube emission spectrum (Ebel model)
#'
#' The photon output of an X-ray tube versus energy: a bremsstrahlung continuum (Ebel 1999) shaped by
#' self-absorption in the anode, plus the anode's characteristic lines, optionally attenuated by a
#' primary-beam filter. Values are relative (an overall constant is dropped) -- what the
#' fundamental-parameters excitation integral needs. Use it to weight excitation across the tube
#' spectrum instead of assuming a single monochromatic energy.
#'
#' The anode self-absorption follows Ebel's depth model: the mean generation depth
#' \eqn{\bar{\rho z}(E) = \bar{\rho z}_{max}\, p_1/p_2} is \emph{energy-dependent} (soft photons are
#' generated over the whole electron range, photons near the endpoint only at the surface, so
#' \eqn{\bar{\rho z} \to 0} as \eqn{E \to E_0}), with
#' \eqn{\bar{\rho z}_{max} = (A/Z)\,[0.787\times10^{-5}\sqrt{J}\,E_0^{1.5} + 0.735\times10^{-6}E_0^2]},
#' \eqn{J = 0.0135\,Z} keV, and \eqn{p_1/p_2} a function of \eqn{\ln(E_0/E)} and the electron
#' backscatter factor \eqn{\eta(Z, E_0)} (Ebel's Eq. 4 constants). The escape factor is
#' \eqn{f = (1 - e^{-\tau})/\tau} with
#' \eqn{\tau = 2\,\mu_{pe}(E)\,\bar{\rho z}\,\sin\alpha_e/\sin\alpha_x}, using the \emph{photoelectric}
#' mass attenuation (as in Ebel; coherent/incoherent scatter merely redirects within the anode).
#' Ebel's stated validity is 5-50 kV tube voltage.
#'
#' \strong{Model the tube exit window as part of \code{filter}.} The Ebel model emits the spectrum at
#' the anode surface; a real tube adds an exit window (commonly 75-300 um Be) that strongly cuts the
#' soft continuum and the anode L series. Include it with any beam filters in \code{tube$filter},
#' e.g. \code{xrf_tube("Rh", 40, filter = "Be 125; Al 100")} -- otherwise the unwindowed L-line /
#' soft-end flux of this model is an overestimate of what reaches the sample.
#'
#' @param tube An \link{xrf_tube} (needs \code{anode} and \code{kv}; \code{filter} optional).
#' @param energy_kev Photon energies to evaluate (keV).
#' @param takeoff_deg X-ray take-off angle \eqn{\alpha_x} (degrees) between the anode surface and the
#'   emerging beam, for the anode self-absorption path.
#' @param electron_incidence_deg Electron-beam incidence angle \eqn{\alpha_e} (degrees) between the
#'   anode surface and the electron beam. The default 90 (normal incidence) reproduces the historical
#'   geometry; side-window tubes are typically ~60-75.
#' @param char_peak_ratio Multiplier on the Ebel characteristic-line intensities (default 1 = Ebel's
#'   absolute K/L line yields, so the line/continuum balance is physical; 0 removes the anode lines,
#'   e.g. when lines and continuum are handled separately). Historical note: this argument previously
#'   set an arbitrary "line height over local continuum" ratio (default 100) -- the Ebel
#'   characteristic model (stopping power, backscatter, depth absorption; see
#'   \code{.xrf_ebel_char_lines}) replaces that heuristic, so values other than 0/1 are now only a
#'   deliberate rescaling of the physical line intensities.
#' @param filter Apply the tube's primary-beam \code{filter} (if any)?
#'
#' @return Photon flux at each \code{energy_kev} (0 at or above the endpoint \code{kv}), on Ebel's
#'   absolute scale: the continuum in photons/(sr mA keV s) and the rendered anode lines carrying
#'   their absolute integrated flux in photons/(sr mA s) (with \code{discrete_lines = TRUE}, that
#'   flux is returned directly). Most consumers use only ratios, but the absolute scale makes the
#'   line/continuum balance -- and any future absolute-intensity work -- physical.
#' @export
#'
#' @references
#' Ebel, H. (1999). X-ray tube spectra. X-Ray Spectrometry 28(4), 255-266.
#'
#' Ebel, H. (2003). X-ray tube spectra: fundamental parameter algorithms for the description of
#' K and L spectra (per-subshell L constants as implemented in XMI-MSIM).
#'
#' @examples
#' e <- seq(1, 40, 0.05)
#' plot(e, xrf_tube_spectrum(xrf_tube("Rh", kv = 40), e), type = "l")
#'
#' @param discrete_lines If TRUE, return the anode characteristic lines as discrete lines rather than rendered Gaussians (default FALSE).
xrf_tube_spectrum <- function(tube, energy_kev, takeoff_deg = 40, electron_incidence_deg = 90,
                              char_peak_ratio = 1,
                              filter = TRUE, discrete_lines = FALSE) {
  stopifnot(inherits(tube, "xrf_tube"))
  anode <- strsplit(trimws(as.character(tube$anode)), "\\s+")[[1]][1]
  if (!(anode %in% all_elements)) stop("tube anode '", tube$anode, "' is not a recognized element.")
  Z <- match(anode, all_elements)
  E0 <- tube$kv
  E <- energy_kev

  # Ebel bremsstrahlung PHOTON-number continuum ~ Z * (E0/E - 1)^x (at x=1 this is exactly the
  # Kramers photon spectrum (E0-E)/E). Do NOT divide by E again -- an extra /E would give an
  # intensity/E^2 shape and over-weight soft photons, and the excitation integral needs photons/keV.
  x <- 1.109 - 0.00435 * Z + 0.00175 * E0

  # Anode self-absorption, full Ebel (1999) depth model (see Details). An earlier version used a
  # single, energy-independent depth 0.787e-5*sqrt(0.0115*Z)*E0^1.5 + 0.735e-6*E0^2 -- i.e. the
  # MAXIMUM depth without the A/Z prefactor, the p1/p2 energy dependence, or Ebel's J = 0.0135*Z --
  # which over-absorbed the soft continuum ~2x below ~6 keV and under-absorbed ~10% at mid energies
  # (verified against the reference implementation). Every polychromatic FP consumer inherits this
  # shape, so the exact formulas matter.
  A_anode <- unname(.xrf_atomic_mass[anode])
  if (!is.finite(A_anode)) A_anode <- 2.3 * Z                       # ~A/Z for heavy elements
  m_bs <- 0.1382 - 0.9211 / sqrt(Z)
  eta <- (0.1904 - 0.2236 * log(Z) + 0.1292 * log(Z)^2 - 0.0149 * log(Z)^3) * E0^m_bs
  rhozmax <- (A_anode / Z) * (0.787e-5 * sqrt(0.0135 * Z) * E0^1.5 + 0.735e-6 * E0^2)  # g/cm^2
  sinfac <- sin(max(electron_incidence_deg, 1e-3) * pi / 180) /
    max(sin(takeoff_deg * pi / 180), 1e-3)
  # f_abs_at(en): Ebel escape factor (1-exp(-tau))/tau at photon energy en, with the energy-dependent
  # mean depth rhoz(en) = rhozmax * p1/p2 (-> 0 at the endpoint: those photons generate at the surface)
  # and the PHOTOELECTRIC anode attenuation (Ebel's tau; scatter only redirects within the anode).
  f_abs_at <- function(en) {
    lu <- log(pmax(E0 / en, 1))
    p1 <- lu * (0.49269 - 1.0987 * eta + 0.78557 * eta^2)
    p2 <- 0.70256 - 1.09865 * eta + 1.0046 * eta^2 + lu
    tau <- 2 * xrf_mass_attenuation(anode, en, type = "photoelectric") * rhozmax * (p1 / p2) * sinfac
    ifelse(is.finite(tau) & tau > 1e-9, (1 - exp(-tau)) / tau, 1)
  }
  # bare continuum at arbitrary energies, with Ebel's absolute constant 1.35e9 photons/(sr mA keV s) --
  # kept (rather than dropped as an arbitrary overall factor) so the continuum and the ABSOLUTE
  # characteristic-line fluxes below share one physical scale (their ratio is the point of the model).
  cont_bare_at <- function(en) {
    v <- 1.35e9 * Z * pmax(E0 / en - 1, 0)^x * f_abs_at(en)
    if (isTRUE(filter)) v <- v * .xrf_filter_transmission(tube$filter, en)
    v
  }
  cont <- cont_bare_at(E)
  cont[!is.finite(cont) | E >= E0 | E <= 0] <- 0

  # Characteristic anode lines with ABSOLUTE Ebel intensities (.xrf_ebel_char_lines): stopping-power,
  # backscatter and depth-absorption factors give each K/L line a flux in photons/sr/mA/s on the same
  # scale as the continuum -- replacing the old "char_peak_ratio x local continuum" heuristic, whose
  # arbitrary line/continuum balance biased every excitation integral toward or away from the anode
  # lines. `flux` (integrated) feeds discrete_lines = TRUE; the default rendering spreads each line as
  # a narrow Gaussian (sigma = s_line) of the same area, because a 0.03-keV line sampled on a ~0.1-keV
  # grid is badly under-integrated (grid-phase jitter) unless taken analytically.
  s_line <- 0.03
  char_lines <- NULL
  if (E0 > 0 && char_peak_ratio > 0) {
    flx <- .xrf_ebel_char_lines(anode, Z, E0, eta, rhozmax, sinfac)
    if (!is.null(flx) && nrow(flx) > 0) {
      flux <- char_peak_ratio * flx$flux
      if (isTRUE(filter)) flux <- flux * .xrf_filter_transmission(tube$filter, flx$energy_kev)
      char_lines <- data.frame(energy_kev = flx$energy_kev,
                               amp = flux / (s_line * sqrt(2 * pi)), flux = flux)
    }
  }
  if (isTRUE(discrete_lines)) {
    return(if (is.null(char_lines)) data.frame(energy_kev = numeric(0), flux = numeric(0))
           else char_lines[, c("energy_kev", "flux")])
  }
  out <- cont
  if (!is.null(char_lines)) {
    for (i in seq_len(nrow(char_lines))) {
      out <- out + char_lines$amp[i] * exp(-0.5 * ((E - char_lines$energy_kev[i]) / s_line)^2)
    }
  }
  out
}
