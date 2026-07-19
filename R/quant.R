
# Standard atomic weights (IUPAC), for oxide stoichiometry.
.xrf_atomic_mass <- c(
  H = 1.008, He = 4.0026, Li = 6.94, Be = 9.0122, B = 10.81, C = 12.011, N = 14.007, O = 15.999,
  F = 18.998, Ne = 20.180, Na = 22.990, Mg = 24.305, Al = 26.982, Si = 28.085, P = 30.974,
  S = 32.06, Cl = 35.45, Ar = 39.948, K = 39.098, Ca = 40.078, Sc = 44.956, Ti = 47.867,
  V = 50.942, Cr = 51.996, Mn = 54.938, Fe = 55.845, Co = 58.933, Ni = 58.693, Cu = 63.546,
  Zn = 65.38, Ga = 69.723, Ge = 72.630, As = 74.922, Se = 78.971, Br = 79.904, Kr = 83.798,
  Rb = 85.468, Sr = 87.62, Y = 88.906, Zr = 91.224, Nb = 92.906, Mo = 95.95, Tc = 98, Ru = 101.07,
  Rh = 102.91, Pd = 106.42, Ag = 107.87, Cd = 112.41, In = 114.82, Sn = 118.71, Sb = 121.76,
  Te = 127.60, I = 126.90, Xe = 131.29, Cs = 132.91, Ba = 137.33, La = 138.91, Ce = 140.12,
  Pr = 140.91, Nd = 144.24, Pm = 145, Sm = 150.36, Eu = 151.96, Gd = 157.25, Tb = 158.93,
  Dy = 162.50, Ho = 164.93, Er = 167.26, Tm = 168.93, Yb = 173.05, Lu = 174.97, Hf = 178.49,
  Ta = 180.95, W = 183.84, Re = 186.21, Os = 190.23, Ir = 192.22, Pt = 195.08, Au = 196.97,
  Hg = 200.59, Tl = 204.38, Pb = 207.2, Bi = 208.98, Po = 209, At = 210, Rn = 222, Fr = 223,
  Ra = 226, Ac = 227, Th = 232.04, Pa = 231.04, U = 238.03
)

# Default oxygen atoms per cation atom for common rock-forming / geological oxides (valence / 2).
# Elements absent here contribute as the metal (0 oxygen) unless overridden.
.xrf_oxide_o_per_cation <- c(
  Na = 0.5, Mg = 1, Al = 1.5, Si = 2, P = 2.5, K = 0.5, Ca = 1, Ti = 2, Cr = 1.5, Mn = 1,
  Fe = 1.5, Ni = 1, Cu = 1, Zn = 1, Ba = 1, Sr = 1, Zr = 2, V = 2.5, Pb = 1, Rb = 0.5,
  La = 1.5, Ce = 2, Nd = 1.5, Y = 1.5, S = 3, Co = 1, Sn = 2, W = 3, Nb = 2.5, Th = 2, U = 3
)

# Representative (approximate, tunable) generic-matrix compositions (mass fractions) for the
# "generalized" self-absorption depth correction in xrf_observed_mass(). These are stand-ins for a
# matrix whose true composition is unknown -- deliberately generic; supply your own named vector for
# anything specific. Fractions need not sum to 1 (they are renormalized).
.xrf_generic_matrix <- list(
  silicate  = c(O = 0.47, Si = 0.29, Al = 0.08, Fe = 0.05, Ca = 0.03, K = 0.026, Na = 0.023,
                Mg = 0.015, Ti = 0.004),
  basalt    = c(O = 0.44, Si = 0.23, Al = 0.08, Fe = 0.09, Ca = 0.07, Mg = 0.045, Na = 0.02,
                Ti = 0.01, K = 0.008),
  soil      = c(O = 0.48, Si = 0.30, Al = 0.07, Fe = 0.04, C = 0.03, Ca = 0.02, K = 0.02, Mg = 0.01),
  carbonate = c(Ca = 0.40, O = 0.48, C = 0.12),
  sulfide   = c(Fe = 0.4655, S = 0.5345),
  cement    = c(Ca = 0.45, O = 0.42, Si = 0.10, Al = 0.03, Fe = 0.02),
  organic   = c(C = 0.50, O = 0.40, N = 0.03, H = 0.07)
)

# Resolve a matrix argument (a built-in name or a named composition vector) to a normalized
# mass-fraction vector.
.xrf_resolve_matrix <- function(matrix) {
  comp <- if (is.character(matrix) && length(matrix) == 1) {
    m <- .xrf_generic_matrix[[matrix]]
    if (is.null(m)) stop("Unknown generic matrix '", matrix, "'. Built-ins: ",
                         paste(names(.xrf_generic_matrix), collapse = ", "),
                         "; or pass a named composition vector.")
    m
  } else if (is.numeric(matrix) && !is.null(names(matrix))) {
    matrix
  } else {
    stop("`matrix` must be a built-in name or a named vector of mass fractions.")
  }
  comp[!is.finite(comp) | comp < 0] <- 0
  comp / sum(comp)
}

# Total mass attenuation (cm^2/g) of a composition (named mass fractions) at given energies.
.xrf_matrix_total_mu <- function(comp, energy_kev) {
  mu <- rep(0, length(energy_kev))
  for (el in names(comp)) {
    if (!(el %in% all_elements)) next
    m <- .xrf_element_total_mu(el, energy_kev)
    m[!is.finite(m)] <- 0
    mu <- mu + comp[[el]] * m
  }
  mu
}

# Warn if a calibration vector was fitted with different settings than are in use. Every stamped setting
# is compared, because the per-element sensitivity S_i (hence the absolute calibration scale) depends not
# only on self_absorption/matrix/normalization but also on beam_energy_kev, detector_type, efficiency,
# excitation(_weighting), coster_kronig and the tube -- a mismatch in any of those silently corrupts the
# absolute values. `current` is a named list of the settings in use.
.xrf_check_calibration_settings <- function(calibration, current) {
  st <- attr(calibration, "settings")
  if (is.null(st)) return(invisible())
  bad <- character(0)
  for (k in intersect(names(st), names(current))) {
    # `matrix` only affects the result under generalized self-absorption
    if (k == "matrix" && !identical(current$self_absorption, "generalized")) next
    same <- isTRUE(tryCatch(all.equal(st[[k]], current[[k]]), error = function(e) FALSE)) ||
      identical(st[[k]], current[[k]])
    if (!same) bad <- c(bad, k)
  }
  if (length(bad) > 0) {
    warning("calibration was fitted with different ", paste(bad, collapse = ", "),
            " settings than used here; the absolute values will be wrong -- match them.", call. = FALSE)
  }
  invisible()
}

# Compact, comparable snapshot of the settings that determine the FP sensitivity scale (stamped on a
# calibration vector and rebuilt on the unknowns for .xrf_check_calibration_settings).
.xrf_calibration_settings <- function(self_absorption, matrix, normalization, beam_energy_kev,
                                      detector_type, efficiency, excitation, excitation_weighting,
                                      coster_kronig, tube) {
  list(self_absorption = self_absorption, matrix = matrix, normalization = normalization,
       beam_energy_kev = beam_energy_kev, detector_type = detector_type, efficiency = efficiency,
       excitation = excitation, excitation_weighting = excitation_weighting, coster_kronig = coster_kronig,
       tube = if (is.null(tube)) NULL else list(anode = tube$anode, kv = tube$kv))
}

# Resolve a sensitivity relative-uncertainty spec (a scalar applied to all elements, or a named per-element
# numeric vector) to a per-element vector aligned with `elements`. Unnamed / missing entries and non-finite or
# negative values become 0 (no systematic contribution). Used to fold the FP-sensitivity SYSTEMATIC error into
# concentration_se / observed_mass_se, in quadrature with the statistical peak-area error.
.xrf_resolve_rel_se <- function(x, elements) {
  if (is.null(x)) return(rep(0, length(elements)))
  v <- if (!is.null(names(x))) suppressWarnings(as.numeric(x[elements]))
       else rep(suppressWarnings(as.numeric(x))[1], length(elements))
  v[!is.finite(v) | v < 0] <- 0
  v
}

#' Compton-normalize deconvolution peak areas
#'
#' Divide each element's net peak area by the fitted Compton (incoherent) scatter area. The Compton
#' peak scales with the average matrix scattering power, so the ratio cancels much of the matrix
#' (absorption) effect for mid-Z analytes -- the standard, robust normalizer for handheld/portable
#' XRF. Requires a deconvolution run with scatter templates (supply a \code{tube} to
#' \link{xrf_deconvolute_gaussian_least_squares}); by default it normalizes to the Compton area, or
#' set \code{by = "rayleigh"} or \code{by = "scatter"} (Compton + Rayleigh).
#'
#' @param object A \code{deconvolution_fit} (from \link{xrf_deconvolute_gaussian_least_squares}) or
#'   its \code{peaks} data frame.
#' @param by Which scatter component to normalize by: "compton" (default), "rayleigh", or "scatter".
#'
#' @return The peaks data frame with added columns \code{peak_area_norm} and \code{peak_area_norm_se}.
#' @export
#'
xrf_compton_normalize <- function(object, by = c("compton", "rayleigh", "scatter")) {
  by <- match.arg(by)
  peaks <- if (inherits(object, "deconvolution_fit")) object$peaks else object
  stopifnot("element" %in% colnames(peaks), "peak_area" %in% colnames(peaks))

  # match BOTH the anode-line scatter (scatter_compton) AND the scattered-continuum pseudo-elements
  # (scatter_compton_continuum) via a prefix pattern -- the continuum carries most of the incoherent
  # scattering power, so omitting it (as an exact "scatter_compton" match did) understates the normalizer
  # whenever scatter_continuum = TRUE.
  pat <- switch(
    by,
    compton  = "^scatter_compton",
    rayleigh = "^scatter_rayleigh",
    scatter  = "^scatter_(compton|rayleigh)"
  )
  denom <- sum(peaks$peak_area[grepl(pat, peaks$element)], na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) {
    stop("No positive ", by, " scatter area found; run the deconvolution with a `tube` so scatter ",
         "templates are fitted.")
  }
  peaks$peak_area_norm <- peaks$peak_area / denom
  if ("peak_area_se" %in% colnames(peaks)) {
    scat <- grepl(pat, peaks$element)
    # Ratio-error propagation for X/Y = area_i / denom, denom = sum of the scatter areas:
    #   Var(X/Y) = Var(X)/Y^2 + X^2 Var(Y)/Y^4 - 2 X Cov(X,Y)/Y^3.
    # When the fit's peak-area COVARIANCE is available (a deconvolution_fit's area_cov, from the Laplace
    # posterior) the denominator variance uses the true correlated sum and the element/denominator covariance
    # Cov(area_i, denom) is included; otherwise denom_se falls back to a quadrature sum (independence) and the
    # Cov(X,Y) term is dropped -- the standard marginal-SE approximation.
    acov <- if (inherits(object, "deconvolution_fit")) object$area_cov else NULL
    if (!is.null(acov) && all(peaks$element %in% rownames(acov))) {
      els <- peaks$element; sc <- els[scat]
      denom_var <- sum(acov[sc, sc, drop = FALSE])
      cov_i_denom <- rowSums(acov[els, sc, drop = FALSE])
      var_norm <- (peaks$peak_area_se / denom)^2 +
        peaks$peak_area^2 * denom_var / denom^4 -
        2 * peaks$peak_area * cov_i_denom / denom^3
      peaks$peak_area_norm_se <- sqrt(pmax(var_norm, 0))
    } else {
      denom_se <- sqrt(sum(peaks$peak_area_se[scat]^2, na.rm = TRUE))
      peaks$peak_area_norm_se <- sqrt((peaks$peak_area_se / denom)^2 +
                                        (peaks$peak_area * denom_se / denom^2)^2)
    }
  }
  peaks
}

#' Fundamental-parameters sensitivity per element
#'
#' The theoretical detected intensity of each element's primary line per unit mass fraction (before
#' matrix absorption): the subshell production (photo- or electron-ionization cross section at the
#' beam energy, times fluorescence yield and radiative branching) times the detector efficiency at
#' the line energy. Values are relative (a single geometry/source constant is dropped), which is
#' what \link{xrf_quantify} needs.
#'
#' @param elements Element symbols.
#' @param beam_energy_kev Excitation energy (keV).
#' @param detector_type,be_window_um,dead_layer_um,active_thickness_um Detector parameters for the efficiency term
#'   (see \link{xrf_detector_efficiency}); \code{efficiency = FALSE} omits it.
#' @param efficiency Include the detector-efficiency term?
#' @param excitation,excitation_weighting,coster_kronig Passed to \link{xrf_relative_peak_intensity}.
#' @param tube Optional \link{xrf_tube}. When supplied (photon mode) the subshell excitation is the
#'   integral of the photoionization cross section over the \link{xrf_tube_spectrum} (polychromatic)
#'   rather than the monochromatic \code{beam_energy_kev}.
#'
#' @return A tibble with \code{element}, \code{primary_energy_kev}, and \code{sensitivity}.
#' @export
#'
xrf_fp_sensitivity <- function(elements, beam_energy_kev, detector_type = NULL,
                               be_window_um = NULL, dead_layer_um = NULL, active_thickness_um = NULL,
                               air_path_cm = NULL, atmosphere = "Air", window = NULL,
                               efficiency = TRUE,
                               excitation = c("photon", "electron"),
                               excitation_weighting = c("cross_section", "jump"),
                               coster_kronig = TRUE, tube = NULL) {
  excitation <- match.arg(excitation)
  excitation_weighting <- match.arg(excitation_weighting)
  xe <- xrftools::x_ray_xrf_energies

  # Polychromatic excitation: when a tube is supplied (photon mode), weight each subshell's
  # production by the integral of its photoionization cross section over the Ebel tube spectrum
  # (continuum + anode lines) instead of a single monochromatic energy.
  poly <- !is.null(tube) && inherits(tube, "xrf_tube") && excitation == "photon"
  if (poly) {
    grid <- seq(0.1, tube$kv * 0.999, length.out = 400L)
    de <- grid[2] - grid[1]
    n_cont <- xrf_tube_spectrum(tube, grid, char_peak_ratio = 0)      # smooth continuum (grid integration)
    char <- xrf_tube_spectrum(tube, grid, discrete_lines = TRUE)      # anode lines (integrated flux, discrete)
  }

  purrr::map_dfr(elements, function(el) {
    lines <- xe[xe$element == el & xe$edge_kev <= beam_energy_kev, , drop = FALSE]
    if (nrow(lines) == 0) {
      return(tibble::tibble(element = el, primary_energy_kev = NA_real_, sensitivity = NA_real_))
    }
    rel <- if (poly) {
      # excitation integral = continuum on the grid + the anode characteristic lines added ANALYTICALLY at
      # their exact energies (sigma(E_line) * line_flux), instead of sampling 0.03-keV lines on the ~0.1-keV
      # grid -- which under-integrated them with up to ~33% grid-phase jitter (worse after the B2/B3 harder
      # continuum). Removes the ~10-14% FP-sensitivity bias for elements excited mainly by an anode line.
      ex_shell <- function(sh) {
        s <- xrf_photoionization_cross_section(el, sh, grid); s[!is.finite(s)] <- 0
        val <- sum(s * n_cont) * de
        if (nrow(char)) {
          sl <- xrf_photoionization_cross_section(el, sh, char$energy_kev); sl[!is.finite(sl)] <- 0
          val <- val + sum(sl * char$flux)
        }
        val
      }
      excit <- vapply(lines$edge, ex_shell, numeric(1))
      # Coster-Kronig L-subshell cascade (mirrors the mono xrf_relative_peak_intensity path, which the
      # poly branch previously bypassed -- so the coster_kronig argument had no effect under a tube and
      # heavy-element L-line sensitivities were under-counted). Redistribute the tube-integrated L1/L2
      # vacancy production down to L3 before radiative decay; non-L rows are unchanged.
      if (coster_kronig) {
        nrw <- length(lines$edge)
        excit <- .xrf_ck_cascade(rep(el, nrw), lines$edge, excit,
                                 rep(ex_shell("L1"), nrw), rep(ex_shell("L2"), nrw),
                                 rep(ex_shell("L3"), nrw))
      }
      excit * xrf_transition_probability(el, lines$trans)   # omega * branching (from emission prob)
    } else {
      xrf_relative_peak_intensity(el, lines$edge, lines$trans, beam_energy_kev,
                                  excitation_weighting = excitation_weighting,
                                  coster_kronig = coster_kronig, excitation = excitation)
    }
    eff <- if (efficiency) {
      # full_energy_peak = TRUE: the sensitivity must predict the measured PHOTOPEAK area, so remove the
      # detector-escape fraction (1 - f_escape). Matters for Ge/CdTe near their K edges (~10-14%).
      xrf_detector_efficiency(lines$energy_kev, detector_type = detector_type,
                              active_thickness_um = active_thickness_um,
                              be_window_um = be_window_um, dead_layer_um = dead_layer_um,
                              air_path_cm = air_path_cm, atmosphere = atmosphere, window = window,
                              full_energy_peak = TRUE)
    } else {
      rep(1, nrow(lines))
    }
    i <- which.max(rel)   # production-primary line (matches the deconvolution's primary_energy_kev)
    tibble::tibble(element = el, primary_energy_kev = lines$energy_kev[i],
                   primary_edge_kev = lines$edge_kev[i], primary_shell = lines$edge[i],
                   sensitivity = (rel * eff)[i])
  })
}

#' First-order fundamental-parameters quantification
#'
#' Convert deconvolution peak areas to mass fractions using \link{xrf_fp_sensitivity} and, optionally,
#' an iterative self-absorption (matrix) correction for an infinitely-thick homogeneous sample:
#' \deqn{C_i \propto A_i \Big/ \big(S_i \cdot M_i\big), \qquad
#'       M_i = 1 / \big(\mu_s(E_0)/\sin\psi_{in} + \mu_s(E_i)/\sin\psi_{out}\big),}
#' with \eqn{\mu_s} the sample mass attenuation (mass-weighted over the current composition plus any
#' \code{dark_matrix}). Concentrations are renormalized to sum to \code{total}. With
#' \code{secondary_fluorescence = TRUE} it also adds first-order secondary (enhancement) fluorescence
#' (Shiraiwa-Fujino): element i is additionally excited by the lines of every element j whose line
#' energy exceeds i's absorption edge (and, with \code{tertiary_fluorescence = TRUE}, the third-order
#' k -> j -> i chain). Excitation is monochromatic at \code{beam_energy_kev} unless a
#' \code{tube} is supplied, in which case the primary excitation is integrated over the polychromatic
#' \link{xrf_tube_spectrum} (Ebel continuum + anode lines). Unmeasured light matrix (e.g. oxygen) can
#' be supplied via \code{dark_matrix}, or computed by stoichiometry with \code{oxide}. Scatter/escape
#' pseudo-elements are ignored.
#'
#' @param object A \code{deconvolution_fit} or its \code{peaks} data frame (needs \code{element},
#'   \code{peak_area}).
#' @param beam_energy_kev Excitation energy (keV).
#' @param detector_type,be_window_um,dead_layer_um,active_thickness_um,efficiency,excitation,excitation_weighting,coster_kronig,tube
#'   Passed to \link{xrf_fp_sensitivity}. Supplying \code{tube = }\link{xrf_tube} integrates the
#'   primary excitation over the polychromatic tube spectrum instead of a single energy.
#' @param self_absorption Apply the iterative self-absorption correction?
#' @param secondary_fluorescence Add secondary (enhancement) fluorescence?
#' @param tertiary_fluorescence Add third-order (k -> j -> i chain) enhancement fluorescence
#'   (implies \code{secondary_fluorescence}). Usually a small (<~1-2\%) correction.
#' @param incidence_deg,takeoff_deg Beam incidence / detector take-off angles (degrees) for the
#'   self-absorption path lengths. These are passed as plain scalars -- \code{xrf_quantify} does not read an
#'   \link{xrf_geometry} object (whose \code{scatter_angle_deg}, used for the Compton shift in the
#'   deconvolution, is a separate quantity); set them to match your instrument, consistently with the
#'   \code{scatter_angle_deg} you used for scatter (a common convention is
#'   \eqn{scatter\_angle \approx 180 - incidence - takeoff}).
#' @param dark_matrix Optional named vector of mass fractions of unmeasured matrix elements (e.g.
#'   \code{c(O = 0.47)}) included in the absorption but not quantified.
#' @param oxide Oxide stoichiometry for oxygen-by-difference. \code{FALSE} (default) makes no oxide
#'   assumption -- elements are reported on a metal basis and no oxygen is added (XRF cannot measure
#'   oxygen directly, so oxide estimates are a modelling choice, not a measurement). \code{TRUE} uses
#'   the built-in common geological oxides; a named vector (e.g. \code{c(Fe = 1, Si = 2)}) overrides
#'   the oxygen-atoms-per-cation ratios. When on, oxygen is computed from the cation concentrations,
#'   added to the absorbing matrix, and reported as an \code{"O"} row.
#' @param total Concentrations are scaled to sum to this (1 for mass fraction, 100 for percent).
#' @param iterations Self-absorption iterations.
#' @param areal_density Sample areal density \eqn{\rho d} (g/cm^2) for a finite-thickness (thin-film /
#'   coating) self-absorption correction. Default \code{Inf} is an infinitely-thick homogeneous sample
#'   (the emitted intensity saturates, kernel \eqn{1/A}); a finite value uses \eqn{(1-e^{-A\rho d})/A}
#'   with \eqn{A=\mu_s(E_0)/\sin\psi_{in}+\mu_s(E_i)/\sin\psi_{out}}, so a thin layer self-absorbs less
#'   (\eqn{\to \rho d} as \eqn{\rho d \to 0}). Applies to both the monochromatic and polychromatic paths.
#'
#' @param sensitivity_rel_se Optional relative 1-sigma uncertainty on the FP sensitivity \eqn{S_i} (the
#'   SYSTEMATIC error from the atomic constants -- fluorescence yields, cross sections, Coster-Kronig; e.g.
#'   PTB reference-free fundamental-parameter values). A scalar applies to all elements; a named numeric
#'   vector is per element (unnamed get 0). Folded in quadrature with the statistical peak-area error into
#'   \code{concentration_se} / \code{observed_mass_se}. Default 0 (statistical error only).
#' @return A tibble with \code{element}, \code{primary_energy_kev}, \code{peak_area},
#'   \code{sensitivity}, \code{concentration} (normalized, closed to \code{total}),
#'   \code{observed_mass} (un-normalized, matrix-corrected \eqn{A_i/(S_i\,\cdot\,\mathrm{selfabs}\,\cdot\,(1+\mathrm{enh}))}),
#'   and their standard errors \code{concentration_se}/\code{observed_mass_se} (first-order, propagated from
#'   the fit's \code{peak_area_se}: the relative area error times the value; \code{NA} when the fit could not
#'   estimate the area error, and treating the FP model factors + the closure sum as known).
#' @export
#'
xrf_quantify <- function(object, beam_energy_kev, detector_type = NULL,
                         be_window_um = NULL, dead_layer_um = NULL, active_thickness_um = NULL,
                         air_path_cm = NULL, atmosphere = "Air", window = NULL,
                         efficiency = TRUE,
                         excitation = c("photon", "electron"),
                         excitation_weighting = c("cross_section", "jump"), coster_kronig = TRUE,
                         tube = NULL, self_absorption = TRUE, secondary_fluorescence = FALSE,
                         tertiary_fluorescence = FALSE, incidence_deg = 45, takeoff_deg = 45,
                         dark_matrix = NULL, oxide = FALSE, total = 1, iterations = 12,
                         areal_density = Inf, sensitivity_rel_se = 0) {
  if (isTRUE(tertiary_fluorescence)) secondary_fluorescence <- TRUE   # tertiary implies secondary
  peaks <- if (inherits(object, "deconvolution_fit")) object$peaks else object
  stopifnot("element" %in% colnames(peaks), "peak_area" %in% colnames(peaks))

  # drop scatter / escape pseudo-elements; keep only real elements with positive area
  keep <- !grepl("^scatter_", peaks$element) & peaks$element %in% all_elements &
    is.finite(peaks$peak_area) & peaks$peak_area > 0
  q <- peaks[keep, , drop = FALSE]
  if (nrow(q) == 0) stop("No real element peaks with positive area to quantify.")

  s <- xrf_fp_sensitivity(q$element, beam_energy_kev, detector_type = detector_type,
                          be_window_um = be_window_um, dead_layer_um = dead_layer_um,
                          active_thickness_um = active_thickness_um,
                          air_path_cm = air_path_cm, atmosphere = atmosphere, window = window,
                          efficiency = efficiency, excitation = match.arg(excitation),
                          excitation_weighting = match.arg(excitation_weighting),
                          coster_kronig = coster_kronig, tube = tube)
  q$primary_energy_kev <- s$primary_energy_kev
  q$primary_edge_kev <- s$primary_edge_kev
  q$primary_shell <- s$primary_shell
  q$sensitivity <- s$sensitivity
  ok <- is.finite(q$sensitivity) & q$sensitivity > 0
  q <- q[ok, , drop = FALSE]
  if (nrow(q) == 0) stop("No elements have a usable fundamental-parameters sensitivity.")

  # oxide stoichiometry (L3): oxygen mass per unit cation mass, r_i = (O per cation) * A_O / A_i.
  # With this on, oxygen is computed from the cation concentrations, included in the absorbing matrix,
  # and reported -- so `dark_matrix` need not be supplied by hand for geological samples.
  els <- q$element
  A_O <- .xrf_atomic_mass[["O"]]
  oxide_on <- !is.null(oxide) && !isFALSE(oxide)
  r <- rep(0, length(els))
  if (oxide_on) {
    opc <- if (is.numeric(oxide)) oxide else .xrf_oxide_o_per_cation   # named-vector override allowed
    o_i <- unname(opc[els]); o_i[is.na(o_i)] <- 0
    A_i <- unname(.xrf_atomic_mass[els])
    r <- o_i * A_O / A_i
    r[!is.finite(r)] <- 0
  }
  # normalize the cations (+ derived oxygen) to sum to `total`
  renorm <- function(cc) cc * (total / (sum(cc) + sum(cc * r)))

  raw <- q$peak_area / q$sensitivity          # un-normalized observed mass; the FP loop refines it in place
  conc <- renorm(raw)

  if (isTRUE(self_absorption) || isTRUE(secondary_fluorescence)) {
    if (match.arg(excitation) == "electron") {
      warning("xrf_quantify: the self-absorption and secondary-fluorescence corrections use the PHOTON ",
              "matrix model (X-ray path-length self-absorption + Shiraiwa-Fujino enhancement). For electron ",
              "excitation this is NOT a ZAF / phi-rho-z correction (no atomic-number/stopping-power/",
              "backscatter terms), so electron-mode concentrations are approximate.", call. = FALSE)
    }
    sin_in <- sin(incidence_deg * pi / 180)
    sin_out <- sin(takeoff_deg * pi / 180)
    Ei <- q$primary_energy_kev
    dm_el <- names(dark_matrix); dm_frac <- as.numeric(dark_matrix)
    # Sample self-absorption uses TOTAL mass attenuation (every interaction -- photoelectric AND scatter --
    # removes a photon from the beam/emission path). This is distinct from the photoelectric mu used below
    # for fluorescence PRODUCTION (tau_E0 / tau_ij), where only photoelectric absorption creates a vacancy.
    # Precompute per-element TOTAL mu at every energy the matrix is evaluated at (the beam energy + each
    # element's line energy). These are composition-INDEPENDENT, so doing the (expensive) mass-attenuation
    # lookups once here -- rather than on every iteration -- makes mu_at() a fast weighted sum.
    mu_energies <- c(beam_energy_kev, Ei)                    # index 1 = beam, index i+1 = line i
    mu_tab <- vapply(mu_energies, function(en) vapply(els, function(el) .xrf_element_total_mu(el, en), numeric(1)),
                     numeric(length(els)))
    if (!is.matrix(mu_tab)) mu_tab <- matrix(mu_tab, nrow = length(els))            # length(els)==1 guard
    muO_tab  <- if (any(r > 0)) vapply(mu_energies, function(en) .xrf_element_total_mu("O", en), numeric(1)) else NULL
    muDM_tab <- if (length(dm_el)) vapply(mu_energies, function(en) vapply(dm_el, function(el) .xrf_element_total_mu(el, en), numeric(1)),
                                          numeric(length(dm_el))) else NULL
    if (!is.null(muDM_tab) && !is.matrix(muDM_tab)) muDM_tab <- matrix(muDM_tab, nrow = length(dm_el))
    mu_at <- function(k) {   # k indexes mu_energies; matrix attenuation over the current mix
      m <- sum(conc * mu_tab[, k])
      o_c <- sum(conc * r); if (o_c > 0 && !is.null(muO_tab)) m <- m + o_c * muO_tab[k]      # derived oxygen
      if (!is.null(muDM_tab)) m <- m + sum(dm_frac * muDM_tab[, k])
      m
    }

    # ---- W1/W10: polychromatic excitation grid (tube spectrum + per-element production) -----------
    # For a tube (polychromatic) source both the primary self-absorption (W1) and the secondary-
    # fluorescence exciter production (W10) belong under an excitation integral over the tube spectrum
    # N(E), not evaluated once at the endpoint kV. Precompute (composition-independent) the tube
    # spectrum, each element's primary-shell production on the grid, and the bare production integral
    # `excit_bare` = integral N(E) sigma_i(E) dE (reduces to sigma_i(E0) for a monochromatic source).
    poly <- !is.null(tube) && inherits(tube, "xrf_tube") && match.arg(excitation) == "photon"
    poly_selfabs <- isTRUE(self_absorption) && poly
    need_grid <- poly && (isTRUE(self_absorption) || isTRUE(secondary_fluorescence))
    if (need_grid) {
      grid <- seq(0.1, tube$kv * 0.999, length.out = 400L)
      dgrid <- grid[2] - grid[1]
      N_grid <- xrf_tube_spectrum(tube, grid)
      # per-element primary-shell production on the grid (composition-independent; same CK basis as S_i)
      P_grid <- vapply(seq_along(els), function(i)
        .xrf_shell_production(els[i], q$primary_shell[i], grid, coster_kronig = isTRUE(coster_kronig)),
        numeric(length(grid)))                                   # [n_grid x n_els]
      if (!is.matrix(P_grid)) P_grid <- matrix(P_grid, ncol = length(els))
      wgt_grid <- P_grid * N_grid                                # excitation weight N(E) * sigma_i(E)
      excit_bare <- colSums(wgt_grid) * dgrid                    # bare (no-absorption) integral per element
      excit_bare[!is.finite(excit_bare) | excit_bare <= 0] <- NA_real_
    }
    if (poly_selfabs) {
      # per-element TOTAL mu on the grid (composition-independent): mu_s(grid) = mu_grid_el %*% conc (+ O + dark)
      mu_grid_el <- vapply(els, function(el) .xrf_element_total_mu(el, grid), numeric(length(grid)))  # [n_grid x n_els]
      if (!is.matrix(mu_grid_el)) mu_grid_el <- matrix(mu_grid_el, ncol = length(els))
      muO_grid  <- if (any(r > 0)) .xrf_element_total_mu("O", grid) else NULL
      muDM_grid <- if (length(dm_el)) vapply(dm_el, function(el) .xrf_element_total_mu(el, grid), numeric(length(grid))) else NULL
      if (!is.null(muDM_grid) && !is.matrix(muDM_grid)) muDM_grid <- matrix(muDM_grid, ncol = length(dm_el))
    }

    # precompute the (concentration-independent) atomic factors for secondary fluorescence
    if (isTRUE(secondary_fluorescence)) {
      omega_j <- xrf_fluorescence_yield(els, q$primary_shell)
      edge_i <- q$primary_edge_kev
      # Exciter/absorber production uses the SHELL-PARTIAL photoionization cross section (with the CK
      # L-cascade) -- the SAME production basis as the primary sensitivity S_i -- not the total
      # photoelectric mu. Total mu counts vacancies in shells that do not emit the observed line, which
      # over-counts the enhancement (~13% for K exciters, several-fold for L-line/heavy exciters).
      # secondary fluorescence is a photon (photoionization) process regardless of the primary
      # excitation mode, so the CK L-cascade follows coster_kronig directly.
      ck <- isTRUE(coster_kronig)
      tau_E0 <- .xrf_shell_production(els, q$primary_shell, beam_energy_kev, coster_kronig = ck)
      # tau_ij[i, j] = shell-partial production of element i (its primary shell) at element j's line energy
      tau_ij <- vapply(seq_along(els), function(j)
        .xrf_shell_production(els, q$primary_shell, Ei[j], coster_kronig = ck), numeric(length(els)))
      if (!is.matrix(tau_ij)) tau_ij <- matrix(tau_ij, nrow = length(els))    # length(els)==1 guard
    }

    # E4: finite-thickness self-absorption kernel. For a homogeneous layer of areal density
    # m = rho*d (g/cm^2), the emitted intensity saturates as (1 - exp(-A*m))/A instead of the thick-sample
    # 1/A, where A = mu_s(E)/sin_in + mu_s(E_i)/sin_out. areal_density = Inf (default) recovers 1/A exactly
    # (thick sample); a finite value models a thin film / coating (-> m as m -> 0, no self-absorption).
    fin_abs <- function(A) {
      if (is.finite(areal_density)) (1 - exp(-A * areal_density)) / A else 1 / A
    }
    for (iter in seq_len(iterations)) {
      mu_e0 <- mu_at(1L)
      mu_ei <- vapply(seq_along(Ei), function(i) mu_at(i + 1L), numeric(1))
      if (!isTRUE(self_absorption)) {
        m_corr <- 1
      } else if (poly_selfabs) {
        # incident-leg absorption folded across the tube spectrum (W1); emergent leg fixed at E_i
        mu_s_grid <- as.vector(mu_grid_el %*% conc)
        o_c <- sum(conc * r); if (o_c > 0 && !is.null(muO_grid)) mu_s_grid <- mu_s_grid + o_c * muO_grid
        if (!is.null(muDM_grid)) mu_s_grid <- mu_s_grid + as.vector(muDM_grid %*% dm_frac)
        m_corr <- vapply(seq_along(els), function(i) {
          if (!is.finite(excit_bare[i])) return(fin_abs(mu_e0 / sin_in + mu_ei[i] / sin_out))  # fallback
          sum(wgt_grid[, i] * fin_abs(mu_s_grid / sin_in + mu_ei[i] / sin_out)) * dgrid / excit_bare[i]
        }, numeric(1))
      } else {
        m_corr <- fin_abs(mu_e0 / sin_in + mu_ei / sin_out)                                 # monochromatic (exact)
      }

      # Secondary/tertiary (enhancement) fluorescence, Shiraiwa-Fujino. B[i, j] is the enhancement of
      # element i per unit mass of element j (element j's line must clear i's edge). Secondary
      # enhancement is B %*% C; tertiary uses j's own enhanced intensity, B %*% (C*(1 + B%*%C)),
      # capturing the k -> j -> i chain.
      enh <- rep(0, length(els))
      if (isTRUE(secondary_fluorescence)) {
        B <- matrix(0, length(els), length(els))
        # W12: polychromatic g-factor INCIDENT leg. The incident term (primary beam penetrating to create
        # exciter j) is folded across the tube spectrum, weighted by exciter j's OWN excitation
        # N(E)*sigma_j(E) = wgt_grid[, j], instead of the single endpoint mu_e0 -- consistent with the W1
        # self-absorption integral above. The endpoint energy has the SMALLEST mu -> deepest penetration ->
        # LARGEST g, so the monochromatic mu_e0 systematically OVER-estimates the enhancement (measured
        # ~132% Fe->Cr for dilute Cr in near-pure Fe vs a ~20-60% expectation / <=22% empirical bound); the
        # excitation-weighted effective mu is larger (mid-spectrum-dominated), giving the physical, shallower
        # enhancement. Precomputed per exciter j; falls back to the endpoint form for a monochromatic source
        # or a missing grid.
        poly_g <- isTRUE(poly) && !is.null(wgt_grid) && is.numeric(dgrid)
        if (poly_g) {
          mu_s_grid2 <- as.vector(mu_grid_el %*% conc)
          o_c2 <- sum(conc * r); if (o_c2 > 0 && !is.null(muO_grid)) mu_s_grid2 <- mu_s_grid2 + o_c2 * muO_grid
          if (!is.null(muDM_grid)) mu_s_grid2 <- mu_s_grid2 + as.vector(muDM_grid %*% dm_frac)
        }
        inc_leg <- vapply(seq_along(els), function(j) {
          muj <- mu_ei[j]
          if (poly_g && is.finite(excit_bare[j]) && excit_bare[j] > 0)
            sum(wgt_grid[, j] * (sin_in / mu_s_grid2) * log(1 + mu_s_grid2 / (sin_in * muj))) * dgrid / excit_bare[j]
          else (sin_in / mu_e0) * log(1 + mu_e0 / (sin_in * muj))
        }, numeric(1))
        for (i in seq_along(els)) {
          for (j in seq_along(els)) {
            if (j == i || !is.finite(Ei[j]) || Ei[j] <= edge_i[i]) next
            muj <- mu_ei[j]                       # sample attenuation at j's line energy
            g <- 0.5 * muj * (inc_leg[j] +
                                (sin_out / mu_ei[i]) * log(1 + mu_ei[i] / (sin_out * muj)))
            # W10: under a tube, weight the enhancer/primary production ratio by the polychromatic
            # excitation integral (integral N(E) sigma(E) dE) instead of the single endpoint energy, so a
            # low-Z exciter excited across the whole spectrum is not undercounted relative to a heavy
            # exciter only reachable near the endpoint. Reduces to the monochromatic sigma_j(E0)/sigma_i(E0)
            # when N -> delta. (The g-factor incident term still uses the endpoint mu_e0 -- a bounded,
            # second-order residual, since (sin_in/mu_e0) ln(1 + mu_e0/(sin_in muj)) -> 1/muj as mu_e0 -> 0.)
            prod_ratio <- if (poly && is.finite(excit_bare[j]) && is.finite(excit_bare[i]) &&
                              excit_bare[i] > 0) {
              excit_bare[j] / excit_bare[i]
            } else {
              tau_E0[j] / tau_E0[i]
            }
            B[i, j] <- omega_j[j] * prod_ratio * (tau_ij[i, j] / muj) * g
          }
        }
        enh2 <- as.vector(B %*% conc)
        enh <- if (isTRUE(tertiary_fluorescence)) as.vector(B %*% (conc * (1 + enh2))) else enh2
      }

      raw <- q$peak_area / (q$sensitivity * m_corr * (1 + enh))   # matrix-corrected, still un-normalized
      conc <- renorm(raw)
    }
  }

  q$concentration <- conc
  # Un-normalized, matrix-corrected observed mass: A_i / (S_i * self_absorption * (1 + enhancement)). Unlike
  # `concentration` it is NOT closed to a fixed sum, so a single element's value is not forced up or down by
  # the others' total -- it moves only through the real matrix corrections (self-absorption + secondary/
  # tertiary fluorescence). This is what the CloudCal "full-fidelity" $Mass reports.
  q$observed_mass <- raw
  # Fit-uncertainty propagation (first order, area-dominated): the FP model factors (S_i, self-absorption
  # m_corr, enhancement) are treated as known, so the RELATIVE error of the fitted peak area carries through
  # to both observed_mass and concentration. concentration_se additionally ignores the covariance introduced
  # by the sum-to-`total` closure (a second-order term when the sum is dominated by well-determined majors).
  rel_se <- if ("peak_area_se" %in% names(q)) q$peak_area_se / q$peak_area else rep(NA_real_, nrow(q))
  # Fold the SYSTEMATIC FP-sensitivity relative uncertainty (sensitivity_rel_se: the 1-sigma on S_i from the
  # atomic constants -- yields / cross sections / CK, e.g. PTB reference-free FP values; 0 = off) in quadrature
  # with the statistical area error, so concentration_se / observed_mass_se report the total uncertainty.
  rel_se <- sqrt(rel_se^2 + .xrf_resolve_rel_se(sensitivity_rel_se, q$element)^2)
  q$observed_mass_se <- raw * rel_se
  q$concentration_se <- conc * rel_se
  out <- tibble::as_tibble(q[, c("element", "primary_energy_kev", "peak_area", "sensitivity",
                                 "concentration", "concentration_se", "observed_mass", "observed_mass_se")])
  if (oxide_on && any(r > 0)) {   # report the stoichiometric oxygen
    out <- dplyr::bind_rows(out, tibble::tibble(
      element = "O", primary_energy_kev = NA_real_, peak_area = NA_real_,
      sensitivity = NA_real_, concentration = sum(conc * r), concentration_se = NA_real_,
      observed_mass = NA_real_, observed_mass_se = NA_real_
    ))
  }
  out
}

#' Per-element observed mass (un-normalized, no forced closure)
#'
#' An alternative to \link{xrf_quantify} that does \strong{not} normalize to a closed sum and can be
#' computed \strong{without knowing the sample matrix} (in particular, no oxygen / oxide-state
#' assumption). Its key property is that the elements are \strong{decoupled}: adding or removing one
#' element does not change any other's result -- unlike a normalized fit, where a wrong oxide guess or
#' a missing light element smears across the whole assay. It reports, per element,
#' \deqn{observed\_mass = A_i / S_i,}
#' where \eqn{S_i} is the fundamental-parameters fluorescence sensitivity (\link{xrf_fp_sensitivity}:
#' excitation x fluorescence yield x branching x detector efficiency).
#'
#' \strong{Interpretation and its limits.} From the thick-sample relation
#' \eqn{A_i = k\,C_i\,S_i / (\mu_s(E_0)/\sin\psi_1 + \mu_s(E_i)/\sin\psi_2)}, we get
#' \eqn{A_i/S_i = k\,C_i\,\tau_{eff,i}}: the areal mass of element i \emph{within its own information
#' depth} \eqn{\tau_{eff,i}}. The value does NOT need \eqn{\mu_s} to be \emph{computed}, but it does
#' \emph{depend} on the matrix through \eqn{\tau_{eff,i}}, which is also different for every element
#' (via \eqn{E_i}). So raw \code{observed_mass} is an areal mass, \strong{not} a concentration: it is
#' not comparable across elements (line-energy/depth differences of ~40x within one silicate) nor
#' across samples of differing matrix (~9x for the same element between an organic and a sulfide
#' matrix). To make it comparable, use \code{self_absorption = "generalized"} (below) and/or
#' \code{normalization} -- and even then only as well as the generic matrix approximates the real one.
#'
#' Optionally divide by the fitted incoherent/coherent \strong{scatter} intensity
#' (\code{normalization = "compton"}), the standard matrix-robust normalizer for portable/geological
#' XRF (the scatter peak measures the light matrix's scattering and penetration), and/or supply a
#' per-element \code{calibration} (from \link{xrf_calibrate}) for absolute units.
#'
#' \strong{Caveats.} (1) The closed form and the generalized correction are exact only for
#' \emph{monochromatic} excitation; with a \code{tube} they are approximations (the absorption sits
#' inside the polychromatic excitation integral, leaving an energy-dependent residual). (2) This is a
#' \emph{primary}-fluorescence quantity: secondary (enhancement) fluorescence is NOT modelled, so
#' elements just below a strong exciter (e.g. Cr/V/Ti/Ca under strong Fe K\eqn{\alpha} in Fe-rich
#' rocks) read high and track the exciter -- use \link{xrf_quantify} when enhancement matters. (3) It
#' assumes a thick, homogeneous sample (thin films / coatings / heterogeneous samples violate this).
#' (4) For a counts-per-channel spectrum \code{peak_area} (and hence \code{observed_mass}) scales with the
#' energy-axis channel width (\code{observed_mass} \eqn{\propto} rate \eqn{\times} \eqn{\Delta E}), so it is
#' \strong{not comparable across MCAs of different gain} -- calibrate (\link{xrf_calibrate}) at the same
#' channel width/gain you will apply, exactly as for the other settings. (This cancels in \link{xrf_quantify}'s
#' normalized \code{concentration}, which closes to a fixed sum, so only the un-normalized value is affected.)
#'
#' @param object A \code{deconvolution_fit} or its \code{peaks} data frame (needs \code{element},
#'   \code{peak_area}; scatter pseudo-elements are used for \code{normalization} and then dropped).
#' @param beam_energy_kev,detector_type,be_window_um,dead_layer_um,active_thickness_um,efficiency,excitation,excitation_weighting,coster_kronig,tube
#'   Passed to \link{xrf_fp_sensitivity} to compute \eqn{S_i}.
#' @param self_absorption "none" (default) returns \eqn{A_i/S_i} (the observed areal mass, where the
#'   real matrix has cancelled into the per-element sampling depth). "generalized" divides out that
#'   per-element depth using a \emph{generic} (documented, tunable) matrix -- multiplying by
#'   \eqn{\mu_g(E_0)/\sin\psi_1 + \mu_g(E_i)/\sin\psi_2} -- to put the elements on a common,
#'   approximately-bulk concentration basis, without pretending to know the true composition or
#'   forcing a closed sum.
#' @param matrix The generic matrix used when \code{self_absorption = "generalized"}: a built-in name
#'   ("silicate", "basalt", "soil", "carbonate", "sulfide", "cement", "organic") or a named vector of
#'   mass fractions.
#' @param incidence_deg,takeoff_deg Beam incidence / detector take-off angles (degrees), used by the
#'   generalized self-absorption.
#' @param normalization "none" (default) applies no scatter normalization; "compton"/"rayleigh"/
#'   "scatter" divides by the corresponding fitted scatter area (matrix/depth normalization without
#'   knowing the matrix). Composes with \code{self_absorption}.
#' @param calibration Optional named vector of per-element response factors (e.g. from
#'   \link{xrf_calibrate}) to scale the relative values to absolute mass. Must have been fitted with
#'   the same \code{self_absorption}/\code{matrix}/\code{normalization} settings.
#'
#' @param sensitivity_rel_se Optional relative 1-sigma uncertainty on the FP sensitivity \eqn{S_i} (SYSTEMATIC,
#'   from the atomic constants; e.g. PTB reference-free fundamental-parameter values). A scalar applies to all
#'   elements; a named numeric vector is per element (unnamed get 0). Folded in quadrature with the statistical
#'   peak-area error into \code{observed_mass_se}. Default 0 (statistical error only).
#' @return A tibble with \code{element}, \code{primary_energy_kev}, \code{peak_area},
#'   \code{sensitivity}, \code{observed_mass} (un-normalized; the sum is diagnostic, not forced), and
#'   \code{observed_mass_se} (first-order standard error propagated from the fit's \code{peak_area_se} as the
#'   relative area error times \code{observed_mass}; \code{NA} when unavailable, and treating the sensitivity /
#'   matrix / normalization / calibration factors as known).
#' @export
#'
xrf_observed_mass <- function(object, beam_energy_kev, detector_type = NULL, be_window_um = NULL,
                              dead_layer_um = NULL, active_thickness_um = NULL,
                              air_path_cm = NULL, atmosphere = "Air", window = NULL, efficiency = TRUE,
                              excitation = c("photon", "electron"),
                              excitation_weighting = c("cross_section", "jump"),
                              coster_kronig = TRUE, tube = NULL,
                              self_absorption = c("none", "generalized"), matrix = "silicate",
                              incidence_deg = 45, takeoff_deg = 45,
                              normalization = c("none", "compton", "rayleigh", "scatter"),
                              calibration = NULL, sensitivity_rel_se = 0) {
  normalization <- match.arg(normalization)
  self_absorption <- match.arg(self_absorption)
  peaks <- if (inherits(object, "deconvolution_fit")) object$peaks else object
  stopifnot("element" %in% colnames(peaks), "peak_area" %in% colnames(peaks))

  keep <- !grepl("^(scatter|sum|pileup|escape)_", peaks$element) & peaks$element %in% all_elements &
    is.finite(peaks$peak_area) & peaks$peak_area > 0
  q <- peaks[keep, , drop = FALSE]
  if (nrow(q) == 0) stop("No real element peaks with positive area.")

  s <- xrf_fp_sensitivity(q$element, beam_energy_kev, detector_type = detector_type,
                          be_window_um = be_window_um, dead_layer_um = dead_layer_um,
                          active_thickness_um = active_thickness_um,
                          air_path_cm = air_path_cm, atmosphere = atmosphere, window = window,
                          efficiency = efficiency, excitation = match.arg(excitation),
                          excitation_weighting = match.arg(excitation_weighting),
                          coster_kronig = coster_kronig, tube = tube)
  q$primary_energy_kev <- s$primary_energy_kev
  q$sensitivity <- s$sensitivity

  obs <- q$peak_area / q$sensitivity                 # matrix-free observed areal mass (relative)
  # Fit-uncertainty propagation (first order, area-dominated): the FP sensitivity and the matrix /
  # normalization / calibration factors are treated as known, so the RELATIVE error of the fitted peak
  # area carries straight through to observed_mass. peak_area_se comes from the deconvolution (NA where the
  # fit could not estimate it, e.g. NNLS-clamped or aliased columns). The scatter-normalizer denominator's
  # uncertainty IS folded in below (when normalization is used); the calibration-factor uncertainty is not
  # (a plain calibration vector carries no SE).
  rel_se <- if ("peak_area_se" %in% names(q)) q$peak_area_se / q$peak_area else rep(NA_real_, nrow(q))

  if (self_absorption == "generalized") {            # remove per-element depth via a generic matrix
    comp <- .xrf_resolve_matrix(matrix)
    sin_in <- sin(incidence_deg * pi / 180); sin_out <- sin(takeoff_deg * pi / 180)
    mu_e0 <- .xrf_matrix_total_mu(comp, beam_energy_kev)
    mu_ei <- vapply(q$primary_energy_kev, function(e) .xrf_matrix_total_mu(comp, e), numeric(1))
    obs <- obs * (mu_e0 / sin_in + mu_ei / sin_out)
  }

  if (normalization != "none") {
    # include the scattered-continuum pseudo-elements (prefix match), not just the anode-line scatter
    pat <- switch(normalization,
                  compton = "^scatter_compton", rayleigh = "^scatter_rayleigh",
                  scatter = "^scatter_(compton|rayleigh)")
    denom <- sum(peaks$peak_area[grepl(pat, peaks$element)], na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) {
      stop("normalization = '", normalization, "' needs a fitted scatter peak; run the ",
           "deconvolution with a `tube`.")
    }
    obs <- obs / denom
    if ("peak_area_se" %in% names(peaks)) {   # ratio error: add the scatter-normalizer's relative SE
      denom_se <- sqrt(sum(peaks$peak_area_se[grepl(pat, peaks$element)]^2, na.rm = TRUE))
      rel_se <- sqrt(rel_se^2 + (denom_se / denom)^2)
    }
  }
  if (!is.null(calibration)) {                       # optional absolute calibration
    if (is.null(names(calibration))) {
      stop("`calibration` must be a named numeric vector keyed by element symbol.")
    }
    .xrf_check_calibration_settings(calibration, .xrf_calibration_settings(
      self_absorption, matrix, normalization, beam_energy_kev, detector_type, efficiency,
      match.arg(excitation), match.arg(excitation_weighting), coster_kronig, tube))
    k <- unname(calibration[q$element])
    miss <- q$element[!is.finite(k)]
    if (length(miss) > 0) {                          # NA (not 1): never mix relative + absolute scales
      warning("No calibration factor for: ", paste(miss, collapse = ", "),
              "; their observed_mass is returned as NA rather than in uncalibrated (relative) units.",
              call. = FALSE)
    }
    k[!is.finite(k)] <- NA_real_
    obs <- obs * k
  }

  # fold the systematic FP-sensitivity relative uncertainty in quadrature (see xrf_quantify)
  rel_se <- sqrt(rel_se^2 + .xrf_resolve_rel_se(sensitivity_rel_se, q$element)^2)
  tibble::as_tibble(data.frame(
    element = q$element, primary_energy_kev = q$primary_energy_kev, peak_area = q$peak_area,
    sensitivity = q$sensitivity, observed_mass = obs, observed_mass_se = obs * rel_se,
    stringsAsFactors = FALSE
  ))
}

#' Fit per-element calibration factors from reference materials
#'
#' Turns the relative \link{xrf_observed_mass} values into absolute mass units by regressing known
#' reference values against the observed values across one or more reference materials. For each
#' element, the factor \eqn{k_i} is the (through-origin) slope of \code{known} vs \code{observed}, so
#' that \code{observed_mass * k_i} recovers the reference quantity. The reference values can be in any
#' consistent unit (mass fraction, wt\%, ppm, areal density) -- \eqn{k_i} carries the units.
#'
#' The slope is fit through the origin (calibration is applied multiplicatively as
#' \code{observed_mass * k_i}, so an intercept could not be reproduced). Crucially, calibrate with the
#' \emph{same} \code{self_absorption} / \code{matrix} / \code{normalization} settings you will use on
#' unknowns -- these are stamped on the result and \link{xrf_observed_mass} warns on a mismatch. Pass
#' the returned vector back as \code{calibration =}; elements you did not calibrate come back as
#' \code{NA} (never silently mixed into the calibrated column). This is a per-element empirical
#' calibration and does not close the sum to 100\%.
#'
#' @param standards A list of reference measurements: each a \code{deconvolution_fit} or a peaks data
#'   frame (as accepted by \link{xrf_observed_mass}).
#' @param values The known reference values, either a list of named numeric vectors (element ->
#'   value), one per standard, or a data frame with an \code{element} column plus one column per
#'   standard, or one row per standard with element columns.
#' @param beam_energy_kev,detector_type,be_window_um,dead_layer_um,active_thickness_um,efficiency,excitation,excitation_weighting,coster_kronig,tube,self_absorption,matrix,incidence_deg,takeoff_deg,normalization
#'   Passed to \link{xrf_observed_mass} to compute the observed values for each standard (use the
#'   settings you will apply to unknowns).
#' @param weighting Least-squares weighting of the per-element through-origin slope across standards.
#'   \code{"none"} (default) is the ordinary slope \eqn{\sum o t / \sum o^2} (homoscedastic; dominated by
#'   the highest-level standard); \code{"relative"} uses inverse-variance weights \eqn{1/o^2} (proportional
#'   error), i.e. \eqn{k = \mathrm{mean}(t/o)}, so standards spanning orders of magnitude contribute evenly.
#'
#' @return A named numeric vector of per-element calibration factors \eqn{k_i}, with an
#'   \code{"n_standards"} attribute (points used per element) and a \code{"settings"} attribute
#'   (the observed-mass settings used). Pass it as \code{calibration =} to \link{xrf_observed_mass}.
#' @export
#'
xrf_calibrate <- function(standards, values, beam_energy_kev, detector_type = NULL,
                          be_window_um = NULL, dead_layer_um = NULL, active_thickness_um = NULL,
                          air_path_cm = NULL, atmosphere = "Air", window = NULL,
                          efficiency = TRUE,
                          excitation = c("photon", "electron"),
                          excitation_weighting = c("cross_section", "jump"), coster_kronig = TRUE,
                          tube = NULL, self_absorption = c("none", "generalized"),
                          matrix = "silicate", incidence_deg = 45, takeoff_deg = 45,
                          normalization = c("none", "compton", "rayleigh", "scatter"),
                          weighting = c("none", "relative")) {
  excitation <- match.arg(excitation)
  excitation_weighting <- match.arg(excitation_weighting)
  self_absorption <- match.arg(self_absorption)
  normalization <- match.arg(normalization)
  weighting <- match.arg(weighting)
  if (inherits(standards, "deconvolution_fit") || is.data.frame(standards)) standards <- list(standards)
  # coerce `values` to a list of named numeric vectors, one per standard
  if (is.numeric(values) && !is.null(names(values))) values <- list(values)   # lone standard
  if (is.data.frame(values)) {
    if ("element" %in% colnames(values)) {                 # element column + one column per standard
      vcols <- setdiff(colnames(values), "element")
      values <- lapply(vcols, function(cc) stats::setNames(values[[cc]], values$element))
    } else {                                               # one row per standard, element columns
      values <- lapply(seq_len(nrow(values)), function(i) unlist(values[i, , drop = TRUE]))
    }
  }
  stopifnot(length(standards) == length(values))

  obs_list <- lapply(standards, function(st)
    xrf_observed_mass(st, beam_energy_kev, detector_type = detector_type, be_window_um = be_window_um,
                      dead_layer_um = dead_layer_um, active_thickness_um = active_thickness_um,
                      air_path_cm = air_path_cm, atmosphere = atmosphere, window = window,
                      efficiency = efficiency, excitation = excitation,
                      excitation_weighting = excitation_weighting, coster_kronig = coster_kronig,
                      tube = tube, self_absorption = self_absorption, matrix = matrix,
                      incidence_deg = incidence_deg, takeoff_deg = takeoff_deg,
                      normalization = normalization))

  els <- unique(unlist(lapply(values, names)))
  k <- stats::setNames(rep(NA_real_, length(els)), els)
  npts <- stats::setNames(rep(0L, length(els)), els)
  for (el in els) {
    o <- c(); t <- c()
    for (j in seq_along(standards)) {
      om <- obs_list[[j]]
      oi <- om$observed_mass[om$element == el]
      ti <- suppressWarnings(as.numeric(values[[j]][el]))
      if (length(oi) == 1 && is.finite(oi) && oi > 0 && length(ti) == 1 && is.finite(ti)) {
        o <- c(o, oi); t <- c(t, ti)
      }
    }
    npts[el] <- length(o)
    # through-origin slope: calibration is applied multiplicatively (observed_mass * k), so an intercept
    # could not be reproduced. "none" is the ordinary (homoscedastic) through-origin slope sum(o t)/sum(o^2),
    # which effectively weights each standard by o^2 and is dominated by the highest-level standard;
    # "relative" uses inverse-variance weights w = 1/o^2 (proportional error), giving k = mean(t/o) so
    # standards spanning orders of magnitude in level contribute evenly.
    if (length(o) >= 1) {
      k[el] <- if (weighting == "relative") sum(t / o) / length(o) else sum(o * t) / sum(o * o)
    }
  }
  attr(k, "n_standards") <- npts
  # record ALL the sensitivity-scale settings so xrf_observed_mass can warn if they are not matched on the unknowns
  attr(k, "settings") <- .xrf_calibration_settings(self_absorption, matrix, normalization, beam_energy_kev,
                                                   detector_type, efficiency, excitation, excitation_weighting,
                                                   coster_kronig, tube)
  k
}
