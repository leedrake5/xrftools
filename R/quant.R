
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

# Warn if a calibration vector was fitted with different observed-mass settings than are in use.
.xrf_check_calibration_settings <- function(calibration, self_absorption, matrix, normalization) {
  st <- attr(calibration, "settings")
  if (is.null(st)) return(invisible())
  bad <- character(0)
  if (!identical(st$self_absorption, self_absorption)) bad <- c(bad, "self_absorption")
  if (!identical(st$normalization, normalization)) bad <- c(bad, "normalization")
  if (self_absorption == "generalized" && !isTRUE(all.equal(st$matrix, matrix))) bad <- c(bad, "matrix")
  if (length(bad) > 0) {
    warning("calibration was fitted with different ", paste(bad, collapse = ", "),
            " settings than used here; the absolute values will be wrong -- match them.", call. = FALSE)
  }
  invisible()
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

  comp <- switch(
    by,
    compton  = "scatter_compton",
    rayleigh = "scatter_rayleigh",
    scatter  = c("scatter_compton", "scatter_rayleigh")
  )
  denom <- sum(peaks$peak_area[peaks$element %in% comp], na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) {
    stop("No positive ", by, " scatter area found; run the deconvolution with a `tube` so scatter ",
         "templates are fitted.")
  }
  peaks$peak_area_norm <- peaks$peak_area / denom
  if ("peak_area_se" %in% colnames(peaks)) peaks$peak_area_norm_se <- peaks$peak_area_se / denom
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
#' @param detector_type,be_window_um,dead_layer_um Detector parameters for the efficiency term
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
                               be_window_um = NULL, dead_layer_um = NULL, efficiency = TRUE,
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
    n_tube <- xrf_tube_spectrum(tube, grid)
  }

  purrr::map_dfr(elements, function(el) {
    lines <- xe[xe$element == el & xe$edge_kev <= beam_energy_kev, , drop = FALSE]
    if (nrow(lines) == 0) {
      return(tibble::tibble(element = el, primary_energy_kev = NA_real_, sensitivity = NA_real_))
    }
    rel <- if (poly) {
      excit <- vapply(lines$edge, function(sh) {
        s <- xrf_photoionization_cross_section(el, sh, grid); s[!is.finite(s)] <- 0
        sum(s * n_tube) * de
      }, numeric(1))
      excit * xrf_transition_probability(el, lines$trans)   # omega * branching (from emission prob)
    } else {
      xrf_relative_peak_intensity(el, lines$edge, lines$trans, beam_energy_kev,
                                  excitation_weighting = excitation_weighting,
                                  coster_kronig = coster_kronig, excitation = excitation)
    }
    eff <- if (efficiency) {
      xrf_detector_efficiency(lines$energy_kev, detector_type = detector_type,
                              be_window_um = be_window_um, dead_layer_um = dead_layer_um)
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
#' @param detector_type,be_window_um,dead_layer_um,efficiency,excitation,excitation_weighting,coster_kronig,tube
#'   Passed to \link{xrf_fp_sensitivity}. Supplying \code{tube = }\link{xrf_tube} integrates the
#'   primary excitation over the polychromatic tube spectrum instead of a single energy.
#' @param self_absorption Apply the iterative self-absorption correction?
#' @param secondary_fluorescence Add secondary (enhancement) fluorescence?
#' @param tertiary_fluorescence Add third-order (k -> j -> i chain) enhancement fluorescence
#'   (implies \code{secondary_fluorescence}). Usually a small (<~1-2\%) correction.
#' @param incidence_deg,takeoff_deg Beam incidence / detector take-off angles (degrees).
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
#'
#' @return A tibble with \code{element}, \code{primary_energy_kev}, \code{peak_area},
#'   \code{sensitivity}, and \code{concentration}.
#' @export
#'
xrf_quantify <- function(object, beam_energy_kev, detector_type = NULL,
                         be_window_um = NULL, dead_layer_um = NULL, efficiency = TRUE,
                         excitation = c("photon", "electron"),
                         excitation_weighting = c("cross_section", "jump"), coster_kronig = TRUE,
                         tube = NULL, self_absorption = TRUE, secondary_fluorescence = FALSE,
                         tertiary_fluorescence = FALSE, incidence_deg = 45, takeoff_deg = 45,
                         dark_matrix = NULL, oxide = FALSE, total = 1, iterations = 12) {
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

  conc <- renorm(q$peak_area / q$sensitivity)

  if (isTRUE(self_absorption) || isTRUE(secondary_fluorescence)) {
    sin_in <- sin(incidence_deg * pi / 180)
    sin_out <- sin(takeoff_deg * pi / 180)
    Ei <- q$primary_energy_kev
    dm_el <- names(dark_matrix); dm_frac <- as.numeric(dark_matrix)
    mu_at <- function(energy) {   # sample mass attenuation at one energy over the current mix
      m <- sum(conc * vapply(els, function(el) .xrf_element_pe_mu(el, energy), numeric(1)))
      o_c <- sum(conc * r)                                  # derived oxygen
      if (o_c > 0) m <- m + o_c * .xrf_element_pe_mu("O", energy)
      if (length(dm_el)) m <- m + sum(dm_frac * vapply(dm_el, function(el) .xrf_element_pe_mu(el, energy), numeric(1)))
      m
    }
    # precompute the (concentration-independent) atomic factors for secondary fluorescence
    if (isTRUE(secondary_fluorescence)) {
      omega_j <- xrf_fluorescence_yield(els, q$primary_shell)
      edge_i <- q$primary_edge_kev
      tau_E0 <- vapply(els, function(el) .xrf_element_pe_mu(el, beam_energy_kev), numeric(1))
      # tau_ij[i, j] = photoelectric absorption of element i at element j's line energy
      tau_ij <- outer(seq_along(els), seq_along(els),
                      Vectorize(function(a, b) .xrf_element_pe_mu(els[a], Ei[b])))
    }

    for (iter in seq_len(iterations)) {
      mu_e0 <- mu_at(beam_energy_kev)
      mu_ei <- vapply(Ei, mu_at, numeric(1))
      m_corr <- if (isTRUE(self_absorption)) 1 / (mu_e0 / sin_in + mu_ei / sin_out) else 1

      # Secondary/tertiary (enhancement) fluorescence, Shiraiwa-Fujino. B[i, j] is the enhancement of
      # element i per unit mass of element j (element j's line must clear i's edge). Secondary
      # enhancement is B %*% C; tertiary uses j's own enhanced intensity, B %*% (C*(1 + B%*%C)),
      # capturing the k -> j -> i chain.
      enh <- rep(0, length(els))
      if (isTRUE(secondary_fluorescence)) {
        B <- matrix(0, length(els), length(els))
        for (i in seq_along(els)) {
          for (j in seq_along(els)) {
            if (j == i || !is.finite(Ei[j]) || Ei[j] <= edge_i[i]) next
            muj <- mu_ei[j]                       # sample attenuation at j's line energy
            g <- 0.5 * muj * ((sin_in / mu_e0) * log(1 + mu_e0 / (sin_in * muj)) +
                                (sin_out / mu_ei[i]) * log(1 + mu_ei[i] / (sin_out * muj)))
            B[i, j] <- omega_j[j] * (tau_E0[j] / tau_E0[i]) * (tau_ij[i, j] / muj) * g
          }
        }
        enh2 <- as.vector(B %*% conc)
        enh <- if (isTRUE(tertiary_fluorescence)) as.vector(B %*% (conc * (1 + enh2))) else enh2
      }

      conc <- renorm(q$peak_area / (q$sensitivity * m_corr * (1 + enh)))
    }
  }

  q$concentration <- conc
  out <- tibble::as_tibble(q[, c("element", "primary_energy_kev", "peak_area", "sensitivity", "concentration")])
  if (oxide_on && any(r > 0)) {   # report the stoichiometric oxygen
    out <- dplyr::bind_rows(out, tibble::tibble(
      element = "O", primary_energy_kev = NA_real_, peak_area = NA_real_,
      sensitivity = NA_real_, concentration = sum(conc * r)
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
#'
#' @param object A \code{deconvolution_fit} or its \code{peaks} data frame (needs \code{element},
#'   \code{peak_area}; scatter pseudo-elements are used for \code{normalization} and then dropped).
#' @param beam_energy_kev,detector_type,be_window_um,dead_layer_um,efficiency,excitation,excitation_weighting,coster_kronig,tube
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
#' @return A tibble with \code{element}, \code{primary_energy_kev}, \code{peak_area},
#'   \code{sensitivity} and \code{observed_mass} (un-normalized; the sum is diagnostic, not forced).
#' @export
#'
xrf_observed_mass <- function(object, beam_energy_kev, detector_type = NULL, be_window_um = NULL,
                              dead_layer_um = NULL, efficiency = TRUE,
                              excitation = c("photon", "electron"),
                              excitation_weighting = c("cross_section", "jump"),
                              coster_kronig = TRUE, tube = NULL,
                              self_absorption = c("none", "generalized"), matrix = "silicate",
                              incidence_deg = 45, takeoff_deg = 45,
                              normalization = c("none", "compton", "rayleigh", "scatter"),
                              calibration = NULL) {
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
                          efficiency = efficiency, excitation = match.arg(excitation),
                          excitation_weighting = match.arg(excitation_weighting),
                          coster_kronig = coster_kronig, tube = tube)
  q$primary_energy_kev <- s$primary_energy_kev
  q$sensitivity <- s$sensitivity

  obs <- q$peak_area / q$sensitivity                 # matrix-free observed areal mass (relative)

  if (self_absorption == "generalized") {            # remove per-element depth via a generic matrix
    comp <- .xrf_resolve_matrix(matrix)
    sin_in <- sin(incidence_deg * pi / 180); sin_out <- sin(takeoff_deg * pi / 180)
    mu_e0 <- .xrf_matrix_total_mu(comp, beam_energy_kev)
    mu_ei <- vapply(q$primary_energy_kev, function(e) .xrf_matrix_total_mu(comp, e), numeric(1))
    obs <- obs * (mu_e0 / sin_in + mu_ei / sin_out)
  }

  if (normalization != "none") {
    comp <- switch(normalization,
                   compton = "scatter_compton", rayleigh = "scatter_rayleigh",
                   scatter = c("scatter_compton", "scatter_rayleigh"))
    denom <- sum(peaks$peak_area[peaks$element %in% comp], na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) {
      stop("normalization = '", normalization, "' needs a fitted scatter peak; run the ",
           "deconvolution with a `tube`.")
    }
    obs <- obs / denom
  }
  if (!is.null(calibration)) {                       # optional absolute calibration
    if (is.null(names(calibration))) {
      stop("`calibration` must be a named numeric vector keyed by element symbol.")
    }
    .xrf_check_calibration_settings(calibration, self_absorption, matrix, normalization)
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

  tibble::as_tibble(data.frame(
    element = q$element, primary_energy_kev = q$primary_energy_kev, peak_area = q$peak_area,
    sensitivity = q$sensitivity, observed_mass = obs, stringsAsFactors = FALSE
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
#' @param beam_energy_kev,detector_type,be_window_um,dead_layer_um,efficiency,excitation,excitation_weighting,coster_kronig,tube,self_absorption,matrix,incidence_deg,takeoff_deg,normalization
#'   Passed to \link{xrf_observed_mass} to compute the observed values for each standard (use the
#'   settings you will apply to unknowns).
#'
#' @return A named numeric vector of per-element calibration factors \eqn{k_i}, with an
#'   \code{"n_standards"} attribute (points used per element) and a \code{"settings"} attribute
#'   (the observed-mass settings used). Pass it as \code{calibration =} to \link{xrf_observed_mass}.
#' @export
#'
xrf_calibrate <- function(standards, values, beam_energy_kev, detector_type = NULL,
                          be_window_um = NULL, dead_layer_um = NULL, efficiency = TRUE,
                          excitation = c("photon", "electron"),
                          excitation_weighting = c("cross_section", "jump"), coster_kronig = TRUE,
                          tube = NULL, self_absorption = c("none", "generalized"),
                          matrix = "silicate", incidence_deg = 45, takeoff_deg = 45,
                          normalization = c("none", "compton", "rayleigh", "scatter")) {
  excitation <- match.arg(excitation)
  excitation_weighting <- match.arg(excitation_weighting)
  self_absorption <- match.arg(self_absorption)
  normalization <- match.arg(normalization)
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
                      dead_layer_um = dead_layer_um, efficiency = efficiency, excitation = excitation,
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
    # through-origin least-squares slope: calibration is applied multiplicatively
    # (observed_mass * k), so an intercept could not be reproduced anyway.
    if (length(o) >= 1) k[el] <- sum(o * t) / sum(o * o)
  }
  attr(k, "n_standards") <- npts
  # record the settings so xrf_observed_mass can warn if they are not matched on the unknowns
  attr(k, "settings") <- list(self_absorption = self_absorption, matrix = matrix,
                              normalization = normalization)
  k
}
