
#' Physical parameters of electron shell transitions
#'
#' See Details
#'
#' Absoroption jump ratios: vectorized over both element and shell.
#' The absorption jump (often written \code{r} sub [K, L1, L2, L3])
#' can be used to calculate the probability that an electron from \code{shell} will be
#' ejected (J sub [K, L1, L2, L3]; the absorption jump ratio).
#'
#' Transition probabilities are the probability that a specified transition takes
#' place rather than another in that shell. Generally given the letter "g" sub [KL2, etc.].
#' See also \link{x_ray_emission_probabilities}.
#'
#' Fluorescence yields are the probability of emission of x-radiation rather than
#' an Auger electron. Generally given the letter omega sub [K, L1, L2, L3]. For
#' L shells, there are also Coster-Kronig (1935) transitions.
#'
#' xrf_relative_peak_intensity attempts to combine these functions to estimate
#' the relative intensities of various peaks. This is an important input for
#' deconvolution, and needs to be tweaked based on sample size.
#'
#'
#' @param element An element symbol or atomic number
#' @param shell A valid electron shell (K, L1, L2, L3). Invalid shells are silently ignored.
#' @param trans A transition (e.g., KL2).
#' @param coster_kronig_trans A sub-transition according to Coster-Kronig (1935).
#' @param beam_energy_kev The beam energy used to excite the electrons. Currently unused, but
#'   should be incorporated into the calculation.
#' @param method Only the Viegele (1973) method is currently supported for xrf_absorption_jump,
#'   and EADL97 tables are used for other parameters.
#'
#' @return A vector of numerics
#' @export
#'
#' @references
#' D.E. Cullen, et al., "Tables and Graphs of Atomic Subshell and
#' Relaxation Data Derived from the LLNL Evaluated Atomic Data Library
#' (EADL), Z = 1 - 100," Lawrence Livermore National Laboratory, UCRL-50400,
#' Vol. 30, October 1991. \url{http://www-nds.iaea.org/epdl97/libsall.htm}.
#'
#' Coster, D. and Kronig De L. (1935). New type of auger effect and its influence on the x-ray spectrum.
#' Pysica 2:13-24. \url{https://doi.org/10.1016/S0031-8914(35)90060-X}.
#'
#' Veigele, W.J. (1973) Atomic and Nuclear Data Tables 5:51-111. p 54 and 55.
#' \url{https://doi.org/10.1016/S0092-640X(73)80015-4}.
#'
#' @examples
#' # absorption jump ratio
#' xrf_absorption_jump("Pb", c("K", "L1", "L2", "L3"))
#'
#' # probability that a K shell electron will be ejected
#' rK <- xrf_absorption_jump("Pb", c("K", "L1", "L2", "L3"))
#' JK <- (rK - 1) / rK
#' xrf_absorption_jump_ratio("Pb", c("K", "L1", "L2", "L3"))
#'
xrf_absorption_jump <- function(element, shell, method = "viegele") {
  method <- match.arg(method)

  if(is.character(element)) {
    z <- match(element, all_elements)
  } else {
    z <- as.integer(element)
  }

  stopifnot(
    all(z > 0, z < 100),
    (length(element) == length(shell)) || (length(element) == 1) || (length(shell) == 1)
  )

  # constant for L1, L2, a / z + c for K, L3
  a <- c("K" = 125, "L1" = 1, "L2" = 1, "L3" = 80.0)[shell]
  b <- c("K" = 1, "L1" = Inf, "L2" = Inf, "L3" = 1)[shell]
  c <- c("K" = 3.5, "L1" = 1.2, "L2" = 1.4, "L3" = 1.5)[shell]

  unname((a / (b * z)) + c)
}

#' @rdname xrf_absorption_jump
#' @export
xrf_absorption_jump_ratio <- function(element, shell, method = "viegele") {
  method <- match.arg(method)
  rK <- xrf_absorption_jump(element = element, shell = shell, method = method)
  (rK - 1) / rK
}

#' @rdname xrf_absorption_jump
#' @export
xrf_transition_probability <- function(element, trans, method = "EADL97") {
  method <- match.arg(method)

  if(is.numeric(element)) {
    element <- all_elements[element]
  }
  stopifnot(is.character(trans))

  # Initialize cache on first use
  .init_physics_cache()

  # Vectorized named vector lookup
  keys <- paste(element, trans, sep = ":")
  unname(.xrf_cache$emission_prob_vec[keys])
}

#' @rdname xrf_absorption_jump
#' @export
xrf_fluorescence_yield <- function(element, shell, method = "EADL97") {
  method <- match.arg(method)

  if(is.numeric(element)) {
    element <- all_elements[element]
  }
  stopifnot(is.character(shell))

  # Initialize cache on first use
  .init_physics_cache()

  # Vectorized named vector lookup
  keys <- paste(element, shell, sep = ":")
  unname(.xrf_cache$fluorescence_yield_vec[keys])
}

#' @rdname xrf_absorption_jump
#' @export
xrf_coster_kronig_probability <- function(element, shell, coster_kronig_trans, method = "EADL97") {
  method <- match.arg(method)
  if(is.numeric(element)) {
    element <- all_elements[element]
  }
  stopifnot(is.character(shell))

  # Initialize cache on first use
  .init_physics_cache()

  # Vectorized named vector lookup
  keys <- paste(element, shell, coster_kronig_trans, sep = ":")
  unname(.xrf_cache$coster_kronig_vec[keys])
}

#' Partial photoionization cross section at a given energy
#'
#' Interpolates the per-subshell photoelectric cross section from \link{x_ray_cross_sections}
#' (K, L1, L2, L3, M1-M3; EPDL97; 5-100 keV grid) at an arbitrary incident energy, using log-log
#' interpolation. Below a subshell's absorption edge the tabulated cross section is \code{NA}, so
#' this returns \code{NA} for incident energies below the edge (the subshell cannot be ionized).
#' Above the tabulated range the cross section is extrapolated as a power law (photoabsorption falls
#' roughly as \eqn{E^{-3}} above an edge). Returns \code{NA} for element/shell combinations absent
#' from the table (e.g. a K edge above the 100 keV grid, such as U).
#'
#' @param element An element symbol or atomic number (vectorized).
#' @param shell A subshell ("K", "L1", "L2", "L3", "M1", "M2", "M3").
#' @param energy_kev The incident (beam) energy in keV.
#'
#' @return A numeric vector of cross sections (same arbitrary units as \link{x_ray_cross_sections};
#'   only ratios within an element are used downstream, so the unit cancels).
#' @export
#'
#' @examples
#' # Pb K vs L3 excitation strength grows apart with energy
#' xrf_photoionization_cross_section("Pb", c("K", "L3"), 90)
#'
xrf_photoionization_cross_section <- function(element, shell, energy_kev) {
  if (is.numeric(element)) element <- all_elements[element]
  .init_physics_cache()

  n <- max(length(element), length(shell), length(energy_kev))
  element <- rep_len(element, n)
  shell <- rep_len(shell, n)
  energy_kev <- rep_len(energy_kev, n)

  key <- paste(element, shell, sep = ":")
  out <- rep(NA_real_, n)
  for (k in unique(key)) {
    g <- .xrf_cache$cross_section_split[[k]]
    if (is.null(g) || nrow(g) == 0) next
    idx <- which(key == k)
    e0 <- energy_kev[idx]
    lo <- min(g$energy_kev)
    hi <- max(g$energy_kev)
    v <- rep(NA_real_, length(e0))           # below the edge (below lo) stays NA

    within <- is.finite(e0) & e0 >= lo & e0 <= hi
    if (any(within)) {
      v[within] <- exp(suppressWarnings(   # edge-doubled EPDL grid points are collapsed by approx()
        stats::approx(log(g$energy_kev), log(g$cross_section), xout = log(e0[within]))$y))
    }
    above <- is.finite(e0) & e0 > hi
    if (any(above) && nrow(g) >= 2) {
      m <- nrow(g)                            # power-law (log-log) extrapolation beyond the grid
      slope <- (log(g$cross_section[m]) - log(g$cross_section[m - 1])) /
        (log(g$energy_kev[m]) - log(g$energy_kev[m - 1]))
      v[above] <- exp(log(g$cross_section[m]) + slope * (log(e0[above]) - log(g$energy_kev[m])))
    }
    out[idx] <- v
  }
  out
}

# electron-shell occupancies, for the (Gryzinski/Bethe) electron-impact ionization cross section
.xrf_shell_occupancy <- c(K = 2, L1 = 2, L2 = 2, L3 = 4, M1 = 2, M2 = 2, M3 = 4, M4 = 4, M5 = 6)

# subshell index used by the Bote-Salvat tables: 1->K, 2->L1, ... 9->M5
.xrf_bs_subshell <- c(K = 1, L1 = 2, L2 = 3, L3 = 4, M1 = 5, M2 = 6, M3 = 7, M4 = 8, M5 = 9)

# Bote & Salvat (2008) electron-impact ionization cross section (cm^2), ported verbatim from the
# public-domain NIST BoteSalvatICX.jl. z, subshell (1..9) and energy_ev are vectorized; the edge
# energy is the Bote-Salvat tabulated value. Returns 0 below threshold or where data are unavailable.
.xrf_bote_salvat_icx <- function(z, subshell, energy_ev) {
  .init_physics_cache()
  bs <- .xrf_cache$bote_salvat
  idx <- match(paste(z, subshell, sep = ":"), .xrf_cache$bote_salvat_key)
  out <- rep(0, length(idx))
  ok <- !is.na(idx) & is.finite(energy_ev)
  if (!any(ok)) return(out)
  i <- idx[ok]
  E <- energy_ev[ok]
  ov <- E / bs$edge_ev[i]
  REV <- 5.10998918e5          # electron rest energy (eV)
  a0cm <- 5.291772108e-9       # Bohr radius (cm)
  xione <- rep(0, length(i))

  lo <- ov > 1 & ov <= 16      # low-energy analytical fit (A coefficients)
  if (any(lo)) {
    li <- i[lo]; o <- ov[lo]; opu <- 1 / (1 + o)
    ffitlo <- bs$a1[li] + bs$a2[li] * o + opu * (bs$a3[li] + opu^2 * (bs$a4[li] + opu^2 * bs$a5[li]))
    xione[lo] <- (o - 1) * (ffitlo / o)^2
  }
  hi <- ov > 16                # relativistic high-energy fit (G coefficients)
  if (any(hi)) {
    hii <- i[hi]; o <- ov[hi]; Eh <- E[hi]
    beta2 <- (Eh * (Eh + 2 * REV)) / ((Eh + REV)^2)
    xx <- sqrt(Eh * (Eh + 2 * REV)) / REV
    ffitup <- ((2 * log(xx) - beta2) * (1 + bs$g1[hii] / xx)) + bs$g2[hii] +
      bs$g3[hii] * sqrt(REV / (Eh + REV)) + bs$g4[hii] / xx
    xione[hi] <- ((bs$anlj[hii] / beta2 * o) / (o + bs$be[hii])) * ffitup
  }
  out[ok] <- 4 * pi * a0cm^2 * pmax(xione, 0)
  out
}

#' Electron-impact ionization cross section (relative)
#'
#' A Bethe-type inner-shell electron-impact ionization cross section, used to weight subshell vacancy
#' production for electron excitation (SEM-EDS) rather than the photoionization cross section used
#' for photon XRF:
#' \deqn{Q_{shell}(E_0) \propto z_{shell}\, \ln(U) / (U\, E_{edge}^2), \quad U = E_0 / E_{edge},}
#' with \eqn{z_{shell}} the subshell electron occupancy and \eqn{U} the overvoltage. \eqn{Q} is zero
#' at threshold (\eqn{U=1}), peaks near \eqn{U\approx e}, then declines -- the qualitative behaviour
#' that makes low-kV beams excite the M series of heavy elements far more strongly than their L or K
#' series. Only relative values matter here (they multiply the emission probabilities and are
#' renormalized per element), so the constant prefactor is dropped. This is an approximate model; for
#' absolute electron-probe work the Bote-Salvat (2008) or Casnati (1982) cross sections are preferred.
#'
#' @param element An element symbol or atomic number.
#' @param shell A subshell ("K", "L1"-"L3", "M1"-"M5").
#' @param beam_energy_kev The electron beam energy (accelerating voltage, keV).
#'
#' @param method "bote-salvat" (default) uses the Bote & Salvat (2008) distorted-wave / plane-wave
#'   Born cross sections -- the state-of-the-art model (as in PENELOPE / DTSA-II), with the
#'   coefficient tables in \link{x_ray_bote_salvat}; it is intrinsically relativistic. "gryzinski"
#'   uses the Gryzinski (1965) classical binary-encounter g-function; "bethe" the simpler
#'   \eqn{\ln(U)/U} form.
#' @param relativistic Applies only to \code{method \%in\% c("gryzinski", "bethe")}: if TRUE
#'   (default) apply a relativistic velocity correction \eqn{2E_0/(\beta^2 m_e c^2)} (\eqn{\to 1}
#'   non-relativistically, growing at high beam energy). Ignored for "bote-salvat" (already
#'   relativistic).
#'
#' @return A vector of cross sections (0 below threshold): cm^2 for "bote-salvat", relative units
#'   for "gryzinski"/"bethe".
#' @export
#'
#' @references
#' Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the
#' distorted-wave and plane-wave Born approximations. Phys. Rev. A 77, 042701.
#'
#' Gryzinski, M. (1965). Classical theory of atomic collisions. I. Theory of inelastic collisions.
#' Physical Review 138, A336. \url{https://doi.org/10.1103/PhysRev.138.A336}.
#'
#' @examples
#' # at 20 kV the M shells of Au have far higher overvoltage than the L shells
#' xrf_electron_ionization_cross_section("Au", c("L3", "M5"), 20)
#'
xrf_electron_ionization_cross_section <- function(element, shell, beam_energy_kev,
                                                  method = c("bote-salvat", "gryzinski", "bethe"),
                                                  relativistic = TRUE) {
  method <- match.arg(method)
  if (is.numeric(element)) element <- all_elements[element]
  .init_physics_cache()
  n <- max(length(element), length(shell), length(beam_energy_kev))
  element <- rep_len(element, n)
  shell <- rep_len(shell, n)
  beam_energy_kev <- rep_len(beam_energy_kev, n)

  if (method == "bote-salvat") {
    q <- .xrf_bote_salvat_icx(match(element, all_elements), unname(.xrf_bs_subshell[shell]),
                              beam_energy_kev * 1000)
    q[!is.finite(q)] <- 0
    return(q)
  }

  e_edge <- unname(.xrf_cache$edge_kev_vec[paste(element, shell, sep = ":")])
  z <- unname(.xrf_shell_occupancy[shell])
  u <- beam_energy_kev / e_edge
  q <- if (method == "gryzinski") {
    # Gryzinski (1965) binary-encounter g-function: rises from 0 at U=1, peaks ~U=3-4, then decays.
    (z / e_edge^2) * (1 / u) * ((u - 1) / (u + 1))^1.5 *
      (1 + (2 / 3) * (1 - 1 / (2 * u)) * log(2.7 + sqrt(pmax(u - 1, 0))))
  } else {
    z * log(u) / (u * e_edge^2)                       # simple Bethe form
  }
  if (relativistic) {
    mc2 <- 510.999                                    # electron rest energy (keV)
    beta2 <- 1 - (mc2 / (beam_energy_kev + mc2))^2
    q <- q * (2 * beam_energy_kev) / (beta2 * mc2)    # -> 1 as E0 -> 0; relativistic rise at high E0
  }
  q[!is.finite(q) | !is.finite(u) | u <= 1] <- 0
  q
}

# Redistribute primary L-subshell vacancies via Coster-Kronig transfer before radiative decay.
# Given per-row primary vacancy factors for L1/L2/L3 (nL1/nL2/nL3), returns the *effective* vacancy
# factor for each row's shell (L rows get the cascaded value; non-L rows are returned unchanged).
# The tabulated L-subshell fluorescence yields already fold in the Coster-Kronig loss channel, so
# L1 is NOT depleted here; the transferred vacancies are simply *added* to L2/L3:
#   V1_eff = nL1;  V2_eff = nL2 + f12*nL1;  V3_eff = nL3 + f23*V2_eff + f13*nL1.
.xrf_ck_cascade <- function(element, shell, n_row, nL1, nL2, nL3) {
  f12 <- dplyr::coalesce(xrf_coster_kronig_probability(element, "L1", "f12"), 0)
  f13 <- dplyr::coalesce(xrf_coster_kronig_probability(element, "L1", "f13"), 0)
  f23 <- dplyr::coalesce(xrf_coster_kronig_probability(element, "L2", "f23"), 0)
  V1 <- nL1
  V2 <- nL2 + f12 * nL1
  V3 <- nL3 + f23 * V2 + f13 * nL1
  out <- n_row
  out[shell == "L1"] <- V1[shell == "L1"]
  out[shell == "L2"] <- V2[shell == "L2"]
  out[shell == "L3"] <- V3[shell == "L3"]
  out
}

#' @rdname xrf_absorption_jump
#' @param excitation_weighting How to weight the primary vacancy production per subshell.
#'   "cross_section" (default) uses the energy-dependent partial photoionization cross section
#'   \eqn{\sigma_{shell}(E_0)} from \link{xrf_photoionization_cross_section} at
#'   \code{beam_energy_kev}, so the K:L branching correctly depends on the beam energy; it falls
#'   back to the (energy-independent) absorption-jump ratio only where the cross section is
#'   unavailable. "jump" uses the absorption-jump ratio everywhere (the classic, beam-independent
#'   behaviour).
#' @param coster_kronig If TRUE (default) redistribute L1/L2 vacancies down the L subshells via
#'   Coster-Kronig transfer (\link{x_ray_coster_kronig_probabilities}) before radiative decay,
#'   correcting the relative intensities of the L family for heavier elements.
#' @param excitation "photon" (default, tube/synchrotron XRF) weights subshell production by the
#'   photoionization cross section (or jump ratio); "electron" (SEM-EDS) uses the electron-impact
#'   ionization cross section \link{xrf_electron_ionization_cross_section} instead, so the M series
#'   of heavy elements is correctly enhanced at low beam energy.
#' @export
xrf_relative_peak_intensity <- function(element, shell, trans, beam_energy_kev = 50,
                                        excitation_weighting = c("cross_section", "jump"),
                                        coster_kronig = TRUE,
                                        excitation = c("photon", "electron")) {
  # The relative intensity of a primary-fluorescence line is
  #   (subshell vacancy production at E0) x (omega * radiative branching).
  # The EADL97 emission probability (via xrf_transition_probability) ALREADY includes the
  # fluorescence yield omega (per shell, sum_trans(emission_probability) == omega), so we must NOT
  # multiply by xrf_fluorescence_yield() again (that historic double-count distorted every
  # cross-shell ratio by omega_shell). The subshell vacancy production is either the partial
  # photoionization cross section sigma_shell(E0) (energy-dependent; #8) or, classically, the
  # absorption-jump ratio; L-subshell vacancies are optionally redistributed by Coster-Kronig (#9).
  excitation_weighting <- match.arg(excitation_weighting)
  excitation <- match.arg(excitation)
  if (is.numeric(element)) element <- all_elements[element]

  n <- max(length(element), length(shell), length(trans))
  element <- rep_len(element, n)
  shell <- rep_len(shell, n)
  trans <- rep_len(trans, n)

  if (excitation == "electron") {
    # Electron-impact excitation: subshell production follows the (overvoltage-dependent) electron
    # ionization cross section rather than photoionization. Relaxation (Coster-Kronig + emission
    # probabilities) is unchanged.
    n_shell <- xrf_electron_ionization_cross_section(element, shell, beam_energy_kev)
    if (coster_kronig) {
      v_eff <- .xrf_ck_cascade(
        element, shell, n_shell,
        xrf_electron_ionization_cross_section(element, "L1", beam_energy_kev),
        xrf_electron_ionization_cross_section(element, "L2", beam_energy_kev),
        xrf_electron_ionization_cross_section(element, "L3", beam_energy_kev)
      )
      is_l <- shell %in% c("L1", "L2", "L3")
      n_shell <- ifelse(is_l, v_eff, n_shell)
    }
    return(n_shell * xrf_transition_probability(element, trans))
  }

  jr <- xrf_absorption_jump_ratio(element, shell)

  if (excitation_weighting == "jump") {
    n_shell <- jr
    if (coster_kronig) {
      n_shell <- .xrf_ck_cascade(
        element, shell, n_shell,
        xrf_absorption_jump_ratio(element, "L1"),
        xrf_absorption_jump_ratio(element, "L2"),
        xrf_absorption_jump_ratio(element, "L3")
      )
    }
  } else {
    # K/M rows: fall back to the jump ratio only where sigma is unavailable (e.g. a K edge above
    # the 100 keV table). Below-edge rows are removed downstream by the edge cut in xrf_energies().
    n_shell <- dplyr::coalesce(
      xrf_photoionization_cross_section(element, shell, beam_energy_kev), jr
    )
    if (coster_kronig) {
      # L primary vacancies: sigma with a 0 (not jump-ratio) fallback, so an L subshell below the
      # beam energy contributes no Coster-Kronig transfer.
      v_eff <- .xrf_ck_cascade(
        element, shell, n_shell,
        dplyr::coalesce(xrf_photoionization_cross_section(element, "L1", beam_energy_kev), 0),
        dplyr::coalesce(xrf_photoionization_cross_section(element, "L2", beam_energy_kev), 0),
        dplyr::coalesce(xrf_photoionization_cross_section(element, "L3", beam_energy_kev), 0)
      )
      is_l <- shell %in% c("L1", "L2", "L3")
      n_shell <- ifelse(is_l, v_eff, n_shell)
    }
  }

  n_shell * xrf_transition_probability(element, trans)
}

#' Get XRF energies for selected elements
#'
#' @param elements Elements or element lists (passed to \link{xrf_element_list}).
#' @param beam_energy_kev Beam energy (keV). For \code{excitation = "photon"} this is the tube
#'   endpoint used to select which lines can be excited (a line is kept when its absorption edge is
#'   at or below \code{beam_energy_kev}). For \code{excitation = "electron"} it is the accelerating
#'   voltage used with \code{overvoltage_min}.
#' @param excitation Excitation mode. "photon" (default, tube/synchrotron XRF) keeps lines whose
#'   edge is below the beam energy and weights production by the photoionization cross section.
#'   "electron" (SEM-EDS) keeps lines whose overvoltage \eqn{U = beam / edge} is at least
#'   \code{overvoltage_min} and weights production by the electron-impact ionization cross section
#'   (\link{xrf_electron_ionization_cross_section}), so the M series of heavy elements is correctly
#'   enhanced at low beam energy.
#' @param overvoltage_min Minimum overvoltage \eqn{U = beam / edge} for a line to be kept in
#'   \code{excitation = "electron"} mode. Useful excitation typically needs \eqn{U \gtrsim 1.2-1.5}.
#' @param excitation_weighting,coster_kronig Passed to \link{xrf_relative_peak_intensity}: how to
#'   weight subshell vacancy production ("cross_section", the energy-dependent default, or "jump"),
#'   and whether to apply the Coster-Kronig L-subshell cascade.
#' @param min_relative_intensity Smallest peak to include, relative to the maximum peak for a given
#'   element.
#' @param ... Used to further \link[dplyr]{filter} the result.
#'
#' @return A subset of \link{x_ray_xrf_energies} with some additional information based on the beam energy.
#' @export
#'
#' @examples
#' xrf_energies("major", 25)
#' xrf_energies("major", 20, excitation = "electron", overvoltage_min = 1.5)
#'
xrf_energies <- function(elements = "everything", beam_energy_kev = 50, ...,
                         excitation = c("photon", "electron"), overvoltage_min = 1,
                         excitation_weighting = c("cross_section", "jump"),
                         coster_kronig = TRUE,
                         min_relative_intensity = 0.01) {
  excitation <- match.arg(excitation)
  excitation_weighting <- match.arg(excitation_weighting)
  elements <- xrf_element_list(elements)

  # line-selection cut: photon keeps edge <= beam; electron keeps overvoltage U = beam/edge >= U_min
  edge_cut <- if (excitation == "electron") beam_energy_kev / overvoltage_min else beam_energy_kev

  xrftools::x_ray_xrf_energies %>%
    dplyr::mutate(
      relative_peak_intensity = xrf_relative_peak_intensity(
        .data$element, .data$edge, .data$trans, !!beam_energy_kev,
        excitation_weighting = !!excitation_weighting, coster_kronig = !!coster_kronig,
        excitation = !!excitation
      )
    ) %>%
    filter(!is.na(.data$relative_peak_intensity)) %>%
    dplyr::group_by(.data$element) %>%
    dplyr::mutate(relative_peak_intensity = .data$relative_peak_intensity / max(.data$relative_peak_intensity)) %>%
    dplyr::filter(
      .data$element %in% !!elements,
      .data$edge_kev <= !!edge_cut,
      .data$relative_peak_intensity >= min_relative_intensity,
      ...
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$z, dplyr::desc(.data$relative_peak_intensity)) %>%
    dplyr::select("element", "trans", "trans_siegbahn", "energy_kev", "relative_peak_intensity", dplyr::everything())
}
