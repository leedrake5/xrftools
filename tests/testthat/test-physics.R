context("physics")

test_that("jump ratios are calculated correctly", {

  expect_identical(xrf_absorption_jump(1:99, "K"), (125.0 / 1:99) + 3.5)
  expect_identical(xrf_absorption_jump(1:99, "L1"), rep(1.2, 99))
  expect_identical(xrf_absorption_jump(1:99, "L2"), rep(1.4, 99))
  expect_identical(xrf_absorption_jump(1:99, "L3"), (80 / 1:99) + 1.5)

  # element symbol or Z as input
  expect_identical(
    xrf_absorption_jump(1:5, "K"),
    xrf_absorption_jump(c("H", "He", "Li", "Be", "B"), "K")
  )

  # vectorized over shell
  expect_identical(
    xrf_absorption_jump("Pb", c("K", "L1", "L2", "L3")),
    c((125.0 / 82) + 3.5, 1.2, 1.4, (80 / 82) + 1.5)
  )

  # vectorized over both
  expect_identical(
    xrf_absorption_jump(c("H", "Pb", "H", "Pb"), c("K", "L1", "L2", "L3")),
    c((125.0 / 1) + 3.5, 1.2, 1.4, (80 / 82) + 1.5)
  )

  # problems are caught
  expect_error(xrf_absorption_jump("Fm", "K"), "z < 100")
  expect_error(xrf_absorption_jump(c("H", "Be"), c("K", "L1", "L2")), "is not TRUE")
  expect_identical(xrf_absorption_jump("H", "not a shell"), NA_real_)

  rK <- xrf_absorption_jump("Pb", c("K", "L1", "L2", "L3"))
  expect_identical(
    xrf_absorption_jump_ratio("Pb", c("K", "L1", "L2", "L3")),
    (rK - 1) / rK
  )
})

test_that("transition probabilities are correct", {
  expect_identical(
    xrf_transition_probability("C", "KL2"),
    x_ray_emission_probabilities %>%
      dplyr::filter(element == "C") %>%
      dplyr::filter(trans == "KL2") %>%
      dplyr::pull(emission_probability)
  )

  expect_identical(
    xrf_transition_probability(c("C", "N", "O"), "KL2"),
    xrf_transition_probability(c(6, 7, 8), "KL2")
  )
})

test_that("fluorescence yields are correct", {
  expect_identical(
    xrf_fluorescence_yield("C", "K"),
    x_ray_fluorescence_yields %>%
      dplyr::filter(element == "C") %>%
      dplyr::filter(shell == "K") %>%
      dplyr::pull(fluorescence_yield)
  )

  expect_identical(
    xrf_fluorescence_yield(c("C", "N", "O"), "K"),
    xrf_fluorescence_yield(c(6, 7, 8), "K")
  )
})

test_that("relative peak intensity does not double-count the fluorescence yield", {
  # EADL97 emission probabilities already fold in omega, so summing the relative intensity over a
  # shell's transitions (divided by that shell's jump ratio) must recover the fluorescence yield.
  # This guards against the historical bug of multiplying by xrf_fluorescence_yield() a second time.
  # pin the classic weighting (jump ratio, no Coster-Kronig) so the invariant isolates the omega
  # question; cross-section weighting and CK deliberately reshape these ratios (tested separately).
  check_shell_invariant <- function(element, shell, tol) {
    trans_prefix <- if (shell == "K") "^K[LMNOP]" else paste0("^", shell, "[LMNOP]")
    trans <- x_ray_emission_probabilities %>%
      dplyr::filter(element == !!element, grepl(trans_prefix, trans), trans != "TOTAL") %>%
      dplyr::pull(trans)
    rel <- xrf_relative_peak_intensity(element, shell, trans,
                                       excitation_weighting = "jump", coster_kronig = FALSE)
    jump <- xrf_absorption_jump_ratio(element, shell)
    omega <- xrf_fluorescence_yield(element, shell)
    expect_equal(sum(rel) / jump, omega, tolerance = tol)
  }

  # K recovers omega tightly for all Z. L2/L3 recover omega tightly only for heavier elements
  # (for light Z the tabulated L radiative branches are a tiny, incompletely-summed fraction of
  # omega). L1 is excluded because Coster-Kronig redistribution means the tabulated L1 radiative
  # emission does not sum to omega_L1.
  for (el in c("Fe", "Mo", "Pb")) check_shell_invariant(el, "K", tol = 0.03)
  for (el in c("Mo", "Pb")) {
    check_shell_invariant(el, "L2", tol = 0.05)
    check_shell_invariant(el, "L3", tol = 0.05)
  }

  # And the fluorescence yield must NOT appear as an explicit factor.
  expect_equal(
    xrf_relative_peak_intensity("Pb", "K", "KL3", excitation_weighting = "jump",
                                coster_kronig = FALSE),
    xrf_absorption_jump_ratio("Pb", "K") * xrf_transition_probability("Pb", "KL3")
  )
})

test_that("photoionization cross section interpolates, is NA below the edge, real to 200 keV", {
  # NA below the Pb K edge (~88 keV), finite above
  expect_true(is.na(xrf_photoionization_cross_section("Pb", "K", 50)))
  expect_true(is.finite(xrf_photoionization_cross_section("Pb", "K", 95)))

  # falls with energy above the edge
  cs <- xrf_photoionization_cross_section("Pb", "K", c(90, 100, 120))
  expect_true(all(diff(cs) < 0))

  # matches a direct log-log interpolation of the (fine EPDL) table
  g <- x_ray_cross_sections[x_ray_cross_sections$element == "Pb" &
                              x_ray_cross_sections$shell == "K", ]
  expect_equal(xrf_photoionization_cross_section("Pb", "K", 95),
               exp(suppressWarnings(approx(log(g$energy_kev), log(g$cross_section), log(95))$y)),
               tolerance = 1e-6)

  # U K edge (115.6 keV) is now WITHIN the EPDL table (to 200 keV) -> a real value, not extrapolated
  expect_true(is.finite(xrf_photoionization_cross_section("U", "K", 130)))
})

test_that("cross-section weighting makes K:L branching depend on beam energy (#8)", {
  kla <- function(E, w) {
    ka <- xrf_relative_peak_intensity("Pb", "K", "KL3", beam_energy_kev = E,
                                      excitation_weighting = w, coster_kronig = FALSE)
    la <- xrf_relative_peak_intensity("Pb", "L3", "L3M5", beam_energy_kev = E,
                                      excitation_weighting = w, coster_kronig = FALSE)
    ka / la
  }
  # jump weighting is frozen; cross-section weighting moves with energy
  expect_equal(kla(90, "jump"), kla(120, "jump"))
  expect_false(isTRUE(all.equal(kla(90, "cross_section"), kla(120, "cross_section"))))
})

test_that("Coster-Kronig redistributes L vacancies toward L3, leaving L1 undepleted (#9)", {
  tr <- c("L1M3", "L2M4", "L3M5")
  no_ck <- xrf_relative_peak_intensity("Au", substr(tr, 1, 2), tr, beam_energy_kev = 60,
                                       excitation_weighting = "cross_section", coster_kronig = FALSE)
  ck <- xrf_relative_peak_intensity("Au", substr(tr, 1, 2), tr, beam_energy_kev = 60,
                                    excitation_weighting = "cross_section", coster_kronig = TRUE)
  # L1-fed line unchanged (yield already accounts for CK loss); L2/L3-fed lines boosted
  expect_equal(ck[1], no_ck[1])
  expect_gt(ck[3], no_ck[3])
  expect_gte(ck[2], no_ck[2])
})

test_that("Gryzinski electron cross section peaks at mid overvoltage; relativistic rise at high E (E1)", {
  Ec <- 7.11   # Fe K edge
  # pure (non-relativistic) Gryzinski g-function peaks around U = 3
  q <- xrf_electron_ionization_cross_section("Fe", "K", Ec * c(0.9, 1.5, 3, 5, 10),
                                             method = "gryzinski", relativistic = FALSE)
  expect_equal(q[1], 0)                       # below threshold
  expect_true(all(q[-1] > 0))
  expect_equal(which.max(q), 3)               # peaks near U = 3
  expect_false(isTRUE(all.equal(
    xrf_electron_ionization_cross_section("Fe", "K", 20, "gryzinski", relativistic = FALSE),
    xrf_electron_ionization_cross_section("Fe", "K", 20, "bethe", relativistic = FALSE))))

  # relativistic correction (gryzinski/bethe) -> ~1 at low beam energy, a rise at high beam energy
  rf <- function(E) xrf_electron_ionization_cross_section("Fe", "K", E, "gryzinski", relativistic = TRUE) /
    xrf_electron_ionization_cross_section("Fe", "K", E, "gryzinski", relativistic = FALSE)
  expect_lt(abs(rf(20) - 1), 0.1)
  expect_gt(rf(150), 1.3)
})

test_that("Bote-Salvat electron cross section is the default and matches the reference (E1/BS)", {
  # default method is bote-salvat, returning an absolute cross section (cm^2)
  q <- xrf_electron_ionization_cross_section("Fe", "K", c(6, 10, 20, 100))   # Fe K edge ~7.08 keV
  expect_equal(q[1], 0)                                    # below threshold
  expect_true(all(q[-1] > 0))
  expect_equal(which.max(q), 3)                            # peaks near 20 keV
  # magnitude is physically sensible (hundreds of barns)
  expect_gt(q[3] / 1e-24, 100); expect_lt(q[3] / 1e-24, 2000)
  # data cover Z 1-99, all subshells
  expect_setequal(sort(unique(x_ray_bote_salvat$subshell)), 1:9)
  expect_equal(max(x_ray_bote_salvat$z), 99)
})

test_that("EPDL cross-section rebuild covers M4/M5 and >100 keV (E4)", {
  # M4/M5 are real now (no M3 proxy)
  expect_true("M5" %in% x_ray_cross_sections$shell)
  expect_true(is.finite(xrf_photoionization_cross_section("Au", "M5", 20)))
  # heavy-element K lines up to ~200 keV are tabulated (U K edge 115.6 keV)
  expect_true(is.finite(xrf_photoionization_cross_section("U", "K", 130)))
  expect_gt(max(x_ray_cross_sections$energy_kev), 150)
  # the new total mass-attenuation table exists with the expected components
  expect_true(all(c("photoelectric", "rayleigh", "compton", "total") %in%
                    colnames(x_ray_mass_attenuation)))
  expect_gt(xrf_mass_attenuation("Si", 10, type = "total"),
            xrf_mass_attenuation("Si", 10, type = "photoelectric") * 0.99)
})

test_that("light elements and excitation mode work in xrf_energies", {
  # C/N/O/F K-lines are now representable (were absent below Ne)
  le <- xrf_energies(c("C", "N", "O", "F"), beam_energy_kev = 5)
  expect_setequal(unique(le$element), c("C", "N", "O", "F"))
  expect_true(all(le$energy_kev < 1))

  # electron mode applies an overvoltage cut: nothing above beam/overvoltage_min survives
  ee <- xrf_energies("everything", 30, excitation = "electron", overvoltage_min = 2)
  expect_true(all(ee$edge_kev <= 30 / 2))
  # a stricter overvoltage keeps the same or fewer lines than a looser one (identical intensities)
  n2 <- nrow(ee)
  n1 <- nrow(xrf_energies("everything", 30, excitation = "electron", overvoltage_min = 1))
  expect_lte(n2, n1)
})

test_that("M lines are available and electron excitation enhances them for heavy elements (#13/#17)", {
  # heavy-element M lines exist (from EADL97) and are prominent under electron (SEM-EDS) excitation;
  # for photon excitation at 20 keV the Au L lines dominate, so the M lines fall below threshold.
  au <- xrf_energies("Au", beam_energy_kev = 20, excitation = "electron")
  expect_true(any(au$edge %in% c("M3", "M4", "M5")))
  # Au Malpha1 ~ 2.12 keV
  ma <- au$energy_kev[au$trans_siegbahn == "Malpha1"]
  expect_equal(ma, 2.12, tolerance = 0.05)

  # U at 15 kV cannot reach the L3 edge (17.2 keV): only M lines survive
  u <- xrf_energies("U", beam_energy_kev = 15)
  expect_true(all(u$edge %in% c("M3", "M4", "M5")))

  # electron excitation raises the M-vs-L family ratio for Au relative to photon (higher M overvoltage)
  mratio <- function(mode) {
    x <- xrf_energies("Au", beam_energy_kev = 20, excitation = mode)
    sum(x$relative_peak_intensity[x$edge %in% c("M3","M4","M5")]) /
      sum(x$relative_peak_intensity[grepl("^L", x$edge)])
  }
  expect_gt(mratio("electron"), mratio("photon"))

  # electron ionization cross section is zero below threshold and positive above
  expect_equal(xrf_electron_ionization_cross_section("Au", "K", 20), 0)   # U < 1 (K edge 80 keV)
  expect_gt(xrf_electron_ionization_cross_section("Au", "M5", 20), 0)
})

test_that("Coster Kronig transitions are correct", {
  expect_identical(
    xrf_coster_kronig_probability("Pb", "L1", "f12"),
    x_ray_coster_kronig_probabilities %>%
      dplyr::filter(element == "Pb", shell == "L1", coster_kronig_trans == "f12") %>%
      dplyr::pull(coster_kronig_prob)
  )

  expect_identical(
    xrf_coster_kronig_probability(c("Pb", "Zn", "U"), "L1", "f12"),
    xrf_coster_kronig_probability(c(82, 30, 92), "L1", "f12")
  )
})

test_that("Coster-Kronig source is switchable to the modern Elam values (item 1)", {
  on.exit(xrf_set_ck_source("EADL97"), add = TRUE)      # always restore the default source
  eadl <- xrf_coster_kronig_probability("Pb", "L1", "f12")
  old <- xrf_set_ck_source("Elam")
  expect_identical(old, "EADL97")                        # setter returns the previous source
  elam <- xrf_coster_kronig_probability("Pb", "L1", "f12")
  expect_equal(elam, 0.040, tolerance = 1e-6)            # the modern Elam value
  expect_false(isTRUE(all.equal(eadl, elam)))            # differs from EADL97
  # the switch propagates through the L cascade to the reported line intensities, then restores cleanly
  ratio <- function() {
    x <- xrf_energies("Pb", 40)
    sum(x$relative_peak_intensity[grepl("^L3", x$edge)]) /
      sum(x$relative_peak_intensity[grepl("^L1|^L2", x$edge)])
  }
  r_elam <- ratio(); xrf_set_ck_source("EADL97"); r_eadl <- ratio()
  expect_false(isTRUE(all.equal(r_eadl, r_elam)))
  expect_equal(xrf_coster_kronig_probability("Pb", "L1", "f12"), eadl)  # back to the default value
})
