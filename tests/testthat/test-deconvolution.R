context("test-deconvolution.R")

test_that("test deconvolution", {

  test_df <- tibble::tibble(
    energy_kev = seq(1, 30, 0.1),
    response = xrftools:::gaussian_fun(energy_kev, mu = 5, sigma = 1, height = 4) +
      xrftools:::gaussian_fun(energy_kev, mu = 10, sigma = 2, height = 2) +
      xrftools:::gaussian_fun(energy_kev, mu = 15, sigma = 0.5, height = 1) +
      rnorm(length(energy_kev), mean = 0, sd = 0.05)
  )

  deconv <- xrf_deconvolute_gaussian_least_squares(
    test_df$energy_kev, test_df$response,
    peaks = tibble::tibble(
      element = c("a", "b", "b"),
      energy_kev = c(5, 10, 15),
      sigma = c(1, 2, 0.5),
      relative_peak_intensity = c(sigma * c(4, 2, 1))
    )
  )

  # make sure peaks line up and are the right height
  expect_equal(deconv$peaks$element, c("a", "b"))
  expect_true(all(abs(deconv$peaks$height - c(4, 2)) < 0.05))

  plot(test_df, type = "l", col = "blue")
  lines(deconv$response$energy_kev, deconv$response$response_fit)
  lines(
    deconv$components$energy_kev[deconv$components$element == "a"],
    deconv$components$response_fit[deconv$components$element == "a"],
    col = "purple"
  )
  lines(
    deconv$components$energy_kev[deconv$components$element == "b"],
    deconv$components$response_fit[deconv$components$element == "b"],
    col = "red"
  )
})

test_that("test base deconvolution function", {

  spec <- read_xrf_example(.which = 7) %>%
    xrf_add_baseline_snip(.values = .spectra$cps, iterations = 20) %>%
    xrf_add_smooth_filter(filter = xrf_filter_gaussian(width = 7, alpha = 2.5), .iter = 5)

  spec_simple <- spec %>%
    dplyr::pull(.spectra) %>%
    dplyr::first() %>%
    dplyr::filter(energy_kev <= 7.5)

  deconv <- xrf_deconvolute_gaussian_least_squares(
    spec_simple$energy_kev,
    spec_simple$smooth - spec_simple$baseline,
    peaks = xrf_energies("lake_sediment")
  )

  # types of output
  expect_is(deconv, "deconvolution_fit")
  expect_setequal(names(deconv), c("fit", "response", "components", "peaks", "coef_cov", "area_cov"))
  expect_is(deconv$fit, "data.frame")
  expect_is(deconv$response, "data.frame")
  expect_is(deconv$components, "data.frame")
  expect_is(deconv$peaks, "data.frame")
  expect_true(is.matrix(deconv$coef_cov) && is.matrix(deconv$area_cov))

  # rows of output
  expect_equal(nrow(deconv$peaks), dplyr::n_distinct(xrf_energies("lake_sediment")$element))
  expect_equal(nrow(deconv$fit), 1)
  expect_equal(nrow(deconv$response), nrow(spec_simple))
  expect_equal(nrow(deconv$components), nrow(spec_simple) * nrow(deconv$peaks))
})

test_that("deconvolution works on spectra objects", {
  spec <- read_xrf_example(.which = 7) %>%
    dplyr::select(SampleIdent, ConditionSet, kV) %>%
    xrf_add_baseline_snip(.values = .spectra$cps, iterations = 20) %>%
    xrf_add_smooth_filter(filter = xrf_filter_gaussian(width = 7, alpha = 2.5), .iter = 5)

  spec_simple <- spec %>%
    dplyr::pull(.spectra) %>%
    dplyr::first()

  deconv <- xrf_deconvolute_gaussian_least_squares(
    spec_simple$energy_kev,
    spec_simple$cps - spec_simple$baseline,
    peaks = xrf_energies("lake_sediment"),
    energy_max_kev = 7.5
  )

  spec_deconv <- spec %>%
    xrf_add_deconvolution_gls(energy_max_kev = 7.5, peaks = xrf_energies("lake_sediment"))

  expect_identical(spec_deconv$.deconvolution_components[[1]], deconv$components)
})

test_that("nonneg (NNLS) clamps absent elements to zero instead of going negative", {
  e <- seq(1, 20, 0.05)
  # noise-free two-peak spectrum so the clamp is deterministic; "ghost" has no signal at 14 keV
  y <- xrftools:::gaussian_fun(e, 5, 0.1, 4000) + xrftools:::gaussian_fun(e, 10, 0.1, 2000)
  peaks <- tibble::tibble(
    element = c("a", "b", "ghost"),
    energy_kev = c(5, 10, 14),
    relative_peak_intensity = 1
  )

  d_nn <- xrf_deconvolute_gaussian_least_squares(e, y, peaks = peaks, default_sigma = 0.1,
                                                 nonneg = TRUE)
  gh <- d_nn$peaks$element == "ghost"
  expect_true(all(d_nn$peaks$height >= 0))
  expect_equivalent(d_nn$peaks$height[gh], 0)             # clamped to exactly zero
  expect_true(is.na(d_nn$peaks$height_se[gh]))            # clamped -> NA (not a confident) SE
  # real peaks are still recovered
  expect_true(all(abs(d_nn$peaks$height[!gh] - c(4000, 2000)) < 1))

  # and with negatively-pushing noise, unconstrained LS would go negative but NNLS stays >= 0
  set.seed(42)
  y2 <- y + rnorm(length(e), 0, 5)
  d_free <- xrf_deconvolute_gaussian_least_squares(e, y2, peaks = peaks, default_sigma = 0.1,
                                                   nonneg = FALSE)
  d_nn2 <- xrf_deconvolute_gaussian_least_squares(e, y2, peaks = peaks, default_sigma = 0.1,
                                                  nonneg = TRUE)
  expect_true(all(d_nn2$peaks$height >= 0))
  expect_true(min(d_free$peaks$height) <= min(d_nn2$peaks$height))
})

test_that("Poisson weighting runs and chi_sq is finite even with zero-count channels", {
  set.seed(7)
  e <- seq(1, 20, 0.05)
  lt <- 15
  truth <- xrftools:::gaussian_fun(e, 5, 0.12, 300) + xrftools:::gaussian_fun(e, 12, 0.12, 120)
  counts <- rpois(length(e), truth * lt)   # many channels are exactly zero
  cps <- counts / lt
  peaks <- tibble::tibble(element = c("a", "b"), energy_kev = c(5, 12), relative_peak_intensity = 1)

  # historic chi_sq = sum(resid^2 / response) would be NaN/Inf here (division by zero-count cps)
  d <- xrf_deconvolute_gaussian_least_squares(e, cps, peaks = peaks, default_sigma = 0.12,
                                              weighting = "poisson", counts = counts, livetime = lt)
  expect_identical(d$fit$weighting, "poisson")
  expect_true(is.finite(d$fit$chi_sq))
  expect_true(all(d$peaks$height >= 0))
  expect_true(all(is.finite(d$peaks$height_se)))

  # asking for poisson without counts warns and falls back
  expect_warning(
    xrf_deconvolute_gaussian_least_squares(e, cps, peaks = peaks, weighting = "poisson"),
    "requires"
  )
})

test_that("detector model gives energy-dependent peak widths", {
  peaks <- tibble::tibble(element = c("lo", "hi"), energy_kev = c(2, 40),
                          relative_peak_intensity = 1)
  e <- seq(0.5, 45, 0.05)
  y <- rep(0, length(e))
  d <- xrf_deconvolute_gaussian_least_squares(e, y, peaks = peaks, detector_type = "SDD")
  # the 40 keV line must be assigned a wider sigma than the 2 keV line
  s_lo <- d$peaks$primary_sigma[d$peaks$element == "lo"]
  s_hi <- d$peaks$primary_sigma[d$peaks$element == "hi"]
  expect_gt(s_hi, s_lo)
  expect_equal(s_lo, xrf_detector_sigma_kev(2, "SDD"))
  expect_equal(s_hi, xrf_detector_sigma_kev(40, "SDD"))
})

test_that("deconvolution function errors are caught", {
  spec <- read_xrf_example(.which = 7)
  expect_error(xrf_add_deconvolution_fun(spec, function(...) "not a list obj"), "not a list")
  expect_error(xrf_add_deconvolution_fun(spec, function(...) list(1)), "empty or missing names")
  expect_error(xrf_add_deconvolution_fun(spec, function(...) list(1, a = "fish")), "empty or missing names")
  expect_silent(xrf_add_deconvolution_fun(spec, function(...) list(a = 1, b = 2)))
})

test_that("active_thickness_um is threaded through the xrf_add_deconvolution_gls WRAPPER", {
  # CloudCal drives the tidy wrapper (not the worker directly), so a dropped `!!active_thickness_um`
  # splice in xrf_add_deconvolution_gls would pass every worker-level test yet silently break production.
  # A wide grid spanning Ba L (~4.5) and K (~32 keV) lets the active-layer thickness reshape the fit.
  e <- seq(2, 40, 0.02)
  sig <- function(E) xrf_detector_sigma_kev(E, "SDD")
  y <- xrftools:::gaussian_fun(e, 4.47, sig(4.47), 500) +
    xrftools:::gaussian_fun(e, 32.19, sig(32.19), 500)
  spec <- tibble::tibble(spectrum_id = "s1",
                         .spectra = list(data.frame(energy_kev = e, smooth = y, baseline = 0)))
  wrap_h <- function(t) {
    o <- spec %>% xrf_add_deconvolution_gls(.spectra$energy_kev, .spectra$smooth - .spectra$baseline,
           peaks = xrf_energies("Ba", 60), detector_type = "SDD", efficiency = TRUE,
           active_thickness_um = t, cache_templates = FALSE)
    p <- o$.deconvolution_peaks[[1]]
    p$height[p$element == "Ba"]
  }
  h_thin <- wrap_h(50); h_thick <- wrap_h(3000)
  expect_gt(h_thin, 0)
  expect_gt(h_thick, 0)
  expect_gt(abs(h_thin - h_thick), 0.1 * h_thick)   # the wrapper genuinely forwards the argument
  expect_equal(wrap_h(NULL), wrap_h(450))           # NULL -> SDD preset, which is now the 450 um SDD
})

test_that("abundance prior accepts a user matrix override (F3)", {
  els <- c("Fe", "Cr", "Ni", "Pb")
  base <- xrftools:::.xrf_abundance_penalty_factor(els, ref_ppm = 100)
  # a steel matrix: Fe/Cr/Ni are majors -> their penalty factors drop toward 0 (not "rare"); Pb, unnamed,
  # keeps its crustal default.
  steel <- xrftools:::.xrf_abundance_penalty_factor(
    els, ref_ppm = 100, abundance = c(Fe = 7e5, Cr = 1.8e5, Ni = 8e4))
  expect_lt(steel[2], base[2])          # Cr no longer penalised as a trace
  expect_lt(steel[3], base[3])          # Ni likewise
  expect_equal(steel[4], base[4])       # Pb (unnamed) keeps its crustal default
  # end-to-end: the plumbed argument is accepted by the deconvolution without error
  e <- seq(2, 15, 0.05); sig <- xrf_detector_sigma_kev(6.4, "SDD")
  y <- xrftools:::gaussian_fun(e, 6.40, sig, 500) + xrftools:::gaussian_fun(e, 5.41, sig, 200) + 5
  pk <- xrf_energies(c("Fe", "Cr", "Mn"), beam_energy_kev = 20)
  d <- xrf_deconvolute_gaussian_least_squares(
    e, y, peaks = pk, detector_type = "SDD",
    abundance_prior = 0.2, abundance_ppm = c(Fe = 7e5, Cr = 1.8e5))
  expect_true(is.data.frame(d$peaks))
})

test_that("clean-line exemption spares elements with an isolated line (F2)", {
  energy <- seq(1, 15, 0.02)
  g <- function(mu, s, h) h * exp(-0.5 * ((energy - mu) / s) ^ 2)
  # As: Ka 10.53 (overlaps Pb La 10.55) + a clean Kb 11.73. Pb: La 10.55 + clean Lb 12.61.
  # Tb: fully isolated pair (control). Ov: a single line buried under the As/Pb overlap (no clean line).
  X <- cbind(As = g(10.53, 0.1, 1) + g(11.73, 0.1, 0.15),
             Pb = g(10.55, 0.1, 1) + g(12.61, 0.1, 0.6),
             Tb = g(6.27, 0.1, 1) + g(6.98, 0.1, 0.5),
             Ov = g(10.54, 0.1, 0.8))
  els <- colnames(X)
  pk <- rbind(
    data.frame(element = "As", energy_kev = c(10.53, 11.73), relative_peak_intensity = c(1, 0.15)),
    data.frame(element = "Pb", energy_kev = c(10.55, 12.61), relative_peak_intensity = c(1, 0.6)),
    data.frame(element = "Tb", energy_kev = c(6.27, 6.98),   relative_peak_intensity = c(1, 0.5)),
    data.frame(element = "Ov", energy_kev = 10.54,           relative_peak_intensity = 1))
  ex <- xrftools:::.xrf_clean_line_exempt(pk, X, els, rep(TRUE, 4), energy)
  names(ex) <- els
  expect_true(ex[["As"]])    # identifiable via its clean Kbeta at 11.73
  expect_true(ex[["Tb"]])    # fully isolated
  expect_true(ex[["Pb"]])    # La overlaps As Ka, but Lb (12.61) is clean
  expect_false(ex[["Ov"]])   # only an overlapped line -> correctly still penalisable
})

test_that("lasso / elastic-net prior gives exact sparsity where least squares does not (item 2)", {
  set.seed(42)
  n <- 400; g <- function(mu) exp(-0.5 * ((seq_len(n) - mu) / 8) ^ 2)
  X <- cbind(A = g(150), B = g(260), Ph = g(158))    # Ph overlaps A; truly absent from y
  rw <- rep(1, n)
  y <- 5 * g(150) + 3 * g(260) + rnorm(n, 0, 0.05)
  scale <- mean(colSums((X * rw) ^ 2)); pf <- c(0, 0, 50)   # penalise only the confounded phantom
  # unconstrained: least squares leaks noise into Ph; ridge shrinks it; lasso (IRL2) zeros it exactly
  ols   <- qr.coef(qr(X * rw), y * rw)
  ridge <- qr.coef(qr(rbind(X * rw, diag(sqrt(0.3 * scale * pf), 3, 3))), c(y * rw, rep(0, 3)))
  lasso <- xrftools:::.xrf_irl2_fit(X, y, rw, scale, pf, alpha = 1,   strength = 0.3, nonneg = FALSE, p = 3)
  enet  <- xrftools:::.xrf_irl2_fit(X, y, rw, scale, pf, alpha = 0.5, strength = 0.3, nonneg = FALSE, p = 3)
  expect_gt(abs(ols[[3]]), 1e-4)               # OLS: Ph non-zero (fits noise)
  expect_lt(abs(ridge[[3]]), abs(ols[[3]]))    # ridge: shrunk
  expect_lt(abs(lasso$coef[3]), 1e-9)          # lasso: snapped to exactly zero
  expect_lt(abs(enet$coef[3]), 1e-9)           # elastic-net likewise
  expect_gt(lasso$coef[1], 0)                  # the real column is retained
})

test_that("lasso prior plugs into the deconvolution; ridge stays the default (item 2)", {
  e <- seq(2, 15, 0.05); sig <- xrf_detector_sigma_kev(6.4, "SDD")
  y <- xrftools:::gaussian_fun(e, 6.40, sig, 500) + xrftools:::gaussian_fun(e, 5.41, sig, 200) + 5
  pk <- xrf_energies(c("Fe", "Cr", "Mn", "Ni", "Co"), beam_energy_kev = 20)
  d_l  <- xrf_deconvolute_gaussian_least_squares(e, y, peaks = pk, detector_type = "SDD",
                                                 abundance_prior = 0.3, prior = "lasso")
  d_en <- xrf_deconvolute_gaussian_least_squares(e, y, peaks = pk, detector_type = "SDD",
                                                 abundance_prior = 0.3, prior = "elastic_net")
  expect_true(is.data.frame(d_l$peaks) && is.data.frame(d_en$peaks))
  expect_gt(max(d_l$peaks$height[d_l$peaks$element == "Fe"]), 0)     # present element retained
  # default prior = ridge is byte-identical to omitting the argument
  d_def <- xrf_deconvolute_gaussian_least_squares(e, y, peaks = pk, detector_type = "SDD", abundance_prior = 0.3)
  d_rid <- xrf_deconvolute_gaussian_least_squares(e, y, peaks = pk, detector_type = "SDD",
                                                  abundance_prior = 0.3, prior = "ridge")
  expect_equal(d_def$peaks$peak_area, d_rid$peaks$peak_area)
})

test_that("deconvolution exposes a Laplace posterior covariance (item 4)", {
  e <- seq(2, 16, 0.03); sig <- function(E) xrf_detector_sigma_kev(E, "SDD")
  y <- xrftools:::gaussian_fun(e, 6.40, sig(6.40), 400) +
       xrftools:::gaussian_fun(e, 3.69, sig(3.69), 250) + 10
  f <- xrf_deconvolute_gaussian_least_squares(e, y, peaks = xrf_energies(c("Fe", "Ca"), 20),
                                              detector_type = "SDD")
  expect_true(is.matrix(f$coef_cov) && is.matrix(f$area_cov))
  expect_equal(dim(f$area_cov), c(nrow(f$peaks), nrow(f$peaks)))
  expect_equal(f$area_cov, t(f$area_cov))                         # symmetric
  # internal consistency: the diagonal of area_cov IS peak_area_se^2
  d <- diag(f$area_cov)[f$peaks$element]; se2 <- f$peaks$peak_area_se^2
  ok <- is.finite(se2)
  expect_equal(as.numeric(d[ok]), as.numeric(se2[ok]), tolerance = 1e-8)
})

test_that("the abundance prior shrinks the Laplace SE of a penalised element (item 4)", {
  sig <- function(E) xrf_detector_sigma_kev(E, "SDD")
  e <- seq(8, 16, 0.02)
  rnd <- function(l, a) { s <- sig(l$energy_kev); y <- numeric(length(e))
    for (i in seq_len(nrow(l))) y <- y + a * l$relative_peak_intensity[i] *
      exp(-0.5 * ((e - l$energy_kev[i]) / s[i]) ^ 2); y }
  y <- rnd(xrf_energies("Pb", 40), 300) + rnd(xrf_energies("As", 40), 40)
  pk <- xrf_energies(c("As", "Pb", "Ga", "Se"), 40)
  # window out the clean As Kbeta (11.7): As Kalpha then overlaps Pb Lalpha with no clean line, so it IS
  # penalised and the prior curvature must reduce its posterior SE.
  se_as <- function(ap) {
    f <- xrf_deconvolute_gaussian_least_squares(e, y, peaks = pk, detector_type = "SDD",
           energy_max_kev = 11.2, abundance_prior = ap, abundance_protect = "Pb")
    max(f$peaks$peak_area_se[f$peaks$element == "As"], na.rm = TRUE)
  }
  expect_lt(se_as(0.5), se_as(0))
})

test_that("Compton normalization uses the fit covariance for its ratio SE (item 4)", {
  sig <- function(E) xrf_detector_sigma_kev(E, "SDD")
  e <- seq(2, 25, 0.03)
  y <- xrftools:::gaussian_fun(e, 6.4, sig(6.4), 400) + xrftools:::gaussian_fun(e, 20.2, sig(20.2), 150) + 30
  f <- xrf_deconvolute_gaussian_least_squares(e, y, peaks = xrf_energies(c("Fe", "Ca"), 40),
         detector_type = "SDD", tube = xrf_tube("Rh", kv = 40),
         geometry = xrf_geometry(scatter_angle_deg = 135), scatter_continuum = TRUE, background = 8)
  n_cov <- xrf_compton_normalize(f)          # deconvolution_fit -> uses area_cov (correlated)
  n_ind <- xrf_compton_normalize(f$peaks)    # plain tibble -> independence fallback
  fe_cov <- n_cov$peak_area_norm_se[n_cov$element == "Fe"]
  fe_ind <- n_ind$peak_area_norm_se[n_ind$element == "Fe"]
  expect_true(is.finite(fe_cov) && fe_cov > 0)
  expect_false(isTRUE(all.equal(fe_cov, fe_ind)))
})
