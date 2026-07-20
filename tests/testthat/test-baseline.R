context("test-baseline.R")

test_that("SNIP baseline works as expected", {
  spectra <- read_xrf_example(.which = 1)
  spectra_baseline <- xrf_add_baseline_snip(spectra)
  expect_true("baseline" %in% names(spectra_baseline$.spectra[[1]]))
  expect_identical(
    spectra_baseline$.spectra[[1]]$baseline,
    xrf_snip_background(spectra$.spectra[[1]]$cps)
  )
})

test_that("baseline package functions work as expected", {
  spectra <- read_xrf_example(.which = 1)
  spectra_baseline <- xrf_add_baseline_pkg(spectra, method = "als", p = 0.003, .clamp = -Inf)
  expect_true("baseline" %in% names(spectra_baseline$.spectra[[1]]))
  expect_identical(
    spectra_baseline$.spectra[[1]]$baseline,
    baseline::getBaseline(
      baseline::baseline(
        matrix(spectra$.spectra[[1]]$cps, nrow = 1),
        method = "als", p = 0.003
      )
    )[1, , drop = TRUE]
  )
})

test_that("tidy evaluation works in baseline package functions", {
  spectra <- tidyr::crossing(
    tibble::tibble(p_val = c(0.003, 0.01, 0.05)),
    read_xrf_example(.which = 1)
  )

  spectra_baseline <- xrf_add_baseline_pkg(spectra, method = "als", p = p_val, .clamp = -Inf)

  expect_identical(
    spectra_baseline$.spectra[[1]]$baseline,
    baseline::getBaseline(
      baseline::baseline(
        matrix(spectra$.spectra[[1]]$cps, nrow = 1),
        method = "als", p = 0.003
      )
    )[1, , drop = TRUE]
  )

  expect_identical(
    spectra_baseline$.spectra[[2]]$baseline,
    baseline::getBaseline(
      baseline::baseline(
        matrix(spectra$.spectra[[2]]$cps, nrow = 1),
        method = "als", p = 0.01
      )
    )[1, , drop = TRUE]
  )
})

test_that("arPLS recovers a smooth background without clipping peaks", {
  e <- seq(1, 20, 0.02)
  # a realistic, gently-varying continuum (slow decay + a broad scatter hump), not a steep exponential
  bg_true <- 40 * exp(-e / 12) + 25 * exp(-0.5 * ((e - 13) / 4) ^ 2) + 8
  peaks <- 300 * exp(-0.5 * ((e - 6.4) / 0.1) ^ 2) +
           200 * exp(-0.5 * ((e - 8.0) / 0.1) ^ 2)        # two narrow peaks
  y <- bg_true + peaks
  z <- xrftools:::.xrf_arpls(y, lambda = 1e5)
  # baseline tracks the continuum in peak-free interior regions (excluding Whittaker boundary bias)
  pf <- e > 3 & e < 17 & abs(e - 6.4) > 0.8 & abs(e - 8.0) > 0.8
  expect_lt(sqrt(mean((z[pf] - bg_true[pf]) ^ 2)) / mean(bg_true[pf]), 0.15)
  # peaks survive in the subtracted spectrum (not absorbed into the baseline)
  net <- y - z
  expect_gt(max(net[abs(e - 6.4) < 0.3]), 0.7 * 300)
  expect_gt(sum(pmax(net, 0)), 0.8 * sum(peaks))
  # a stiffer lambda gives a smoother (smaller-curvature) baseline
  z_soft  <- xrftools:::.xrf_arpls(y, lambda = 1e3)
  z_stiff <- xrftools:::.xrf_arpls(y, lambda = 1e7)
  expect_lt(sum(diff(z_stiff, differences = 2) ^ 2), sum(diff(z_soft, differences = 2) ^ 2))
  # non-finite input is tolerated
  y2 <- y; y2[c(5, 500)] <- NA
  expect_true(all(is.finite(xrftools:::.xrf_arpls(y2))))
})

test_that("xrf_add_baseline_arpls plugs into the spectra pipeline", {
  spectra <- read_xrf_example(.which = 1)
  out <- xrf_add_baseline_arpls(spectra, lambda = 1e5)
  bl <- out$.spectra[[1]]$baseline
  cps <- spectra$.spectra[[1]]$cps
  expect_true("baseline" %in% names(out$.spectra[[1]]))
  expect_equal(length(bl), length(cps))
  expect_true(all(is.finite(bl)))
  expect_true(mean(bl) < mean(cps))               # a baseline sits below the data on average
  # tidy-eval of a per-row parameter works (mirrors the SNIP/pkg wrappers)
  multi <- tidyr::crossing(tibble::tibble(lam = c(1e4, 1e6)), read_xrf_example(.which = 1))
  expect_silent(xrf_add_baseline_arpls(multi, lambda = lam))
})
