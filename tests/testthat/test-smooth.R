context("test-smooth.R")

test_that("smoothing works", {
  spec <- read_xrf_example(.which = 3)
  smoothed <- xrf_add_smooth_filter(spec, filter = xrf_filter_gaussian(width = 3, alpha = 2.5), .iter = 1)
  expect_true("smooth" %in% names(smoothed$.spectra[[1]]))
  expect_equal(
    smoothed$.spectra[[1]]$smooth[c(-1, -nrow(smoothed$.spectra[[1]]))],
    stats::filter(spec$.spectra[[1]]$cps, filter = xrf_filter_gaussian(width = 3, alpha = 2.5)) %>%
      as.numeric() %>%
      .[c(-1, -nrow(smoothed$.spectra[[1]]))]
  )
})

test_that("filter generating functions work", {
  expect_equal(sum(xrf_filter_gaussian()), 1)
  expect_equal(sum(xrf_filter_pyramid()), 1)
})

test_that("xrf_filter_gaussian coerces a non-odd-integer width to the nearest odd integer, with a warning", {
  # a valid odd width passes through unchanged and silently
  expect_silent(f7 <- xrf_filter_gaussian(width = 7))
  expect_equal(length(f7), 7)

  # a non-integer width (e.g. from a continuous optimizer) is snapped to the nearest odd integer and warns
  expect_warning(f <- xrf_filter_gaussian(width = 20.327), "odd integer")
  expect_equal(length(f), 21)        # 20.327 -> 21
  expect_equal(sum(f), 1)            # still a normalized kernel

  # an even integer is coerced too
  expect_warning(f20 <- xrf_filter_gaussian(width = 20), "coercing")
  expect_equal(length(f20) %% 2, 1)  # odd number of taps

  # widths below the minimum are floored to the smallest valid odd width (3)
  expect_warning(f1 <- xrf_filter_gaussian(width = 1))
  expect_equal(length(f1), 3)
})
