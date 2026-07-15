
#' Calculate baselines using SNIP
#'
#' Based on \link[Peaks]{SpectrumBackground} (Peaks package). The most useful parameter to change is the
#' \code{iterations} argument (higher value will result in a more conservative baseline estimation).
#'
#' ... Passed to \link[Peaks]{SpectrumBackground}
#' @inheritParams xrf_add_baseline_pkg
#'
#' @export
#'
xrf_add_baseline_snip <- function(.spectra, .values = .data$.spectra$cps, ...,
                                  energy_aware = FALSE, .energy_kev = .data$.spectra$energy_kev,
                                  .clamp = -Inf, .env = parent.frame()) {
  # have to load the Peaks DLLs for some reason
  .First.lib <- NULL; rm(.First.lib)
  withr::with_namespace("Peaks", .First.lib(dirname(find.package("Peaks")), "Peaks"))

  .values <- enquo(.values)

  if (!energy_aware) {
    return(xrf_add_baseline(
      .spectra,
      Peaks::SpectrumBackground,
      !!.values,
      ...,
      .clamp = .clamp,
      .env = .env
    ))
  }

  # Energy-aware SNIP: the clipping window is constant in CHANNELS, but the peak FWHM grows with
  # energy (FWHM^2 = noise^2 + F*eps*E), so a window right at one energy is wrong elsewhere. We
  # resample onto a sqrt(E) abscissa -- where the peak width is approximately constant per channel --
  # run SNIP there with a fixed window, and map the background back. This keeps the window from
  # climbing into the wings of the broad high-energy K-lines / Compton hump.
  .energy_kev <- enquo(.energy_kev)
  dots <- quos(...)
  .spectra$.spectra <- purrr::map(purrr::transpose(.spectra), function(spectrum) {
    x <- spectrum$.spectra
    vals <- rlang::eval_tidy(.values, data = spectrum, env = .env)
    energy <- rlang::eval_tidy(.energy_kev, data = spectrum, env = .env)
    args <- purrr::map(dots, rlang::eval_tidy, data = spectrum, env = .env)
    bg <- do.call(snip_energy_aware, c(list(values = vals, energy_kev = energy), args))
    x$baseline <- pmax(bg, .clamp)
    x
  })

  .spectra
}

# SNIP on a sqrt(E)-warped, uniformly-resampled axis, mapped back to the original energy grid.
snip_energy_aware <- function(values, energy_kev, ...) {
  stopifnot(length(values) == length(energy_kev))
  finite <- is.finite(energy_kev) & is.finite(values)
  if (sum(finite) < 3) return(rep(NA_real_, length(values)))
  u <- sqrt(pmax(energy_kev, 0))
  u_uniform <- seq(min(u[finite]), max(u[finite]), length.out = length(u))
  # ties at u = 0 (from clamped non-positive energies) are collapsed by approx() -- harmless here
  v_u <- suppressWarnings(stats::approx(u, values, xout = u_uniform, rule = 2))$y
  bg_u <- Peaks::SpectrumBackground(v_u, ...)
  suppressWarnings(stats::approx(u_uniform, bg_u, xout = u, rule = 2))$y
}

#' Caclulate baselines using the baseline package
#'
#' Wraps a call to the \link[baseline]{baseline}(). See documentation in that package for
#' references for baseline correction.
#'
#' @param .values The value column in each spectra object
#' @param ... Passed to \link[baseline]{baseline}. Evaluated rowwise using objects in \code{spectra}.
#' @inheritParams xrf_add_baseline
#'
#' @references
#' See \link[baseline]{baseline}
#'
#' @export
#'
xrf_add_baseline_pkg <- function(.spectra, .values = .data$.spectra$cps, ..., .clamp = -Inf, .env = parent.frame()) {
  .values <- enquo(.values)
  xrf_add_baseline(
    .spectra,
    function(x, ...) {
      obj <- baseline::baseline(matrix(x, nrow = 1), ...)
      baseline::getBaseline(obj)[1, , drop = TRUE]
    },
    !!.values,
    ...,
    .clamp = .clamp,
    .env = .env
  )
}

#' Add a baseline estimate from a function
#'
#' @param .spectra A spectra_df
#' @param .fun A function that receives a vector of values and outputs a vector of values
#'   of the same length
#' @param ... Passed to \code{.fun}, evaluated within .spectra
#' @param .clamp The minimum allowable value for background
#' @param .env The calling environment
#'
#' @return .spectra with a modified .spectra column
#' @export
#'
#' @importFrom rlang enquo quos !! !!!
#'
xrf_add_baseline <- function(.spectra, .fun, ..., .clamp = -Inf, .env = parent.frame()) {
  stopifnot(
    is.numeric(.clamp), length(.clamp) == 1
  )

  dots <- quos(...)
  .spectra$.spectra <- purrr::map(purrr::transpose(.spectra), function(spectrum) {
    x <- spectrum$.spectra

    # args evaluated within each row
    args <- purrr::map(dots, rlang::eval_tidy, data = spectrum, env = .env)

    x$baseline <- do.call(.fun, args)
    x$baseline <- pmax(x$baseline, .clamp)

    x
  })

  .spectra
}
