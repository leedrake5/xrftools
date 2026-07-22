#' SNIP background (native, order-2 clipping)
#'
#' Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP) baseline. A native, dependency-free
#' replacement for the (CRAN-archived) \code{Peaks::SpectrumBackground} default (order 2, growing window, no
#' log-log-sqrt transform): each channel is iteratively clipped to the mean of its two neighbours at a
#' half-window \eqn{p = 1 \ldots \code{iterations}}, leaving the smooth continuum under the peaks. Reproduces
#' the reference implementation's output exactly (max abs difference 0 on real spectra).
#'
#' @param y Numeric intensity vector.
#' @param iterations Largest clipping half-window; higher = more conservative (flatter) baseline. Default 100.
#' @param decreasing If \code{TRUE} the half-window shrinks (\code{iterations..1}) instead of growing. Default FALSE.
#' @param ... Ignored; tolerates the legacy \code{Peaks} arguments (\code{order}/\code{window}/\code{smoothing}/
#'   \code{compton}), which the deconvolution never set away from their order-2, no-smoothing defaults.
#' @return The estimated baseline, same length as \code{y}.
#' @export
xrf_snip_background <- function(y, iterations = 100, decreasing = FALSE, ...) {
  z <- as.numeric(y)
  n <- length(z)
  if (n < 3L) return(z)
  # interpolate any non-finite so the clip does not propagate NA (the reference required finite input)
  if (anyNA(z) || any(!is.finite(z))) {
    bad <- !is.finite(z); ok <- which(!bad)
    if (length(ok) >= 2L) z[bad] <- suppressWarnings(stats::approx(ok, z[ok], xout = which(bad), rule = 2)$y)
    z[!is.finite(z)] <- 0
  }
  iterations <- max(1L, as.integer(iterations))
  ps <- if (isTRUE(decreasing)) seq.int(iterations, 1L) else seq.int(1L, iterations)
  for (p in ps) {
    if (p < 1L || p >= n) next
    idx <- (p + 1L):(n - p)
    # vectorised assignment evaluates the RHS against the pre-step z (double-buffered), matching the reference
    z[idx] <- pmin(z[idx], 0.5 * (z[idx - p] + z[idx + p]))
  }
  z
}

#' Calculate baselines using SNIP
#'
#' A native SNIP baseline (\link{xrf_snip_background}). The most useful parameter to change is the
#' \code{iterations} argument (higher value will result in a more conservative baseline estimation).
#'
#' With \code{energy_aware = TRUE} the SNIP window is warped onto a constant-resolution abscissa built from
#' the detector FWHM model; \code{detector_type} / \code{fano} / \code{epsilon_ev} / \code{noise_fwhm_ev}
#' select that model and \strong{default to an SDD (silicon)} when not supplied -- pass the detector you are
#' actually using (matching the deconvolution) so the warp uses the right resolution curve.
#'
#' ... Passed to \link{xrf_snip_background} (\code{iterations}, \code{decreasing}).
#' @inheritParams xrf_add_baseline_pkg
#' @param energy_aware Warp onto the constant-resolution abscissa before SNIP (see Details).
#' @param .energy_kev Channel energies (keV); used only when \code{energy_aware = TRUE} to build the warp.
#' @param detector_type,fano,epsilon_ev,noise_fwhm_ev Detector resolution model for the energy-aware warp
#'   (see \link{xrf_detector_fwhm_ev}); default an SDD.
#'
#' @export
#'
xrf_add_baseline_snip <- function(.spectra, .values = .data$.spectra$cps, ...,
                                  energy_aware = FALSE, .energy_kev = .data$.spectra$energy_kev,
                                  detector_type = NULL, fano = NULL, epsilon_ev = NULL,
                                  noise_fwhm_ev = NULL,
                                  .clamp = -Inf, .env = parent.frame()) {
  .values <- enquo(.values)

  if (!energy_aware) {
    return(xrf_add_baseline(
      .spectra,
      xrf_snip_background,
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
    bg <- do.call(snip_energy_aware, c(list(values = vals, energy_kev = energy,
                                            detector_type = detector_type, fano = fano,
                                            epsilon_ev = epsilon_ev, noise_fwhm_ev = noise_fwhm_ev), args))
    x$baseline <- pmax(bg, .clamp)
    x
  })

  .spectra
}

# SNIP on a resolution-warped, uniformly-resampled axis, mapped back to the original energy grid.
# The clipping window is constant in CHANNELS, but the peak FWHM varies with energy, so a fixed window is
# wrong on the raw axis. The abscissa that makes the peak width constant per channel is
#   u(E) = integral dE / FWHM(E),   FWHM(E) = sqrt(noise^2 + (2.3548)^2 F eps E)  (the detector model),
# on which one FWHM spans a constant step in u. This reduces to ~sqrt(E) only where the Fano term dominates;
# a bare sqrt(E) warp over-compresses the LOW-energy axis (where the electronic-noise floor makes FWHM ~
# constant), which lets SNIP strip part of the smooth light-element continuum. Using the real resolution
# curve stays near-linear there and removes that light-element baseline bias.
snip_energy_aware <- function(values, energy_kev, ..., detector_type = NULL, fano = NULL,
                              epsilon_ev = NULL, noise_fwhm_ev = NULL) {
  stopifnot(length(values) == length(energy_kev))
  finite <- is.finite(energy_kev) & is.finite(values)
  if (sum(finite) < 3) return(rep(NA_real_, length(values)))
  fwhm <- xrf_detector_fwhm_ev(energy_kev, detector_type = detector_type, fano = fano,
                               epsilon_ev = epsilon_ev, noise_fwhm_ev = noise_fwhm_ev) / 1000  # keV
  good <- is.finite(fwhm) & fwhm > 0
  if (!any(good)) return(rep(NA_real_, length(values)))
  fwhm[!good] <- min(fwhm[good])
  # constant-resolution abscissa u(E) = cumulative integral of dE / FWHM(E) (energies sorted ascending,
  # then mapped back to the original channel order).
  o <- order(energy_kev)
  e_s <- energy_kev[o]; f_s <- fwhm[o]
  du <- c(0, diff(e_s) / (0.5 * (f_s[-length(f_s)] + f_s[-1])))
  du[!is.finite(du) | du < 0] <- 0
  u <- numeric(length(energy_kev)); u[o] <- cumsum(du)
  u_uniform <- seq(min(u[finite]), max(u[finite]), length.out = length(u))
  # ties (from any zero-width intervals) are collapsed by approx() -- harmless here
  v_u <- suppressWarnings(stats::approx(u, values, xout = u_uniform, rule = 2))$y
  bg_u <- xrf_snip_background(v_u, ...)
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

#' Calculate baselines using asymmetrically reweighted penalized least squares (arPLS)
#'
#' A modern Whittaker-smoother baseline (Baek et al. 2015): it fits a smooth curve penalised by its second
#' differences (smoothness \code{lambda}) and iteratively down-weights channels that sit \strong{above} the
#' current estimate (peaks) while keeping those at or below it (background), via an asymmetric logistic
#' reweighting driven by the noise level of the below-baseline residuals. Relative to SNIP it is tunable
#' through a single smoothness parameter, follows a smooth continuum (e.g. a scatter / Compton hump) more
#' gently, and does not depend on a channel-window width. Larger \code{lambda} gives a stiffer (smoother)
#' baseline.
#'
#' @inheritParams xrf_add_baseline_pkg
#' @param ... arPLS solver controls, evaluated rowwise and passed to the fit: \code{lambda} (smoothness
#'   penalty on the baseline's second differences; larger = stiffer / smoother, default 1e5), \code{ratio}
#'   (convergence tolerance on the relative change of the weight vector between iterations, default 0.01), and
#'   \code{max_iter} (maximum reweighting iterations, default 50).
#'
#' @references
#' Baek, S.-J., Park, A., Ahn, Y.-J., Choo, J. (2015). Baseline correction using asymmetrically reweighted
#' penalized least squares smoothing. \emph{Analyst} 140(1), 250-257. \doi{10.1039/C4AN01061B}
#'
#' @export
#'
xrf_add_baseline_arpls <- function(.spectra, .values = .data$.spectra$cps, ...,
                                   .clamp = -Inf, .env = parent.frame()) {
  # lambda / ratio / max_iter flow through `...` (evaluated rowwise, like xrf_add_baseline_pkg's method/p),
  # so a per-row parameter column works; the inner function supplies the defaults when they are omitted.
  .values <- enquo(.values)
  xrf_add_baseline(
    .spectra,
    function(x, lambda = 1e5, ratio = 0.01, max_iter = 50L) {
      .xrf_arpls(as.numeric(x), lambda = lambda, ratio = ratio, max_iter = max_iter)
    },
    !!.values,
    ...,
    .clamp = .clamp,
    .env = .env
  )
}

# arPLS core (Baek et al. 2015): solve (W + lambda D2'D2) z = W y for the baseline z, then reweight
#   w_i = 1 / (1 + exp(2 (d_i - (2 s - m)) / s)),   d = y - z,   (m, s) = mean/SD of the BELOW-baseline
# residuals d < 0. Peaks (d >> 0) get w -> 0 and are ignored; background (d <= 0) gets w -> 1. D2 is the
# second-difference operator, so the penalised system is sparse pentadiagonal and is solved with a sparse
# Cholesky (Matrix). Non-finite inputs are linearly interpolated before the solve.
.xrf_arpls <- function(y, lambda = 1e5, ratio = 0.01, max_iter = 50L) {
  y <- as.numeric(y); n <- length(y)
  if (n < 3) return(rep(NA_real_, n))
  fin <- is.finite(y)
  if (!all(fin)) {
    if (sum(fin) < 3) return(rep(NA_real_, n))
    y <- stats::approx(which(fin), y[fin], xout = seq_len(n), rule = 2)$y
  }
  lambda <- max(as.numeric(lambda), 0)
  k <- seq_len(n - 2)
  D <- Matrix::sparseMatrix(i = rep(k, 3), j = c(k, k + 1, k + 2),
                            x = rep(c(1, -2, 1), each = n - 2), dims = c(n - 2, n))
  DtD <- Matrix::crossprod(D)
  w <- rep(1, n); z <- y
  for (i in seq_len(max(1L, as.integer(max_iter)))) {
    z <- as.numeric(Matrix::solve(Matrix::Diagonal(x = w) + lambda * DtD, w * y))
    d <- y - z
    neg <- d[d < 0]
    if (length(neg) < 2) break
    m <- mean(neg); s <- stats::sd(neg)
    if (!is.finite(s) || s <= 0) break
    wt <- 1 / (1 + exp(2 * (d - (2 * s - m)) / s))
    denom <- sqrt(sum(w * w))
    if (is.finite(denom) && denom > 0 && sqrt(sum((wt - w)^2)) / denom < ratio) break
    w <- wt
  }
  z
}

#' Extract the modelled baseline from a deconvolution fit
#'
#' The deconvolution's non-element templates -- the fitted \code{background} basis, the
#' Rayleigh/Compton scatter lines, the scattered-continuum templates (including the second-order
#' Compton continuum) and any free-amplitude sum peaks -- together ARE a physics-informed model of the
#' spectrum's baseline. This extracts that modelled baseline on the fit's energy grid, so it can be
#' used the way a SNIP/arPLS estimate would be: net-peak extraction, or a normalization denominator
#' (e.g. a backscatter reference channel for a Lucas-Tooth calibration) that is free of analyte peaks
#' by construction. Unlike a clipping estimator it cannot eat slow peak shoulders, and unlike SNIP it
#' reproduces sharp physical structure (the Compton edge, the scatter peaks themselves if included).
#'
#' Escape satellites cannot be separated out (they are summed into their parent element's template),
#' and the fixed-amplitude two-pass pile-up model (\code{pileup_tau}) contributes to
#' \code{response_fit} but not to \code{components}, so it is not included; the free-amplitude
#' \code{sum_peaks} templates are.
#'
#' @param fit A \code{deconvolution_fit} (from \link{xrf_deconvolute_gaussian_least_squares}), fitted
#'   against an \strong{un-baselined} response with a fitted \code{background} (and usually
#'   \code{scatter_continuum = TRUE}); see the "Working recipe" in \link{xrf_add_deconvolution_gls}.
#' @param include Which template classes count as baseline: any of \code{"background"} (the smooth
#'   fitted basis), \code{"scatter"} (anode-line Rayleigh/Compton + scattered continuum), and
#'   \code{"pileup"} (free-amplitude \code{sum_} coincidence templates). Default: all three. For a
#'   Compton-normalization denominator you may want \code{include = "scatter"} only.
#'
#' @return A tibble with \code{energy_kev} and \code{baseline} (same units as the fitted response),
#'   on the fit's (possibly windowed) energy grid.
#' @export
#'
xrf_fit_baseline <- function(fit, include = c("background", "scatter", "pileup")) {
  stopifnot(inherits(fit, "deconvolution_fit"))
  pats <- c(background = "^background_", scatter = "^scatter_", pileup = "^(pileup_|sum_)")
  include <- match.arg(include, names(pats), several.ok = TRUE)
  pat <- paste(pats[include], collapse = "|")
  comp <- fit$components
  energy <- fit$response$energy_kev
  sel <- grepl(pat, comp$element)
  baseline <- if (any(sel)) {
    rowSums(matrix(comp$response_fit[sel], nrow = length(energy)))
  } else {
    rep(0, length(energy))
  }
  tibble::tibble(energy_kev = energy, baseline = baseline)
}

#' Model the baseline with the physics-informed deconvolution
#'
#' A drop-in alternative to \link{xrf_add_baseline_snip} / \link{xrf_add_baseline_arpls} that models
#' the baseline \emph{physically} instead of clipping it: the full deconvolution is run against the
#' \strong{raw} spectrum (\code{.values = cps}) with a fitted smooth \code{background} basis and --
#' when a \code{tube} is supplied -- the Rayleigh/Compton scatter templates and scattered-continuum
#' shapes (first- and second-order Compton). The fitted non-element part (\link{xrf_fit_baseline})
#' is then written into each spectrum's \code{baseline} column, so downstream
#' \code{cps - baseline} workflows (net peaks, Lucas-Tooth normalization points, backscatter
#' references) work unchanged -- but the baseline now contains the Compton edge, scatter humps and
#' peak shelves that a clipping estimator smears, and is guaranteed peak-free because the element
#' peaks were modelled out, not clipped around.
#'
#' Convenience defaults injected when not supplied: \code{.values = .spectra$cps} (the raw spectrum
#' -- required for baseline modelling), \code{background = 8}, and \code{scatter_continuum = TRUE}
#' whenever a \code{tube} is given. All other arguments pass through to
#' \link{xrf_add_deconvolution_gls} (consider \code{background_smooth} for a stiffer background and
#' \code{tailing = TRUE} to include peak shelves in the model). The \code{.deconvolution_*} columns
#' are kept on the result (the peak areas came for free). Channels outside a restricted fit window
#' (\code{energy_min_kev}/\code{energy_max_kev}) get \code{NA} baseline.
#'
#' @inheritParams xrf_add_baseline
#' @param ... Passed to \link{xrf_add_deconvolution_gls} (e.g. \code{peaks}, \code{detector_type},
#'   \code{tube}, \code{background}, \code{background_smooth}, \code{tailing}).
#' @param baseline_include Template classes counted as baseline (see \link{xrf_fit_baseline}).
#'
#' @return \code{.spectra} with a \code{baseline} column added to each spectrum (plus the
#'   \code{.deconvolution_*} columns from the underlying fit).
#' @export
#'
xrf_add_baseline_deconvolution <- function(.spectra, ..., .clamp = -Inf,
                                           baseline_include = c("background", "scatter", "pileup")) {
  dots <- rlang::enquos(...)
  if (!(".values" %in% names(dots))) dots$.values <- rlang::quo(.data$.spectra$cps)
  if (!("background" %in% names(dots))) dots$background <- 8L
  if (("tube" %in% names(dots)) && !("scatter_continuum" %in% names(dots))) {
    dots$scatter_continuum <- TRUE
  }
  fitted <- rlang::inject(xrf_add_deconvolution_gls(.spectra, !!!dots, .env = parent.frame()))

  fitted$.spectra <- purrr::map2(fitted$.spectra, fitted$.deconvolution_components, function(x, comp) {
    energy_fit <- comp$energy_kev[comp$element == comp$element[1]]
    pats <- c(background = "^background_", scatter = "^scatter_", pileup = "^(pileup_|sum_)")
    pat <- paste(pats[match.arg(baseline_include, names(pats), several.ok = TRUE)], collapse = "|")
    sel <- grepl(pat, comp$element)
    bl <- if (any(sel)) rowSums(matrix(comp$response_fit[sel], nrow = length(energy_fit)))
          else rep(0, length(energy_fit))
    x$baseline <- pmax(bl[match(x$energy_kev, energy_fit)], .clamp)
    x
  })
  fitted
}
