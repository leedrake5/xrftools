
#' Deconvolute a spectrum using gaussian least squares
#'
#' @inheritParams xrf_add_deconvolution_fun
#' @param .values,response Values to use for the deconvolution
#' @param energy_kev,.energy_kev Energy values corresponding to response/.values
#' @param energy_min_kev,energy_max_kev Constrain the energies that should be deconvoluted
#' @param peaks A data frame with (at least) columns "element" and "energy_kev". Optional columns
#'   "sigma" and "relative_peak_intensity" are used when present.
#' @param default_sigma The default standard deviation of the peaks (keV), used when neither a
#'   \code{sigma} column nor a detector model is supplied.
#' @param detector_type,fano,epsilon_ev,noise_fwhm_ev Detector resolution parameters. When any of
#'   these is supplied (and \code{peaks} has no \code{sigma} column) the per-line peak width is
#'   computed from the energy-dependent resolution model \link{xrf_detector_sigma_kev} instead of a
#'   constant \code{default_sigma}. \code{detector_type} is one of \link{xrf_detector_presets}
#'   (e.g. "SDD", "HPGe", "CdTe"); the scalars override the preset.
#' @param efficiency If TRUE, reweight each line by the energy-dependent detector efficiency
#'   (\link{xrf_detector_efficiency}); corrects line shape and high-energy sensitivity. Uses the
#'   \code{detector_type} geometry (override window/dead layer with \code{be_window_um} /
#'   \code{dead_layer_um}).
#' @param escape If TRUE, add detector escape-peak satellites (\link{xrf_escape_peaks}), tied to each
#'   parent line's amplitude. Important for Ge/CdTe at high energy.
#' @param be_window_um,dead_layer_um Detector window / dead-layer thickness (microns) for the
#'   efficiency model; \code{NULL} uses the \code{detector_type} preset.
#' @param active_thickness_um Active-volume thickness (microns) for the efficiency model; \code{NULL}
#'   uses the \code{detector_type} preset (SDD = 450 um, a typical modern silicon-drift detector).
#'   Governs high-energy transmission (thin Si becomes transparent to heavy-element K lines), so it
#'   only matters when \code{efficiency = TRUE}.
#' @param nonneg If TRUE (default) coefficients are constrained to be non-negative via
#'   \link[nnls]{nnls} (peak areas / concentrations cannot be negative). Set FALSE for the classic
#'   unconstrained least-squares behaviour.
#' @param weighting Counting-statistics weighting. "poisson" performs an iteratively-reweighted fit
#'   with weights 1/variance derived from the modelled counts (requires \code{counts}/\code{.counts}
#'   and \code{livetime}/\code{.livetime}); "none" is ordinary (unweighted) least squares; "auto"
#'   (default) uses "poisson" when counts and live time are available and "none" otherwise.
#' @param counts,.counts,livetime,.livetime Raw gross counts per channel and the acquisition live
#'   time (seconds), used for Poisson weighting and for the count-space goodness-of-fit. Optional.
#' @param tube,geometry Optional \link{xrf_tube} / \link{xrf_geometry} objects. When a \code{tube}
#'   is supplied, Rayleigh and Compton tube-scatter templates are added to the fit automatically
#'   (see \code{scatter} and \link{xrf_scatter_peaks}); \code{geometry} sets the Compton scatter
#'   angle. (Not otherwise consumed yet -- reserved for fundamental-parameters work.)
#' @param scatter Whether to include Rayleigh/Compton scatter templates. \code{NULL} (default)
#'   enables them when a \code{tube} is supplied; \code{TRUE}/\code{FALSE} force on/off.
#' @param compton_broadening Multiplier on the detector width for the Doppler-broadened Compton
#'   peak (passed to \link{xrf_scatter_peaks}).
#' @param tail,step,beta Hypermet line-shape parameters (passed to \link{xrf_lineshape}): relative
#'   low-energy tail amplitude, shelf/step amplitude, and tail decay length. Defaults (0, 0, NULL)
#'   give a pure Gaussian. The reported \code{peak_area} includes the tail area when \code{tail > 0}.
#' @param cache_templates If TRUE (default) reuse the design matrix across successive calls that
#'   share the same energy grid, peaks and line-shape (e.g. a batch from one instrument), skipping
#'   the Gaussian re-evaluation.
#' @param refine_calibration If TRUE, refine a 2-parameter energy calibration (zero + gain) of the
#'   template centroids by variable projection before the final fit -- robust to small gain/zero
#'   drift that would otherwise bias amplitudes. Disables template caching for that call.
#' @param sum_peaks If TRUE, add coincidence (pile-up) sum-peak templates at E_i + E_j for the
#'   strongest lines, so pile-up is not misassigned to real elements.
#' @param count_rate Optional total count rate (cps); when supplied and 0, the sum-peak step is
#'   skipped. Reserved for scaling the pile-up contribution.
#' @param pileup_tau Optional detector pulse-pair resolution time (seconds). When supplied with
#'   \code{counts}/\code{.counts} and \code{livetime}/\code{.livetime}, a physically-constrained
#'   two-pass pile-up correction is applied (coincidence peak areas fixed to \eqn{2\tau/T\,A_iA_j});
#'   this supersedes the free-amplitude \code{sum_peaks} model.
#' @param use_qr If TRUE (default) unconstrained fits (\code{nonneg = FALSE}) use fast QR
#'   decomposition; set FALSE for fork-safe \code{lm()} under multicore parallelism. Ignored when
#'   \code{nonneg = TRUE}.
#'
#' @return A modified version of .spectra
#' @export
#'
xrf_add_deconvolution_gls <- function(.spectra, .energy_kev = .data$.spectra$energy_kev,
                                      .values = .data$.spectra$cps - .data$.spectra$baseline,
                                      energy_min_kev = -Inf, energy_max_kev = Inf,
                                      peaks = xrf_energies(), default_sigma = 0.07,
                                      detector_type = NULL, fano = NULL, epsilon_ev = NULL,
                                      noise_fwhm_ev = NULL,
                                      efficiency = FALSE, escape = FALSE,
                                      be_window_um = NULL, dead_layer_um = NULL,
                                      active_thickness_um = NULL,
                                      air_path_cm = NULL, atmosphere = "Air", window = NULL,
                                      nonneg = TRUE, weighting = c("auto", "poisson", "none"),
                                      .counts = NULL, .livetime = NULL,
                                      tube = NULL, geometry = NULL,
                                      scatter = NULL, compton_broadening = 2,
                                      tail = 0, step = 0, beta = NULL,
                                      cache_templates = TRUE, refine_calibration = FALSE,
                                      sum_peaks = FALSE, count_rate = NULL, pileup_tau = NULL,
                                      use_qr = TRUE, abundance_prior = 0, abundance_ref_ppm = 100,
                                      .env = parent.frame()) {
  .energy_kev <- enquo(.energy_kev)
  .values <- enquo(.values)
  peaks <- enquo(peaks)
  default_sigma <- enquo(default_sigma)
  energy_min_kev <- enquo(energy_min_kev)
  energy_max_kev <- enquo(energy_max_kev)
  detector_type <- enquo(detector_type)
  fano <- enquo(fano)
  epsilon_ev <- enquo(epsilon_ev)
  noise_fwhm_ev <- enquo(noise_fwhm_ev)
  efficiency <- enquo(efficiency)
  escape <- enquo(escape)
  be_window_um <- enquo(be_window_um)
  dead_layer_um <- enquo(dead_layer_um)
  active_thickness_um <- enquo(active_thickness_um)
  air_path_cm <- enquo(air_path_cm)
  atmosphere <- enquo(atmosphere)
  window <- enquo(window)
  nonneg <- enquo(nonneg)
  weighting <- enquo(weighting)
  .counts <- enquo(.counts)
  .livetime <- enquo(.livetime)
  tube <- enquo(tube)
  geometry <- enquo(geometry)
  scatter <- enquo(scatter)
  compton_broadening <- enquo(compton_broadening)
  tail <- enquo(tail)
  step <- enquo(step)
  beta <- enquo(beta)
  cache_templates <- enquo(cache_templates)
  refine_calibration <- enquo(refine_calibration)
  sum_peaks <- enquo(sum_peaks)
  count_rate <- enquo(count_rate)
  pileup_tau <- enquo(pileup_tau)
  use_qr <- enquo(use_qr)
  abundance_prior <- enquo(abundance_prior)
  abundance_ref_ppm <- enquo(abundance_ref_ppm)

  xrf_add_deconvolution_fun(
    .spectra,
    xrf_deconvolute_gaussian_least_squares,
    !!.energy_kev,
    !!.values,
    peaks = !!peaks,
    default_sigma = !!default_sigma,
    energy_min_kev = !!energy_min_kev,
    energy_max_kev = !!energy_max_kev,
    detector_type = !!detector_type,
    fano = !!fano,
    epsilon_ev = !!epsilon_ev,
    noise_fwhm_ev = !!noise_fwhm_ev,
    efficiency = !!efficiency,
    escape = !!escape,
    be_window_um = !!be_window_um,
    dead_layer_um = !!dead_layer_um,
    active_thickness_um = !!active_thickness_um,
    air_path_cm = !!air_path_cm, atmosphere = !!atmosphere, window = !!window,
    nonneg = !!nonneg,
    weighting = !!weighting,
    counts = !!.counts,
    livetime = !!.livetime,
    tube = !!tube,
    geometry = !!geometry,
    scatter = !!scatter,
    compton_broadening = !!compton_broadening,
    tail = !!tail,
    step = !!step,
    beta = !!beta,
    cache_templates = !!cache_templates,
    refine_calibration = !!refine_calibration,
    sum_peaks = !!sum_peaks,
    count_rate = !!count_rate,
    pileup_tau = !!pileup_tau,
    use_qr = !!use_qr,
    abundance_prior = !!abundance_prior,
    abundance_ref_ppm = !!abundance_ref_ppm,
    .env = .env
  )
}

#' @rdname xrf_add_deconvolution_gls
#' @export
xrf_deconvolute_gaussian_least_squares <- function(energy_kev, response, peaks = xrf_energies(),
                                                   default_sigma = 0.07,
                                                   energy_min_kev = -Inf, energy_max_kev = Inf,
                                                   detector_type = NULL, fano = NULL,
                                                   epsilon_ev = NULL, noise_fwhm_ev = NULL,
                                                   efficiency = FALSE, escape = FALSE,
                                                   be_window_um = NULL, dead_layer_um = NULL,
                                                   active_thickness_um = NULL,
                                                   air_path_cm = NULL, atmosphere = "Air", window = NULL,
                                                   nonneg = TRUE,
                                                   weighting = c("auto", "poisson", "none"),
                                                   counts = NULL, livetime = NULL,
                                                   tube = NULL, geometry = NULL,
                                                   scatter = NULL, compton_broadening = 2,
                                                   tail = 0, step = 0, beta = NULL,
                                                   cache_templates = TRUE, refine_calibration = FALSE,
                                                   sum_peaks = FALSE, count_rate = NULL,
                                                   pileup_tau = NULL,
                                                   use_qr = TRUE,
                                                   abundance_prior = 0, abundance_ref_ppm = 100) {

  # ---- true two-pass pile-up correction (#E3) ----------------------------------------------------
  # When a pulse-pair resolution time `pileup_tau` (s) is supplied with counts + livetime, model the
  # coincidence peaks with PHYSICALLY-FIXED amplitudes (proportional to the product of the parent
  # line areas), rather than the free-amplitude single-pass `sum_peaks`: fit the elements (pass 1),
  # subtract the predicted pile-up spectrum, then refit the elements (pass 2).
  if (!is.null(pileup_tau) && is.finite(pileup_tau) && pileup_tau > 0 &&
      !is.null(counts) && !is.null(livetime) && is.finite(livetime) && livetime > 0) {
    args <- mget(names(formals()))
    pass1 <- do.call(xrf_deconvolute_gaussian_least_squares,
                     utils::modifyList(args, list(sum_peaks = FALSE, pileup_tau = NULL)))
    pu <- .xrf_pileup_model(pass1$peaks, energy_kev, pileup_tau, livetime,
                            energy_min_kev, energy_max_kev)
    pass2 <- do.call(xrf_deconvolute_gaussian_least_squares,
                     utils::modifyList(args, list(response = response - pu$model,
                                                  sum_peaks = FALSE, pileup_tau = NULL)))
    if (nrow(pu$peaks) > 0) {
      pass2$peaks <- dplyr::bind_rows(pass2$peaks, pu$peaks)
      pass2$response$response_fit <- pass2$response$response_fit + pu$model
    }
    return(pass2)
  }

  weighting <- match.arg(weighting)

  stopifnot(
    "energy_kev" %in% colnames(peaks),
    "element" %in% colnames(peaks),
    is.numeric(energy_kev), is.numeric(response),
    all(is.finite(energy_kev)), all(is.finite(response)),
    length(energy_kev) == length(response)
  )

  if(!("relative_peak_intensity" %in% colnames(peaks))) {
    peaks$relative_peak_intensity <- 1
  } else {
    stopifnot(is.numeric(peaks$relative_peak_intensity))
  }

  # Optional detector efficiency: reweight each line by the fraction actually recorded in the
  # photopeak (window/dead-layer transmission x active-volume absorption). Corrects intra-element
  # line shape and, for quantification, absolute sensitivity; largest at high energy where thin Si
  # detectors are transparent to heavy-element K-lines.
  if (isTRUE(efficiency)) {
    peaks$relative_peak_intensity <- peaks$relative_peak_intensity *
      xrf_detector_efficiency(peaks$energy_kev, detector_type = detector_type,
                              active_thickness_um = active_thickness_um,
                              be_window_um = be_window_um, dead_layer_um = dead_layer_um,
                              air_path_cm = air_path_cm, atmosphere = atmosphere, window = window)
  }

  # Optional escape peaks: satellites at (parent energy - detector K X-ray), tied to the parent's
  # amplitude. Built from the (efficiency-corrected) parent intensities, so append after efficiency.
  if (isTRUE(escape)) {
    peaks <- dplyr::bind_rows(peaks, xrf_escape_peaks(peaks, detector_type = detector_type))
  }

  # Optional Rayleigh/Compton tube-scatter templates: appended as two extra "elements" so they are
  # fit (non-negatively) jointly with the element lines instead of being left in the baseline.
  # Enabled automatically when a tube is supplied; force with scatter = TRUE / FALSE.
  include_scatter <- if (is.null(scatter)) {
    !is.null(tube) && inherits(tube, "xrf_tube")
  } else {
    isTRUE(scatter)
  }
  if (include_scatter) {
    if (is.null(tube) || !inherits(tube, "xrf_tube")) {
      stop("scatter = TRUE requires tube = xrf_tube(...).")
    }
    scatter_peaks <- xrf_scatter_peaks(
      tube, if (is.null(geometry)) xrf_geometry() else geometry,
      detector_type = detector_type, fano = fano, epsilon_ev = epsilon_ev,
      noise_fwhm_ev = noise_fwhm_ev, default_sigma = default_sigma,
      compton_broadening = compton_broadening
    )
    peaks <- dplyr::bind_rows(peaks, scatter_peaks)
  }

  # Per-line peak width: an explicit `sigma` column wins; otherwise, if any detector-resolution
  # parameter is given, compute the energy-dependent sigma(E); otherwise fall back to the constant
  # default_sigma. A single constant width is only correct near ~8-10 keV, so a detector model is
  # strongly preferred for wide-range work.
  detector_given <- !is.null(detector_type) || !is.null(fano) ||
    !is.null(epsilon_ev) || !is.null(noise_fwhm_ev)
  model_sigma <- function(e) {
    xrf_detector_sigma_kev(e, detector_type = detector_type, fano = fano,
                           epsilon_ev = epsilon_ev, noise_fwhm_ev = noise_fwhm_ev)
  }
  if(!("sigma" %in% colnames(peaks))) {
    peaks$sigma <- if (detector_given) model_sigma(peaks$energy_kev) else default_sigma
  } else {
    stopifnot(is.numeric(peaks$sigma))
    fallback <- if (detector_given) model_sigma(peaks$energy_kev) else default_sigma
    peaks$sigma <- dplyr::coalesce(peaks$sigma, fallback)
  }

  # Optional coincidence (pile-up) sum peaks at E_i + E_j (#21): at high count rate two photons are
  # recorded as one event at the summed energy. Giving the strongest lines their own sum templates
  # (fit freely) keeps pile-up from being misassigned to real trace elements. Simplified single-pass
  # model; the pile-up rate scales with `count_rate` (used only to skip the step when clearly low).
  if (isTRUE(sum_peaks) && (is.null(count_rate) || count_rate > 0)) {
    lc <- peaks[is.finite(peaks$energy_kev) & is.finite(peaks$sigma) &
                  peaks$relative_peak_intensity > 0, , drop = FALSE]
    lc <- lc[order(-lc$relative_peak_intensity), , drop = FALSE]
    max_lines <- min(nrow(lc), 5L)
    if (max_lines >= 1) {
      lc <- lc[seq_len(max_lines), ]
      pairs <- utils::combn(max_lines, 2, simplify = FALSE)
      pairs <- c(pairs, lapply(seq_len(max_lines), function(i) c(i, i)))   # + 2x self-sums
      sum_rows <- purrr::map_dfr(pairs, function(ij) {
        tibble::tibble(
          element = paste0("sum_", round(lc$energy_kev[ij[1]] + lc$energy_kev[ij[2]], 2)),
          energy_kev = lc$energy_kev[ij[1]] + lc$energy_kev[ij[2]],
          sigma = sqrt(lc$sigma[ij[1]]^2 + lc$sigma[ij[2]]^2),
          relative_peak_intensity = 1
        )
      })
      sum_rows <- sum_rows[!duplicated(round(sum_rows$energy_kev, 3)) &
                             sum_rows$energy_kev >= energy_min_kev &
                             sum_rows$energy_kev <= energy_max_kev, , drop = FALSE]
      if (nrow(sum_rows) > 0) peaks <- dplyr::bind_rows(peaks, sum_rows)
    }
  }

  # filter energies/responses (and counts, if supplied) to the requested window
  within_range <- (energy_kev >= energy_min_kev) & (energy_kev <= energy_max_kev)
  energy_kev <- energy_kev[within_range]
  response <- response[within_range]
  if(!is.null(counts)) {
    stopifnot(length(counts) == length(within_range))
    counts <- counts[within_range]
  }

  # check finite-ness of inputs
  stopifnot(
    all(is.finite(peaks$sigma)),
    all(is.finite(peaks$relative_peak_intensity))
  )

  # filter, normalize peak heights (scale to sigma)
  # the peak with a relative height of 1 should be the peak with the greatest intensity
  peaks <- peaks %>%
    dplyr::filter(
      .data$energy_kev >= !!energy_min_kev,
      .data$energy_kev <= !!energy_max_kev,
      .data$relative_peak_intensity > 0
    ) %>%
    dplyr::group_by(.data$element) %>%
    dplyr::mutate(
      relative_peak_intensity = .data$relative_peak_intensity / max(.data$relative_peak_intensity),
      relative_peak_height = .data$relative_peak_intensity / .data$sigma
    ) %>%
    dplyr::mutate(
      relative_peak_height = .data$relative_peak_height /
        .data$relative_peak_height[which.max(.data$relative_peak_intensity)]
    ) %>%
    dplyr::ungroup()

  # Build the per-element template design matrix. Reused across a batch (#15): if the energy grid,
  # the (normalized) peaks and the line-shape are identical to the previous call, the cached matrix
  # is returned instead of re-evaluating every Gaussian.
  peaks_tmpl <- peaks %>%
    dplyr::select("element", "energy_kev", "sigma", "relative_peak_height", "relative_peak_intensity")

  # per-line response: a windowed pure Gaussian (evaluate exp() only within +/-6 sigma, #19) unless a
  # Hypermet tail/shelf is requested, in which case the full line shape is evaluated (#20).
  build_line <- function(mu, sigma, height) {
    if (tail == 0 && step == 0) {
      out <- numeric(length(energy_kev))
      idx <- which(energy_kev >= mu - 6 * sigma & energy_kev <= mu + 6 * sigma)
      if (length(idx)) out[idx] <- height * exp(-0.5 * ((energy_kev[idx] - mu) / sigma) ^ 2)
      out
    } else {
      xrf_lineshape(energy_kev, mu = mu, sigma = sigma, height = height,
                    tail = tail, step = step, beta = beta)
    }
  }
  # assemble the per-element template matrix + metadata from a peaks table
  build_X <- function(pk) {
    built <- pk %>%
      dplyr::mutate(
        response_element = purrr::pmap(
          list(.data$energy_kev, .data$sigma, .data$relative_peak_height), build_line
        )
      ) %>%
      dplyr::group_by(.data$element) %>%
      dplyr::summarise(
        primary_energy_kev = .data$energy_kev[which.max(.data$relative_peak_intensity)],
        primary_sigma = .data$sigma[which.max(.data$relative_peak_intensity)],
        response_element = list(purrr::reduce(.data$response_element, `+`))
      ) %>%
      dplyr::ungroup()
    Xb <- do.call(cbind, built$response_element)
    colnames(Xb) <- built$element
    list(X = Xb, meta = built[, c("element", "primary_energy_kev", "primary_sigma")])
  }

  # ---- optional energy-calibration refinement (#16) ----------------------------------------------
  # A small zero/gain miscalibration shifts every centroid, and because templates are non-linear in
  # the centroid this injects large derivative-shaped residuals. Refine the 2-parameter map
  # E' = zero + gain * E of the template centroids by variable projection: for each (zero, gain) the
  # amplitudes are solved (fast NNLS) and the residual sum of squares minimized over (zero, gain).
  if (isTRUE(refine_calibration)) {
    rss_for <- function(par) {
      pk <- peaks_tmpl
      pk$energy_kev <- par[1] + par[2] * pk$energy_kev
      Xc <- build_X(pk)$X
      cf <- nnls::nnls(Xc, response)$x
      sum((response - as.vector(Xc %*% cf)) ^ 2)
    }
    opt <- tryCatch(
      stats::optim(c(0, 1), rss_for, method = "Nelder-Mead",
                   control = list(reltol = 1e-9, maxit = 200)),
      error = function(e) NULL
    )
    if (!is.null(opt) && opt$convergence == 0) {
      peaks_tmpl$energy_kev <- opt$par[1] + opt$par[2] * peaks_tmpl$energy_kev
    }
  }

  # ---- design matrix (cached across a batch, #15) ------------------------------------------------
  tmpl_key <- list(energy_kev = energy_kev, peaks = peaks_tmpl, tail = tail, step = step, beta = beta)
  if (cache_templates && !isTRUE(refine_calibration) && !is.null(.xrf_template_cache$key) &&
      identical(.xrf_template_cache$key, tmpl_key)) {
    responses <- .xrf_template_cache$responses
    X <- .xrf_template_cache$X
  } else {
    bx <- build_X(peaks_tmpl)
    X <- bx$X
    responses <- bx$meta
    if (cache_templates && !isTRUE(refine_calibration)) {
      .xrf_template_cache$key <- tmpl_key
      .xrf_template_cache$responses <- responses
      .xrf_template_cache$X <- X
    }
  }
  elements <- responses$element
  n <- length(response)
  p <- ncol(X)

  # ---- abundance-informed ridge prior (opt-in; #abundance) ---------------------------------------
  # Per-column penalty factor ~ 1/crustal_abundance (0 for scatter_*/sum_* and off-table columns). The
  # absolute penalty applied in solve_fit is `abundance_prior * data_scale * pen_factor`, so a degenerate
  # overlap splits toward the more abundant element; strength 0 disables it (identical legacy fit).
  use_prior <- is.finite(abundance_prior) && abundance_prior > 0
  # pen_factor = (1/abundance) x collinearity, so only rare AND overlap-confounded columns are penalised.
  pen_factor <- if (use_prior)
      .xrf_abundance_penalty_factor(elements, abundance_ref_ppm) * .xrf_template_collinearity(X)
    else rep(0, p)

  # ---- counting-statistics weighting -------------------------------------------------------------
  # var(response_i [cps]) = gross_counts_i / livetime^2 (Poisson), so w_i = livetime^2 / gross_counts.
  # `response` is net (baseline-subtracted) cps; gross counts are reconstructed from `counts`.
  have_counts <- !is.null(counts) && !is.null(livetime) &&
    length(livetime) == 1 && is.finite(livetime) && livetime > 0 &&
    length(counts) == n
  if (weighting == "auto") weighting <- if (have_counts) "poisson" else "none"
  if (weighting == "poisson" && !have_counts) {
    warning("weighting = 'poisson' requires `counts` and `livetime`; falling back to 'none'.")
    weighting <- "none"
  }
  baseline_cps <- if (have_counts) (counts / livetime) - response else NULL

  # ---- solver ------------------------------------------------------------------------------------
  # Returns the coefficient vector and the "active" set (the estimable / non-clamped columns).
  # For NNLS the active set is the passive set reported by the solver (its exactly-non-zero
  # columns), which is more reliable than a numerical `coef > 0` test on floating-point dust.
  solve_fit <- function(w) {
    rw <- sqrt(w)
    Xw <- X * rw          # row-scale (recycles rw down each column)
    yw <- response * rw
    if (use_prior && any(pen_factor > 0)) {
      # Tikhonov augmentation: append p pseudo-observation rows sqrt(lambda_j) with target 0, so the
      # solve minimises ||W^1/2(y-Xb)||^2 + sum lambda_j b_j^2. lambda scales to the data (mean column
      # self-energy) so `abundance_prior` is a dimensionless strength; penalty rows carry weight 1 (a
      # prior is not a Poisson observation), so this is correct on every IRLS reweighting pass.
      scale <- mean(colSums(Xw^2))
      if (is.finite(scale) && scale > 0) {
        lambda <- abundance_prior * scale * pen_factor
        Xw <- rbind(Xw, diag(sqrt(lambda), p, p))
        yw <- c(yw, rep(0, p))
      }
    }
    if (nonneg) {
      nn <- nnls::nnls(Xw, yw)
      active <- if (nn$nsetp > 0) sort(nn$passive[seq_len(nn$nsetp)]) else integer(0)
      cf <- nn$x
      # Treat coefficients negligible relative to the largest as clamped (drop NNLS numerical
      # dust): they are set exactly to zero and excluded from the active set, so they report an
      # NA (rather than a spurious) standard error.
      if (length(active) > 0) {
        tol <- max(cf[active]) * 1e-8
        active <- active[cf[active] > tol]
      }
      cf[setdiff(seq_len(p), active)] <- 0
      list(coef = cf, active = active, aliased = integer(0))
    } else if (use_qr) {
      cf <- qr.coef(qr(Xw), yw)
      aliased <- which(is.na(cf))     # rank-deficient columns pivoted out by the QR
      cf[aliased] <- 0
      list(coef = cf, active = setdiff(seq_len(p), aliased), aliased = aliased)
    } else if (use_prior && any(pen_factor > 0)) {
      cf <- qr.coef(qr(Xw), yw)      # lm.wfit can't take the penalty rows; ridge via augmented QR
      aliased <- which(is.na(cf))
      cf[aliased] <- 0
      list(coef = cf, active = setdiff(seq_len(p), aliased), aliased = aliased)
    } else {
      cf <- unname(stats::lm.wfit(X, response, w = w)$coefficients)
      aliased <- which(is.na(cf))
      cf[aliased] <- 0
      list(coef = cf, active = setdiff(seq_len(p), aliased), aliased = aliased)
    }
  }

  if (weighting == "poisson") {
    # iteratively reweighted: model-based Poisson weights (avoids the zero-count bias of data weights)
    w <- rep(1, n)
    for (iter in 1:3) {
      coef <- solve_fit(w)$coef
      fit_iter <- as.vector(X %*% coef)
      model_gross_counts <- pmax((fit_iter + baseline_cps) * livetime, 1)
      w <- livetime^2 / model_gross_counts
    }
    final <- solve_fit(w)
  } else {
    w <- rep(1, n)
    final <- solve_fit(w)
  }
  coef <- final$coef
  names(coef) <- elements
  active <- final$active
  aliased <- final$aliased
  response_fit <- as.vector(X %*% coef)
  residuals <- response - response_fit

  # ---- rank / collinearity diagnostics (#22) -----------------------------------------------------
  rw <- sqrt(w)
  Xw <- X * rw
  # condition number of the (weighted) design: large values (>~1e3) flag near-collinear templates
  # whose amplitudes are poorly separated (e.g. overlapping lines of different elements).
  condition_number <- tryCatch(kappa(Xw, exact = FALSE), error = function(e) NA_real_)
  if (length(aliased) > 0) {
    warning("Deconvolution design matrix is rank-deficient; the following element templates are ",
            "aliased (not separable) and are reported as NA rather than 0: ",
            paste(elements[aliased], collapse = ", "), call. = FALSE)
  }

  # ---- coefficient standard errors ---------------------------------------------------------------
  # From the covariance of the estimable (active) coefficient set; clamped/aliased columns get NA.
  rank <- if (nonneg) length(active) else qr(Xw)$rank
  sigma2 <- sum(w * residuals^2) / max(n - length(active), 1)
  coef_se <- rep(NA_real_, p)
  names(coef_se) <- elements
  if (length(active) > 0) {
    Xa <- Xw[, active, drop = FALSE]
    qra <- qr(Xa)
    ra <- qra$rank
    if (ra == ncol(Xa)) {
      Ra <- qr.R(qra)
      coef_se[active] <- sqrt(rowSums(backsolve(Ra, diag(ncol(Xa)))^2) * sigma2)
    } else {
      # rank-deficient active set: report SEs only for the estimable pivots, NA for the aliased
      Ra <- qr.R(qra)[1:ra, 1:ra, drop = FALSE]
      piv <- qra$pivot[1:ra]
      coef_se[active[piv]] <- sqrt(rowSums(backsolve(Ra, diag(ra))^2) * sigma2)
    }
  }

  # ---- goodness of fit ---------------------------------------------------------------------------
  ss_res <- sum(residuals^2)
  ss_tot <- sum((response - mean(response))^2)
  r2 <- 1 - ss_res / ss_tot
  r2_adj <- 1 - (1 - r2) * (n - 1) / (n - rank)

  # Reduced chi-square. In count space (Pearson) when counts are available -- never divide by the
  # measured value (the historic sum(resid^2 / response) is NaN/Inf on the many zero-count channels
  # of a real spectrum). Falls back to a guarded net-space diagnostic otherwise.
  if (have_counts) {
    mu <- pmax((response_fit + baseline_cps) * livetime, .Machine$double.eps)
    chi_sq <- sum((counts - mu)^2 / pmax(mu, 1)) / max(n - rank, 1)
  } else {
    valid <- response > 0
    chi_sq <- if (any(valid)) {
      sum(residuals[valid]^2 / response[valid]) / max(sum(valid) - rank, 1)
    } else {
      NA_real_
    }
  }

  df <- tibble::tibble(
    energy_kev = energy_kev,
    response = response,
    response_fit = response_fit
  )

  n_energy <- length(energy_kev)
  n_elements <- length(elements)
  components <- tibble::tibble(
    energy_kev = rep(energy_kev, n_elements),
    element = rep(elements, each = n_energy),
    response_fit = as.vector(X) * rep(coef, each = n_energy),
    height = rep(coef, each = n_energy)
  )

  # Build peaks output. peak_area is the primary line's Gaussian area, plus the Hypermet tail area
  # (integral of the exponential tail ~ tail * beta) when a tail is modelled.
  responses$height <- coef[elements]
  responses$height_se <- coef_se[elements]
  if (length(aliased) > 0) {           # aliased (rank-deficient) templates: report NA, not 0 (#22)
    responses$height[aliased] <- NA_real_
    responses$height_se[aliased] <- NA_real_
  }
  tail_beta <- if (is.null(beta) || !is.finite(beta) || beta <= 0) responses$primary_sigma else beta
  area_per_height <- responses$primary_sigma * sqrt(2 * pi) + if (tail != 0) tail * tail_beta else 0
  responses$peak_area <- responses$height * area_per_height
  responses$peak_area_se <- responses$height_se * area_per_height

  structure(
    list(
      fit = tibble::tibble(
        r2 = r2,
        r2_adj = r2_adj,
        chi_sq = chi_sq,
        rank = rank,
        condition_number = condition_number,
        nonneg = nonneg,
        weighting = weighting
      ),
      response = df,
      components = components,
      peaks = responses
    ),
    class = "deconvolution_fit"
  )
}

#' Add a deconvolution based on a function
#'
#' @param .spectra A \link{xrf_spectra} tibble.
#' @param .fun A function that
#' @param ... Passed to .fun
#' @param .env Calling environment
#'
#' @return .spectra with a .deconvolution column, and a deconvolution column added to
#'   the .spectra column
#' @export
#'
xrf_add_deconvolution_fun <- function(.spectra, .fun, ..., .env = parent.frame()) {
  dots <- quos(...)
  fun_wrap <- function(spectrum) {
    # args evaluated within each row
    args <- purrr::map(dots, rlang::eval_tidy, data = spectrum, env = .env)
    do.call(.fun, args)
  }

  # calculate deconvolutions
  deconvs <- purrr::map(purrr::transpose(.spectra), fun_wrap)

  # place each element of the output in a column in the data frame

  # check for named lists
  deconv_not_list <- !purrr::map_lgl(deconvs, is.list)
  if(any(deconv_not_list)) {
    stop(
      "Result of deconvolution function was not a list at positions ",
      paste(which(deconv_not_list), collapse = ", ")
    )
  }

  deconv_names <- purrr::map(deconvs, names)
  deconv_no_names <- purrr::map_lgl(
    deconv_names,
    function(x) is.null(x) || ("" %in% x)
  )

  if(any(deconv_no_names)) {
    stop(
      "Result of deconvolution function had empty or missing names at positions ",
      paste(which(deconv_no_names), collapse = ", ")
    )
  }

  for(name in unique(do.call(c, deconv_names))) {
    .spectra[[paste(".deconvolution", name, sep = "_")]] <- purrr::map(deconvs, name)
  }

  .spectra
}

gaussian_fun <- function(energy_kev, mu = 0, sigma = 1, height = 1) {
  height * exp(-0.5 * ((energy_kev - mu) / sigma) ^ 2)
}

#' XRF/EDS detector line shape (Gaussian with optional Hypermet tail and shelf)
#'
#' A real semiconductor-detector peak is not a pure Gaussian: incomplete charge collection adds a
#' low-energy exponential tail and a flat shelf (step). This returns the Gaussian plus optional
#' Hypermet tail/step terms (Campbell):
#' \deqn{h\big[\,e^{-\frac{(E-\mu)^2}{2\sigma^2}} + \tfrac{step}{2}\,\mathrm{erfc}(\tfrac{E-\mu}{\sqrt2\sigma})
#'   + tail\,e^{(E-\mu)/\beta}\,\mathrm{erfc}(\tfrac{E-\mu}{\sqrt2\sigma} + \tfrac{\sigma}{\sqrt2\beta})\,\big].}
#' With \code{tail = step = 0} (the default) it is exactly a Gaussian, so existing behaviour is
#' preserved. The tail carries a few percent to >10\% of the counts (worst at low energy and for
#' CdTe hole-tailing at high energy) and, if unmodelled, leaks into lower-energy neighbours.
#'
#' @param energy_kev Energies at which to evaluate (keV).
#' @param mu Peak centroid (keV).
#' @param sigma Gaussian standard deviation (keV).
#' @param height Peak height.
#' @param tail Relative amplitude of the low-energy exponential tail (0 = none).
#' @param step Relative amplitude of the low-energy shelf/step (0 = none).
#' @param beta Tail decay length (keV); defaults to \code{sigma}.
#'
#' @return A numeric vector the length of \code{energy_kev}.
#' @export
#'
xrf_lineshape <- function(energy_kev, mu = 0, sigma = 1, height = 1, tail = 0, step = 0, beta = NULL) {
  erfc <- function(x) 2 * stats::pnorm(x * sqrt(2), lower.tail = FALSE)
  d <- energy_kev - mu
  out <- exp(-0.5 * (d / sigma) ^ 2)
  if (step != 0) out <- out + step * 0.5 * erfc(d / (sqrt(2) * sigma))
  if (tail != 0) {
    b <- if (is.null(beta) || !is.finite(beta) || beta <= 0) sigma else beta
    out <- out + tail * exp(pmin(d / b, 50)) * erfc(d / (sqrt(2) * sigma) + sigma / (sqrt(2) * b))
  }
  height * out
}

.empty_pileup <- function() {
  tibble::tibble(element = character(0), primary_energy_kev = numeric(0), primary_sigma = numeric(0),
                 height = numeric(0), height_se = numeric(0), peak_area = numeric(0),
                 peak_area_se = numeric(0))
}

# Fixed-amplitude coincidence (pile-up) model from a first-pass element fit. For the strongest lines,
# the sum peak at E_i + E_j has a physically-fixed area: coincidence counts = 2*tau/T*A_i*A_j (i!=j)
# or tau/T*A_i^2 (i==j), with T the live time; its width adds in quadrature (sqrt(sigma_i^2+sigma_j^2)).
.xrf_pileup_model <- function(peaks, energy_kev, pileup_tau, livetime, energy_min_kev, energy_max_kev,
                              max_lines = 6L) {
  model <- numeric(length(energy_kev))
  el <- peaks[!grepl("^(scatter|sum|pileup|escape)_", peaks$element) &
                is.finite(peaks$peak_area) & peaks$peak_area > 0 &
                is.finite(peaks$primary_energy_kev) & is.finite(peaks$primary_sigma), , drop = FALSE]
  if (nrow(el) == 0) return(list(model = model, peaks = .empty_pileup()))
  el <- el[order(-el$peak_area), , drop = FALSE]
  k <- min(nrow(el), max_lines)
  el <- el[seq_len(k), ]

  out <- list()
  for (i in seq_len(k)) for (j in i:k) {
    e_sum <- el$primary_energy_kev[i] + el$primary_energy_kev[j]
    if (!is.finite(e_sum) || e_sum < energy_min_kev || e_sum > energy_max_kev) next
    area <- (if (i == j) 1 else 2) * (pileup_tau / livetime) * el$peak_area[i] * el$peak_area[j]
    if (!is.finite(area) || area <= 0) next
    s_sum <- sqrt(el$primary_sigma[i]^2 + el$primary_sigma[j]^2)
    height <- area / (s_sum * sqrt(2 * pi))
    model <- model + height * exp(-0.5 * ((energy_kev - e_sum) / s_sum)^2)
    out[[length(out) + 1]] <- tibble::tibble(
      element = paste0("pileup_", el$element[i], "_", el$element[j]),
      primary_energy_kev = e_sum, primary_sigma = s_sum,
      height = height, height_se = NA_real_, peak_area = area, peak_area_se = NA_real_
    )
  }
  list(model = model, peaks = if (length(out)) dplyr::bind_rows(out) else .empty_pileup())
}
