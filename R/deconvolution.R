
#' Deconvolute a spectrum using gaussian least squares
#'
#' @inheritParams xrf_add_deconvolution_fun
#' @param .values,response Values to use for the deconvolution
#' @param energy_kev,.energy_kev Energy values corresponding to response/.values
#' @param energy_min_kev,energy_max_kev Constrain the energies that should be deconvoluted
#' @param peaks A data frame of columns with columns "element", "
#' @param default_sigma The default standard deviation of the peaks (can also be passed in as
#'   a columm in \code{peaks} if this value varies between peaks)
#'
#' @return A modified version of .spectra
#' @export
#'
xrf_add_deconvolution_gls <- function(.spectra, .energy_kev = .data$.spectra$energy_kev,
                                      .values = .data$.spectra$cps - .data$.spectra$baseline,
                                      energy_min_kev = -Inf, energy_max_kev = Inf,
                                      peaks = xrf_energies(), default_sigma = 0.07, .env = parent.frame()) {
  .energy_kev <- enquo(.energy_kev)
  .values <- enquo(.values)
  peaks <- enquo(peaks)
  default_sigma <- enquo(default_sigma)
  energy_min_kev <- enquo(energy_min_kev)
  energy_max_kev <- enquo(energy_max_kev)

  xrf_add_deconvolution_fun(
    .spectra,
    xrf_deconvolute_gaussian_least_squares,
    !!.energy_kev,
    !!.values,
    peaks = !!peaks,
    default_sigma = !!default_sigma,
    energy_min_kev = !!energy_min_kev,
    energy_max_kev = !!energy_max_kev,
    .env = .env
  )
}

#' @rdname xrf_add_deconvolution_gls
#' @export
xrf_deconvolute_gaussian_least_squares <- function(energy_kev, response, peaks = xrf_energies(), default_sigma = 0.07,
                                                   energy_min_kev = -Inf, energy_max_kev = Inf) {

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

  if(!("sigma" %in% colnames(peaks))) {
    peaks$sigma <- default_sigma
  } else {
    stopifnot(is.numeric(peaks$sigma))
    peaks$sigma <- dplyr::coalesce(peaks$sigma, default_sigma)
  }

  # filter energies and responses
  within_range <- (energy_kev >= energy_min_kev) & (energy_kev <= energy_max_kev)
  energy_kev <- energy_kev[within_range]
  response <- response[within_range]

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

  # calculate gaussian responses for each row in peaks
  # summarise combined response per element
  responses <- peaks %>%
    dplyr::select("element", "energy_kev", "sigma", "relative_peak_height", "relative_peak_intensity") %>%
    dplyr::mutate(
      response_element = purrr::pmap(
        list(.data$energy_kev, .data$sigma, .data$relative_peak_height),
        function(mu, sigma, height) gaussian_fun(!!energy_kev, mu = mu, sigma = sigma, height = height)
      )
    ) %>%
    # sum all the responses per element
    dplyr::group_by(.data$element) %>%
    dplyr::summarise(
      primary_energy_kev = .data$energy_kev[which.max(.data$relative_peak_intensity)],
      primary_sigma = .data$sigma[which.max(.data$relative_peak_intensity)],
      response_element = list(purrr::reduce(.data$response_element, `+`))
    ) %>%
    dplyr::ungroup()

  # least squares estimation using QR decomposition (faster than lm())
  X <- do.call(cbind, responses$response_element)
  colnames(X) <- responses$element
  elements <- responses$element
  n <- length(response)
  p <- ncol(X)

  # QR decomposition for solving and SE calculation
  qr_X <- qr(X)
  coef <- qr.coef(qr_X, response)
  coef[is.na(coef)] <- 0  # Handle singular columns like lm() does
  names(coef) <- elements

  # Calculate fitted values and residuals
  response_fit <- qr.fitted(qr_X, response)
  residuals <- response - response_fit

  # Calculate standard errors
  rank <- qr_X$rank
  sigma2 <- sum(residuals^2) / max(n - rank, 1)

  # Standard errors from (X'X)^-1 diagonal using chol2inv on R
  coef_se <- rep(NA_real_, p)
  names(coef_se) <- elements
  if (rank > 0 && rank == p) {
    # Full rank case - fast path
    R <- qr.R(qr_X)
    XtX_inv_diag <- rowSums(backsolve(R, diag(p))^2)
    coef_se <- sqrt(XtX_inv_diag * sigma2)
    names(coef_se) <- elements
  } else if (rank > 0) {
    # Rank deficient - use chol2inv on the non-singular part
    R <- qr.R(qr_X)[1:rank, 1:rank, drop = FALSE]
    R_inv <- backsolve(R, diag(rank))
    pivot <- qr_X$pivot[1:rank]
    coef_se[pivot] <- sqrt(rowSums(R_inv^2) * sigma2)
  }

  # Calculate R-squared
  ss_res <- sum(residuals^2)
  ss_tot <- sum((response - mean(response))^2)
  r2 <- 1 - ss_res / ss_tot
  r2_adj <- 1 - (1 - r2) * (n - 1) / (n - p)

  # Build response data frame
  df <- tibble::tibble(
    energy_kev = energy_kev,
    response = response,
    response_fit = response_fit
  )

  # Build peaks output
  responses$response_element <- NULL
  responses$height <- coef[elements]
  responses$height_se <- coef_se[elements]
  responses$peak_area <- responses$height * responses$primary_sigma * sqrt(2 * pi)
  responses$peak_area_se <- responses$height_se * responses$primary_sigma * sqrt(2 * pi)

  # Build components efficiently using base R
  n_energy <- length(energy_kev)
  n_elements <- length(elements)
  components <- tibble::tibble(
    energy_kev = rep(energy_kev, n_elements),
    element = rep(elements, each = n_energy),
    response_fit = as.vector(X) * rep(coef, each = n_energy),
    height = rep(coef, each = n_energy)
  )

  structure(
    list(
      fit = tibble::tibble(
        r2 = r2,
        r2_adj = r2_adj,
        chi_sq = sum(residuals^2 / response)
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
