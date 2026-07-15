
#' Convert detector channels to energy
#'
#' Most multichannel-analyzer vendors (and all SEM-EDS systems) store a spectrum as counts versus
#' channel number together with a linear (occasionally quadratic) energy calibration, rather than a
#' pre-computed energy axis. This converts channel index to energy using
#'
#' \deqn{E = zero + gain \cdot channel + quad \cdot channel^2 .}
#'
#' \code{zero} and \code{gain} must be in the same energy unit as the returned value. To get
#' \code{energy_kev} (the unit used throughout this package), supply \code{zero}/\code{gain} in keV
#' (i.e. keV/channel for \code{gain}). A very common mistake is to mix keV/channel and eV/channel,
#' which produces a 1000x error -- be explicit about units.
#'
#' @param channel Channel index (numeric vector). Channels are usually 0-based; pass
#'   \code{seq_along(counts) - 1} if your MCA counts from zero.
#' @param gain Energy per channel (slope). Same unit as the desired output (keV/channel for keV).
#' @param zero Energy at channel 0 (intercept, same unit as output). Default 0.
#' @param quad Optional quadratic term (energy per channel^2). Default 0 (pure linear).
#'
#' @return A numeric vector of energies (same unit as \code{zero}/\code{gain}).
#' @export
#'
#' @examples
#' # 2048-channel spectrum at 20 eV/channel with a 0 offset -> keV
#' xrf_calibrate_energy(0:2047, gain = 0.02, zero = 0)[1:5]
#'
xrf_calibrate_energy <- function(channel, gain, zero = 0, quad = 0) {
  stopifnot(
    is.numeric(channel), is.numeric(gain), length(gain) == 1,
    is.numeric(zero), length(zero) == 1, is.numeric(quad), length(quad) == 1
  )
  zero + gain * channel + quad * channel^2
}

#' Build a spectra object from a raw counts vector and an energy calibration
#'
#' A vendor-neutral entry point for spectra that arrive as counts-per-channel plus a calibration
#' (the usual SEM-EDS / generic-MCA case), rather than through a vendor-specific reader such as
#' \link{read_xrf_panalytical}. It applies \link{xrf_calibrate_energy} to build the \code{energy_kev}
#' axis and, when \code{livetime} is supplied, a \code{cps} column (\code{counts / livetime}); the
#' raw \code{counts} and the live time are retained so that Poisson-weighted deconvolution
#' (\code{weighting = "poisson"}) can use them.
#'
#' @param counts Integer/numeric vector of counts per channel.
#' @param gain,zero,quad Energy calibration (keV), passed to \link{xrf_calibrate_energy}.
#' @param channel Channel indices for \code{counts}. Default \code{seq_along(counts) - 1} (0-based).
#' @param livetime Acquisition live time in seconds. If finite and > 0, a \code{cps} column is
#'   added and stored alongside as the \code{LiveTime} metadata column; otherwise \code{cps} is set
#'   equal to \code{counts}.
#' @param ... Additional metadata columns (e.g. \code{SampleIdent = "std-1"}), passed to
#'   \link{xrf_spectra}.
#' @param .validate Passed to \link{xrf_spectra}.
#'
#' @return A spectra tibble with one row.
#' @export
#'
#' @examples
#' cts <- rpois(2048, 5)
#' xrf_spectra_from_counts(cts, gain = 0.02, zero = 0, livetime = 60, SampleIdent = "demo")
#'
xrf_spectra_from_counts <- function(counts, gain, zero = 0, quad = 0,
                                    channel = seq_along(counts) - 1L,
                                    livetime = NA_real_, ..., .validate = TRUE) {
  stopifnot(is.numeric(counts), length(counts) >= 1, length(channel) == length(counts))
  energy_kev <- xrf_calibrate_energy(channel, gain = gain, zero = zero, quad = quad)
  has_lt <- length(livetime) == 1 && is.finite(livetime) && livetime > 0
  cps <- if (has_lt) counts / livetime else as.numeric(counts)
  spect <- tibble::tibble(
    channel = channel,
    energy_kev = energy_kev,
    counts = as.numeric(counts),
    cps = cps
  )
  xrf_spectra(spect, ..., LiveTime = if (has_lt) livetime else NA_real_, .validate = .validate)
}
