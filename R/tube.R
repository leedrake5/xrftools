
# Transmission through a primary-beam filter parsed from a tube$filter string like "Al" or "Al 100".
# The number is a thickness in microns (default 100). Returns 1 (no filter) when unparseable.
.xrf_filter_transmission <- function(filter, energy_kev) {
  if (is.null(filter) || length(filter) != 1 || is.na(filter) || !nzchar(trimws(filter))) {
    return(rep(1, length(energy_kev)))
  }
  toks <- strsplit(trimws(as.character(filter)), "\\s+")[[1]]
  mat <- toks[1]
  if (!(mat %in% all_elements) && !(mat %in% c("CdTe"))) return(rep(1, length(energy_kev)))
  thick_um <- suppressWarnings(as.numeric(toks[2])); if (is.na(thick_um)) thick_um <- 100
  rho <- .xrf_material_info(mat)$density
  if (!is.finite(rho)) return(rep(1, length(energy_kev)))   # unknown density -> no filter
  exp(-xrf_mass_attenuation(mat, energy_kev, type = "total") * rho * thick_um * 1e-4)
}

#' X-ray tube emission spectrum (Ebel model)
#'
#' The photon output of an X-ray tube versus energy: a bremsstrahlung continuum (Ebel 2001) shaped by
#' self-absorption in the anode, plus the anode's characteristic lines, optionally attenuated by a
#' primary-beam filter. Values are relative (an overall constant is dropped) -- what the
#' fundamental-parameters excitation integral needs. Use it to weight excitation across the tube
#' spectrum instead of assuming a single monochromatic energy.
#'
#' @param tube An \link{xrf_tube} (needs \code{anode} and \code{kv}; \code{filter} optional).
#' @param energy_kev Photon energies to evaluate (keV).
#' @param takeoff_deg Anode take-off angle (degrees) for the anode self-absorption.
#' @param char_peak_ratio Approximate peak height of the strongest characteristic line relative to
#'   the local continuum (tube characteristic lines tower over the continuum; this sets the
#'   line/continuum balance and is a modelling approximation).
#' @param filter Apply the tube's primary-beam \code{filter} (if any)?
#'
#' @return Relative photon flux at each \code{energy_kev} (0 at or above the endpoint \code{kv}).
#' @export
#'
#' @references
#' Ebel, H. (2001). X-ray tube spectra. X-Ray Spectrometry 28, 255.
#'
#' @examples
#' e <- seq(1, 40, 0.05)
#' plot(e, xrf_tube_spectrum(xrf_tube("Rh", kv = 40), e), type = "l")
#'
xrf_tube_spectrum <- function(tube, energy_kev, takeoff_deg = 40, char_peak_ratio = 100,
                              filter = TRUE) {
  stopifnot(inherits(tube, "xrf_tube"))
  anode <- strsplit(trimws(as.character(tube$anode)), "\\s+")[[1]][1]
  if (!(anode %in% all_elements)) stop("tube anode '", tube$anode, "' is not a recognized element.")
  Z <- match(anode, all_elements)
  E0 <- tube$kv
  E <- energy_kev

  # Ebel bremsstrahlung continuum ~ Z * (E0/E - 1)^x / E, with the Ebel energy exponent
  x <- 1.109 - 0.00435 * Z + 0.00175 * E0
  cont <- Z * pmax(E0 / E - 1, 0)^x / E

  # self-absorption of the continuum escaping the anode (Ebel), using the anode mass attenuation and
  # an approximate mean production mass-depth
  mu_a <- xrf_mass_attenuation(anode, E, type = "total")
  sin_eps <- sin(takeoff_deg * pi / 180)
  m_depth <- 0.787e-5 * sqrt(11.5 * Z) * E0^1.5 + 0.735e-6 * E0^2   # g/cm^2 (approx.)
  tau <- 2 * mu_a * m_depth / max(sin_eps, 1e-3)
  f_a <- ifelse(tau > 1e-9, (1 - exp(-tau)) / tau, 1)
  cont <- cont * f_a
  cont[!is.finite(cont) | E >= E0 | E <= 0] <- 0

  if (isTRUE(filter)) cont <- cont * .xrf_filter_transmission(tube$filter, E)

  # characteristic lines of the anode (as narrow Gaussians towering over the local continuum)
  out <- cont
  if (E0 > 0 && char_peak_ratio > 0) {
    lines <- tryCatch(xrf_energies(anode, beam_energy_kev = E0, min_relative_intensity = 0.02),
                      error = function(e) NULL)
    if (!is.null(lines) && nrow(lines) > 0) {
      s_line <- 0.03
      cont_at <- function(en) {
        v <- Z * pmax(E0 / en - 1, 0)^x / en *
          ifelse(2 * xrf_mass_attenuation(anode, en, "total") * m_depth / max(sin_eps, 1e-3) > 1e-9,
                 (1 - exp(-2 * xrf_mass_attenuation(anode, en, "total") * m_depth / max(sin_eps, 1e-3))) /
                   (2 * xrf_mass_attenuation(anode, en, "total") * m_depth / max(sin_eps, 1e-3)), 1)
        if (isTRUE(filter)) v <- v * .xrf_filter_transmission(tube$filter, en)
        v
      }
      for (i in seq_len(nrow(lines))) {
        ei <- lines$energy_kev[i]
        if (ei >= E0) next
        amp <- char_peak_ratio * lines$relative_peak_intensity[i] * cont_at(ei)
        out <- out + amp * exp(-0.5 * ((E - ei) / s_line)^2)
      }
    }
  }
  out
}

# Polychromatic excitation of a subshell integrated over the tube spectrum:
#   integral of sigma_shell(E) * N_tube(E) dE from the edge up to the endpoint.
.xrf_tube_excitation <- function(tube, element, shell, n_grid = 400L) {
  E0 <- tube$kv
  grid <- seq(0.1, E0 * 0.999, length.out = n_grid)
  de <- grid[2] - grid[1]
  n_tube <- xrf_tube_spectrum(tube, grid)
  sig <- xrf_photoionization_cross_section(element, shell, grid)
  sig[!is.finite(sig)] <- 0                       # 0 below the edge
  sum(sig * n_tube) * de
}
