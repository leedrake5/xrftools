
# Escape lines (detector characteristic X-rays that can leave the crystal) for an absorber element:
# the K edge, the Kalpha/Kbeta energies, and their emission probabilities (omega * branching).
.xrf_escape_lines <- function(absorber) {
  xe <- xrftools::x_ray_xrf_energies
  ep <- xrftools::x_ray_emission_probabilities
  ka_E <- xe$energy_kev[xe$element == absorber & xe$trans_siegbahn == "Kalpha1"][1]
  kb_E <- xe$energy_kev[xe$element == absorber & xe$trans_siegbahn == "Kbeta1"][1]
  k_edge <- xe$edge_kev[xe$element == absorber & xe$edge == "K"][1]
  ka_p <- sum(ep$emission_probability[ep$element == absorber & ep$trans %in% c("KL2", "KL3")])
  kb_p <- sum(ep$emission_probability[ep$element == absorber & ep$trans %in% c("KM2", "KM3")])
  lines <- data.frame(
    energy = c(ka_E, kb_E),
    emission_prob = c(ka_p, kb_p),
    label = c("Ka", "Kb"),
    stringsAsFactors = FALSE
  )
  lines <- lines[is.finite(lines$energy) & is.finite(lines$emission_prob) & lines$emission_prob > 0, ]
  list(k_edge = k_edge, lines = lines)
}

# Total detector escape probability at each incident energy: the fraction of photoelectric events whose
# characteristic X-ray (Si/Ge K, or Cd/Te K for CdTe) escapes the crystal, so the event is NOT recorded in
# the full-energy peak. Same per-line formula as xrf_escape_peaks() (kept in sync so the escaped fraction
# and the escape-satellite amplitudes conserve total photoabsorption), summed over the detector's K lines
# and absorbers. Used by xrf_detector_efficiency(full_energy_peak = TRUE) to turn the photoabsorption
# probability into the full-energy-peak fraction (1 - f_escape).
.xrf_detector_escape_fraction <- function(energy_kev, material) {
  info <- .xrf_material_info(material)
  absorbers <- names(info$composition)
  f_tot <- rep(0, length(energy_kev))
  mu_i_pe  <- xrf_mass_attenuation(material, energy_kev, type = "photoelectric")  # for the absorber share
  mu_i_tot <- xrf_mass_attenuation(material, energy_kev, type = "total")          # for the escape bracket (W13)
  for (A in absorbers) {
    esc <- .xrf_escape_lines(A)
    if (is.na(esc$k_edge) || nrow(esc$lines) == 0) next
    jump_k <- xrf_absorption_jump_ratio(A, "K")
    comp_A <- info$composition[[A]]
    above <- is.finite(energy_kev) & energy_kev > esc$k_edge
    if (!any(above)) next
    w_A <- comp_A * .xrf_element_pe_mu(A, energy_kev) / mu_i_pe  # this absorber's share of the absorption
    for (j in seq_len(nrow(esc$lines))) {
      mu_esc_tot <- xrf_mass_attenuation(material, esc$lines$energy[j], type = "total")
      bracket <- 1 - (mu_esc_tot / mu_i_tot) * log(1 + mu_i_tot / mu_esc_tot)
      f <- 0.5 * jump_k * esc$lines$emission_prob[j] * bracket * w_A
      f[!above | !is.finite(f) | f < 0] <- 0
      f_tot <- f_tot + f
    }
  }
  pmin(pmax(f_tot, 0), 0.99)
}

#' Build detector escape-peak templates
#'
#' When a photon of energy \eqn{E} is photoelectrically absorbed in the detector, the detector's own
#' characteristic X-ray (Si K, Ge K, or Cd/Te K for CdTe) can escape the crystal, so the event is
#' recorded at \eqn{E - E_{esc}} -- an "escape peak". These are small for Si (~1%) but reach
#' 10-15\% just above the Ge K edge and are large for CdTe, where unmodelled they masquerade as
#' spurious elements. This returns escape satellite rows for each parent line in \code{peaks} whose
#' energy is above the detector's K edge; each is tagged with the parent's \code{element} so that,
#' when appended to \code{peaks}, it is tied to the same fitted amplitude (no free parameter).
#'
#' \strong{K-shell escape only.} Only the detector's K X-rays are modelled. For CdTe this means no escape
#' is produced below the Cd/Te K edges (~26.7/31.8 keV) -- i.e. across most of its usual working band -- even
#' though the Cd/Te \strong{L} X-rays (~3-4.6 keV) can escape there; the escape fraction (and hence the
#' \code{full_energy_peak} loss in \link{xrf_detector_efficiency}) is therefore a lower bound for CdTe below
#' its K edges. L-escape is negligible for Si/Ge (their L X-rays are absorbed within a few nm).
#'
#' The escape fraction follows
#' \deqn{f = \tfrac{1}{2}\,\omega_K\,(1-1/r)\,[\,1 - (\mu_{esc}/\mu_i)\ln(1 + \mu_i/\mu_{esc})\,]}
#' with \eqn{(1-1/r)} the detector K jump ratio, \eqn{\omega_K\cdot}branching supplied by the
#' emission probability, and \eqn{\mu_i}, \eqn{\mu_{esc}} the detector attenuation at the parent and
#' escape energies (the bracket tends to 1 just above the edge, where escape is maximal). For CdTe
#' the Cd and Te escapes are each weighted by that element's share of the absorption.
#'
#' @param peaks A peaks data frame (needs \code{element}, \code{energy_kev}; uses
#'   \code{relative_peak_intensity} if present).
#' @param detector_type A preset name (\link{xrf_detector_presets}) giving the detector material.
#' @param active_material Override the detector material directly ("Si", "Ge", "CdTe").
#' @param min_escape_fraction Drop escape peaks weaker than this fraction of their parent.
#'
#' @return A tibble of escape peak templates (\code{element}, \code{trans}, \code{energy_kev},
#'   \code{relative_peak_intensity}), possibly with zero rows.
#' @export
#'
#' @examples
#' pk <- data.frame(element = "As", energy_kev = 10.53, relative_peak_intensity = 1)
#' xrf_escape_peaks(pk, detector_type = "HPGe")   # Ge escape ~0.6 keV below As Ka
#'
xrf_escape_peaks <- function(peaks, detector_type = NULL, active_material = NULL,
                             min_escape_fraction = 1e-4) {
  material <- if (!is.null(active_material)) {
    active_material
  } else {
    .resolve_detector_geometry(detector_type)$active_material
  }
  info <- .xrf_material_info(material)
  absorbers <- names(info$composition)

  if (!("relative_peak_intensity" %in% colnames(peaks))) peaks$relative_peak_intensity <- 1

  empty <- tibble::tibble(
    element = character(0), trans = character(0), energy_kev = numeric(0),
    relative_peak_intensity = numeric(0)
  )
  out <- list()

  for (A in absorbers) {
    esc <- .xrf_escape_lines(A)
    if (is.na(esc$k_edge) || nrow(esc$lines) == 0) next
    jump_k <- xrf_absorption_jump_ratio(A, "K")
    comp_A <- info$composition[[A]]

    parent_ok <- which(is.finite(peaks$energy_kev) & peaks$energy_kev > esc$k_edge)
    if (length(parent_ok) == 0) next
    E_par <- peaks$energy_kev[parent_ok]
    # w_A (this absorber's share of the absorption) follows PHOTOELECTRIC absorption -- escape only happens
    # after a photoelectric event -- so its denominator is the material photoelectric mu.
    mu_i_pe <- xrf_mass_attenuation(material, E_par, type = "photoelectric")
    # The escape bracket models the depth distribution of the absorbed photon and the reabsorption of the
    # escaping detector X-ray, both governed by TOTAL attenuation (W13; below ~30 keV total ~= PE so this is
    # unchanged, but for Ge/CdTe above ~70 keV where Compton makes PE < total it was ~15-30% too low).
    mu_i_tot <- xrf_mass_attenuation(material, E_par, type = "total")
    w_A <- comp_A * .xrf_element_pe_mu(A, E_par) / mu_i_pe

    for (j in seq_len(nrow(esc$lines))) {
      E_esc <- esc$lines$energy[j]
      mu_esc_tot <- xrf_mass_attenuation(material, E_esc, type = "total")
      bracket <- 1 - (mu_esc_tot / mu_i_tot) * log(1 + mu_i_tot / mu_esc_tot)
      f <- 0.5 * jump_k * esc$lines$emission_prob[j] * bracket * w_A
      esc_E <- E_par - E_esc
      keep <- is.finite(f) & f > min_escape_fraction & esc_E > 0
      if (!any(keep)) next
      out[[length(out) + 1]] <- tibble::tibble(
        element = peaks$element[parent_ok][keep],
        trans = paste0("escape_", A, esc$lines$label[j]),
        energy_kev = esc_E[keep],
        relative_peak_intensity = peaks$relative_peak_intensity[parent_ok][keep] * f[keep]
      )
    }
  }

  if (length(out) == 0) return(empty)
  dplyr::bind_rows(out)
}
