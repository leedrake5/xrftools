# Add light-element (Z = 6-9: C, N, O, F) K-line energies.
#
# The fluorescence yields and K-shell emission probabilities for C/N/O/F already ship in the
# package (x_ray_fluorescence_yields, x_ray_emission_probabilities start at Z = 6), but the
# *energies* were missing: x_ray_xrf_energies started at Ne (Z = 10). Without energy rows these
# elements can never appear in xrf_energies(), so light-element K-lines -- the backbone of
# SEM-EDS and low-energy work (O is ubiquitous) -- were unrepresentable. Here we add the K edge
# and the Kalpha (KL2/KL3) transition energies. For these light elements the L2/L3 spin-orbit
# splitting is unresolvable, so KL2 and KL3 share the Kalpha energy (as the existing table already
# does for Na).
#
# Energies (keV) from the X-ray Data Booklet (Thompson et al., LBNL):
#   K binding energy (edge):  C 0.2842  N 0.4099  O 0.5431  F 0.6967
#   Kalpha1 emission:         C 0.2770  N 0.3924  O 0.5249  F 0.6768

library(tibble)
library(dplyr)

setwd(rprojroot::find_root(rprojroot::is_r_package))

light <- tibble::tribble(
  ~element, ~z, ~edge_kev, ~kalpha_kev,
  "C",       6,  0.2842,    0.2770,
  "N",       7,  0.4099,    0.3924,
  "O",       8,  0.5431,    0.5249,
  "F",       9,  0.6967,    0.6768
)

# ---- x_ray_xrf_energies (the table consumed by xrf_energies) ------------------------------------
load("data/x_ray_xrf_energies.rda")

new_xrf <- light %>%
  tidyr::crossing(tibble::tibble(
    trans = c("KL3", "KL2"),
    trans_siegbahn = c("Kalpha1", "Kalpha2")
  )) %>%
  transmute(
    element, z,
    trans, edge = "K", trans_siegbahn,
    energy_kev = kalpha_kev,
    edge_kev
  ) %>%
  select(all_of(colnames(x_ray_xrf_energies)))

# idempotent: drop any prior C/N/O/F rows, then add and re-sort by Z
x_ray_xrf_energies <- x_ray_xrf_energies %>%
  filter(!(element %in% light$element)) %>%
  bind_rows(new_xrf) %>%
  arrange(z, trans_siegbahn)

usethis::use_data(x_ray_xrf_energies, overwrite = TRUE)

# ---- x_ray_energies (NIST-style source table, for provenance / future rebuilds) -----------------
load("data/x_ray_energies.rda")

template_row <- x_ray_energies[0, ]
mk_rows <- function(el, trans, kev) {
  out <- template_row[rep(1, length(trans)), ]
  out$element <- el
  out$trans <- trans
  out$direct_kev <- kev
  out
}

new_src <- dplyr::bind_rows(lapply(seq_len(nrow(light)), function(i) {
  el <- light$element[i]
  dplyr::bind_rows(
    mk_rows(el, "K edge", light$edge_kev[i]),
    mk_rows(el, c("KL2", "KL3"), c(light$kalpha_kev[i], light$kalpha_kev[i]))
  )
}))

x_ray_energies <- x_ray_energies %>%
  filter(!(element %in% light$element)) %>%
  bind_rows(new_src)

usethis::use_data(x_ray_energies, overwrite = TRUE)
