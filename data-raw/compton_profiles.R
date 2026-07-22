# Parse the Biggs/DABAX Compton profiles (see modern_sources/xraylib_DATA_NOTE.md) into
# x_ray_compton_profiles: the TOTAL electron momentum density J(p_z) per element (atomic units).
# J is symmetric in p_z; the table gives p_z >= 0. Normalization: integral of J over all p_z = Z.
library(tibble); library(dplyr)

all_elements <- c(
  "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
  "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb",
  "Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm",
  "Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
  "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
  "Md","No")

txt <- readLines("data-raw/modern_sources/ComptonProfiles_biggs.dat")
starts <- grep("^#S", txt)
blocks <- Map(function(a, b) txt[a:b], starts, c(starts[-1] - 1L, length(txt)))
x_ray_compton_profiles <- bind_rows(lapply(blocks, function(bl) {
  z <- as.integer(strsplit(trimws(sub("^#S", "", bl[1])), "\\s+")[[1]][1])
  rows <- bl[!grepl("^#", bl)]
  rows <- rows[nzchar(trimws(rows))]
  m <- read.table(text = rows)
  tibble(element = all_elements[z], z = z, pz_au = m[[1]], j_total = m[[2]])
}))

# sanity: 2 * integral_0^pmax J dpz ~= Z (the grid to pz = 100 captures ~all electrons)
norms <- x_ray_compton_profiles %>% group_by(element, z) %>%
  summarise(n_el = 2 * sum(diff(pz_au) * (head(j_total, -1) + tail(j_total, -1)) / 2), .groups = "drop")
stopifnot(max(abs(norms$n_el - norms$z) / norms$z) < 0.08)
stopifnot(all(x_ray_compton_profiles$j_total >= 0))
cat(sprintf("compton profiles: %d elements, %d rows; worst norm error %.1f%%\n",
            max(x_ray_compton_profiles$z), nrow(x_ray_compton_profiles),
            100 * max(abs(norms$n_el - norms$z) / norms$z)))
usethis::use_data(x_ray_compton_profiles, overwrite = TRUE)

# ---- per-SHELL profiles with occupancies and binding energies (for the bound-electron threshold /
# Compton-defect treatment): shell i contributes to scattering only where the energy transfer
# exceeds its binding energy B_i. UBIND is in eV; J_shell is per electron (total = sum occ * J).
x_ray_compton_profiles_shells <- bind_rows(lapply(blocks, function(bl) {
  z <- as.integer(strsplit(trimws(sub("^#S", "", bl[1])), "\\s+")[[1]][1])
  occ <- as.numeric(strsplit(trimws(sub("^#UOCCUP", "", bl[grepl("^#UOCCUP", bl)][1])), "\\s+")[[1]])
  bind_ev <- as.numeric(strsplit(trimws(sub("^#UBIND", "", bl[grepl("^#UBIND", bl)][1])), "\\s+")[[1]])
  stopifnot(length(occ) == length(bind_ev))
  rows <- bl[!grepl("^#", bl)]; rows <- rows[nzchar(trimws(rows))]
  m <- read.table(text = rows)
  stopifnot(ncol(m) == 2 + length(occ))
  # DABAX transcription repair: a per-shell value with occ * j > 1.2 * total at the same pz is
  # physically impossible (the total IS the occupancy-weighted shell sum) -- e.g. Eu shell 3 at
  # pz = 30 reads 4.610 where its neighbours are 0.011 / 0.0017 (a dropped exponent). Replace such
  # points by log-interpolation along pz between their valid neighbours.
  for (k in seq_along(occ)) {
    jk <- m[[2 + k]]
    bad <- which(occ[k] * jk > 1.2 * m[[2]] & jk > 0)
    if (length(bad)) {
      okp <- setdiff(which(jk > 0), bad)
      jk[bad] <- exp(approx(m[[1]][okp], log(jk[okp]), xout = m[[1]][bad], rule = 2)$y)
      m[[2 + k]] <- jk
      message(sprintf("repaired %d corrupt point(s): Z=%d shell %d", length(bad), z, k))
    }
  }
  bind_rows(lapply(seq_along(occ), function(k) {
    tibble(element = all_elements[z], z = z, shell = k, occupancy = occ[k],
           binding_kev = bind_ev[k] / 1000, pz_au = m[[1]], j_shell = m[[2 + k]])
  }))
}))
# sanity: occupancies close over Z; the occupancy-weighted shell sum reproduces the total profile
occ_tot <- x_ray_compton_profiles_shells %>% distinct(z, shell, occupancy) %>%
  group_by(z) %>% summarise(n = sum(occupancy), .groups = "drop")
stopifnot(all(occ_tot$n == occ_tot$z))
chk <- x_ray_compton_profiles_shells %>% group_by(z, pz_au) %>%
  summarise(tot = sum(occupancy * j_shell), .groups = "drop") %>%
  left_join(x_ray_compton_profiles, by = c("z", "pz_au")) %>%
  group_by(z) %>% mutate(j0 = max(j_total)) %>% ungroup() %>%
  filter(j_total > 0.02 * j0)                       # judge agreement where J is significant
stopifnot(max(abs(chk$tot - chk$j_total) / chk$j_total) < 0.06)
cat(sprintf("per-shell profiles: %d rows, %d elements\n",
            nrow(x_ray_compton_profiles_shells), max(x_ray_compton_profiles_shells$z)))
usethis::use_data(x_ray_compton_profiles_shells, overwrite = TRUE)
