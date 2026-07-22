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
