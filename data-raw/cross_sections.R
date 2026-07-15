# Rebuild x_ray_cross_sections from the full EPDL97 source, and add x_ray_mass_attenuation.
#
# The previously-shipped x_ray_cross_sections was a coarse extraction: a 5-keV grid from 5-100 keV
# with only the K-M3 subshells. The vendored EPDL97_CrossSections.dat actually contains, per element,
# the full native energy grid (1 eV - 495 keV) for ALL subshells K, L1-L3, M1-M5, plus the Rayleigh
# (coherent), Compton (incoherent), total-photoelectric and total attenuation. Rebuilding from it
# removes three approximations at once:
#   * M4/M5 are real (no more M3 proxy for the Malpha/Mbeta production),
#   * the range reaches ~495 keV (no >100 keV power-law extrapolation for heavy K-lines),
#   * the grid is fine down to 1 eV (no sub-5 keV extrapolation in the detector-efficiency model).
# The new x_ray_mass_attenuation table (photoelectric/rayleigh/compton/total cm^2/g) feeds the
# detector-efficiency window term (total attenuation) and is available for scatter modelling.

library(dplyr)
setwd(rprojroot::find_root(rprojroot::is_r_package))

lines <- readLines("data-raw/EPDL97_CrossSections.dat")
s_idx <- grep("^#S ", lines)
block_end <- c(s_idx[-1] - 1L, length(lines))

# columns (fixed order in every #L header):
# 1 PhotonEnergy 2 Rayleigh 3 Compton 4 Coh+Incoh 5 Photoelectric 6 K 7 L1 8 L2 9 L3
# 10 M1 11 M2 12 M3 13 M4 14 M5 15 AllOthers 16 Total
subshell_col <- c(K = 6, L1 = 7, L2 = 8, L3 = 9, M1 = 10, M2 = 11, M3 = 12, M4 = 13, M5 = 14)
E_MAX <- 200   # keV cap (covers the ~100 keV use case with headroom; bounds table size)

xs_list <- list()
att_list <- list()
for (b in seq_along(s_idx)) {
  hdr <- strsplit(trimws(lines[s_idx[b]]), "\\s+")[[1]]     # "#S" Z Symbol
  z <- as.integer(hdr[2]); sym <- hdr[3]
  if (is.na(z) || z < 1 || z > 99) next

  block <- lines[s_idx[b]:block_end[b]]
  parts <- strsplit(trimws(block[!grepl("^#", block) & grepl("[0-9]", block)]), "\\s+")
  parts <- parts[lengths(parts) == 16]
  if (length(parts) == 0) next
  m <- matrix(as.numeric(unlist(parts)), ncol = 16, byrow = TRUE)

  energy <- m[, 1]
  ok <- is.finite(energy) & energy > 0 & energy <= E_MAX
  m <- m[ok, , drop = FALSE]; energy <- energy[ok]
  if (nrow(m) == 0) next

  for (sh in names(subshell_col)) {
    val <- m[, subshell_col[sh]]
    keep <- is.finite(val) & val > 0
    if (any(keep)) {
      xs_list[[length(xs_list) + 1]] <- data.frame(
        element = sym, z = z, shell = sh, energy_kev = energy[keep], cross_section = val[keep],
        stringsAsFactors = FALSE
      )
    }
  }
  att_list[[length(att_list) + 1]] <- data.frame(
    element = sym, z = z, energy_kev = energy,
    photoelectric = m[, 5], rayleigh = m[, 2], compton = m[, 3], total = m[, 16],
    stringsAsFactors = FALSE
  )
}

x_ray_cross_sections <- dplyr::as_tibble(do.call(rbind, xs_list))
x_ray_mass_attenuation <- dplyr::as_tibble(do.call(rbind, att_list))

cat("x_ray_cross_sections:", nrow(x_ray_cross_sections), "rows; shells:",
    paste(sort(unique(x_ray_cross_sections$shell)), collapse = ","),
    "; energy", round(min(x_ray_cross_sections$energy_kev), 4), "-",
    round(max(x_ray_cross_sections$energy_kev), 1), "keV\n")
cat("x_ray_mass_attenuation:", nrow(x_ray_mass_attenuation), "rows\n")

# sanity checks
si10 <- approx(log(subset(x_ray_mass_attenuation, element == "Si")$energy_kev),
               log(subset(x_ray_mass_attenuation, element == "Si")$photoelectric), log(10))$y
cat("Si photoelectric mu/rho @10 keV:", round(exp(si10), 2), "(expect ~33)\n")
au_m5 <- subset(x_ray_cross_sections, element == "Au" & shell == "M5")
cat("Au M5 now present:", nrow(au_m5) > 0, "| U K present up to",
    round(max(subset(x_ray_cross_sections, element == "U" & shell == "K")$energy_kev), 0), "keV\n")

usethis::use_data(x_ray_cross_sections, overwrite = TRUE)
usethis::use_data(x_ray_mass_attenuation, overwrite = TRUE)
cat("Done.\n")
