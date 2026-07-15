# Parse the Bote & Salvat (2008) electron-impact ionization cross-section coefficients from the
# NIST BoteSalvatICX.jl package (public domain / Unlicense; vendored as data-raw/BoteSalvatICX_*).
#
# Each element (Z = 1..99) carries, per available subshell (1->K, 2->L1, 3->L2, 4->L3, 5->M1 ...
# 9->M5): the edge energy (eV), and the fitting coefficients Be, Anlj, G[4] and A[5] used by the
# Bote-Salvat cross-section formula (see xrf_electron_ionization_cross_section, method="bote-salvat").
#
# Reference: Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact
# with the distorted-wave and plane-wave Born approximations. Phys. Rev. A 77, 042701.

library(dplyr)
setwd(rprojroot::find_root(rprojroot::is_r_package))

txt <- paste(readLines("data-raw/BoteSalvatICX_xione.jl"), collapse = "\n")
# keep only the const data tuple (stop at the first docstring after it)
txt <- sub('.*const BoteSalvatElectron = \\(', "", txt)
txt <- sub('\n"""[[:print:][:space:]]*$', "", txt)

# split into per-element chunks
chunks <- strsplit(txt, "BoteSalvatElementDatum\\(")[[1]]
chunks <- chunks[nzchar(trimws(chunks))]

parse_vec <- function(s) as.numeric(strsplit(gsub("[\\[\\]]", "", s, perl = TRUE), ",")[[1]])
parse_mat <- function(s) {
  rows <- strsplit(gsub("[\\[\\]]", "", s, perl = TRUE), ";")[[1]]
  do.call(rbind, lapply(rows, function(r) as.numeric(strsplit(trimws(r), "\\s+")[[1]])))
}

rows <- list()
for (ch in chunks) {
  z <- as.integer(sub("^\\s*([0-9]+)\\s*,.*", "\\1", ch))
  if (is.na(z)) next
  brk <- regmatches(ch, gregexpr("\\[[^\\]]*\\]", ch, perl = TRUE))[[1]]
  if (length(brk) < 5) next
  be <- parse_vec(brk[1]); anlj <- parse_vec(brk[2]); g <- parse_mat(brk[3])
  edge <- parse_vec(brk[4]); a <- parse_mat(brk[5])
  ne <- length(edge)
  rows[[length(rows) + 1]] <- data.frame(
    z = z, subshell = seq_len(ne), edge_ev = edge, be = be[seq_len(ne)], anlj = anlj[seq_len(ne)],
    g1 = g[, 1], g2 = g[, 2], g3 = g[, 3], g4 = g[, 4],
    a1 = a[, 1], a2 = a[, 2], a3 = a[, 3], a4 = a[, 4], a5 = a[, 5]
  )
}
x_ray_bote_salvat <- dplyr::as_tibble(do.call(rbind, rows))

cat("Parsed", length(unique(x_ray_bote_salvat$z)), "elements,", nrow(x_ray_bote_salvat), "subshell rows\n")
cat("Z range:", paste(range(x_ray_bote_salvat$z), collapse = "-"),
    "| subshells:", paste(range(x_ray_bote_salvat$subshell), collapse = "-"), "\n")
# spot check Fe (Z=26) K-shell (subshell 1): edge ~7112 eV
fek <- subset(x_ray_bote_salvat, z == 26 & subshell == 1)
cat("Fe K edge (eV):", round(fek$edge_ev, 1), " a1..a5:", paste(signif(c(fek$a1,fek$a2,fek$a3,fek$a4,fek$a5),3), collapse=","), "\n")

usethis::use_data(x_ray_bote_salvat, overwrite = TRUE)
cat("Done.\n")
