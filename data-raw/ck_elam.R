# Generate x_ray_coster_kronig_probabilities_elam (data/) from the modern Elam, Ravel & Sieber (2002)
# database (data-raw/modern_sources/elam.dat, vendored from XrayDB). Same schema as the built-in
# x_ray_coster_kronig_probabilities: the DIRECT Coster-Kronig probabilities (the "CK" lines, NOT "CKtotal",
# which fold in intermediate states) as f_ij per source subshell, where i/j are the subshell indices within
# the L or M family (L1->1 .. M5->5), so a transition SHELL_i -> SHELL_j is labelled f{i}{j} (e.g. L1->L3 = f13).
#
# Run with:  Rscript data-raw/ck_elam.R   (from the package root)

ln <- readLines("data-raw/modern_sources/elam.dat"); ln <- ln[!grepl("^//", ln)]
sp <- function(s) strsplit(s, "[[:space:]]+")[[1]]
subidx <- function(sh) suppressWarnings(as.integer(substr(sh, 2, nchar(sh))))

rows <- list(); cur <- NA_character_; z <- NA_integer_; cur_edge <- NA_character_
for (s in trimws(ln)) {
  if (startsWith(s, "Element ")) { f <- sp(s); cur <- f[2]; z <- as.integer(f[3]) }
  else if (startsWith(s, "Edge ")) { cur_edge <- sp(s)[2] }
  else if (startsWith(s, "CK ") && !startsWith(s, "CKtotal")) {
    f <- sp(s)[-1]
    for (k in seq(1, length(f) - 1, 2)) {
      i <- subidx(cur_edge); j <- subidx(f[k])
      if (is.na(i) || is.na(j)) next
      rows[[length(rows) + 1]] <- data.frame(
        element = cur, z = z, shell = cur_edge,
        coster_kronig_trans = sprintf("f%d%d", i, j),
        coster_kronig_prob = as.numeric(f[k + 1]), stringsAsFactors = FALSE)
    }
  }
}
x_ray_coster_kronig_probabilities_elam <- tibble::as_tibble(do.call(rbind, rows))
save(x_ray_coster_kronig_probabilities_elam,
     file = "data/x_ray_coster_kronig_probabilities_elam.rda", compress = "xz")
cat("wrote", nrow(x_ray_coster_kronig_probabilities_elam), "Elam CK rows;",
    "shells:", paste(sort(unique(x_ray_coster_kronig_probabilities_elam$shell)), collapse = " "), "\n")
