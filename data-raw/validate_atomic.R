# Validation harness for the atomic-data modernization (MODERNIZATION_PLAN.md item 1).
#
# Parses the modern Elam, Ravel & Sieber (2002) consolidated database (data-raw/modern_sources/elam.dat,
# vendored from the XrayDB project) and quantifies how the package's current (EADL97 / Veigele-vintage)
# fluorescence yields, absorption jump ratios, edge energies and Coster-Kronig probabilities differ from it.
# This is the prerequisite for any table swap: it says, per element and shell, how large the modernization
# delta actually is, and flags where the two references disagree beyond a tolerance.
#
# Run with:  Rscript data-raw/validate_atomic.R   (from the package root, with the package loadable)

suppressMessages(library(devtools)); suppressMessages(load_all(".", quiet = TRUE))

elam_path <- "data-raw/modern_sources/elam.dat"
stopifnot(file.exists(elam_path))
ln <- readLines(elam_path); ln <- ln[!grepl("^//", ln)]

# --- parse edges (energy / fluorescence yield / jump ratio) and Coster-Kronig probabilities ---------
edges <- list(); cks <- list(); cur <- NA_character_; cur_edge <- NA_character_
sp <- function(s) strsplit(s, "[[:space:]]+")[[1]]
for (s in trimws(ln)) {
  if (startsWith(s, "Element ")) {
    cur <- sp(s)[2]
  } else if (startsWith(s, "Edge ")) {
    f <- sp(s); cur_edge <- f[2]
    edges[[length(edges) + 1]] <- data.frame(
      el = cur, edge = f[2], e_ev = as.numeric(f[3]),
      yield = as.numeric(f[4]), jump = as.numeric(f[5]), stringsAsFactors = FALSE)
  } else if (startsWith(s, "CK ") && !startsWith(s, "CKtotal")) {
    f <- sp(s)[-1]
    for (k in seq(1, length(f) - 1, 2))
      cks[[length(cks) + 1]] <- data.frame(el = cur, from = cur_edge, to = f[k],
                                            prob = as.numeric(f[k + 1]), stringsAsFactors = FALSE)
  }
}
elam <- do.call(rbind, edges); elam_ck <- do.call(rbind, cks)

# --- compare a representative light -> heavy set ----------------------------------------------------
pct <- function(a, b) if (is.na(a) || is.na(b) || b == 0) "NA" else sprintf("%+.0f%%", 100 * (a - b) / b)
cmp <- function(elem, edge) {
  r <- elam[elam$el == elem & elam$edge == edge, ]
  if (!nrow(r)) return(NULL)
  py <- tryCatch(xrf_fluorescence_yield(elem, edge), error = function(e) NA_real_)
  pj <- tryCatch(xrf_absorption_jump_ratio(elem, edge), error = function(e) NA_real_)
  # CONVENTION: the package's xrf_absorption_jump_ratio returns the jump FACTOR (r-1)/r (the fraction of the
  # photoabsorption captured by the subshell), whereas Elam stores the jump RATIO r itself. Compare like for
  # like by converting Elam to the factor.
  ejf <- (r$jump - 1) / r$jump
  data.frame(elem, edge,
             yield_pkg = round(py, 4), yield_elam = round(r$yield, 4), dYield = pct(py, r$yield),
             jumpfac_pkg = round(pj, 3), jumpfac_elam = round(ejf, 3), dJump = pct(pj, ejf),
             stringsAsFactors = FALSE)
}
targets <- list(c("Ca","K"), c("Ti","K"), c("Fe","K"), c("Zn","K"), c("Sr","K"), c("Mo","K"),
                c("Ba","L3"), c("Nd","L3"), c("W","L3"), c("Pb","L3"), c("U","L3"))
cat("== fluorescence yield & jump ratio: package (EADL97 / Veigele) vs Elam 2002 ==\n")
tab <- do.call(rbind, Filter(Negate(is.null), lapply(targets, function(x) cmp(x[1], x[2]))))
print(tab, row.names = FALSE)

# summary of the yield deltas
dy <- suppressWarnings(as.numeric(sub("%", "", tab$dYield)))
dj <- suppressWarnings(as.numeric(sub("%", "", tab$dJump)))
cat(sprintf("\nyield delta: median |%%| = %.1f%%, max |%%| = %.1f%%\n",
            median(abs(dy), na.rm = TRUE), max(abs(dy), na.rm = TRUE)))
cat(sprintf("jump  delta: median |%%| = %.1f%%, max |%%| = %.1f%%  (NA = no package jump ratio, e.g. M/L gaps)\n",
            median(abs(dj), na.rm = TRUE), max(abs(dj), na.rm = TRUE)))

# --- Coster-Kronig spot check (Pb L) ----------------------------------------------------------------
cat("\n== Coster-Kronig f_ij: package vs Elam (Pb L-subshells) ==\n")
for (tr in list(c("L1","L2","f12"), c("L1","L3","f13"), c("L2","L3","f23"))) {
  pk <- tryCatch(xrf_coster_kronig_probability("Pb", tr[1], tr[3]), error = function(e) NA_real_)
  er <- elam_ck[elam_ck$el == "Pb" & elam_ck$from == tr[1] & elam_ck$to == tr[2], "prob"]
  cat(sprintf("  Pb %s->%s (%s): pkg=%.3f  elam=%.3f  %s\n", tr[1], tr[2], tr[3],
              pk, if (length(er)) er[1] else NA_real_, pct(pk, if (length(er)) er[1] else NA_real_)))
}

# --- photoabsorption cross section: package "photoelectric" mass attenuation vs Elam's log-log spline ---
# Elam stores total photoabsorption as a natural cubic spline in (ln E[eV], ln mu[cm^2/g]) with the stored
# second derivatives (Numerical Recipes format). Evaluate it and compare to xrf_mass_attenuation(.,"photoelectric").
photo <- list(); cur <- NA_character_; inpho <- FALSE; acc <- NULL
for (s in trimws(ln)) {
  if (startsWith(s, "Element ")) cur <- sp(s)[2]
  if (inpho) {
    if (grepl("^(Edge|Scatter|EndElement|Lines|CK|Element|End)", s)) {
      if (length(acc)) photo[[cur]] <- do.call(rbind, acc); inpho <- FALSE
    } else if (nchar(s) && grepl("^[0-9.+-]", s)) {
      acc[[length(acc) + 1]] <- as.numeric(sp(s))
    }
  }
  if (startsWith(s, "Photo")) { inpho <- TRUE; acc <- list() }
}
splint <- function(xa, ya, y2a, x) {
  if (x <= xa[1]) return(ya[1]); if (x >= xa[length(xa)]) return(ya[length(ya)])
  klo <- max(which(xa <= x)); khi <- klo + 1; h <- xa[khi] - xa[klo]
  a <- (xa[khi] - x) / h; b <- (x - xa[klo]) / h
  a * ya[klo] + b * ya[khi] + ((a^3 - a) * y2a[klo] + (b^3 - b) * y2a[khi]) * h * h / 6
}
elam_photo <- function(el, EkeV) { m <- photo[[el]]; if (is.null(m)) return(NA_real_); exp(splint(m[, 1], m[, 2], m[, 3], log(EkeV * 1000))) }
cat("\n== photoabsorption mu (cm^2/g): package photoelectric vs Elam, at 1.5x & 3x the K/L3 edge ==\n")
pcmp <- function(elem, edge, mult) {
  ee <- elam[elam$el == elem & elam$edge == edge, "e_ev"] / 1000
  if (!length(ee)) return(NULL)
  E <- ee * mult
  pp <- tryCatch(xrf_mass_attenuation(elem, E, type = "photoelectric"), error = function(e) NA_real_)
  data.frame(elem, edge, E_keV = round(E, 2), mu_pkg = round(pp, 1), mu_elam = round(elam_photo(elem, E), 1),
             d = pct(pp, elam_photo(elem, E)), stringsAsFactors = FALSE)
}
pt <- do.call(rbind, Filter(Negate(is.null), c(
  lapply(list(c("Fe","K"), c("Sr","K"), c("Mo","K")), function(x) pcmp(x[1], x[2], 1.5)),
  lapply(list(c("Fe","K"), c("Sr","K"), c("Mo","K")), function(x) pcmp(x[1], x[2], 3)),
  lapply(list(c("Pb","L3"), c("W","L3")),            function(x) pcmp(x[1], x[2], 1.5)))))
print(pt, row.names = FALSE)
dpp <- suppressWarnings(as.numeric(sub("%", "", pt$d)))
cat(sprintf("photoabsorption delta: median |%%| = %.1f%%, max |%%| = %.1f%%\n",
            median(abs(dpp), na.rm = TRUE), max(abs(dpp), na.rm = TRUE)))
