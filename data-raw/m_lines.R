# Add M-shell emission lines (Malpha, Mbeta, Mgamma, Mzeta) for heavy elements.
#
# The K and L machinery already works for M shells once the data are present; this adds (a) the
# M-line ENERGIES to x_ray_xrf_energies and (b) the M-line EMISSION PROBABILITIES to
# x_ray_emission_probabilities. Both are built from the EADL97 files already vendored in data-raw/:
#   - EADL97_BindingEnergies.dat   -> subshell binding energies (line energy = E_initial - E_final)
#   - EADL97_MShellRadiativeRates.dat -> radiative emission probabilities (already omega-weighted,
#                                        like the K/L emission_probabilities: per-shell TOTAL == omega_M)
#
# These are the analytical lines for low-kV / SEM-EDS work on heavy elements (Au, Pt, Pb, U ...),
# where a <30 keV beam excites the M series (~2-3.5 keV) but not the L or K edges.

library(dplyr)
setwd(rprojroot::find_root(rprojroot::is_r_package))

all_elements <- c(
  "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
  "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y",
  "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce",
  "Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir",
  "Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm",
  "Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc",
  "Lv","Ts","Og"
)

# ---- binding energies: Z + subshells in EADL order --------------------------------------------
be_lines <- readLines("data-raw/EADL97_BindingEnergies.dat")
be_lines <- be_lines[!grepl("^#", be_lines) & grepl("[0-9]", be_lines)]
be <- do.call(rbind, lapply(strsplit(trimws(be_lines), "\\s+"), function(x) as.numeric(x)))
be <- as.data.frame(be)
subshells <- c("K","L1","L2","L3","M1","M2","M3","M4","M5",
               "N1","N2","N3","N4","N5","N6","N7","O1","O2","O3","O4","O5")
colnames(be)[1] <- "Z"
colnames(be)[2:(1 + length(subshells))] <- subshells
BE <- function(z, sh) be[[sh]][match(z, be$Z)]

# ---- M-shell radiative rates: one #S block per initial subshell (M1..M5) ----------------------
read_rate_block <- function(path, initial) {
  ln <- readLines(path)
  starts <- grep("^#S ", ln)
  hdr <- grep(paste0("^#S .* ", initial, " emission"), ln)
  stopifnot(length(hdr) == 1)
  cols_line <- ln[grep("^#L", ln[hdr:length(ln)])[1] + hdr - 1]
  cols <- strsplit(trimws(sub("^#L", "", cols_line)), "\\s+")[[1]]
  # data rows: after the #L line, until the next #S (or end)
  block_end <- c(starts[starts > hdr], length(ln) + 1)[1] - 1
  data_ln <- ln[(grep("^#L", ln)[grep("^#L", ln) > hdr][1] + 1):block_end]
  data_ln <- data_ln[!grepl("^#", data_ln) & grepl("[0-9]", data_ln)]
  m <- do.call(rbind, lapply(strsplit(trimws(data_ln), "\\s+"), as.numeric))
  df <- as.data.frame(m)
  colnames(df) <- cols
  df
}
rates <- lapply(c("M1","M2","M3","M4","M5"), function(s)
  read_rate_block("data-raw/EADL97_MShellRadiativeRates.dat", s))
names(rates) <- c("M1","M2","M3","M4","M5")
RATE <- function(z, initial, trans) {
  df <- rates[[initial]]
  if (!(trans %in% colnames(df))) return(NA_real_)
  df[[trans]][match(z, df$Z)]
}

# ---- main M lines: initial subshell, final subshell, transition + Siegbahn name ----------------
m_defs <- tibble::tribble(
  ~edge, ~final, ~trans,  ~trans_siegbahn,
  "M5",  "N7",   "M5N7",  "Malpha1",
  "M5",  "N6",   "M5N6",  "Malpha2",
  "M4",  "N6",   "M4N6",  "Mbeta",
  "M3",  "N5",   "M3N5",  "Mgamma",
  "M5",  "N3",   "M5N3",  "Mzeta1",
  "M4",  "N2",   "M4N2",  "Mzeta2"
)

Z_range <- 50:92   # Sn..U (M lines land ~0.4-3.5 keV; filtered below)

m_rows <- do.call(rbind, lapply(Z_range, function(z) {
  do.call(rbind, lapply(seq_len(nrow(m_defs)), function(i) {
    d <- m_defs[i, ]
    e_init <- BE(z, d$edge); e_fin <- BE(z, d$final)
    energy <- e_init - e_fin
    rate <- RATE(z, d$edge, d$trans)
    data.frame(element = all_elements[z], z = z, trans = d$trans, edge = d$edge,
               trans_siegbahn = d$trans_siegbahn, energy_kev = energy,
               edge_kev = e_init, emission_probability = rate,
               stringsAsFactors = FALSE)
  }))
}))

# keep physically sensible, measurable lines
m_rows <- m_rows %>%
  filter(is.finite(energy_kev), energy_kev > 0.3, energy_kev < 6,
         is.finite(emission_probability), emission_probability > 0,
         is.finite(edge_kev), edge_kev > 0)

cat("Built", nrow(m_rows), "M lines for Z", min(m_rows$z), "-", max(m_rows$z), "\n")
cat("Au Malpha1 energy (expect ~2.12):",
    round(m_rows$energy_kev[m_rows$element == "Au" & m_rows$trans == "M5N7"], 3), "\n")

# sanity: per-shell M5 TOTAL vs omega_M5 (confirms rates are omega-inclusive like K/L)
load("data/x_ray_fluorescence_yields.rda")
au_m5_sum <- sum(rates$M5[rates$M5$Z == 79, setdiff(colnames(rates$M5), c("Z","TOTAL"))], na.rm = TRUE)
au_m5_omega <- x_ray_fluorescence_yields$fluorescence_yield[
  x_ray_fluorescence_yields$element == "Au" & x_ray_fluorescence_yields$shell == "M5"]
cat("Au M5: sum(rates)=", round(au_m5_sum, 4), " omega_M5=", round(au_m5_omega, 4),
    " (ratio", round(au_m5_sum / au_m5_omega, 3), "=> omega-inclusive if ~1)\n")

# ---- append to x_ray_xrf_energies -------------------------------------------------------------
load("data/x_ray_xrf_energies.rda")
new_xrf <- m_rows %>%
  transmute(element, z, trans, edge, trans_siegbahn, energy_kev, edge_kev) %>%
  select(all_of(colnames(x_ray_xrf_energies)))
x_ray_xrf_energies <- x_ray_xrf_energies %>%
  filter(!(trans %in% unique(m_defs$trans))) %>%   # idempotent
  bind_rows(new_xrf) %>%
  arrange(z, trans_siegbahn)
usethis::use_data(x_ray_xrf_energies, overwrite = TRUE)

# ---- append to x_ray_emission_probabilities ---------------------------------------------------
load("data/x_ray_emission_probabilities.rda")
new_ep <- m_rows %>% transmute(element, z, trans, emission_probability)
x_ray_emission_probabilities <- x_ray_emission_probabilities %>%
  filter(!(trans %in% unique(m_defs$trans))) %>%   # idempotent
  bind_rows(new_ep)
usethis::use_data(x_ray_emission_probabilities, overwrite = TRUE)

cat("Done: x_ray_xrf_energies and x_ray_emission_probabilities updated.\n")
