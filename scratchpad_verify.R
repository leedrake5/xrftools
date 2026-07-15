# load package via devtools if possible, else source needed files
suppressMessages({
  ok <- requireNamespace("devtools", quietly=TRUE)
})
if (ok) {
  suppressMessages(devtools::load_all(".", quiet=TRUE))
} else {
  # fall back: source everything
  for (f in list.files("R", full.names=TRUE)) try(source(f), silent=TRUE)
}

comp_sil  <- xrftools:::.xrf_resolve_matrix("silicate")
comp_carb <- xrftools:::.xrf_resolve_matrix("carbonate")
mu <- xrftools:::.xrf_matrix_total_mu

E0 <- 40
sinp <- sin(45*pi/180)
lines <- c(K=3.314, Ca=3.692, Fe=6.404, Sr=14.165)

tau_eff <- function(comp, Ei) 1/(mu(comp,E0)/sinp + mu(comp,Ei)/sinp)

cat(sprintf("%-4s %10s %10s %14s %14s\n","el","tau_sil","tau_carb","carb/sil","sil/carb"))
for (nm in names(lines)) {
  Ei <- lines[[nm]]
  ts <- tau_eff(comp_sil, Ei)
  tc <- tau_eff(comp_carb, Ei)
  cat(sprintf("%-4s %10.4g %10.4g %14.3f %14.3f\n", nm, ts, tc, tc/ts, ts/tc))
}
cat("\nrecovered/true = tau_eff(unknown=carb)/tau_eff(standard=sil) = carb/sil column\n")
