# Parse the xraylib FF.dat / SF.dat tables (see modern_sources/xraylib_DATA_NOTE.md) into
# x_ray_form_factors (atomic form factor F(q, Z), Rayleigh) and x_ray_incoherent_functions
# (incoherent scattering function S(q, Z), Compton). q in 1/angstrom, = E[keV] sin(theta/2)/12.39842.
library(tibble)
library(dplyr)

all_elements <- c(
  "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
  "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb",
  "Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm",
  "Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
  "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm")

parse_xraylib_grid <- function(path, value_name) {
  txt <- readLines(path)
  i <- 1L; z <- 0L; out <- list()
  while (i <= length(txt)) {
    n <- suppressWarnings(as.integer(trimws(txt[i])))
    stopifnot(is.finite(n), n > 0)
    z <- z + 1L
    block <- read.table(text = txt[(i + 1L):(i + n)], col.names = c("q", "value", "spline2"))
    out[[z]] <- tibble(element = all_elements[z], z = z,
                       q_inv_angstrom = block$q, !!value_name := block$value)
    i <- i + n + 1L
  }
  bind_rows(out)
}

x_ray_form_factors <- parse_xraylib_grid("data-raw/modern_sources/FF.dat", "form_factor")
x_ray_incoherent_functions <- parse_xraylib_grid("data-raw/modern_sources/SF.dat", "incoherent_function")

# sanity: F(0, Z) = Z; S(0, Z) = 0; S saturates toward Z
ff0 <- x_ray_form_factors %>% group_by(z) %>% slice_min(q_inv_angstrom, n = 1)
stopifnot(max(abs(ff0$form_factor - ff0$z)) < 0.51)
sf0 <- x_ray_incoherent_functions %>% group_by(z) %>% slice_min(q_inv_angstrom, n = 1)
stopifnot(all(sf0$incoherent_function < 0.01))
sfmax <- x_ray_incoherent_functions %>% group_by(z) %>% slice_max(q_inv_angstrom, n = 1)
stopifnot(all(abs(sfmax$incoherent_function - sfmax$z) / sfmax$z < 0.05))

cat(sprintf("form factors: %d elements, %d rows; incoherent: %d elements, %d rows\n",
            max(x_ray_form_factors$z), nrow(x_ray_form_factors),
            max(x_ray_incoherent_functions$z), nrow(x_ray_incoherent_functions)))
usethis::use_data(x_ray_form_factors, overwrite = TRUE)
usethis::use_data(x_ray_incoherent_functions, overwrite = TRUE)
