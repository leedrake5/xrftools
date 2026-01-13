
# Package environment for cached lookup tables
.xrf_cache <- new.env(parent = emptyenv())

# Lazy initialization of lookup tables using named vectors
# Named vector lookup with match() is faster than left_join for repeated lookups
.init_physics_cache <- function() {
  if (is.null(.xrf_cache$initialized)) {
    # Build named vector lookup for emission probabilities
    ep <- x_ray_emission_probabilities
    ep_keys <- paste(ep$element, ep$trans, sep = ":")
    .xrf_cache$emission_prob_vec <- setNames(ep$emission_probability, ep_keys)

    # Build named vector lookup for fluorescence yields
    fy <- x_ray_fluorescence_yields
    fy_keys <- paste(fy$element, fy$shell, sep = ":")
    .xrf_cache$fluorescence_yield_vec <- setNames(fy$fluorescence_yield, fy_keys)

    # Build named vector lookup for Coster-Kronig probabilities
    ck <- x_ray_coster_kronig_probabilities
    ck_keys <- paste(ck$element, ck$shell, ck$coster_kronig_trans, sep = ":")
    .xrf_cache$coster_kronig_vec <- setNames(ck$coster_kronig_prob, ck_keys)

    .xrf_cache$initialized <- TRUE
  }
  invisible()
}
