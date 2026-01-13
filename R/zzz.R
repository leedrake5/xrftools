
# Package environment for cached lookup tables
.xrf_cache <- new.env(parent = emptyenv())

# Lazy initialization of lookup tables using fast hash-based environments
# This avoids data.table namespace issues while providing O(1) lookups
.init_physics_cache <- function() {
  if (is.null(.xrf_cache$initialized)) {
    # Build hash-based lookup for emission probabilities (element:trans -> probability)
    ep <- x_ray_emission_probabilities
    .xrf_cache$emission_prob_hash <- new.env(hash = TRUE, parent = emptyenv())
    for (i in seq_len(nrow(ep))) {
      key <- paste(ep$element[i], ep$trans[i], sep = ":")
      .xrf_cache$emission_prob_hash[[key]] <- ep$emission_probability[i]
    }

    # Build hash-based lookup for fluorescence yields (element:shell -> yield)
    fy <- x_ray_fluorescence_yields
    .xrf_cache$fluorescence_yield_hash <- new.env(hash = TRUE, parent = emptyenv())
    for (i in seq_len(nrow(fy))) {
      key <- paste(fy$element[i], fy$shell[i], sep = ":")
      .xrf_cache$fluorescence_yield_hash[[key]] <- fy$fluorescence_yield[i]
    }

    # Build hash-based lookup for Coster-Kronig probabilities
    ck <- x_ray_coster_kronig_probabilities
    .xrf_cache$coster_kronig_hash <- new.env(hash = TRUE, parent = emptyenv())
    for (i in seq_len(nrow(ck))) {
      key <- paste(ck$element[i], ck$shell[i], ck$coster_kronig_trans[i], sep = ":")
      .xrf_cache$coster_kronig_hash[[key]] <- ck$coster_kronig_prob[i]
    }

    .xrf_cache$initialized <- TRUE
  }
  invisible()
}

# Vectorized hash lookup helper
.hash_lookup <- function(hash_env, keys, default = NA_real_) {
  vapply(keys, function(k) {
    val <- hash_env[[k]]
    if (is.null(val)) default else val
  }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
}
