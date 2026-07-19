# Crustal abundance of the elements and the abundance-informed deconvolution prior.
#
# The "fit everything" least-squares deconvolution is ill-posed wherever two element templates overlap
# (e.g. Tb Lalpha+Lbeta on Fe Kalpha+Kbeta, or Ba Lalpha on Ti Kalpha): the data does not determine how the
# shared intensity splits, and an unconstrained NNLS partitions it arbitrarily -- often loading a physically
# implausible, rare element. An abundance-informed Tikhonov (ridge) prior resolves this: each element column
# is shrunk toward zero with a penalty inversely proportional to its crustal abundance, so along a degenerate
# direction the shared intensity flows to the MORE abundant element (the split approaches the abundance ratio,
# since for a flat likelihood the ridge minimises sum lambda_j beta_j^2 -> beta_a/beta_b = lambda_b/lambda_a).
# Well-determined peaks are barely affected (the penalty only bites where the data is uninformative). The
# prior is OPT-IN (strength 0 by default) so it never changes legacy behaviour.

# Approximate continental-crust abundance (ppm by mass); ~0 for synthetic / short-lived nuclides (they get
# an ~infinite bar). These are order-of-magnitude proxies for the abundance PRIOR only -- a mix of upper-
# and bulk-crust reference values, not an internally consistent single reference, and the exact major-
# element values do not matter (their penalty factor -> 0). Only the RELATIVE ordering of the trace values
# drives the prior. (As raised to ~4.8, the upper-crust value, so this common trace is not over-penalised.)
.xrf_crustal_abundance_ppm <- c(
  O=461000, Si=282000, Al=82300, Fe=56300, Ca=41500, Na=23600, Mg=23300, K=20900, Ti=5650, H=1400,
  P=1050, Mn=950, F=585, Ba=425, Sr=370, S=350, C=200, Zr=165, Cl=145, V=120, Cr=102, Rb=90, Ni=84, Zn=70,
  Cu=60, Ce=66.5, Nd=41.5, La=39, Y=33, Co=25, Sc=22, Li=20, Nb=20, Ga=19, Pb=14, B=10, Th=9.6, Pr=9.2,
  Sm=7.05, Gd=6.2, Dy=5.2, As=4.8, Er=3.5, Yb=3.2, Hf=3, Cs=3, Be=2.8, Sn=2.3, U=2.7, Br=2.4, Ta=2, Eu=2,
  Ge=1.5, Ho=1.3, W=1.25, Mo=1.2, Tb=1.2, Tl=0.85, Lu=0.8, Tm=0.52, I=0.45, In=0.25, Sb=0.2,
  Cd=0.15, Hg=0.085, Ag=0.075, Se=0.05, Ar=1.2, Pd=0.015, Bi=0.009, Os=0.0015, Pt=0.005, Au=0.004,
  Te=0.001, Ru=0.001, Rh=0.001, Ir=0.001, Re=0.0007, Kr=1e-4, Xe=3e-5,
  Tc=1e-9, Pm=1e-9, Po=2e-10, At=1e-12, Rn=4e-13, Ra=9e-7, Ac=5e-10, Pa=1.4e-6, Fr=1e-18, Np=1e-12, Pu=1e-12)

# Per-column abundance penalty FACTOR for a set of design-matrix columns: ref_ppm/abundance (capped only to
# avoid numerical blow-up on ~0-abundance nuclides). Non-element columns (scatter_*, sum_*, off-table) get 0.
# On its own this over-penalises rare-but-well-determined traces, so it is multiplied by a per-column
# IDENTIFIABILITY weight (see .xrf_template_vif_weight) before use: the product `pf_j * w_j` penalises an
# element only where it is BOTH rare AND not separately estimable from the data. A rare peak the fit can
# resolve (VIF ~ 1, e.g. an isolated Nb/Y/Th/U line, OR two nearby-but-separable peaks) gets w ~ 0 -> ~no
# penalty; a rare template genuinely confounded with another (VIF >> 1, e.g. Tb aliased onto Fe, Ru on the
# scatter hump) gets w ~ 1 and large pf -> crushed, with the abundant partner (small pf) keeping the intensity.
.xrf_abundance_penalty_factor <- function(elements, ref_ppm = 100, cap = 1e4, abundance = NULL) {
  ab <- suppressWarnings(as.numeric(.xrf_crustal_abundance_ppm[elements]))
  # Optional user matrix prior: a named ppm vector overrides the crustal default for the named elements
  # (e.g. a steel matrix c(Fe=7e5, Cr=1.8e5, Ni=8e4) so its known majors are not treated as rare). Elements
  # absent from `abundance` keep their crustal value; non-positive / non-finite overrides are ignored.
  if (!is.null(abundance) && length(abundance) && !is.null(names(abundance))) {
    user <- suppressWarnings(as.numeric(abundance[elements]))
    hit <- is.finite(user) & user > 0
    ab[hit] <- user[hit]
  }
  pf <- ref_ppm / ab
  pf[!is.finite(pf)] <- 0
  pmin(pf, cap)
}

# Per-column identifiability weight in [0,1), from the variance-inflation factor of the (weighted) design
# actually being solved: VIF_j = G_jj * (G^-1)_jj with G = X'X (VIF_j = 1/(1 - R_j^2), R_j^2 = fraction of
# column j explained by the others). w_j = max(0, 1 - vif_min/VIF_j): a column the fit can estimate
# (VIF <= vif_min, e.g. two peaks >~2 sigma apart) gets EXACTLY 0 -> no penalty regardless of raw spectral
# overlap or how rare the element is; a near-collinear / aliased column (VIF >> vif_min) gets w -> 1. The
# vif_min floor is what makes an estimable-but-rare trace safe: without it, a tiny weight times a large
# 1/abundance factor still shrinks the trace. This targets true NON-identifiability, unlike a raw cosine
# overlap (which fires on any two nearby-but-separable peaks). A tiny ridge stabilises the inverse for
# exactly-aliased columns (VIF -> huge, w -> 1). Cheap: one p x p solve (p ~ tens). Pass the weighted design
# Xw so the gate reflects the problem the solver sees. vif_min = 2 penalises only columns >50% explained by
# the rest (the standard "start worrying" VIF threshold), so estimable traces are left untouched.
.xrf_template_vif_weight <- function(X, vif_min = 2) {
  G <- crossprod(X)                      # p x p Gram of the (weighted) design
  d <- diag(G); p <- ncol(G)
  ok <- is.finite(d) & d > 0
  if (!any(ok)) return(rep(0, p))
  Ginv <- tryCatch(solve(G + diag(1e-8 * mean(d[ok]), p, p)), error = function(e) NULL)
  if (is.null(Ginv)) return(rep(0, p))
  vif <- d * diag(Ginv)
  vif[!is.finite(vif) | vif < 1] <- 1
  w <- 1 - vif_min / vif
  w[!ok] <- 0
  pmin(pmax(w, 0), 1)
}

# Per-element "has a spectrally isolated (clean) line" flag, used to exempt such elements from the abundance
# penalty. The per-element VIF gate reads an element whose DOMINANT lines overlap another template as
# non-identifiable; because the element's design column is the sum of its lines and is dominated by that
# strong line, even a genuinely CLEAN but WEAKER line barely lowers R^2, so VIF stays high and the large
# 1/abundance factor crushes the trace -- although the clean line ALONE identifies the amplitude (finite
# variance, not a degenerate direction). Here, for each element column, we ask whether any of its real lines
# is isolated: element j supplies more than `iso_thresh` of the TOTAL element-template response at that
# line's channel. Uses the fixed template shapes X (compute once). Exempting these is phantom-safe: an
# isolated line means the direction is NOT degenerate, so the data -- not the prior -- sets the amplitude
# (a true phantom has no signal there and NNLS zeros it regardless). Returns a logical vector over columns.
.xrf_clean_line_exempt <- function(peaks_tmpl, X, elements, is_element_col, energy_kev,
                                   iso_thresh = 0.7, min_line = 0.02) {
  p <- ncol(X)
  exempt <- rep(FALSE, p)
  el_cols <- which(is_element_col)
  if (!length(el_cols)) return(exempt)
  tot <- rowSums(X[, el_cols, drop = FALSE])          # total ELEMENT-template response per channel
  ne <- length(energy_kev)
  for (jj in el_cols) {
    lines <- peaks_tmpl[peaks_tmpl$element == elements[jj] &
                          peaks_tmpl$relative_peak_intensity >= min_line, , drop = FALSE]
    if (!nrow(lines)) next
    idx <- findInterval(lines$energy_kev, energy_kev)
    idx <- pmax(pmin(idx, ne), 1L)
    frac <- X[cbind(idx, jj)] / pmax(tot[idx], .Machine$double.eps)
    if (any(is.finite(frac) & frac >= iso_thresh)) exempt[jj] <- TRUE
  }
  exempt
}

# L1 / elastic-net solver for the deconvolution, used when prior = "lasso" / "elastic_net". Minimises
#   ||W^{1/2}(y - X b)||^2 + sum_j base_j [ (1-alpha) b_j^2 + alpha |b_j| ],   base_j = strength * scale * pf_j
# with pf_j the SAME per-column gated abundance penalty the ridge uses (0 for protected / clean-line /
# non-element columns, so they are never shrunk) and `strength` = abundance_prior on the SAME scale as the
# ridge. The |b_j| (L1) term is handled by iteratively reweighted L2 (majorize-minimize): at each pass the
# effective ridge penalty on column j is base_j[(1-alpha) + alpha/(2(|b_j|+eps))], so a shrinking column's
# penalty grows without bound and it is driven to (snapped to) exactly zero. Reuses the augmented NNLS / QR
# solve of the ridge path, so non-negativity and the strength convention carry over unchanged -- no extra
# package dependency. NOTE: under non-negativity NNLS already zeros clearly-absent elements, so this differs
# from the ridge mainly in the unconstrained fit and on confounded rare columns.
.xrf_irl2_fit <- function(X, y, rw, scale, pen_factor, alpha, strength, nonneg, p, iters = 12L) {
  Xrw <- X * rw; yrw <- y * rw
  base <- strength * scale * pmax(pen_factor, 0)
  solve_aug <- function(lam) {
    cf <- if (nonneg) nnls::nnls(rbind(Xrw, diag(sqrt(lam), p, p)), c(yrw, rep(0, p)))$x
          else qr.coef(qr(rbind(Xrw, diag(sqrt(lam), p, p))), c(yrw, rep(0, p)))
    cf[!is.finite(cf)] <- 0
    cf
  }
  b <- solve_aug(base)                                   # ridge-penalised initialisation
  for (t in seq_len(iters)) {
    eps <- 1e-6 * max(abs(b), 1)
    b <- solve_aug(base * ((1 - alpha) + alpha / (2 * (abs(b) + eps))))
  }
  b[abs(b) < 1e-6 * max(abs(b), 1)] <- 0                 # snap the MM residual to an exact zero
  # effective L2 penalty at the solution (the local curvature the Laplace SE conditions on): for active
  # columns base_j[(1-alpha)+alpha/(2|b_j|)]; approximate for L1 (the posterior is not Gaussian at |b|=0).
  eps <- 1e-6 * max(abs(b), 1)
  pen_final <- base * ((1 - alpha) + alpha / (2 * (abs(b) + eps)))
  list(coef = b, active = which(abs(b) > 0), aliased = integer(0), penalty = pen_final)
}
