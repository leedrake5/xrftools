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

# ppm by mass, upper continental crust; ~0 for synthetic / short-lived nuclides (they get an ~infinite bar).
.xrf_crustal_abundance_ppm <- c(
  O=461000, Si=282000, Al=82300, Fe=56300, Ca=41500, Na=23600, Mg=23300, K=20900, Ti=5650, H=1400,
  P=1050, Mn=950, F=585, Ba=425, Sr=370, S=350, C=200, Zr=165, Cl=145, V=120, Cr=102, Ni=84, Zn=70,
  Cu=60, Ce=66.5, Nd=41.5, La=39, Y=33, Co=25, Sc=22, Li=20, Nb=20, Ga=19, Pb=14, B=10, Th=9.6, Pr=9.2,
  Sm=7.05, Gd=6.2, Dy=5.2, Er=3.5, Yb=3.2, Hf=3, Cs=3, Be=2.8, Sn=2.3, U=2.7, Br=2.4, Ta=2, Eu=2,
  As=1.8, Ge=1.5, Ho=1.3, W=1.25, Mo=1.2, Tb=1.2, Tl=0.85, Lu=0.8, Tm=0.52, I=0.45, In=0.25, Sb=0.2,
  Cd=0.15, Hg=0.085, Ag=0.075, Se=0.05, Ar=1.2, Pd=0.015, Bi=0.009, Os=0.0015, Pt=0.005, Au=0.004,
  Te=0.001, Ru=0.001, Rh=0.001, Ir=0.001, Re=0.0007, Kr=1e-4, Xe=3e-5,
  Tc=1e-9, Pm=1e-9, Po=2e-10, At=1e-12, Rn=4e-13, Ra=9e-7, Ac=5e-10, Pa=1.4e-6, Fr=1e-18, Np=1e-12, Pu=1e-12)

# Per-column abundance penalty FACTOR for a set of design-matrix columns: ref_ppm/abundance (capped only to
# avoid numerical blow-up on ~0-abundance nuclides). Non-element columns (scatter_*, sum_*, off-table) get 0.
# On its own this over-penalises rare-but-well-determined traces, so it is multiplied by a per-column
# COLLINEARITY weight (see .xrf_template_collinearity) before use: the product `pf_j * collin_j` penalises an
# element only where it is BOTH rare AND explained by another template. A rare isolated peak (Nb, Y, Th, U)
# has collin ~ 0 -> ~no penalty; a rare template overlapping an abundant one (Tb on Fe, Ru on the scatter
# hump) has collin ~ 1 and large pf -> crushed, with the abundant partner (small pf) keeping the intensity.
.xrf_abundance_penalty_factor <- function(elements, ref_ppm = 100, cap = 1e4) {
  ab <- suppressWarnings(as.numeric(.xrf_crustal_abundance_ppm[elements]))
  pf <- ref_ppm / ab
  pf[!is.finite(pf)] <- 0
  pmin(pf, cap)
}

# Per-column collinearity weight in [0,1]: the maximum cosine overlap of each design column with any OTHER
# column. Disjoint peaks -> 0 (dot product of non-overlapping supports); overlapping templates -> ~1. Used to
# confine the abundance ridge to columns the data cannot separate. Cheap: one Gram matrix (p ~ tens).
.xrf_template_collinearity <- function(X) {
  g <- crossprod(X)                      # p x p Gram (columns are non-negative templates)
  nrm <- sqrt(diag(g)); nrm[nrm <= 0] <- Inf
  cs <- g / outer(nrm, nrm)              # cosine similarity
  diag(cs) <- 0
  w <- apply(cs, 2, max)
  w[!is.finite(w)] <- 0
  pmin(pmax(w, 0), 1)
}
