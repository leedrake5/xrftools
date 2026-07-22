
#' Deconvolute a spectrum using gaussian least squares
#'
#' @inheritParams xrf_add_deconvolution_fun
#' @param .values,response Values to use for the deconvolution
#' @param energy_kev,.energy_kev Energy values corresponding to response/.values
#' @param energy_min_kev,energy_max_kev Constrain the energies that should be deconvoluted
#' @param peaks A data frame with (at least) columns "element" and "energy_kev". Optional columns
#'   "sigma" and "relative_peak_intensity" are used when present.
#' @param default_sigma The default standard deviation of the peaks (keV), used when neither a
#'   \code{sigma} column nor a detector model is supplied.
#' @param detector_type,fano,epsilon_ev,noise_fwhm_ev Detector resolution parameters. When any of
#'   these is supplied (and \code{peaks} has no \code{sigma} column) the per-line peak width is
#'   computed from the energy-dependent resolution model \link{xrf_detector_sigma_kev} instead of a
#'   constant \code{default_sigma}. \code{detector_type} is one of \link{xrf_detector_presets}
#'   (e.g. "SDD", "HPGe", "CdTe"); the scalars override the preset.
#' @param efficiency If TRUE, reweight each line by the energy-dependent detector efficiency
#'   (\link{xrf_detector_efficiency}); corrects line shape and high-energy sensitivity. Uses the
#'   \code{detector_type} geometry (override window/dead layer with \code{be_window_um} /
#'   \code{dead_layer_um}).
#' @param escape If TRUE, add detector escape-peak satellites (\link{xrf_escape_peaks}), tied to each
#'   parent line's amplitude. Important for Ge/CdTe at high energy.
#' @param be_window_um,dead_layer_um Detector window / dead-layer thickness (microns) for the
#'   efficiency model; \code{NULL} uses the \code{detector_type} preset.
#' @param active_thickness_um Active-volume thickness (microns) for the efficiency model; \code{NULL}
#'   uses the \code{detector_type} preset (SDD = 450 um, a typical modern silicon-drift detector).
#'   Governs high-energy transmission (thin Si becomes transparent to heavy-element K lines), so it
#'   only matters when \code{efficiency = TRUE}.
#' @param air_path_cm,atmosphere,window Measurement-path terms for the efficiency model: air/He path length
#'   (cm), atmosphere (\code{"Air"}/\code{"He"}/\code{"Vacuum"}) and an optional snout window; \code{NULL}/
#'   \code{"Air"} reproduce the historical no-path behaviour.
#' @param nonneg If TRUE (default) coefficients are constrained to be non-negative via
#'   \link[nnls]{nnls} (peak areas / concentrations cannot be negative). Set FALSE for the classic
#'   unconstrained least-squares behaviour.
#' @param weighting Counting-statistics weighting. "poisson" performs an iteratively-reweighted fit
#'   with weights 1/variance derived from the modelled counts (requires \code{counts}/\code{.counts}
#'   and \code{livetime}/\code{.livetime}); "none" is ordinary (unweighted) least squares; "auto"
#'   (default) uses "poisson" when counts and live time are available and "none" otherwise.
#' @param counts,.counts,livetime,.livetime Raw gross counts per channel and the acquisition live
#'   time (seconds), used for Poisson weighting and for the count-space goodness-of-fit. Optional.
#' @param tube,geometry Optional \link{xrf_tube} / \link{xrf_geometry} objects. When a \code{tube}
#'   is supplied, Rayleigh and Compton tube-scatter templates are added to the fit automatically
#'   (see \code{scatter} and \link{xrf_scatter_peaks}); \code{geometry} sets the Compton scatter
#'   angle. (Not otherwise consumed yet -- reserved for fundamental-parameters work.)
#' @param scatter Whether to include Rayleigh/Compton scatter templates. \code{NULL} (default)
#'   enables them when a \code{tube} is supplied; \code{TRUE}/\code{FALSE} force on/off.
#' @param compton_broadening Multiplier on the detector width for the Doppler-broadened Compton
#'   peak (passed to \link{xrf_scatter_peaks}).
#' @param scatter_continuum If TRUE (needs a \code{tube}), also add scattered-\emph{continuum} templates
#'   (Rayleigh + Compton of the tube bremsstrahlung, \link{xrf_scatter_continuum}) -- the broad scatter
#'   background / Compton hump that \link{xrf_scatter_peaks} (anode lines only) misses. Default FALSE.
#'   \strong{Working recipe:} fit against an UN-baselined response (\code{.values = cps}) together with a
#'   fitted \code{background} (6-12), so the continuum + background REPLACE the SNIP baseline. This recovers
#'   trace areas cleanly (r^2 = 1 on synthetics; SNIP-only over-estimates traces 100-300\% under a Compton
#'   hump). Enabling it on the default baseline-subtracted \code{.values = cps - baseline}, or without a
#'   fitted \code{background}, double-counts the baseline and biases real peaks -- a warning then fires.
#' @param scatter_scatterer Representative sample scatterer weighting both the anode-line scatter
#'   templates (\link{xrf_scatter_peaks}) and the scattered-continuum shape
#'   (\link{xrf_scatter_continuum}); default "Si".
#' @param background Number of broad, non-negative Gaussian background basis functions to fit jointly with
#'   the peaks (default 0 = off). When \eqn{\ge 1}, the deconvolution models a smooth background \emph{as
#'   part of the fit}, so you can feed it an \strong{un-baselined} response (\code{.values = cps}) instead
#'   of pre-subtracting a SNIP baseline. This is the intended companion to \code{scatter_continuum}: with a
#'   fitted background there is no subtracted baseline for the broad scatter templates to double-count, so
#'   real trace areas are preserved. Excluded from quantification. Typical values ~6-12.
#' @param background_smooth Roughness penalty (dimensionless, default 0 = off) on the fitted
#'   \code{background} basis: penalises the second differences of the bump coefficients (a P-spline),
#'   so a \emph{generous} background (e.g. \code{background = 15-25}) stays smooth instead of chasing
#'   peak residuals -- more baseline flexibility without peak-stealing. Strength ~1 is a gentle
#'   default that leaves well-resolved smooth features essentially untouched; tens-to-hundreds give a
#'   stiff spline. Applied with \code{prior = "ridge"} (the default) only.
#' @param tail,step,beta Hypermet line-shape parameters (passed to \link{xrf_lineshape}): relative
#'   low-energy tail amplitude, shelf/step amplitude, and tail decay length. Defaults (0, 0, NULL)
#'   give a pure Gaussian. The reported \code{peak_area} includes the tail area when \code{tail > 0}.
#'   These are \emph{global} defaults: per-line \code{tail}/\code{step}/\code{beta} \emph{columns} on
#'   \code{peaks} override them row-wise (the reported \code{peak_area} then uses each element's
#'   primary-line values).
#' @param tailing If TRUE, fill the per-line tail/step/beta of every line that has no explicit value
#'   from the energy-dependent incomplete-charge-collection model \link{xrf_detector_tailing} -- Si
#'   tailing grows toward low energy, Ge/CdTe hole-tailing toward high energy. The shelves of strong
#'   peaks (above all the big scatter features) are a real part of the measured continuum, so this
#'   mainly improves BASELINE fidelity. Default FALSE (pure Gaussians / global scalars).
#' @param compton_tail Relative amplitude of a low-energy exponential tail on the anode-line Compton
#'   scatter templates (passed to \link{xrf_scatter_peaks}; \code{compton_shape = "gaussian"} mode
#'   only). Typical values ~0.1-0.5; default 0 (off).
#' @param compton_shape Shape of the anode-line Compton scatter features (passed to
#'   \link{xrf_scatter_peaks}): "profile" (default) uses the exact impulse-approximation Doppler
#'   shape from the Biggs Compton profiles -- the physical asymmetric hump -- while "gaussian" is the
#'   legacy broadened Gaussian. Compton normalizations calibrated in one mode must not be mixed with
#'   the other (the fitted template's \code{peak_area} is its strongest single component).
#' @param cache_templates If TRUE (default) reuse the design matrix across successive calls that
#'   share the same energy grid, peaks and line-shape (e.g. a batch from one instrument), skipping
#'   the Gaussian re-evaluation.
#' @param refine_calibration If TRUE, refine a 2-parameter energy calibration (zero + gain) of the
#'   template centroids by variable projection before the final fit -- robust to small gain/zero
#'   drift that would otherwise bias amplitudes. Disables template caching for that call.
#' @param sum_peaks If TRUE, add coincidence (pile-up) sum-peak templates at E_i + E_j for the
#'   strongest lines, so pile-up is not misassigned to real elements.
#' @param count_rate Optional total count rate (cps); when supplied and 0, the sum-peak step is
#'   skipped. Reserved for scaling the pile-up contribution.
#' @param pileup_tau Optional detector pulse-pair resolution time (seconds). When supplied with
#'   \code{counts}/\code{.counts} and \code{livetime}/\code{.livetime}, a physically-constrained
#'   two-pass pile-up correction is applied: coincidence peaks with areas fixed to
#'   \eqn{2\tau A_i A_j/\Delta E} (\eqn{\Delta E} the channel width; live time cancels in rate space)
#'   are subtracted and the elements refitted, and the reported \code{height}/\code{peak_area} (and
#'   SEs/covariances) of all non-pileup templates are then multiplied by \eqn{e^{2\tau R_{tot}}}
#'   (\eqn{R_{tot}} the gross count rate) to restore the counts each line lost \emph{to} pile-up --
#'   a uniform factor, so ratios are unchanged but absolute quantities
#'   (\link{xrf_observed_mass}, calibrated results) are corrected. \code{response}/\code{components}
#'   stay in measured (piled-up) space. This supersedes the free-amplitude \code{sum_peaks} model.
#' @param use_qr If TRUE (default) unconstrained fits (\code{nonneg = FALSE}) use fast QR
#'   decomposition; set FALSE for fork-safe \code{lm()} under multicore parallelism. Ignored when
#'   \code{nonneg = TRUE}.
#' @param abundance_prior Strength of an optional crustal-abundance ridge prior (default 0 = off, an
#'   identical fit). When > 0, each element template is shrunk with a penalty proportional to
#'   \code{abundance_ref_ppm}/(crustal abundance), \strong{gated} by a variance-inflation-factor
#'   identifiability weight so only columns that are \emph{both} rare \emph{and} not separately estimable
#'   are penalised (\link{x_ray_xrf_energies}; see \code{R/abundance.R}). This resolves degenerate
#'   overlaps (e.g. a rare element's lines coinciding with an abundant one's) toward the more abundant
#'   element. Typical useful strengths are ~0.1-0.5. \strong{Caveat:} it is a biasing prior -- a
#'   genuinely-present but truly-unresolvable trace that overlaps a common line can be shrunk low, so leave
#'   it at 0 for the most defensible trace quantification and raise it only to suppress phantom peaks.
#' @param abundance_ref_ppm Pivot abundance (ppm) at which the abundance penalty factor equals 1 (default
#'   100); elements rarer than this get a larger penalty factor, more abundant ones a smaller one. Only
#'   used when \code{abundance_prior > 0}.
#' @param abundance_protect Character vector of element symbols to \strong{exempt} from the abundance prior
#'   (never penalised). Use it to protect elements you know are present -- in particular heavy elements read
#'   via L-lines that overlap an abundant K-line or the scatter peak (Pb L\eqn{\alpha} under As K\eqn{\alpha},
#'   U L\eqn{\gamma} under the Rayleigh line), which the VIF gate would otherwise read as unidentifiable and
#'   crush. Protecting them makes the tiebreaker explicit: at each coincidence the prior penalises only the
#'   rare unprotected competitor, so the shared intensity flows to the protected element while phantoms are
#'   still removed. Default none. Only used when \code{abundance_prior > 0}.
#' @param abundance_ppm Optional named numeric vector of per-element abundances (ppm by mass) that
#'   \strong{override} the built-in crustal reference for the named elements -- e.g. a known matrix
#'   \code{c(Fe = 7e5, Cr = 1.8e5, Ni = 8e4)} so its major/known elements are not treated as "rare" and
#'   penalised. Elements not named keep their crustal default. Only used when \code{abundance_prior > 0}.
#' @param prior Form of the abundance penalty: \code{"ridge"} (default, L2 Tikhonov, legacy behaviour),
#'   \code{"lasso"} (L1) or \code{"elastic_net"} (an L1/L2 mix, \code{enet_alpha}). All three reuse the same
#'   VIF / clean-line / \code{abundance_protect} gating and the same \code{abundance_prior} strength; the L1
#'   forms are realised by iteratively reweighted L2 (a majorize-minimize scheme, no extra dependency) and
#'   snap penalised-and-shrunken columns to \strong{exactly} zero. \strong{Note:} under the default
#'   non-negativity, NNLS already zeros clearly-absent elements, so the L1 forms differ from the ridge mainly
#'   (a) in the unconstrained \code{nonneg = FALSE} fit (where they zero what least squares would leave
#'   non-zero) and (b) by removing confounded rare columns more aggressively.
#' @param enet_alpha Elastic-net mixing in [0,1] for \code{prior = "elastic_net"} (1 = lasso, 0 = ridge-like);
#'   default 0.5. Ignored otherwise.
#'
#' @return A modified version of .spectra
#' @export
#'
xrf_add_deconvolution_gls <- function(.spectra, .energy_kev = .data$.spectra$energy_kev,
                                      .values = .data$.spectra$cps - .data$.spectra$baseline,
                                      energy_min_kev = -Inf, energy_max_kev = Inf,
                                      peaks = xrf_energies(), default_sigma = 0.07,
                                      detector_type = NULL, fano = NULL, epsilon_ev = NULL,
                                      noise_fwhm_ev = NULL,
                                      efficiency = FALSE, escape = FALSE,
                                      be_window_um = NULL, dead_layer_um = NULL,
                                      active_thickness_um = NULL,
                                      air_path_cm = NULL, atmosphere = "Air", window = NULL,
                                      nonneg = TRUE, weighting = c("auto", "poisson", "none"),
                                      .counts = NULL, .livetime = NULL,
                                      tube = NULL, geometry = NULL,
                                      scatter = NULL, compton_broadening = 2,
                                      scatter_continuum = FALSE, scatter_scatterer = "Si",
                                      background = 0L, background_smooth = 0,
                                      tail = 0, step = 0, beta = NULL, tailing = FALSE,
                                      compton_tail = 0, compton_shape = c("profile", "gaussian"),
                                      cache_templates = TRUE, refine_calibration = FALSE,
                                      sum_peaks = FALSE, count_rate = NULL, pileup_tau = NULL,
                                      use_qr = TRUE, abundance_prior = 0, abundance_ref_ppm = 100,
                                      abundance_protect = character(0), abundance_ppm = NULL,
                                      prior = c("ridge", "lasso", "elastic_net"), enet_alpha = 0.5,
                                      .cores = 1L,
                                      .env = parent.frame()) {
  .energy_kev <- enquo(.energy_kev)
  .values <- enquo(.values)
  peaks <- enquo(peaks)
  default_sigma <- enquo(default_sigma)
  energy_min_kev <- enquo(energy_min_kev)
  energy_max_kev <- enquo(energy_max_kev)
  detector_type <- enquo(detector_type)
  fano <- enquo(fano)
  epsilon_ev <- enquo(epsilon_ev)
  noise_fwhm_ev <- enquo(noise_fwhm_ev)
  efficiency <- enquo(efficiency)
  escape <- enquo(escape)
  be_window_um <- enquo(be_window_um)
  dead_layer_um <- enquo(dead_layer_um)
  active_thickness_um <- enquo(active_thickness_um)
  air_path_cm <- enquo(air_path_cm)
  atmosphere <- enquo(atmosphere)
  window <- enquo(window)
  nonneg <- enquo(nonneg)
  weighting <- enquo(weighting)
  .counts <- enquo(.counts)
  .livetime <- enquo(.livetime)
  tube <- enquo(tube)
  geometry <- enquo(geometry)
  scatter <- enquo(scatter)
  compton_broadening <- enquo(compton_broadening)
  scatter_continuum <- enquo(scatter_continuum)
  scatter_scatterer <- enquo(scatter_scatterer)
  background <- enquo(background)
  background_smooth <- enquo(background_smooth)
  tail <- enquo(tail)
  step <- enquo(step)
  beta <- enquo(beta)
  tailing <- enquo(tailing)
  compton_tail <- enquo(compton_tail)
  compton_shape <- enquo(compton_shape)
  cache_templates <- enquo(cache_templates)
  refine_calibration <- enquo(refine_calibration)
  sum_peaks <- enquo(sum_peaks)
  count_rate <- enquo(count_rate)
  pileup_tau <- enquo(pileup_tau)
  use_qr <- enquo(use_qr)
  abundance_prior <- enquo(abundance_prior)
  abundance_ref_ppm <- enquo(abundance_ref_ppm)
  abundance_protect <- enquo(abundance_protect)
  abundance_ppm <- enquo(abundance_ppm)
  prior <- enquo(prior)
  enet_alpha <- enquo(enet_alpha)

  xrf_add_deconvolution_fun(
    .spectra,
    xrf_deconvolute_gaussian_least_squares,
    !!.energy_kev,
    !!.values,
    peaks = !!peaks,
    default_sigma = !!default_sigma,
    energy_min_kev = !!energy_min_kev,
    energy_max_kev = !!energy_max_kev,
    detector_type = !!detector_type,
    fano = !!fano,
    epsilon_ev = !!epsilon_ev,
    noise_fwhm_ev = !!noise_fwhm_ev,
    efficiency = !!efficiency,
    escape = !!escape,
    be_window_um = !!be_window_um,
    dead_layer_um = !!dead_layer_um,
    active_thickness_um = !!active_thickness_um,
    air_path_cm = !!air_path_cm, atmosphere = !!atmosphere, window = !!window,
    nonneg = !!nonneg,
    weighting = !!weighting,
    counts = !!.counts,
    livetime = !!.livetime,
    tube = !!tube,
    geometry = !!geometry,
    scatter = !!scatter,
    compton_broadening = !!compton_broadening,
    scatter_continuum = !!scatter_continuum,
    scatter_scatterer = !!scatter_scatterer,
    background = !!background,
    background_smooth = !!background_smooth,
    tail = !!tail,
    step = !!step,
    beta = !!beta,
    tailing = !!tailing,
    compton_tail = !!compton_tail,
    compton_shape = !!compton_shape,
    cache_templates = !!cache_templates,
    refine_calibration = !!refine_calibration,
    sum_peaks = !!sum_peaks,
    count_rate = !!count_rate,
    pileup_tau = !!pileup_tau,
    use_qr = !!use_qr,
    abundance_prior = !!abundance_prior,
    abundance_ref_ppm = !!abundance_ref_ppm,
    abundance_protect = !!abundance_protect,
    abundance_ppm = !!abundance_ppm,
    prior = !!prior,
    enet_alpha = !!enet_alpha,
    .cores = .cores,
    .env = .env
  )
}

#' @rdname xrf_add_deconvolution_gls
#' @export
xrf_deconvolute_gaussian_least_squares <- function(energy_kev, response, peaks = xrf_energies(),
                                                   default_sigma = 0.07,
                                                   energy_min_kev = -Inf, energy_max_kev = Inf,
                                                   detector_type = NULL, fano = NULL,
                                                   epsilon_ev = NULL, noise_fwhm_ev = NULL,
                                                   efficiency = FALSE, escape = FALSE,
                                                   be_window_um = NULL, dead_layer_um = NULL,
                                                   active_thickness_um = NULL,
                                                   air_path_cm = NULL, atmosphere = "Air", window = NULL,
                                                   nonneg = TRUE,
                                                   weighting = c("auto", "poisson", "none"),
                                                   counts = NULL, livetime = NULL,
                                                   tube = NULL, geometry = NULL,
                                                   scatter = NULL, compton_broadening = 2,
                                                   scatter_continuum = FALSE, scatter_scatterer = "Si",
                                                   background = 0L, background_smooth = 0,
                                                   tail = 0, step = 0, beta = NULL, tailing = FALSE,
                                                   compton_tail = 0, compton_shape = c("profile", "gaussian"),
                                                   cache_templates = TRUE, refine_calibration = FALSE,
                                                   sum_peaks = FALSE, count_rate = NULL,
                                                   pileup_tau = NULL,
                                                   use_qr = TRUE,
                                                   abundance_prior = 0, abundance_ref_ppm = 100,
                                                   abundance_protect = character(0), abundance_ppm = NULL,
                                                   prior = c("ridge", "lasso", "elastic_net"), enet_alpha = 0.5) {

  # ---- true two-pass pile-up correction (#E3) ----------------------------------------------------
  # When a pulse-pair resolution time `pileup_tau` (s) is supplied with counts + livetime, model the
  # coincidence peaks with PHYSICALLY-FIXED amplitudes (proportional to the product of the parent
  # line areas), rather than the free-amplitude single-pass `sum_peaks`: fit the elements (pass 1),
  # subtract the predicted pile-up spectrum, then refit the elements (pass 2).
  if (!is.null(pileup_tau) && is.finite(pileup_tau) && pileup_tau > 0 &&
      !is.null(counts) && !is.null(livetime) && is.finite(livetime) && livetime > 0) {
    args <- mget(names(formals()))
    pass1 <- do.call(xrf_deconvolute_gaussian_least_squares,
                     utils::modifyList(args, list(sum_peaks = FALSE, pileup_tau = NULL)))
    pu <- .xrf_pileup_model(pass1$peaks, energy_kev, pileup_tau, livetime,
                            energy_min_kev, energy_max_kev)
    pass2 <- do.call(xrf_deconvolute_gaussian_least_squares,
                     utils::modifyList(args, list(response = response - pu$model,
                                                  sum_peaks = FALSE, pileup_tau = NULL)))
    if (nrow(pu$peaks) > 0) {
      pass2$peaks <- dplyr::bind_rows(pass2$peaks, pu$peaks)
      # pu$model is built on the FULL input grid, but pass2 filters to the fit window, so its
      # response_fit is only the in-window channels. Subset pu$model to the same window before
      # adding (otherwise the lengths differ -> hard error whenever energy_min/max restrict the grid).
      pu_within <- (energy_kev >= energy_min_kev) & (energy_kev <= energy_max_kev)
      pass2$response$response_fit <- pass2$response$response_fit + pu$model[pu_within]
    }
    # Restore the parents' pile-up LOSSES (the two passes above only remove the sum-peak artifacts;
    # the counts that migrated into them -- and into out-of-window sums -- came out of the parent
    # lines). Under Poisson arrivals a pulse survives un-piled with probability exp(-2 tau R_tot),
    # R_tot the gross count rate, so the true line intensities exceed the fitted ones by
    # exp(+2 tau R_tot). The factor is UNIFORM (any pulse piles with any other): ratios, Compton
    # normalization and closed concentrations are unchanged; absolute quantities (xrf_observed_mass,
    # calibrated masses) are corrected. Covariances scale by the factor squared. The pileup_ rows
    # keep their measured coincidence areas, and response/components stay in measured space.
    R_tot <- sum(counts, na.rm = TRUE) / livetime
    if (is.finite(R_tot) && R_tot > 0) {
      k_loss <- exp(2 * pileup_tau * R_tot)
      not_pu <- !grepl("^pileup_", pass2$peaks$element)
      for (cc in intersect(c("height", "height_se", "peak_area", "peak_area_se"), names(pass2$peaks))) {
        pass2$peaks[[cc]][not_pu] <- pass2$peaks[[cc]][not_pu] * k_loss
      }
      if (!is.null(pass2$coef_cov)) pass2$coef_cov <- pass2$coef_cov * k_loss^2
      if (!is.null(pass2$area_cov)) pass2$area_cov <- pass2$area_cov * k_loss^2
    }
    return(pass2)
  }

  weighting <- match.arg(weighting)
  prior <- match.arg(prior)
  if (is.finite(background_smooth) && background_smooth > 0 && prior != "ridge") {
    warning("background_smooth is only applied with prior = 'ridge'; ignored for the L1/elastic-net path.",
            call. = FALSE)
  }

  stopifnot(
    "energy_kev" %in% colnames(peaks),
    "element" %in% colnames(peaks),
    is.numeric(energy_kev), is.numeric(response),
    all(is.finite(energy_kev)), all(is.finite(response)),
    length(energy_kev) == length(response)
  )

  if(!("relative_peak_intensity" %in% colnames(peaks))) {
    peaks$relative_peak_intensity <- 1
  } else {
    stopifnot(is.numeric(peaks$relative_peak_intensity))
  }

  # Optional Rayleigh/Compton tube-scatter templates: appended as two extra "elements" so they are
  # fit (non-negatively) jointly with the element lines instead of being left in the baseline.
  # Enabled automatically when a tube is supplied; force with scatter = TRUE / FALSE. Appended
  # BEFORE the efficiency step below so a multi-line scatter template's INTERNAL shape (e.g. the
  # anode L components at ~3 keV vs K at ~20 keV) gets the same detector-efficiency weighting as
  # the element lines -- unweighted, the shared-amplitude template misfits both regions and pushes
  # residuals into whatever elements sit there.
  include_scatter <- if (is.null(scatter)) {
    !is.null(tube) && inherits(tube, "xrf_tube")
  } else {
    isTRUE(scatter)
  }
  if (include_scatter) {
    if (is.null(tube) || !inherits(tube, "xrf_tube")) {
      stop("scatter = TRUE requires tube = xrf_tube(...).")
    }
    scatter_peaks <- xrf_scatter_peaks(
      tube, if (is.null(geometry)) xrf_geometry() else geometry,
      detector_type = detector_type, fano = fano, epsilon_ev = epsilon_ev,
      noise_fwhm_ev = noise_fwhm_ev, default_sigma = default_sigma,
      compton_broadening = compton_broadening, scatterer = scatter_scatterer,
      compton_tail = compton_tail, compton_shape = compton_shape
    )
    peaks <- dplyr::bind_rows(peaks, scatter_peaks)
  }

  # Optional detector efficiency: reweight each line by the fraction actually recorded in the
  # photopeak (window/dead-layer transmission x active-volume absorption). Corrects intra-element
  # line shape and, for quantification, absolute sensitivity; largest at high energy where thin Si
  # detectors are transparent to heavy-element K-lines. Applies to the element AND the anode-line
  # scatter templates (the scattered-continuum templates below weight themselves with the same
  # settings).
  if (isTRUE(efficiency)) {
    peaks$relative_peak_intensity <- peaks$relative_peak_intensity *
      xrf_detector_efficiency(peaks$energy_kev, detector_type = detector_type,
                              active_thickness_um = active_thickness_um,
                              be_window_um = be_window_um, dead_layer_um = dead_layer_um,
                              air_path_cm = air_path_cm, atmosphere = atmosphere, window = window)
  }

  # Optional escape peaks: satellites at (parent energy - detector K X-ray), tied to the parent's
  # amplitude. Built from the (efficiency-corrected) parent intensities, so append after efficiency.
  # Scatter templates get escape satellites too (a 20 keV Rayleigh line has a real Si-escape peak at
  # ~18.3 keV), tied to the scatter amplitudes.
  if (isTRUE(escape)) {
    peaks <- dplyr::bind_rows(peaks, xrf_escape_peaks(peaks, detector_type = detector_type))
  }

  # Optional scattered bremsstrahlung CONTINUUM templates (Rayleigh + Compton of the tube continuum, #E1),
  # opt-in via scatter_continuum = TRUE (needs a tube). Two extra fixed-shape free-amplitude "elements".
  # IMPORTANT: these model the smooth scatter BACKGROUND, so they must be fit against a spectrum whose
  # background is NOT already removed, together with a fitted background (`background >= 1`). The working
  # recipe is: `.values = cps` (UN-baselined) + `scatter_continuum = TRUE` + `background = 6..12`, which
  # recovers trace areas cleanly (validated to r2 = 1 on synthetics, incl. a case where SNIP-only over-
  # estimates traces by 100-300% under a Compton hump). If instead `response` is baseline-subtracted (the
  # default `.values = cps - baseline`), SNIP has already removed most of this continuum and the broad
  # templates then bleed real peak area -- hence the warning when no fitted background accompanies it.
  if (isTRUE(scatter_continuum)) {
    if (is.null(tube) || !inherits(tube, "xrf_tube")) {
      stop("scatter_continuum = TRUE requires tube = xrf_tube(...).")
    }
    if (!(is.numeric(background) && length(background) == 1 && isTRUE(background >= 1))) {
      warning("scatter_continuum = TRUE models the scatter BACKGROUND and needs a companion fitted ",
              "background: use `.values = cps` (UN-baselined) with `background = 6..12`. Without a fitted ",
              "background (or on a SNIP-baseline-subtracted response) the broad continuum templates ",
              "double-count the baseline and bias real peak areas.", call. = FALSE)
    }
    cont_peaks <- xrf_scatter_continuum(
      tube, if (is.null(geometry)) xrf_geometry() else geometry,
      detector_type = detector_type, fano = fano, epsilon_ev = epsilon_ev,
      noise_fwhm_ev = noise_fwhm_ev, default_sigma = default_sigma,
      compton_broadening = compton_broadening, scatterer = scatter_scatterer,
      active_thickness_um = active_thickness_um, be_window_um = be_window_um,
      dead_layer_um = dead_layer_um, air_path_cm = air_path_cm, atmosphere = atmosphere,
      window = window
    )
    if (nrow(cont_peaks) > 0) peaks <- dplyr::bind_rows(peaks, cont_peaks)
  }

  # Optional fittable smooth background basis (#E1 baseline reconciliation): `background` broad, non-negative
  # Gaussian bumps spanning the fit window, fit jointly with the peaks so the deconvolution models the
  # continuum/background ITSELF -- fed an UN-baselined spectrum (.values = cps) -- instead of pre-subtracting
  # a SNIP baseline. This is what lets scatter_continuum work without double-counting: with a fitted
  # background (+ scattered-continuum templates) there is no subtracted baseline to bleed real trace area
  # into. Each bump is a free-amplitude "element" (background_k), excluded from quantification. Default 0 = off.
  if (is.numeric(background) && length(background) == 1 && isTRUE(background >= 1)) {
    nb <- as.integer(background)
    elo <- max(min(energy_kev), energy_min_kev); ehi <- min(max(energy_kev), energy_max_kev)
    if (is.finite(elo) && is.finite(ehi) && ehi > elo) {
      # nb == 1: one broad bump CENTERED on the window (not at the low edge, which gave a lopsided ramp);
      # nb >= 2: evenly spaced, overlapping (sigma >> detector line width).
      bg_E <- if (nb == 1) (elo + ehi) / 2 else seq(elo, ehi, length.out = nb)
      bg_sigma <- if (nb == 1) (ehi - elo) / 2 else (ehi - elo) / (nb - 1)
      peaks <- dplyr::bind_rows(peaks, tibble::tibble(
        element = paste0("background_", seq_len(nb)),
        trans = paste0("bg_", seq_len(nb)),
        energy_kev = bg_E, sigma = bg_sigma, relative_peak_intensity = 1
      ))
    }
  }

  # Per-line peak width: an explicit `sigma` column wins; otherwise, if any detector-resolution
  # parameter is given, compute the energy-dependent sigma(E); otherwise fall back to the constant
  # default_sigma. A single constant width is only correct near ~8-10 keV, so a detector model is
  # strongly preferred for wide-range work.
  detector_given <- !is.null(detector_type) || !is.null(fano) ||
    !is.null(epsilon_ev) || !is.null(noise_fwhm_ev)
  model_sigma <- function(e) {
    xrf_detector_sigma_kev(e, detector_type = detector_type, fano = fano,
                           epsilon_ev = epsilon_ev, noise_fwhm_ev = noise_fwhm_ev)
  }
  if(!("sigma" %in% colnames(peaks))) {
    peaks$sigma <- if (detector_given) model_sigma(peaks$energy_kev) else default_sigma
  } else {
    stopifnot(is.numeric(peaks$sigma))
    fallback <- if (detector_given) model_sigma(peaks$energy_kev) else default_sigma
    peaks$sigma <- dplyr::coalesce(peaks$sigma, fallback)
  }

  # Optional coincidence (pile-up) sum peaks at E_i + E_j (#21): at high count rate two photons are
  # recorded as one event at the summed energy. Giving the strongest lines their own sum templates
  # (fit freely) keeps pile-up from being misassigned to real trace elements. Simplified single-pass
  # model; the pile-up rate scales with `count_rate` (used only to skip the step when clearly low).
  if (isTRUE(sum_peaks) && (is.null(count_rate) || count_rate > 0)) {
    lc <- peaks[is.finite(peaks$energy_kev) & is.finite(peaks$sigma) &
                  peaks$relative_peak_intensity > 0 &
                  !grepl("^(scatter|background|sum|pileup|escape)_", peaks$element), , drop = FALSE]
    # Rank candidate pile-up parents by their actual amplitude IN THE DATA (pile-up scales with the line
    # count rate), not the raw per-element template `relative_peak_intensity`, which is not comparable
    # across elements (every element's primary line is ~1). Use the peak response within +/- sigma of each
    # line as a cheap, cross-element-comparable amplitude proxy so the sum templates land on the sums of the
    # genuinely dominant peaks. (Spurious templates are harmless -- NNLS zeros them -- but a MISSED true
    # pile-up location is not corrected, which is what the old template-intensity ranking risked.)
    amp_proxy <- vapply(seq_len(nrow(lc)), function(k) {
      inw <- energy_kev >= lc$energy_kev[k] - lc$sigma[k] & energy_kev <= lc$energy_kev[k] + lc$sigma[k]
      if (any(inw)) max(response[inw]) else 0
    }, numeric(1))
    lc <- lc[order(-amp_proxy), , drop = FALSE]
    max_lines <- min(nrow(lc), 5L)
    if (max_lines >= 1) {
      lc <- lc[seq_len(max_lines), ]
      pairs <- utils::combn(max_lines, 2, simplify = FALSE)
      pairs <- c(pairs, lapply(seq_len(max_lines), function(i) c(i, i)))   # + 2x self-sums
      sum_rows <- purrr::map_dfr(pairs, function(ij) {
        tibble::tibble(
          element = paste0("sum_", round(lc$energy_kev[ij[1]] + lc$energy_kev[ij[2]], 2)),
          energy_kev = lc$energy_kev[ij[1]] + lc$energy_kev[ij[2]],
          sigma = sqrt(lc$sigma[ij[1]]^2 + lc$sigma[ij[2]]^2),
          relative_peak_intensity = 1
        )
      })
      sum_rows <- sum_rows[!duplicated(round(sum_rows$energy_kev, 3)) &
                             sum_rows$energy_kev >= energy_min_kev &
                             sum_rows$energy_kev <= energy_max_kev, , drop = FALSE]
      if (nrow(sum_rows) > 0) peaks <- dplyr::bind_rows(peaks, sum_rows)
    }
  }

  # filter energies/responses (and counts, if supplied) to the requested window
  within_range <- (energy_kev >= energy_min_kev) & (energy_kev <= energy_max_kev)
  energy_kev <- energy_kev[within_range]
  response <- response[within_range]
  if(!is.null(counts)) {
    stopifnot(length(counts) == length(within_range))
    counts <- counts[within_range]
  }

  # check finite-ness of inputs
  stopifnot(
    all(is.finite(peaks$sigma)),
    all(is.finite(peaks$relative_peak_intensity))
  )

  # filter, normalize peak heights (scale to sigma)
  # the peak with a relative height of 1 should be the peak with the greatest intensity
  peaks <- peaks %>%
    dplyr::filter(
      .data$energy_kev >= !!energy_min_kev,
      .data$energy_kev <= !!energy_max_kev,
      .data$relative_peak_intensity > 0
    ) %>%
    dplyr::group_by(.data$element) %>%
    dplyr::mutate(
      relative_peak_intensity = .data$relative_peak_intensity / max(.data$relative_peak_intensity),
      relative_peak_height = .data$relative_peak_intensity / .data$sigma
    ) %>%
    dplyr::mutate(
      relative_peak_height = .data$relative_peak_height /
        .data$relative_peak_height[which.max(.data$relative_peak_intensity)]
    ) %>%
    dplyr::ungroup()

  # ---- per-line Hypermet columns --------------------------------------------------------------
  # The `tail`/`step`/`beta` scalars stay as GLOBAL defaults; optional per-line `tail`/`step`/`beta`
  # COLUMNS on `peaks` override them row-wise (e.g. the asymmetric Compton tail set by
  # xrf_scatter_peaks(compton_tail = ...)), and `tailing = TRUE` fills the remaining rows from the
  # energy-dependent incomplete-charge-collection model xrf_detector_tailing(). Per-line shelves
  # matter for BASELINE fidelity: the step of every strong peak -- above all the big scatter
  # features -- is a genuine part of the measured continuum, and a single global scalar cannot
  # describe Si ICC (grows at low E) and the high-energy scatter shelves at once.
  user_tail <- if ("tail" %in% colnames(peaks)) peaks$tail else rep(NA_real_, nrow(peaks))
  user_step <- if ("step" %in% colnames(peaks)) peaks$step else rep(NA_real_, nrow(peaks))
  user_beta <- if ("beta" %in% colnames(peaks)) peaks$beta else rep(NA_real_, nrow(peaks))
  if (isTRUE(tailing)) {
    dtl <- xrf_detector_tailing(peaks$energy_kev, detector_type = detector_type, fano = fano,
                                epsilon_ev = epsilon_ev, noise_fwhm_ev = noise_fwhm_ev)
    peaks$tail <- dplyr::coalesce(user_tail, dtl$tail)
    peaks$step <- dplyr::coalesce(user_step, dtl$step)
    peaks$beta <- dplyr::coalesce(user_beta, dtl$beta)
  } else {
    peaks$tail <- dplyr::coalesce(user_tail, tail)
    peaks$step <- dplyr::coalesce(user_step, step)
    peaks$beta <- dplyr::coalesce(user_beta, if (is.null(beta)) NA_real_ else beta)
  }

  # Build the per-element template design matrix. Reused across a batch (#15): if the energy grid,
  # the (normalized) peaks and the line-shape are identical to the previous call, the cached matrix
  # is returned instead of re-evaluating every Gaussian.
  peaks_tmpl <- peaks %>%
    dplyr::select("element", "energy_kev", "sigma", "relative_peak_height", "relative_peak_intensity",
                  "tail", "step", "beta")

  # per-line response: a windowed pure Gaussian (evaluate exp() only within +/-6 sigma, #19) unless a
  # Hypermet tail/shelf is requested for THAT line, in which case the full shape is evaluated (#20).
  build_line <- function(mu, sigma, height, tl, st, be) {
    if (tl == 0 && st == 0) {
      out <- numeric(length(energy_kev))
      idx <- which(energy_kev >= mu - 6 * sigma & energy_kev <= mu + 6 * sigma)
      if (length(idx)) out[idx] <- height * exp(-0.5 * ((energy_kev[idx] - mu) / sigma) ^ 2)
      out
    } else {
      xrf_lineshape(energy_kev, mu = mu, sigma = sigma, height = height,
                    tail = tl, step = st, beta = be)
    }
  }
  # assemble the per-element template matrix + metadata from a peaks table
  build_X <- function(pk) {
    built <- pk %>%
      dplyr::mutate(
        response_element = purrr::pmap(
          list(.data$energy_kev, .data$sigma, .data$relative_peak_height,
               .data$tail, .data$step, .data$beta), build_line
        )
      ) %>%
      dplyr::group_by(.data$element) %>%
      dplyr::summarise(
        primary_energy_kev = .data$energy_kev[which.max(.data$relative_peak_intensity)],
        primary_sigma = .data$sigma[which.max(.data$relative_peak_intensity)],
        primary_tail = .data$tail[which.max(.data$relative_peak_intensity)],
        primary_beta = .data$beta[which.max(.data$relative_peak_intensity)],
        response_element = list(purrr::reduce(.data$response_element, `+`))
      ) %>%
      dplyr::ungroup()
    Xb <- do.call(cbind, built$response_element)
    colnames(Xb) <- built$element
    list(X = Xb, meta = built[, c("element", "primary_energy_kev", "primary_sigma",
                                  "primary_tail", "primary_beta")])
  }

  # ---- optional energy-calibration refinement (#16) ----------------------------------------------
  # A small zero/gain miscalibration shifts every centroid, and because templates are non-linear in
  # the centroid this injects large derivative-shaped residuals. Refine the 2-parameter map
  # E' = zero + gain * E of the template centroids by variable projection: for each (zero, gain) the
  # amplitudes are solved (fast NNLS) and the residual sum of squares minimized over (zero, gain).
  if (isTRUE(refine_calibration)) {
    # Weight the calibration objective (and its inner solve) by counting statistics -- Poisson data weights
    # w = livetime^2 / gross_counts when counts/livetime are available, else a variance proxy 1/max(response, .)
    # -- so the (zero, gain) alignment is not dominated by the few tallest peaks and weak low-energy lines get
    # their statistical say (consistent with the weighted final fit rather than the old unweighted RSS).
    w_ref <- if (!is.null(counts) && !is.null(livetime) && length(livetime) == 1 &&
                 is.finite(livetime) && livetime > 0 && length(counts) == length(response)) {
      livetime^2 / pmax(counts, 1)
    } else {
      1 / pmax(response, max(response, na.rm = TRUE) * 1e-3, 1)
    }
    rw_ref <- sqrt(w_ref)
    rss_for <- function(par) {
      pk <- peaks_tmpl
      pk$energy_kev <- par[1] + par[2] * pk$energy_kev
      Xc <- build_X(pk)$X
      cf <- nnls::nnls(Xc * rw_ref, response * rw_ref)$x
      sum((response - as.vector(Xc %*% cf)) ^ 2 * w_ref)
    }
    opt <- tryCatch(
      stats::optim(c(0, 1), rss_for, method = "Nelder-Mead",
                   control = list(reltol = 1e-9, maxit = 200)),
      error = function(e) NULL
    )
    if (!is.null(opt) && opt$convergence == 0) {
      peaks_tmpl$energy_kev <- opt$par[1] + opt$par[2] * peaks_tmpl$energy_kev
    }
  }

  # ---- design matrix (cached across a batch, #15) ------------------------------------------------
  tmpl_key <- list(energy_kev = energy_kev, peaks = peaks_tmpl, tail = tail, step = step, beta = beta)
  cached <- if (cache_templates && !isTRUE(refine_calibration)) .xrf_kv_get(.xrf_template_cache, tmpl_key)
  if (!is.null(cached)) {
    responses <- cached$responses
    X <- cached$X
  } else {
    bx <- build_X(peaks_tmpl)
    X <- bx$X
    responses <- bx$meta
    if (cache_templates && !isTRUE(refine_calibration)) {
      .xrf_kv_put(.xrf_template_cache, tmpl_key, list(responses = responses, X = X))
    }
  }
  elements <- responses$element
  n <- length(response)
  p <- ncol(X)

  # ---- abundance-informed ridge prior (opt-in; #abundance) ---------------------------------------
  # Per-column penalty factor ~ 1/crustal_abundance (0 for scatter_*/sum_* and off-table columns). The
  # absolute penalty applied in solve_fit is `abundance_prior * data_scale * pen_factor`, so a degenerate
  # overlap splits toward the more abundant element; strength 0 disables it (identical legacy fit).
  use_prior <- is.finite(abundance_prior) && abundance_prior > 0
  # Composition-independent abundance factor ~ 1/crustal_abundance, computed once. It is gated per solve
  # by a VIF identifiability weight on the WEIGHTED design (inside solve_fit, so it updates each IRLS pass
  # and reflects the problem the solver sees), so only rare AND non-identifiable columns are penalised.
  pf_abund <- if (use_prior) .xrf_abundance_penalty_factor(elements, abundance_ref_ppm, abundance = abundance_ppm) else rep(0, p)
  # Protect user-specified elements from the abundance prior (pen_factor 0 -> never shrunk). This is the
  # abundance-tiebreaker made explicit: an element you KNOW is present (e.g. Pb/U read via L-lines that
  # overlap an abundant K-line or the scatter peak) keeps its intensity, so at the coincidence the prior
  # penalises only the rare, unprotected competitor and the shared area flows to the protected element.
  # (The VIF gate otherwise reads a heavy-L element whose STRONG lines overlap as unidentifiable -- even
  # when it has a weaker clean line -- and crushes it.)
  if (use_prior && length(abundance_protect)) {
    pf_abund[elements %in% as.character(abundance_protect)] <- 0
  }
  # Element columns only (exclude scatter/background/continuum/sum/pileup pseudo-elements): the abundance
  # ridge's data-scale is computed over these so that enabling the broad `background` / scatter-continuum
  # templates -- whose column self-energies are ~10-30x an element line's -- does not silently inflate the
  # effective `abundance_prior` strength.
  is_element_col <- !grepl("^(background|scatter|sum|pileup|escape)_", elements)
  if (!any(is_element_col)) is_element_col <- rep(TRUE, p)
  # F2: elements with a spectrally isolated (clean) line are identifiable even if their STRONG lines overlap
  # another template, so exempt them from the abundance penalty -- an imprecise-but-estimable trace must not
  # be crushed (e.g. As via its clean Kbeta when As Kalpha overlaps Pb Lalpha). Computed once on the fixed
  # template shapes; only meaningful when the prior is active.
  clean_exempt <- if (use_prior && any(pf_abund > 0))
    .xrf_clean_line_exempt(peaks_tmpl, X, elements, is_element_col, energy_kev) else rep(FALSE, p)

  # ---- counting-statistics weighting -------------------------------------------------------------
  # var(response_i [cps]) = gross_counts_i / livetime^2 (Poisson), so w_i = livetime^2 / gross_counts.
  # `response` is net (baseline-subtracted) cps; gross counts are reconstructed from `counts`.
  have_counts <- !is.null(counts) && !is.null(livetime) &&
    length(livetime) == 1 && is.finite(livetime) && livetime > 0 &&
    length(counts) == n
  if (weighting == "auto") weighting <- if (have_counts) "poisson" else "none"
  if (weighting == "poisson" && !have_counts) {
    warning("weighting = 'poisson' requires `counts` and `livetime`; falling back to 'none'.")
    weighting <- "none"
  }
  baseline_cps <- if (have_counts) (counts / livetime) - response else NULL

  # ---- solver ------------------------------------------------------------------------------------
  # Returns the coefficient vector and the "active" set (the estimable / non-clamped columns).
  # For NNLS the active set is the passive set reported by the solver (its exactly-non-zero
  # columns), which is more reliable than a numerical `coef > 0` test on floating-point dust.
  solve_fit <- function(w) {
    rw <- sqrt(w)
    Xw <- X * rw          # row-scale (recycles rw down each column)
    yw <- response * rw
    ridge_applied <- FALSE
    penalty <- rep(0, p)          # effective per-column L2 penalty applied (0 = unpenalised); for the Laplace SE
    if (use_prior && any(pf_abund > 0)) {
      # VIF identifiability gate on the WEIGHTED design actually being solved: an estimable column
      # (VIF ~ 1) is left unpenalised even if it spectrally overlaps another, so only genuinely
      # non-identifiable rare columns are shrunk.
      pen_factor <- pf_abund * .xrf_template_vif_weight(Xw)
      pen_factor[clean_exempt] <- 0   # F2: never penalise an element that has a clean (isolated) line
      # Tikhonov augmentation: append p pseudo-observation rows sqrt(lambda_j) with target 0, so the
      # solve minimises ||W^1/2(y-Xb)||^2 + sum lambda_j b_j^2. lambda scales to the data (mean column
      # self-energy) so `abundance_prior` is a dimensionless strength; penalty rows carry weight 1 (a
      # prior is not a Poisson observation), so this is correct on every IRLS reweighting pass.
      scale <- mean(colSums(Xw[, is_element_col, drop = FALSE]^2))
      if (any(pen_factor > 0) && is.finite(scale) && scale > 0) {
        # L1 / elastic-net path (IRL2, .xrf_irl2_fit): penalise |b_j| (and, for elastic_net, also b_j^2) with the SAME
        # per-column gated penalty factor as the ridge, so only rare AND non-identifiable columns are shrunk
        # -- but lasso drives them to EXACTLY zero (sparse element selection) rather than merely shrinking.
        if (prior != "ridge") {
          return(.xrf_irl2_fit(X, response, rw, scale, pen_factor,
                               alpha = if (prior == "lasso") 1 else enet_alpha,
                               strength = abundance_prior, nonneg = nonneg, p = p))
        }
        lambda <- abundance_prior * scale * pen_factor
        penalty <- lambda
        Xw <- rbind(Xw, diag(sqrt(lambda), p, p))
        yw <- c(yw, rep(0, p))
        ridge_applied <- TRUE
      }
    }
    # ---- background roughness penalty (P-spline style) -------------------------------------------
    # `background_smooth` > 0 penalises the SECOND DIFFERENCES of the fitted background-bump
    # coefficients (lambda ||D2 b_bg||^2 via augmented rows), so a generous background basis stays
    # smooth instead of chasing peak residuals -- more baseline flexibility WITHOUT peak-stealing.
    # Scaled to the data like the abundance ridge (dimensionless strength); element columns are
    # untouched; the diagonal of lambda D2'D2 feeds the Laplace SE penalty (a documented diagonal
    # approximation). Skipped under the L1/elastic-net prior path (which solves separately).
    if (is.finite(background_smooth) && background_smooth > 0 && prior == "ridge") {
      idx_bg <- which(grepl("^background_", elements))
      if (length(idx_bg) >= 3) {
        idx_bg <- idx_bg[order(as.integer(sub("^background_", "", elements[idx_bg])))]
        # scaled by the mean column self-energy over n_bg^2, so strength ~1 is a GENTLE default that
        # barely perturbs a well-resolved smooth feature; tens-to-hundreds give a stiff spline
        sc_bg <- mean(colSums(Xw[seq_len(n), idx_bg, drop = FALSE]^2)) / length(idx_bg)^2
        if (is.finite(sc_bg) && sc_bg > 0) {
          D2 <- diff(diag(length(idx_bg)), differences = 2L)
          Arows <- matrix(0, nrow(D2), p)
          Arows[, idx_bg] <- sqrt(background_smooth * sc_bg) * D2
          Xw <- rbind(Xw, Arows)
          yw <- c(yw, rep(0, nrow(D2)))
          penalty[idx_bg] <- penalty[idx_bg] + background_smooth * sc_bg * diag(crossprod(D2))
          ridge_applied <- TRUE
        }
      }
    }
    if (nonneg) {
      nn <- nnls::nnls(Xw, yw)
      active <- if (nn$nsetp > 0) sort(nn$passive[seq_len(nn$nsetp)]) else integer(0)
      cf <- nn$x
      # Treat coefficients negligible relative to the largest as clamped (drop NNLS numerical
      # dust): they are set exactly to zero and excluded from the active set, so they report an
      # NA (rather than a spurious) standard error.
      if (length(active) > 0) {
        tol <- max(cf[active]) * 1e-8
        active <- active[cf[active] > tol]
      }
      cf[setdiff(seq_len(p), active)] <- 0
      list(coef = cf, active = active, aliased = integer(0), penalty = penalty)
    } else if (use_qr) {
      cf <- qr.coef(qr(Xw), yw)
      aliased <- which(is.na(cf))     # rank-deficient columns pivoted out by the QR
      cf[aliased] <- 0
      list(coef = cf, active = setdiff(seq_len(p), aliased), aliased = aliased, penalty = penalty)
    } else if (ridge_applied) {
      cf <- qr.coef(qr(Xw), yw)      # lm.wfit can't take the penalty rows; ridge via augmented QR
      aliased <- which(is.na(cf))
      cf[aliased] <- 0
      list(coef = cf, active = setdiff(seq_len(p), aliased), aliased = aliased, penalty = penalty)
    } else {
      cf <- unname(stats::lm.wfit(X, response, w = w)$coefficients)
      aliased <- which(is.na(cf))
      cf[aliased] <- 0
      list(coef = cf, active = setdiff(seq_len(p), aliased), aliased = aliased, penalty = penalty)
    }
  }

  if (weighting == "poisson") {
    # iteratively reweighted: model-based Poisson weights (avoids the zero-count bias of data weights)
    w <- rep(1, n)
    for (iter in 1:3) {
      coef <- solve_fit(w)$coef
      fit_iter <- as.vector(X %*% coef)
      model_gross_counts <- pmax((fit_iter + baseline_cps) * livetime, 1)
      w <- livetime^2 / model_gross_counts
    }
    final <- solve_fit(w)
  } else {
    w <- rep(1, n)
    final <- solve_fit(w)
  }
  coef <- final$coef
  names(coef) <- elements
  active <- final$active
  aliased <- final$aliased
  response_fit <- as.vector(X %*% coef)
  residuals <- response - response_fit

  # ---- rank / collinearity diagnostics (#22) -----------------------------------------------------
  rw <- sqrt(w)
  Xw <- X * rw
  # condition number of the (weighted) design: large values (>~1e3) flag near-collinear templates
  # whose amplitudes are poorly separated (e.g. overlapping lines of different elements).
  # over the ELEMENT columns only -- the diagnostic is documented as element-template collinearity, and the
  # broad background_/scatter-continuum columns would otherwise dominate kappa and mask the real overlaps.
  condition_number <- tryCatch(kappa(Xw[, is_element_col, drop = FALSE], exact = FALSE),
                               error = function(e) NA_real_)
  if (length(aliased) > 0) {
    warning("Deconvolution design matrix is rank-deficient; the following element templates are ",
            "aliased (not separable) and are reported as NA rather than 0: ",
            paste(elements[aliased], collapse = ", "), call. = FALSE)
  }

  # ---- coefficient standard errors + Laplace posterior covariance --------------------------------
  # SEs and covariance from the estimable (active) coefficient set. The Laplace posterior covariance is
  #   Sigma_a = sigma2 (Xa' W Xa + Lambda_a)^-1,
  # with Lambda_a the per-column penalty Hessian returned by the solver (0 for unpenalised columns, so an
  # unpenalised fit is byte-identical to before). Including the prior curvature makes each SE the shrinkage-
  # consistent posterior SD of the (biased) penalised estimate, and the OFF-DIAGONAL terms carry the
  # coefficient correlations that ratio / closure quantities need downstream (e.g. Compton normalization).
  # Computed from R of the QR of the augmented [Xa; sqrt(Lambda_a)] (stable; the same design the fit solved).
  # Clamped / aliased columns get NA. NOTE (NNLS): with non-negativity the inference conditions on the active,
  # non-boundary set -- a standard approximation that ignores the inequality constraints; for the L1 prior the
  # local L2-equivalent curvature is used (approximate, since the posterior is non-Gaussian at exactly 0).
  rank <- if (nonneg) length(active) else qr(Xw)$rank
  sigma2 <- sum(w * residuals^2) / max(n - length(active), 1)
  coef_se <- rep(NA_real_, p); names(coef_se) <- elements
  coef_cov <- matrix(0, p, p, dimnames = list(elements, elements))
  penalty <- if (!is.null(final$penalty)) final$penalty else rep(0, p)
  if (length(active) > 0) {
    na_ <- length(active)
    Xaug <- rbind(Xw[, active, drop = FALSE], diag(sqrt(pmax(penalty[active], 0)), na_, na_))
    qra <- qr(Xaug); ra <- qra$rank
    if (ra == na_) {
      Rinv <- backsolve(qr.R(qra), diag(na_))
      covP <- tcrossprod(Rinv) * sigma2          # covariance in the QR's (possibly pivoted) column order
      piv <- qra$pivot
      covA <- matrix(0, na_, na_); covA[piv, piv] <- covP   # unpermute to `active` column order
      coef_cov[active, active] <- covA
      coef_se[active] <- sqrt(pmax(diag(covA), 0))
    } else {
      # rank-deficient active set: estimable-pivot SEs only (diagonal), NA for the aliased
      Ra <- qr.R(qra)[1:ra, 1:ra, drop = FALSE]; piv <- qra$pivot[1:ra]
      coef_se[active[piv]] <- sqrt(rowSums(backsolve(Ra, diag(ra))^2) * sigma2)
    }
  }

  # ---- goodness of fit ---------------------------------------------------------------------------
  ss_res <- sum(residuals^2)
  ss_tot <- sum((response - mean(response))^2)
  r2 <- 1 - ss_res / ss_tot
  r2_adj <- 1 - (1 - r2) * (n - 1) / (n - rank)

  # Reduced chi-square. In count space (Pearson) when counts are available -- never divide by the
  # measured value (the historic sum(resid^2 / response) is NaN/Inf on the many zero-count channels
  # of a real spectrum). Falls back to a guarded net-space diagnostic otherwise.
  if (have_counts) {
    mu <- pmax((response_fit + baseline_cps) * livetime, .Machine$double.eps)
    chi_sq <- sum((counts - mu)^2 / pmax(mu, 1)) / max(n - rank, 1)
  } else {
    valid <- response > 0
    chi_sq <- if (any(valid)) {
      sum(residuals[valid]^2 / response[valid]) / max(sum(valid) - rank, 1)
    } else {
      NA_real_
    }
  }

  df <- tibble::tibble(
    energy_kev = energy_kev,
    response = response,
    response_fit = response_fit
  )

  n_energy <- length(energy_kev)
  n_elements <- length(elements)
  components <- tibble::tibble(
    energy_kev = rep(energy_kev, n_elements),
    element = rep(elements, each = n_energy),
    response_fit = as.vector(X) * rep(coef, each = n_energy),
    height = rep(coef, each = n_energy)
  )

  # Build peaks output. peak_area is the primary line's Gaussian area, plus the Hypermet tail area
  # (integral of the exponential tail ~ tail * beta) when a tail is modelled.
  responses$height <- coef[elements]
  responses$height_se <- coef_se[elements]
  if (length(aliased) > 0) {           # aliased (rank-deficient) templates: report NA, not 0 (#22)
    responses$height[aliased] <- NA_real_
    responses$height_se[aliased] <- NA_real_
  }
  tail_beta <- ifelse(is.finite(responses$primary_beta) & responses$primary_beta > 0,
                      responses$primary_beta, responses$primary_sigma)
  # Exact area of the Hypermet exponential tail per unit height: the integral of
  # tail*exp(d/beta)*erfc(d/(sqrt2 sigma) + sigma/(sqrt2 beta)) over d is tail * 2*beta*exp(-sigma^2/(2 beta^2)),
  # NOT tail*beta -- the dropped 2*exp(-sigma^2/(2 beta^2)) factor depends on sigma (grows with line energy),
  # so the old form biased peak_area energy-dependently (~18% for beta~sigma, up to ~2x for long tails).
  # Uses each element's PRIMARY LINE tail/beta (per-line columns), which reduce to the global scalars
  # when no columns were supplied.
  tail_area <- ifelse(responses$primary_tail != 0,
                      responses$primary_tail * 2 * tail_beta *
                        exp(-responses$primary_sigma^2 / (2 * tail_beta^2)), 0)
  area_per_height <- responses$primary_sigma * sqrt(2 * pi) + tail_area
  responses$peak_area <- responses$height * area_per_height
  responses$peak_area_se <- responses$height_se * area_per_height

  # peak-area covariance (a_i a_j Cov(coef_i, coef_j)): the correlated companion of peak_area_se, consumed by
  # ratio / closure error propagation (e.g. xrf_compton_normalize). Element order matches `elements`.
  area_cov <- outer(area_per_height, area_per_height) * coef_cov
  dimnames(area_cov) <- list(elements, elements)

  structure(
    list(
      fit = tibble::tibble(
        r2 = r2,
        r2_adj = r2_adj,
        chi_sq = chi_sq,
        rank = rank,
        condition_number = condition_number,
        nonneg = nonneg,
        weighting = weighting
      ),
      response = df,
      components = components,
      peaks = responses,
      coef_cov = coef_cov,
      area_cov = area_cov
    ),
    class = "deconvolution_fit"
  )
}

#' Add a deconvolution based on a function
#'
#' @param .spectra A \link{xrf_spectra} tibble.
#' @param .fun A function that
#' @param ... Passed to .fun
#' @param .cores Number of processes for the per-spectrum map (default 1 = serial). Values > 1 use
#'   \code{parallel::mclapply} (fork-based; on Windows this silently runs serially). Forked workers
#'   inherit -- but cannot write back to -- the template/line/sensitivity caches, so for best scaling
#'   on a repetitive batch, run one spectrum per beam condition first to warm the caches, then the
#'   rest in parallel.
#' @param .env Calling environment
#'
#' @return .spectra with a .deconvolution column, and a deconvolution column added to
#'   the .spectra column
#' @export
#'
xrf_add_deconvolution_fun <- function(.spectra, .fun, ..., .cores = 1L, .env = parent.frame()) {
  dots <- quos(...)
  fun_wrap <- function(spectrum) {
    # args evaluated within each row
    args <- purrr::map(dots, rlang::eval_tidy, data = spectrum, env = .env)
    do.call(.fun, args)
  }

  # calculate deconvolutions
  deconvs <- if (.cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(purrr::transpose(.spectra), fun_wrap, mc.cores = .cores)
  } else {
    purrr::map(purrr::transpose(.spectra), fun_wrap)
  }

  # place each element of the output in a column in the data frame

  # check for named lists
  deconv_not_list <- !purrr::map_lgl(deconvs, is.list)
  if(any(deconv_not_list)) {
    stop(
      "Result of deconvolution function was not a list at positions ",
      paste(which(deconv_not_list), collapse = ", ")
    )
  }

  deconv_names <- purrr::map(deconvs, names)
  deconv_no_names <- purrr::map_lgl(
    deconv_names,
    function(x) is.null(x) || ("" %in% x)
  )

  if(any(deconv_no_names)) {
    stop(
      "Result of deconvolution function had empty or missing names at positions ",
      paste(which(deconv_no_names), collapse = ", ")
    )
  }

  for(name in unique(do.call(c, deconv_names))) {
    .spectra[[paste(".deconvolution", name, sep = "_")]] <- purrr::map(deconvs, name)
  }

  .spectra
}

gaussian_fun <- function(energy_kev, mu = 0, sigma = 1, height = 1) {
  height * exp(-0.5 * ((energy_kev - mu) / sigma) ^ 2)
}

#' XRF/EDS detector line shape (Gaussian with optional Hypermet tail and shelf)
#'
#' A real semiconductor-detector peak is not a pure Gaussian: incomplete charge collection adds a
#' low-energy exponential tail and a flat shelf (step). This returns the Gaussian plus optional
#' Hypermet tail/step terms (Campbell):
#' \deqn{h\big[\,e^{-\frac{(E-\mu)^2}{2\sigma^2}} + \tfrac{step}{2}\,\mathrm{erfc}(\tfrac{E-\mu}{\sqrt2\sigma})
#'   + tail\,e^{(E-\mu)/\beta}\,\mathrm{erfc}(\tfrac{E-\mu}{\sqrt2\sigma} + \tfrac{\sigma}{\sqrt2\beta})\,\big].}
#' With \code{tail = step = 0} (the default) it is exactly a Gaussian, so existing behaviour is
#' preserved. The tail carries a few percent to >10\% of the counts (worst at low energy and for
#' CdTe hole-tailing at high energy) and, if unmodelled, leaks into lower-energy neighbours.
#'
#' @param energy_kev Energies at which to evaluate (keV).
#' @param mu Peak centroid (keV).
#' @param sigma Gaussian standard deviation (keV).
#' @param height Peak height.
#' @param tail Relative amplitude of the low-energy exponential tail (0 = none).
#' @param step Relative amplitude of the low-energy shelf/step (0 = none).
#' @param beta Tail decay length (keV); defaults to \code{sigma}.
#'
#' @return A numeric vector the length of \code{energy_kev}.
#' @export
#'
xrf_lineshape <- function(energy_kev, mu = 0, sigma = 1, height = 1, tail = 0, step = 0, beta = NULL) {
  erfc <- function(x) 2 * stats::pnorm(x * sqrt(2), lower.tail = FALSE)
  d <- energy_kev - mu
  out <- exp(-0.5 * (d / sigma) ^ 2)
  if (step != 0) out <- out + step * 0.5 * erfc(d / (sqrt(2) * sigma))
  if (tail != 0) {
    b <- if (is.null(beta) || !is.finite(beta) || beta <= 0) sigma else beta
    out <- out + tail * exp(pmin(d / b, 50)) * erfc(d / (sqrt(2) * sigma) + sigma / (sqrt(2) * b))
  }
  height * out
}

.empty_pileup <- function() {
  tibble::tibble(element = character(0), primary_energy_kev = numeric(0), primary_sigma = numeric(0),
                 height = numeric(0), height_se = numeric(0), peak_area = numeric(0),
                 peak_area_se = numeric(0))
}

# Fixed-amplitude coincidence (pile-up) model from a first-pass element fit. For the strongest lines,
# the coincidence RATE is r_sum = 2*tau*R_i*R_j (i!=j) or tau*R_i^2 (i==j), with R the line count
# RATE (cps). The response here is a cps spectrum, so peak_area = height*sigma*sqrt(2*pi) has units
# cps*keV = R*dE (dE = channel width), i.e. R = peak_area/dE. The sum peak's own area (cps*keV) is
# therefore r_sum*dE = 2*tau*A_i*A_j/dE (self: tau*A_i^2/dE) -- live time CANCELS in rate space; the
# surviving normalizer is the energy channel width dE, not the live time. Width adds in quadrature.
.xrf_pileup_model <- function(peaks, energy_kev, pileup_tau, livetime, energy_min_kev, energy_max_kev,
                              max_lines = 6L) {
  model <- numeric(length(energy_kev))
  dE <- abs(stats::median(diff(energy_kev)))            # energy channel width (keV); abs() so a descending
                                                        # grid does not give a negative dE (which silently
                                                        # zeroed the pile-up); median approximates a non-uniform grid.
  if (!is.finite(dE) || dE <= 0) return(list(model = model, peaks = .empty_pileup()))
  el <- peaks[!grepl("^(scatter|sum|pileup|escape|background)_", peaks$element) &
                is.finite(peaks$peak_area) & peaks$peak_area > 0 &
                is.finite(peaks$primary_energy_kev) & is.finite(peaks$primary_sigma), , drop = FALSE]
  if (nrow(el) == 0) return(list(model = model, peaks = .empty_pileup()))
  el <- el[order(-el$peak_area), , drop = FALSE]
  k <- min(nrow(el), max_lines)
  el <- el[seq_len(k), ]

  out <- list()
  for (i in seq_len(k)) for (j in i:k) {
    e_sum <- el$primary_energy_kev[i] + el$primary_energy_kev[j]
    if (!is.finite(e_sum) || e_sum < energy_min_kev || e_sum > energy_max_kev) next
    area <- (if (i == j) 1 else 2) * pileup_tau * el$peak_area[i] * el$peak_area[j] / dE
    if (!is.finite(area) || area <= 0) next
    s_sum <- sqrt(el$primary_sigma[i]^2 + el$primary_sigma[j]^2)
    height <- area / (s_sum * sqrt(2 * pi))
    model <- model + height * exp(-0.5 * ((energy_kev - e_sum) / s_sum)^2)
    out[[length(out) + 1]] <- tibble::tibble(
      element = paste0("pileup_", el$element[i], "_", el$element[j]),
      primary_energy_kev = e_sum, primary_sigma = s_sum,
      height = height, height_se = NA_real_, peak_area = area, peak_area_se = NA_real_
    )
  }
  list(model = model, peaks = if (length(out)) dplyr::bind_rows(out) else .empty_pileup())
}
