
#' Elemental X-ray energies and intensities
#'
#' Transition energies that may be found in X-Ray spectra and that may be useful in the analysis
#' of XRF data. This dataset (from NIST 2018) includes most transition energies
#' for all elements. References in the \code{ref} column can be found at the
#' \href{https://physics.nist.gov/PhysRefData/XrayTrans/Html/refs.html}{NIST X-Ray Transitions Database References}
#' page. Relative intensities are not available for most transitions, but where they are they are provided
#' according to Salem et al. (1974). Similar information is summarised in Kaye and Laby (1995).
#'
#' @references
#' National Institute of Standards and Technology (NIST): X-Ray Transition Energies Database.
#' Retrieved August 2018. \url{https://physics.nist.gov/PhysRefData/XrayTrans/Html/search.html}.
#'
#' Kaye, G. W. C., and T. H. Laby. Tables of Physical and Chemical Constants and Some Mathematical Functions.
#' 16th edition. Essex, England; New York: Longman Sc & Tech, 1995. Table 4.2.1:
#' Kaye and Laby Tables of Physical and Chemical Constants (NPL).
#'
#' Salem, S.I., S.L. Panossian, and R.A. Krause. "Experimental K and L Relative X-Ray Emission Rates."
#' Atomic Data and Nuclear Data Tables 14, no. 2 (August 1974): 91–109.
#' \doi{10.1016/S0092-640X(74)80017-3}.
#'
"x_ray_energies"

#' X-Ray transitions
#'
#' A data frame keeping track of the formal and siegbahn names of various transitions
#'
"x_ray_transitions"


#' EADL97 Atomic data
#'
#' Data from the EADL97 database.
#'
#' @source
#' The PyMca5 python package (\url{https://pypi.org/project/PyMca5/#files}),
#' which obtained the data from the EADL page (\url{https://nds.iaea.org/epdl97/libsall.htm})
#'
#' @references
#' D.E. Cullen, et al., "Tables and Graphs of Atomic Subshell and
#' Relaxation Data Derived from the LLNL Evaluated Atomic Data Library
#' (EADL), Z = 1 - 100," Lawrence Livermore National Laboratory, UCRL-50400,
#' Vol. 30, October 1991
#'
"x_ray_cross_sections"

#' Photon mass attenuation coefficients
#'
#' Per-element mass attenuation coefficients \eqn{\mu/\rho} (cm^2/g) from the EPDL97 library, on the
#' native energy grid (~5 eV to 200 keV): \code{photoelectric}, \code{rayleigh} (coherent),
#' \code{compton} (incoherent) and \code{total}. Used by \link{xrf_mass_attenuation},
#' \link{xrf_detector_efficiency} and \link{xrf_escape_peaks}.
#'
#' @source
#' EPDL97 (\url{https://nds.iaea.org/epdl97/libsall.htm}), via the specfile conversion vendored
#' in \code{data-raw/EPDL97_CrossSections.dat}.
#'
"x_ray_mass_attenuation"

#' Bote-Salvat electron-impact ionization coefficients
#'
#' Fitting coefficients for the Bote & Salvat (2008) inner-shell electron-impact ionization cross
#' sections, per element (Z = 1-99) and subshell (1 = K, 2 = L1, ... 9 = M5): the tabulated edge
#' energy (\code{edge_ev}) and the coefficients \code{be}, \code{anlj}, \code{g1}-\code{g4} and
#' \code{a1}-\code{a5}. Consumed by \link{xrf_electron_ionization_cross_section} with
#' \code{method = "bote-salvat"}.
#'
#' @source
#' NIST BoteSalvatICX.jl (\url{https://github.com/usnistgov/BoteSalvatICX.jl}), public domain.
#'
#' @references
#' Bote, D. & Salvat, F. (2008). Calculations of inner-shell ionization by electron impact with the
#' distorted-wave and plane-wave Born approximations. Phys. Rev. A 77, 042701.
#'
"x_ray_bote_salvat"

#' @rdname x_ray_cross_sections
"x_ray_fluorescence_yields"

#' @rdname x_ray_cross_sections
"x_ray_emission_probabilities"

#' @rdname x_ray_cross_sections
"x_ray_coster_kronig_probabilities"

#' Coster-Kronig transition probabilities (modern Elam 2002 source)
#'
#' A drop-in modern alternative to \link{x_ray_coster_kronig_probabilities}, parsed from the Elam, Ravel &
#' Sieber (2002) database (via the XrayDB project); identical schema. Opt in with
#' \code{xrf_set_ck_source("Elam")} (or \code{options(xrftools.ck_source = "Elam")}) --
#' \link{xrf_coster_kronig_probability} and the L/M cascades then read these values. Differs from EADL97 most
#' for the L1 transitions (up to ~40 percent for heavy elements), reshaping heavy-element L-line ratios.
#'
#' @source Elam, W.T., Ravel, B.D., Sieber, J.R. (2002). \emph{Radiation Physics and Chemistry} 63, 121-128.
#' @seealso \link{x_ray_coster_kronig_probabilities}
"x_ray_coster_kronig_probabilities_elam"

#' XRF quantification energies
#'
#' A form of \link{x_ray_energies} in a more suitable form for XRF quantification.
#' Best accessed via \link{xrf_energies}.
#'
#' @source
#' NIST (2018), EADL97. See \link{x_ray_energies} and \link{x_ray_emission_probabilities}.
#'
#' @references
#' National Institute of Standards and Technology (NIST): X-Ray Transition Energies Database.
#' Retrieved August 2018. \url{https://physics.nist.gov/PhysRefData/XrayTrans/Html/search.html}.
#'
"x_ray_xrf_energies"

#' Read XRF example spectra
#'
#' @param .dir The subdirectory from which to read
#' @param ... Used to \link[dplyr]{filter} spectra
#' @param .which Used to subset files before they are read. Use TRUE for all.
#'
#' @return A spectra tibble
#' @export
#'
#' @examples
#' read_xrf_example(.which = 1:10)
#' read_xrf_example(SampleIdent == "oreas 22d")
#'
read_xrf_example <- function(..., .dir = c("Panalytical"), .which = TRUE) {
  .dir <- match.arg(.dir)
  example_dir <- system.file("spectra_files", package = "xrftools")
  if(.dir == "Panalytical") {
    subdir <- file.path(example_dir, .dir)
    spectra <- read_xrf_panalytical(
      list.files(subdir, "\\.mp2$", full.names = TRUE)[.which]
    )
  }

  dplyr::filter(spectra, ...)
}

#' Atomic form factors and incoherent scattering functions
#'
#' Per-element grids of the atomic form factor \eqn{F(q, Z)} (coherent/Rayleigh scattering) and the
#' incoherent scattering function \eqn{S(q, Z)} (Compton), versus the momentum-transfer variable
#' \eqn{q = \sin(\theta/2)/\lambda = E[\mathrm{keV}]\,\sin(\theta/2)/12.39842} (1/angstrom).
#' \eqn{F(0, Z) = Z} and falls steeply with \eqn{q} (high energy / large angle cannot scatter
#' coherently); \eqn{S(0, Z) = 0} (binding suppression) and saturates to \eqn{Z}. These turn the
#' angle-\emph{integrated} scatter cross sections into exact \emph{fixed-angle} differentials:
#' \deqn{d\sigma_R/d\Omega \propto (1 + \cos^2\theta)\, F^2(q, Z), \qquad
#'       d\sigma_C/d\Omega \propto KN(E, \theta)\, S(q, Z),}
#' which shape the Rayleigh/Compton templates of \link{xrf_scatter_peaks} and
#' \link{xrf_scatter_continuum} at the instrument's actual scattering angle.
#'
#' @source Vendored from the xraylib project (Schoonjans et al.; BSD-3-Clause,
#'   \url{https://github.com/tschoonj/xraylib}), which derives them from EPDL97.
#' @references
#' Hubbell, J.H., Veigele, W.J., Briggs, E.A., Brown, R.T., Cromer, D.T., Howerton, R.J. (1975).
#' Atomic form factors, incoherent scattering functions, and photon scattering cross sections.
#' \emph{Journal of Physical and Chemical Reference Data} 4, 471-538. \doi{10.1063/1.555523}
#'
#' Schoonjans, T., Brunetti, A., Golosio, B., Sanchez del Rio, M., Sole, V.A., Ferrero, C.,
#' Vincze, L. (2011). The xraylib library for X-ray-matter interactions. \emph{Spectrochimica Acta
#' Part B} 66, 776-784. \doi{10.1016/j.sab.2011.09.011}
"x_ray_form_factors"

#' @rdname x_ray_form_factors
"x_ray_incoherent_functions"

#' Compton profiles of the elements
#'
#' The total electron momentum density \eqn{J(p_z)} per element (Biggs, Mendelsohn & Mann 1975,
#' Hartree-Fock; \eqn{p_z} in atomic units, \eqn{J} symmetric with \eqn{\int J\,dp_z = Z}). In the
#' impulse approximation the Compton-scattered energy distribution of a line \eqn{E_0} at angle
#' \eqn{\theta} is \eqn{\rho(E_2) \propto J(|p_z(E_2)|)\,|dp_z/dE_2|} with the Ribberfors relation
#' \deqn{p_z = 137.036\,\frac{E_0 - E_2 - E_0 E_2 (1-\cos\theta)/m c^2}
#'   {\sqrt{E_0^2 + E_2^2 - 2 E_0 E_2 \cos\theta}},}
#' which is the physical, asymmetric Doppler broadening of the Compton peak -- used by
#' \link{xrf_scatter_peaks} (\code{compton_shape = "profile"}) in place of a broadened Gaussian.
#'
#' @source DABAX copy vendored via the xraylib project (BSD-3-Clause,
#'   \url{https://github.com/tschoonj/xraylib}).
#' @references
#' Biggs, F., Mendelsohn, L.B., Mann, J.B. (1975). Hartree-Fock Compton profiles for the elements.
#' \emph{Atomic Data and Nuclear Data Tables} 16, 201-309. \doi{10.1016/0092-640X(75)90030-3}
#'
#' Ribberfors, R. (1975). Relationship of the relativistic Compton cross section to the momentum
#' distribution of bound electron states. \emph{Physical Review B} 12, 2067-2074.
#' \doi{10.1103/PhysRevB.12.2067}
"x_ray_compton_profiles"

#' @rdname x_ray_compton_profiles
#' @details \code{x_ray_compton_profiles_shells} carries the per-SUBSHELL profiles (per electron;
#'   the occupancy-weighted sum reproduces the total) together with each subshell's occupancy and
#'   binding energy. A subshell contributes to Compton scattering only where the energy transfer
#'   exceeds its binding energy -- the constraint that displaces core-shell contributions to the
#'   low-energy side of the Compton line and produces the Compton defect. One transcription error in
#'   the DABAX source (Eu subshell 3 at \eqn{p_z = 30}) is repaired at build time by log-interpolation
#'   (see \code{data-raw/compton_profiles.R}).
"x_ray_compton_profiles_shells"
