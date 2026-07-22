FF.dat (atomic form factors F(q,Z)) and SF.dat (incoherent scattering functions S(q,Z)) are vendored
verbatim from the xraylib project (https://github.com/tschoonj/xraylib, data/ directory), BSD-3-Clause,
Copyright (c) Tom Schoonjans and contributors. Underlying data: EPDL97 / Hubbell et al. (1975),
J. Phys. Chem. Ref. Data 4, 471-538. Format per element: one line with the grid length N, then N rows of
"q[1/angstrom]  value  spline-second-derivative" (the third column is xraylib's spline helper, unused here).
q is the momentum-transfer variable sin(theta/2)/lambda = E[keV] * sin(theta/2) / 12.39842.
Parsed by data-raw/scatter_factors.R into x_ray_form_factors / x_ray_incoherent_functions.
ComptonProfiles.dat (Biggs J(p_z)) exists in the same xraylib directory -- the natural NEXT vendoring
step for exact Doppler-broadened Compton line shapes.

ComptonProfiles_biggs.dat is the DABAX/xraylib copy (data/comptonprofiles/ComptonProfiles.dat) of the
Biggs, Mendelsohn & Mann (1975) Hartree-Fock Compton profiles (At. Data Nucl. Data Tables 16, 201-309),
Z = 1-102: per element a #S block with columns pz [atomic units], total J(pz), then per-subshell J's.
Parsed (total only) by data-raw/compton_profiles.R into x_ray_compton_profiles.
