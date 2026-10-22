# Packed-bed correlations for the parameters of a fixed bed.

function compute_permeability(porosity, d_p)
    return 4 / 150 * ((porosity / (1 - porosity))^2) * (d_p / 2)^2 * porosity
end

"""
    compute_dispersion(D_m, v, d_p)

Axial dispersion coefficient [m²/s] of a packed bed, `0.7 D_m + 0.5 v d_p`,
for molecular diffusivity `D_m` [m²/s], interstitial velocity `v` [m/s] and
particle diameter `d_p` [m].
"""
function compute_dispersion(D_m, V0_inter, d_p)
    return 0.7 * D_m + 0.5 * V0_inter * d_p
end

"""
    compute_radial_dispersion(D_m, v, d_p)

Radial dispersion coefficient [m²/s] of a packed bed, `0.7 D_m + 0.1 v d_p`,
for molecular diffusivity `D_m` [m²/s], interstitial velocity `v` [m/s] and
particle diameter `d_p` [m]: the form of [`compute_dispersion`](@ref) with a
radial Péclet number `v d_p / D` of 10 at high velocity, where along the axis
it is 2.
"""
function compute_radial_dispersion(D_m, v, d_p)
    return 0.7 * D_m + 0.1 * v * d_p
end
