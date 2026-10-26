"""
    AbstractAxialDispersion

Abstract base type for models of the axial dispersion of the gas in a bed.

- [`ConstantDispersion`](@ref): the unit's `AxialDispersion` parameter, the
  same for every component, everywhere and at all times;
- [`WakaoFunazkri`](@ref): from the molecular diffusivity of each component
  and the gas velocity at each face.
"""
abstract type AbstractAxialDispersion end

"""
    ConstantDispersion() <: AbstractAxialDispersion

Axial dispersion given by the unit's `AxialDispersion` parameter [m²/s].
"""
struct ConstantDispersion <: AbstractAxialDispersion end

"""
    WakaoFunazkri(; particle_diameter, porosity, diffusivities) <: AbstractAxialDispersion

Axial dispersion by the correlation of Wakao and Funazkri (1978),

    ε·D_ax,i / D_m,i = 20 + 0.5·Sc_i·Re_p,   Re_p = ρ u d_p / μ,   Sc_i = μ / (ρ D_m,i),

that is `D_ax,i = (20 D_m,i + 0.5 |u| d_p) / ε`, with `u` the superficial
velocity [m/s], `d_p` the `particle_diameter` [m] and `ε` the bed `porosity`.
`D_m,i` is the molecular diffusivity of component `i`: `diffusivities` [m²/s]
at the feed, scaled by the unit's transport model (see
[`TransportReference`](@ref)). The unit's `DispersionScale` parameter
multiplies it (1 by default).
"""
struct WakaoFunazkri{T, D} <: AbstractAxialDispersion
    d_p::T              # Particle diameter [m]
    ε::T                # Bed porosity
    diffusivities::D    # Molecular diffusivity of each component at the feed [m²/s]
end

WakaoFunazkri(; particle_diameter, porosity, diffusivities) =
    WakaoFunazkri(float(particle_diameter), float(porosity), SVector{length(diffusivities)}(float.(diffusivities)))

"""
    dispersion_coefficient(model::WakaoFunazkri, D, u)

Axial dispersion coefficient [m²/s] of a component of molecular diffusivity
`D` [m²/s] at superficial velocity `u` [m/s].
"""
dispersion_coefficient(model::WakaoFunazkri, D, u) = (20 * D + abs(u) * model.d_p / 2) / model.ε
