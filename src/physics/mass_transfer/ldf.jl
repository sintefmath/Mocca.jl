"""
    LinearDrivingForce{T, D, M} <: AbstractMassTransfer

Linear Driving Force (LDF) mass transfer model with single lumped resistance.

The rate equation is:
```
∂q/∂t = k * (q* - q)
```
where `k = 15 * ε_p * D_p / r_p² * C / q*` and `D_p = D / τ`. A component the
sorbent does not take up (`q* = 0` at any concentration, such as one without
isotherm data) has no mass transfer.

`D` is the molecular diffusivity `D_m`, or, with a macropore radius
`r_pore`, molecular and Knudsen diffusion in series (Bosanquet):
`1/D = 1/D_m + 1/D_K`, with `D_K = 2/3·r_pore·√(8RT/(πM))` for a component
of molar mass `M` at temperature `T`.
"""
struct LinearDrivingForce{T, D, M} <: AbstractMassTransfer
    D_m::D
    τ::T
    ϵ_p::T
    d_p::T
    r_pore::T
    molar_masses::M
end

"""
    LinearDrivingForce(D_m, τ, ϵ_p, d_p; macropore_radius = 0.0, molar_masses = nothing)
    LinearDrivingForce(; D_m, τ, ϵ_p, d_p, macropore_radius = 0.0, molar_masses = nothing)

Construct a Linear Driving Force mass transfer model.

- `D_m`: molecular diffusivity [m²/s], one value for all components or one per component
- `τ`: tortuosity [-]
- `ϵ_p`: particle porosity [-]
- `d_p`: particle diameter [m]
- `macropore_radius` [m]: with a positive radius, Knudsen diffusion in the
  macropores adds to the resistance; it needs `molar_masses` [kg/mol], one per
  component. 0 leaves molecular diffusion alone.
"""
function LinearDrivingForce(D_m, τ, ϵ_p, d_p; macropore_radius = 0.0, molar_masses = nothing)
    T = promote_type(typeof(τ), typeof(ϵ_p), typeof(d_p), typeof(macropore_radius))
    macropore_radius >= 0 || error("macropore_radius must not be negative, got $macropore_radius")
    M = nothing
    if macropore_radius > 0
        isnothing(molar_masses) && error("Knudsen diffusion (macropore_radius > 0) needs the molar_masses of the components")
        D_m isa Number || length(D_m) == length(molar_masses) ||
            error("molar_masses has $(length(molar_masses)) values for $(length(D_m)) diffusivities")
        M = SVector{length(molar_masses), Float64}(molar_masses)
    end
    return LinearDrivingForce{T, typeof(D_m), typeof(M)}(D_m, τ, ϵ_p, d_p, macropore_radius, M)
end
LinearDrivingForce(; D_m, τ, ϵ_p, d_p, kwarg...) = LinearDrivingForce(D_m, τ, ϵ_p, d_p; kwarg...)

"Whether `mt` includes Knudsen diffusion in the macropores."
has_knudsen(mt::LinearDrivingForce) = !isnothing(mt.molar_masses)

"""
    pore_diffusivity(mt::LinearDrivingForce, T, diffusivity_factors = 1) → D

The diffusivity of each component in the macropores at temperature `T` [K],
before dividing by the tortuosity: the molecular diffusivity, times
`diffusivity_factors` where it follows the state of the gas, in series with
Knudsen diffusion if `mt` has a macropore radius.
"""
function pore_diffusivity(mt::LinearDrivingForce, T, diffusivity_factors = 1)
    D_m = mt.D_m .* diffusivity_factors
    return has_knudsen(mt) ? _bosanquet.(D_m, _knudsen_diffusivity.(mt.r_pore, mt.molar_masses, T)) : D_m
end

_knudsen_diffusivity(r, M, T) = 2 / 3 * r * sqrt(8 * GAS_CONSTANT * T / (π * M))
_bosanquet(D_m, D_K) = 1 / (1 / D_m + 1 / D_K)

function compute_mass_transfer_rate(mt::LinearDrivingForce, C, q, qstar, T)
    return _ldf_rate(mt, C, q, qstar, pore_diffusivity(mt, T))
end

function compute_mass_transfer_rate(mt::LinearDrivingForce, C, q, qstar)
    has_knudsen(mt) && error("Knudsen diffusion depends on temperature; call compute_mass_transfer_rate(mt, C, q, qstar, T)")
    return _ldf_rate(mt, C, q, qstar, mt.D_m)
end

function _ldf_rate(mt::LinearDrivingForce, C, q, qstar, D)
    D_p = D ./ mt.τ
    r_p = mt.d_p / 2.0
    k = _capacity_ratio.(C, qstar) .* 15 .* mt.ϵ_p .* D_p ./ r_p^2
    return k .* (qstar .- q)
end

# Molecular diffusion through the pores follows the state of the gas; Knudsen
# diffusion, if any, only the temperature
compute_mass_transfer_rate(mt::LinearDrivingForce, C, q, qstar, T, diffusivity_factors) =
    _ldf_rate(mt, C, q, qstar, pore_diffusivity(mt, T, diffusivity_factors))

# C/q*, or zero where the sorbent takes up none of a component, where C/q*
# would be C/0. ifelse evaluates both sides, so the derivatives keep the same
# pattern whichever applies.
@inline _capacity_ratio(c, qstar) = ifelse(qstar > 0, c / qstar, zero(c / one(qstar)))
