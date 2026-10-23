"""
    LinearDrivingForce{T, D} <: AbstractMassTransfer

Linear Driving Force (LDF) mass transfer model with single lumped resistance.

The rate equation is:
```
∂q/∂t = k * (q* - q)
```
where `k = 15 * ε_p * D_p / r_p² * C / q*` and `D_p = D_m / τ`. A component the
sorbent does not take up (`q* = 0` at any concentration, such as one without
isotherm data) has no mass transfer.
"""
struct LinearDrivingForce{T, D} <: AbstractMassTransfer
    D_m::D
    τ::T
    ϵ_p::T
    d_p::T
end

"""
    LinearDrivingForce(; D_m, τ, ϵ_p, d_p)

Construct a Linear Driving Force mass transfer model from named parameters.

- `D_m`: molecular diffusivity [m²/s], one value for all components or one per component
- `τ`: tortuosity [-]
- `ϵ_p`: particle porosity [-]
- `d_p`: particle diameter [m]
"""
LinearDrivingForce(; D_m, τ, ϵ_p, d_p) = LinearDrivingForce(D_m, τ, ϵ_p, d_p)


function compute_mass_transfer_rate(mt::LinearDrivingForce, C, q, qstar)
    D_p = mt.D_m ./ mt.τ
    r_p = mt.d_p / 2.0
    k = _capacity_ratio.(C, qstar) .* 15 .* mt.ϵ_p .* D_p ./ r_p^2
    return k .* (qstar .- q)
end

# Macropore diffusion is molecular diffusion through the pores
compute_mass_transfer_rate(mt::LinearDrivingForce, C, q, qstar, diffusivity_factor) =
    compute_mass_transfer_rate(LinearDrivingForce(mt.D_m .* diffusivity_factor, mt.τ, mt.ϵ_p, mt.d_p), C, q, qstar)

# C/q*, or zero where the sorbent takes up none of a component, where C/q*
# would be C/0. ifelse evaluates both sides, so the derivatives keep the same
# pattern whichever applies.
@inline _capacity_ratio(c, qstar) = ifelse(qstar > 0, c / qstar, zero(c / one(qstar)))
