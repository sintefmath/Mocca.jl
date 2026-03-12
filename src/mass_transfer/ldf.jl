"""
    LinearDrivingForce{T} <: AbstractMassTransfer

Linear Driving Force (LDF) mass transfer model with single lumped resistance.

The rate equation is:
```
∂q/∂t = k * (q* - q)
```
where `k = 15 * ε_p * D_p / r_p² * C / q*` and `D_p = D_m / τ`.
"""
struct LinearDrivingForce{T} <: AbstractMassTransfer
    D_m::T
    τ::T
    ϵ_p::T
    d_p::T
end


function compute_mass_transfer_rate(mt::LinearDrivingForce, C, q, qstar)
    D_p = mt.D_m / mt.τ
    r_p = mt.d_p / 2.0
    k = C ./ qstar .* 15 .* mt.ϵ_p .* D_p ./ r_p^2
    return k .* (qstar .- q)
end
