"""
    AbstractMassTransfer

Abstract base type for mass transfer models.

Implementations provide `compute_mass_transfer_rate(mt, C, q, qstar, T)`, or,
for a model that does not depend on temperature,
`compute_mass_transfer_rate(mt, C, q, qstar)`; and, if the rate depends on
the molecular diffusivity of the gas,
`compute_mass_transfer_rate(mt, C, q, qstar, T, diffusivity_factors)`, with
the diffusivity of each component scaled by `diffusivity_factors` (see
[`diffusivity_factors`](@ref)).
"""
abstract type AbstractMassTransfer end

"""
    compute_mass_transfer_rate(mt::AbstractMassTransfer, C, q, qstar, T) → rate

Rate of change of the loading of each component [mol/(m³ adsorbent·s)], for
gas concentrations `C` [mol/m³], loadings `q` and equilibrium loadings `qstar`
[mol/m³ adsorbent], at temperature `T` [K]. The total pressure is
`sum(C)·R·T`.

A model whose rate does not depend on temperature may define the method
without `T`, which this one falls back to.
"""
function compute_mass_transfer_rate end

compute_mass_transfer_rate(mt::AbstractMassTransfer, C, q, qstar, T) = compute_mass_transfer_rate(mt, C, q, qstar)

"""
    compute_mass_transfer_rate(mt, C, q, qstar, T, diffusivity_factors) → rate

The rate with the molecular diffusivity of each component of the gas its
`diffusivity_factors` times that of `mt`. By default the rate does not depend on it.
"""
compute_mass_transfer_rate(mt::AbstractMassTransfer, C, q, qstar, T, diffusivity_factors) = compute_mass_transfer_rate(mt, C, q, qstar, T)
