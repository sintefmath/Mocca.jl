"""
    AbstractMassTransfer

Abstract base type for mass transfer models.

Implementations must provide:
- `compute_mass_transfer_rate(mt, C, q, qstar)`: Calculate mass transfer rate

and, if the rate depends on the molecular diffusivity of the gas,
- `compute_mass_transfer_rate(mt, C, q, qstar, diffusivity_factors)`, with the
  diffusivity of each component scaled by `diffusivity_factors` (see
  [`diffusivity_factors`](@ref)).
"""
abstract type AbstractMassTransfer end

"""
    compute_mass_transfer_rate(mt::AbstractMassTransfer, C, q, qstar) → rate

Compute mass transfer rate [mol/(m³ adsorbent·s)] for the adsorption process.
"""
function compute_mass_transfer_rate end


"""
    compute_mass_transfer_rate(mt, C, q, qstar, diffusivity_factors) → rate

The rate with the molecular diffusivity of each component of the gas its
`diffusivity_factors` times that of `mt`. By default the rate does not depend on it.
"""
compute_mass_transfer_rate(mt::AbstractMassTransfer, C, q, qstar, diffusivity_factors) = compute_mass_transfer_rate(mt, C, q, qstar)
