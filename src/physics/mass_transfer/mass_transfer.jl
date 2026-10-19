"""
    AbstractMassTransfer

Abstract base type for mass transfer models.

Implementations must provide:
- `compute_mass_transfer_rate(mt, C, q, qstar)`: Calculate mass transfer rate
"""
abstract type AbstractMassTransfer end

"""
    compute_mass_transfer_rate(mt::AbstractMassTransfer, C, q, qstar) → rate

Compute mass transfer rate [mol/(kg·s)] for the adsorption process.
"""
function compute_mass_transfer_rate end

include("ldf.jl")
