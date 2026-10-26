"""
    AbstractGasTransport

Abstract base type for models of the transport properties of a gas mixture:
its viscosity and the molecular diffusivity of each component, as they change
with temperature, pressure and composition.

A unit holds such a model with its values at the state at which the unit's
viscosity and diffusivities are given (see [`TransportReference`](@ref)), and
scales those by the ratio of the model's value at the state of a cell to its
value there. The given values keep their meaning, and stay what design
variables and gradients act on.

Implementations must provide:
- `mixture_viscosity(model, y, T)` [Pa·s];
- `molecular_diffusivities(model, y, T, P)` [m²/s]: the diffusivity of each
  component in the mixture, as an `SVector`.

As for [`AbstractGasThermodynamics`](@ref), each must work for
automatic-differentiation numbers and not allocate.
"""
abstract type AbstractGasTransport end

"""
    ConstantTransport() <: AbstractGasTransport

Transport properties that do not change: each unit uses its viscosity and
diffusivities as given. The default.
"""
struct ConstantTransport <: AbstractGasTransport end

"""
    mixture_viscosity(model::AbstractGasTransport, y, T)

Viscosity [Pa·s] of the gas mixture of mole fractions `y` at temperature `T` [K].
"""
function mixture_viscosity end

"""
    molecular_diffusivities(model::AbstractGasTransport, y, T, P) → SVector

Molecular diffusivity [m²/s] of each component in the gas mixture of mole
fractions `y` at temperature `T` [K] and pressure `P` [Pa].
"""
function molecular_diffusivities end

"""
    TransportReference(model, y, T, P)
    TransportReference(ConstantTransport())

A transport model with its values at the reference state, mole fractions `y`,
temperature `T` [K] and pressure `P` [Pa]: the state at which a unit's
viscosity and diffusivities are given, the feed for the beds the interface
builds. [`viscosity_factor`](@ref) and [`diffusivity_factors`](@ref) give the
ratio of the model's values at another state to those at the reference.
"""
struct TransportReference{M<:AbstractGasTransport, T<:Real, D}
    model::M
    viscosity::T            # At the reference [Pa·s]
    diffusivities::D        # Of each component at the reference [m²/s]
    pressure::T             # Of the reference [Pa]
end

function TransportReference(model::AbstractGasTransport, y::AbstractVector, T, P)
    D = molecular_diffusivities(model, y, T, P)
    return TransportReference(model, float(mixture_viscosity(model, y, T)), float.(D), float(P))
end
TransportReference(model::ConstantTransport, y::AbstractVector, T, P) = TransportReference(model)
TransportReference(model::ConstantTransport) = TransportReference(model, 1.0, nothing, 1.0)

"Whether the transport properties of `r` stay at their reference values."
is_constant_transport(r::TransportReference) = r.model isa ConstantTransport

"""
    viscosity_factor(r::TransportReference, y, T)

Viscosity at mole fractions `y` and temperature `T` [K] relative to that at
the reference of `r`.
"""
viscosity_factor(r::TransportReference, y, T) = mixture_viscosity(r.model, y, T) / r.viscosity
viscosity_factor(::TransportReference{ConstantTransport}, y, T) = one(T)

"""
    diffusivity_factors(r::TransportReference, y, T, P) → SVector

Molecular diffusivity of each component at mole fractions `y`, temperature
`T` [K] and pressure `P` [Pa] relative to that at the reference of `r`.
"""
diffusivity_factors(r::TransportReference, y, T, P) = molecular_diffusivities(r.model, y, T, P) ./ r.diffusivities
