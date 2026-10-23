"""
    AbstractGasTransport

Abstract base type for models of the transport properties of a gas mixture:
its viscosity, thermal conductivity and molecular diffusivity, as they change
with temperature, pressure and composition.

A unit holds such a model with the state at which its viscosity, thermal
conductivity and diffusivity parameters hold (see
[`TransportReference`](@ref)), and scales each parameter by the ratio of the
model's value at the state of a cell to its value at that reference. The
parameters keep their meaning, and stay what design variables and gradients
act on.

Implementations must provide:
- `mixture_viscosity(model, y, T)` [Pa·s];
- `mixture_thermal_conductivity(model, y, T)` [W/(m·K)];
- `diffusivity_scale(model, T, P)`: proportional to the molecular diffusivity
  at temperature `T` [K] and pressure `P` [Pa], in any normalisation, as only
  its ratios are used.

As for [`AbstractGasThermodynamics`](@ref), each must work for
automatic-differentiation numbers and not allocate.
"""
abstract type AbstractGasTransport end

"""
    ConstantTransport() <: AbstractGasTransport

Transport properties that do not change: each unit uses its viscosity,
thermal conductivity and diffusivity parameters as they are. The default.
"""
struct ConstantTransport <: AbstractGasTransport end

"""
    mixture_viscosity(model::AbstractGasTransport, y, T)

Viscosity [Pa·s] of the gas mixture of mole fractions `y` at temperature `T` [K].
"""
function mixture_viscosity end

"""
    mixture_thermal_conductivity(model::AbstractGasTransport, y, T)

Thermal conductivity [W/(m·K)] of the gas mixture of mole fractions `y` at
temperature `T` [K].
"""
function mixture_thermal_conductivity end

"""
    diffusivity_scale(model::AbstractGasTransport, T, P)

A number proportional to the molecular diffusivity at temperature `T` [K] and
pressure `P` [Pa]; only its ratios are used.
"""
function diffusivity_scale end

"""
    TransportReference(model, y, T, P)
    TransportReference(ConstantTransport())

A transport model with its values at the reference state, mole fractions `y`,
temperature `T` [K] and pressure `P` [Pa]: the state at which a unit's
viscosity, thermal conductivity and diffusivity parameters hold, the feed for
the beds the interface builds. [`viscosity_factor`](@ref),
[`conductivity_factor`](@ref) and [`diffusivity_factor`](@ref) give the ratio
of the model's value at another state to that at the reference.
"""
struct TransportReference{M<:AbstractGasTransport, T<:Real}
    model::M
    viscosity::T                # At the reference [Pa·s]
    thermal_conductivity::T     # At the reference [W/(m·K)]
    diffusivity_scale::T        # diffusivity_scale at the reference
end

TransportReference(model::AbstractGasTransport, y::AbstractVector, T, P) = TransportReference(model,
    float(mixture_viscosity(model, y, T)), float(mixture_thermal_conductivity(model, y, T)), float(diffusivity_scale(model, T, P)))
TransportReference(model::ConstantTransport, y::AbstractVector, T, P) = TransportReference(model)
TransportReference(model::ConstantTransport) = TransportReference(model, 1.0, 1.0, 1.0)

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
    conductivity_factor(r::TransportReference, y, T)

Gas thermal conductivity at mole fractions `y` and temperature `T` [K]
relative to that at the reference of `r`.
"""
conductivity_factor(r::TransportReference, y, T) = mixture_thermal_conductivity(r.model, y, T) / r.thermal_conductivity
conductivity_factor(::TransportReference{ConstantTransport}, y, T) = one(T)

"""
    diffusivity_factor(r::TransportReference, T, P)

Molecular diffusivity at temperature `T` [K] and pressure `P` [Pa] relative to
that at the reference of `r`.
"""
diffusivity_factor(r::TransportReference, T, P) = diffusivity_scale(r.model, T, P) / r.diffusivity_scale
diffusivity_factor(::TransportReference{ConstantTransport}, T, P) = one(T / P)
