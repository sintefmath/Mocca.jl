"""
    AbstractSourceTerm

A contribution to one of the conservation laws of a distributed unit that is
not a flux: sorption uptake, reaction, heat of sorption, wall exchange, ...

Each unit type lists the terms it needs, and the generic residuals sum them:

    R_i = ∂M_i/∂t + ∇·F_i − Σ_s source(s, model, state, state0, cell, i, Δt)

The term tuples are known at compile time, so an empty tuple costs nothing.

A term implements whichever of these it contributes to:
- [`mass_source`](@ref) for the `:ComponentMasses` law
- [`energy_source`](@ref) for the `:ColumnConservedEnergy` law
- [`wall_source`](@ref) for the `:WallConservedEnergy` law
- [`heat_capacity`](@ref) for the effective heat capacity multiplying `∂T/∂t`
"""
abstract type AbstractSourceTerm end

"""
    mass_source_terms(sys) → Tuple

Source terms entering the per-component gas mass balance of `sys`.
"""
mass_source_terms(::MoccaSystem) = ()

"""
    energy_source_terms(sys) → Tuple

Source terms entering the gas-phase energy balance of `sys`.
"""
energy_source_terms(::MoccaSystem) = ()

"""
    wall_source_terms(sys) → Tuple

Source terms entering the wall energy balance of `sys`.
"""
wall_source_terms(::MoccaSystem) = ()

"""
    heat_capacity_terms(sys) → Tuple

Terms whose [`heat_capacity`](@ref) values add up to the effective heat
capacity [J/K] multiplying `∂T/∂t` in the gas-phase energy balance, apart from
the gas itself.
"""
heat_capacity_terms(::MoccaSystem) = ()

"""
    mass_source(term, model, state, state0, cell, component, Δt)

Molar source [mol/s] of `component` in `cell`. Positive adds gas.
"""
function mass_source end

"""
    energy_source(term, model, state, state0, cell, Δt)

Heat source [W] in `cell` for the gas-phase energy balance. Positive heats the gas.
"""
function energy_source end

"""
    wall_source(term, model, state, state0, cell, Δt)

Heat source [W] in wall `cell`. Positive heats the wall.
"""
function wall_source end

"""
    heat_capacity(term, model, state, cell)

Heat capacity [J/K] contributed by `term` in `cell`.
"""
function heat_capacity end

# Sum `f(term, args...)` over a tuple of terms. Recursion over the tuple keeps
# this type-stable; the empty sum is `false`, the additive identity for any
# number type (including AD types).
@inline sum_terms(f, ::Tuple{}, args...) = false
@inline sum_terms(f, terms::Tuple, args...) = f(first(terms), args...) + sum_terms(f, Base.tail(terms), args...)
