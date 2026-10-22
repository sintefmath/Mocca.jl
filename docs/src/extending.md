```@meta
CurrentModule = Mocca
```

# Extending Mocca

Most studies need only the [Reference](@ref) functions. This page is for adding
your own kinds of step and unit, and for working with the simulation model
that a [`Process`](@ref) builds. That model is a
[Jutul](https://github.com/sintefmath/Jutul.jl) `MultiModel`, with one model
per bed and device, and the steps of the cycle become its forces.

## Adding a kind of step

The set of steps is open: define a subtype of [`AbstractStep`](@ref) and a
method of [`step_settings`](@ref).

```@docs
AbstractStep
StepContext
compile_cycle
step_settings
end_pressure
step_beds
default_step_name
register_step!
find_device
valve
boundary
start_pressure
ramp
```

## Adding a thermodynamics model

The heat capacities and the equation of state of the gas come from an
[`AbstractGasThermodynamics`](@ref) model, held by the [`Gas`](@ref) and by
every unit with a gas phase (the `thermodynamics` field of its system; a bed
also has `adsorbed_thermodynamics` for the adsorbed phase). Mocca has
[`IdealGasConstantCp`](@ref) and [`IdealGasShomate`](@ref). Another model,
such as NASA polynomials, is a subtype with one method,
[`component_heat_capacities`](@ref): the molar heat capacity of each
component [J/(mol·K)] at a temperature, as an `SVector`.

```julia
using Mocca, StaticArrays

# NASA 7-coefficient polynomials, one temperature range: cp/R = a1 + a2 T + ... + a5 T⁴
struct NASA7{N} <: Mocca.AbstractGasThermodynamics
    a::NTuple{5, SVector{N, Float64}}
end

function Mocca.component_heat_capacities(m::NASA7, T)
    a = m.a
    return Mocca.GAS_CONSTANT .* (a[1] .+ T .* (a[2] .+ T .* (a[3] .+ T .* (a[4] .+ T .* a[5]))))
end

gas = Gas(["CO2", "N2"]; thermodynamics = NASA7((SVector(2.36, 3.53), SVector(8.98e-3, -1.24e-4),
    SVector(-7.12e-6, -5.03e-7), SVector(2.46e-9, 2.44e-9), SVector(-1.44e-13, -1.41e-12))))
```

The rest follows:
[`component_enthalpies`](@ref), relative to 298.15 K, by Gauss–Legendre
quadrature of the heat capacities (add a method if the integral is known),
and from these [`molar_heat_capacity`](@ref), [`molar_enthalpy`](@ref) and
[`heat_capacity_ratio`](@ref) of a mixture. [`molar_concentration`](@ref) is
that of the ideal gas; the energy balance of the column also assumes an ideal
gas, so a model with another equation of state needs more than a method for
it.

The functions are evaluated in every cell at every Newton iteration with
automatic-differentiation numbers for `T`, so they must accept any `Real`, and
should not allocate: keep the parameters in `SVector`s or tuples, and the
struct `isbits`, as [`DualSiteLangmuir`](@ref) does.

The model's heat capacities enter the energy balance through the gas and
adsorbed-phase heat capacity variables `C_pg` and `C_pa`, the sensible heat of
sorption and the enthalpy a stream brings into a bed (see [Model
equations](@ref)). A model that is not one of Mocca's cannot be written to a
process file.

## The layer underneath

The objects of the [Reference](@ref) build these, which can also be used
directly. A
[`Flowsheet`](@ref) holds a [`UnitNetwork`](@ref): the unit models and the
connections between their ports, from which the simulation model is built.

```@docs
fixed_bed_model
MoccaSystem
DistributedUnit
LumpedUnit
FlowElement
component_names
AbstractFlowsheet
UnitNetwork
network
add_unit!
connect!
set_boundary!
PortRef
AbstractConnection
connection_ports
check_connection
add_cross_terms!
StreamConnection
holdup_units
```

## Adding a kind of unit

A new kind of unit is a `MoccaSystem` with these methods.

```@docs
stream_equations
setup_unit_state
setup_unit_parameters
holds_sorbent
AxialDispersion
```

## Devices, schedules and steady states

```@docs
FlowDevice
setup_device_model
AbstractDeviceLaw
device_flow
LinearValve
VolumetricFlow
PortStateCT
StreamCT
Stage
setup_schedule
parallel_stages
PortCondition
ExponentialRamp
port_conductance
default_conductance
setup_flowsheet_model
setup_flowsheet_state
setup_flowsheet_parameters
DualSiteLangmuir
AbstractIsotherm
AbstractMassTransfer
LinearDrivingForce
simulate_to_cyclic_steady_state
newton_cyclic_steady_state
simulate_until_mass_balance
CycleProgress
cycle_change
cyclic_steady_state_gradient
force_gradients
StreamObjective
VacuumPumpObjective
CompressorObjective
SimulationFailure
```
