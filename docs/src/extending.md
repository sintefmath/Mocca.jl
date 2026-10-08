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
cycle_change
cyclic_steady_state_gradient
force_gradients
StreamObjective
VacuumPumpObjective
CompressorObjective
SimulationFailure
```
