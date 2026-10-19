```@meta
CurrentModule = Mocca
```

# Setting up processes

A process is described in layers, from the physics to the operation. Each
layer is a plain object that can be edited, saved and reused on its own:

| Layer | Objects |
|---|---|
| Gas | [`Gas`](@ref), [`Feed`](@ref), [`Component`](@ref) |
| Materials | [`Adsorbent`](@ref) with an isotherm by component ([`dual_site_langmuir`](@ref), [`langmuir`](@ref)), [`Column`](@ref), [`Wall`](@ref) |
| Flowsheet | [`Flowsheet`](@ref), beds ([`add_bed!`](@ref)) and devices with roles ([`add_device!`](@ref)) |
| Operation | a [`Cycle`](@ref) of steps ([`Pressurise`](@ref), [`Adsorb`](@ref), [`Blowdown`](@ref), [`Evacuate`](@ref), [`Purge`](@ref), [`Equalise`](@ref), [`Idle`](@ref), [`Concurrent`](@ref), or your own), [`TimeStepping`](@ref) |
| Process | [`Process`](@ref): [`simulate`](@ref), [`cyclic_steady_state`](@ref), [`performance`](@ref) |
| Design | [`add_variable!`](@ref), [`performance_gradient`](@ref), [`optimise!`](@ref) |
| Files | [`read_process`](@ref), [`write_process`](@ref), [`load_adsorbent`](@ref); see [Process files](@ref) |

```julia
using Mocca
using Mocca.Units

gas = Gas(["CO2", "N2"])
feed = Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), pressure = 1bar)
fs = Flowsheet(gas, feed)
add_bed!(fs, :Bed, load_adsorbent("zeolite_13x_haghpanah_2013"),
    Column(length = 1.0, radius = 0.1445, outer_radius = 0.162, porosity = 0.37))
add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
add_device!(fs, :V_product, :product; at = (:Bed, :top))
add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
cycle = Cycle([Pressurise(15.0), Adsorb(15.0; velocity = 1.0),
    Blowdown(30.0; to = 0.2bar), Evacuate(40.0; to = 0.1bar)])
process = Process(fs, cycle)
performance(cyclic_steady_state(process))
```

## Units

Everything is in plain SI units: Pa, K, m, s, mol, kg and J.

```@docs
Units
Units.celsius
```

## Gas and feed

```@docs
Component
component
register_component!
Gas
composition
Feed
```

## Isotherms, adsorbents and columns

Isotherm parameters are given by component name, so an adsorbent does not
depend on the order of the components of a gas, and components can be added
or removed without editing parallel arrays.

```@docs
AbstractIsothermData
dual_site_langmuir
langmuir
build_isotherm
Adsorbent
equilibrium_loading
isotherm_table
Column
Wall
```

## Flowsheets

Devices have roles. A role sets the device's default boundary condition and
valve size, tells the steps which device to open, and tells
[`performance`](@ref) which streams are feed and product.

```@docs
add_bed!
fixed_bed_model
add_device!
DEVICE_ROLES
default_conductance
print_flowsheet_diagram
```

## Cycles and steps

A step acts on a bed and opens the devices it needs, found by their roles;
every other device stays closed. Pressure carries over from step to step: a
ramp starts from the pressure the step before left the bed at, around the
cycle, unless `from` is given.

```@docs
Cycle
Pressurise
Adsorb
Blowdown
Evacuate
Purge
Equalise
Idle
Concurrent
TimeStepping
compile_cycle
```

### Adding a kind of step

The set of steps is open: define a subtype of [`AbstractStep`](@ref) and a
method of [`step_settings`](@ref).

```@docs
AbstractStep
StepContext
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

## Processes

```@docs
Process
simulate(::Process)
initial_state
stages
cyclic_steady_state
ProcessResult
performance
validate
MoccaInputError
```

## Design variables, gradients and optimisation

```@docs
add_variable!
StepPlace
FeedPlace
ParameterPlace
performance_gradient
optimise!
reset_step_pattern!
```

## Files and the library

```@docs
read_process
write_process
library_processes
load_adsorbent
save_adsorbent
library_adsorbents
register_isotherm!
```

## Ready-made processes

```@docs
four_stage_vsa_process
two_bed_vsa_process
lpp_vsa_process
dual_adsorbent_vsa_process
haghpanah_2013_gas
haghpanah_2013_column
krishnamurthy_2014_gas
krishnamurthy_2014_column
convert_constants
```

## Deprecated

These functions of Mocca 0.1.0 still work, with a deprecation warning when
Julia runs with `--depwarn=yes`:

| Deprecated | Use instead |
|---|---|
| `HaghpanahConstants`, `adsorptionConstants` with `FixedBed(constants)`, `setup_process_model(system, constants)`, `DualSiteLangmuir(constants)` | [`Gas`](@ref), [`Adsorbent`](@ref), [`Column`](@ref) and [`add_bed!`](@ref); [`convert_constants`](@ref) converts the structs |
| `setup_boundary_conditions`, `setup_forces(model, stage_times, bcs)` and the `AdsorptionBC`, `PressurisationBC`, `BlowdownBC`, `EvacuationBC` boundary conditions | a flowsheet with devices and a [`Cycle`](@ref), run as a [`Process`](@ref) |
| `parse_input`, `setup_mocca_case` and their JSON format | [`read_process`](@ref) and [Process files](@ref) |
| `four_stage_vsa_flowsheet`, `two_bed_vsa_flowsheet`, `lpp_vsa_flowsheet`, `dual_adsorbent_vsa_flowsheet` | [`four_stage_vsa_process`](@ref), [`two_bed_vsa_process`](@ref), [`lpp_vsa_process`](@ref), [`dual_adsorbent_vsa_process`](@ref) |
| `wet_flue_gas_isotherm`, `wet_flue_gas_constants`, `wet_flue_gas_system`, `wet_flue_gas_initial_state` | [`load_adsorbent`](@ref), [`krishnamurthy_2014_gas`](@ref), [`krishnamurthy_2014_column`](@ref) |
| `Column()` for the unit entity | `Unit()` |

## The layer underneath

The objects above build these, which can also be used directly.

```@docs
Flowsheet
Stage
setup_schedule
parallel_stages
PortCondition
ExponentialRamp
port_conductance
setup_flowsheet_model
DualSiteLangmuir
AbstractIsotherm
AbstractMassTransfer
LinearDrivingForce
simulate_to_cyclic_steady_state
newton_cyclic_steady_state
simulate_until_mass_balance
cyclic_steady_state_gradient
```
