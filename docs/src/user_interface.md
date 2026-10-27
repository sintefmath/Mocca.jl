```@meta
CurrentModule = Mocca
```

# Reference

The functions and types for setting up, running and measuring a process.
A process is described in layers, from the physics to the operation. Each
layer is a plain object that can be edited, saved and reused on its own:

| Layer | Objects |
|---|---|
| Gas | [`Gas`](@ref), [`Feed`](@ref), [`Component`](@ref) |
| Materials | [`Adsorbent`](@ref) with an isotherm by component ([`dual_site_langmuir`](@ref), [`langmuir`](@ref), [`virial`](@ref), [`toth_gab`](@ref)), [`Column`](@ref), [`Wall`](@ref) |
| Flowsheet | [`Flowsheet`](@ref), beds ([`add_bed!`](@ref)), moving beds ([`add_moving_bed!`](@ref), [`connect_solids!`](@ref)) and devices with roles ([`add_device!`](@ref)) |
| Operation | a [`Cycle`](@ref) of steps ([`Pressurise`](@ref), [`Adsorb`](@ref), [`Blowdown`](@ref), [`Evacuate`](@ref), [`Purge`](@ref), [`Equalise`](@ref), [`Idle`](@ref), [`Concurrent`](@ref), or your own), [`TimeStepping`](@ref); or device settings held steady |
| Process | [`Process`](@ref), a [`CyclicProcess`](@ref) or a [`SteadyProcess`](@ref): [`simulate`](@ref), [`steady_state`](@ref) ([`cyclic_steady_state`](@ref)), [`performance`](@ref) |
| Design | [`add_variable!`](@ref), [`performance_gradient`](@ref), [`optimise!`](@ref), [`pareto_front!`](@ref), [`explore`](@ref); measures: [`add_measure!`](@ref), [`add_integral!`](@ref) |
| Files | [`read_process`](@ref), [`write_process`](@ref), [`load_adsorbent`](@ref); see [Process files](@ref) |
| Results | [`ProcessResult`](@ref) ([`CyclicResult`](@ref), [`SteadyResult`](@ref)), [`write_results`](@ref), plots; see [Results and plots](@ref) |

The [Building a process from its parts](@ref) example uses each layer in turn.
To add your own kinds of step or unit, see [Extending Mocca](@ref).

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

The heat capacities of a gas come from its thermodynamics model: constant by
default, as given by `heat_capacities` per mass, or changing with temperature
with `thermodynamics = :shomate`:

```julia
gas = Gas(["CO2", "N2", "H2O"]; thermodynamics = :shomate)
molar_heat_capacity(gas.thermodynamics, [0.15, 0.82, 0.03], 350.0)  # J/(mol·K)
```

```@docs
AbstractGasThermodynamics
IdealGasConstantCp
IdealGasShomate
molar_heat_capacity
molar_enthalpy
molar_concentration
heat_capacity_ratio
component_heat_capacities
component_enthalpies
THERMODYNAMIC_REFERENCE_TEMPERATURE
```

Its viscosity and diffusivities are those given, the values at the feed,
unless it has a transport model, by which they change from there with the
state of each cell; the axial dispersion of a bed then follows them and the
gas velocity:

```julia
gas = Gas(["CO2", "N2", "H2O"]; thermodynamics = :polynomial, transport = :chapman_enskog)
```

```@docs
IdealGasPolynomial
AbstractGasTransport
ConstantTransport
ChapmanEnskogTransport
chapman_enskog_transport
TransportReference
mixture_viscosity
molecular_diffusivities
viscosity_factor
diffusivity_factors
gas_diffusivities
AbstractAxialDispersion
ConstantDispersion
WakaoFunazkri
DispersionScale
```

## Isotherms, adsorbents and columns

Isotherm parameters are given by component name, so an adsorbent does not
depend on the order of the components of a gas, and components can be added
or removed without editing parallel arrays.

```@docs
AbstractIsothermData
dual_site_langmuir
langmuir
virial
VirialData
toth_gab
lumped_ldf
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
Flowsheet
add_bed!
add_moving_bed!
add_device!
DEVICE_ROLES
print_flowsheet_diagram
```

## Cycles and steps

A step acts on a bed and opens the devices it needs, found by their roles;
every other device stays closed. Each step that changes the pressure starts
where the step before ended, wrapping around the cycle; give `from` to start
elsewhere.

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
```

## Processes

Every process is a [`Process`](@ref): a [`CyclicProcess`](@ref), operated by
a cycle of steps, or a [`SteadyProcess`](@ref), whose devices keep their
settings. Both are run, settled, measured, varied and optimised through the
same functions. `Process(fs, cycle)` builds a cyclic process.

```@docs
Process
CyclicProcess
ProcessResult
simulate(::CyclicProcess)
initial_state
stages
steady_state(::CyclicProcess)
cyclic_steady_state
performance
pump_heat_capacity_ratio
jacket_heat
validate
MoccaInputError
```

## Steady processes

A moving bed process has no cycle of steps: its sorbent circulates through
sections that each do one job, and its devices keep their settings. It runs
as a [`SteadyProcess`](@ref) until nothing changes any more. The
[Moving bed TSA](examples/mbtsa_mondino_2022.md) example builds one.

```@docs
SteadyProcess
SteadyResult
steady_state(::SteadyProcess)
simulate(::SteadyProcess)
operating_stage
```

## Performance measures

Every measure of a process's performance is a function of integrals over its
cycle, or rates at steady state, and of quantities such as the cycle time.
[`add_measure!`](@ref) adds a measure of your own, from those or from an
integral added with [`add_integral!`](@ref), and [`performance`](@ref)
reports it. [`performance_gradient`](@ref) and [`optimise!`](@ref)
differentiate it without derivative code of its own:

```julia
add_integral!(p, :n2_extract, StreamObjective(:V_vacuum; component = 2))
add_measure!(p, :nitrogen_fraction, s -> s.n2_extract / s.total_extract)
performance_gradient(p; measures = [:nitrogen_fraction, :recovery])
```

```@docs
performance_integrals
performance_quantities
performance_measures
performance_sums
add_integral!
add_measure!
```

## Design variables, gradients and optimisation

The [Optimisation](optimisation.md) page of the manual shows how these fit
together.

A design variable sets one or more places of a process: a field of a step
or of the feed, a unit parameter, the geometry of a bed, or for a steady
process a device setting; several places tie them to one value, and a
[`ScaledPlace`](@ref) to a multiple of it. Its gradient needs no code of the
place's own: setting the variable a little up and down recompiles the cycle
or rebuilds a unit, without simulating, and the changes in the device
settings, stage durations, unit parameters and quantities such as the
adsorbent volume go through the adjoints.
[`performance_gradient`](@ref) differentiates the performance at the state
the process settles to. The cyclic steady state is the fixed point ``x =
\Phi(x)`` of the cycle map, and the gradients follow from adjoints over one
cycle and the implicit function theorem (see
[`cyclic_steady_state_gradient`](@ref)). A steady state is the fixed point of
a cycle of one long backward-Euler step: a state that the step returns
satisfies ``F(x) = 0``, whatever the step's length. The same machinery then
gives the gradients of a steady process, from the adjoint of the steady
equations at the steady state, for unit parameters such as the sorbent flow
or the water temperatures and for device settings such as the purge rate.
They are as accurate as the steady state they are taken at.

For a steady process, `performance(result)` and the variables read the
process's current settings and parameters: measure a result before changing
the process to compute finite differences.

```@docs
add_variable!
Place
StepPlace
FeedPlace
ParameterPlace
SettingPlace
GeometryPlace
ScaledPlace
performance_gradient
optimise!
pareto_front!
explore
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
dac_tvsa_process
haghpanah_2013_gas
haghpanah_2013_column
krishnamurthy_2014_gas
krishnamurthy_2014_column
stampi_bombelli_2020_gas
stampi_bombelli_2020_column
```

## Index

```@index
Pages = ["user_interface.md", "results.md"]
```
