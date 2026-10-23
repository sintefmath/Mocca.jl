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
| Materials | [`Adsorbent`](@ref) with an isotherm by component ([`dual_site_langmuir`](@ref), [`langmuir`](@ref)), [`Column`](@ref), [`Wall`](@ref) |
| Flowsheet | [`Flowsheet`](@ref), beds ([`add_bed!`](@ref)) and devices with roles ([`add_device!`](@ref)) |
| Operation | a [`Cycle`](@ref) of steps ([`Pressurise`](@ref), [`Adsorb`](@ref), [`Blowdown`](@ref), [`Evacuate`](@ref), [`Purge`](@ref), [`Equalise`](@ref), [`Idle`](@ref), [`Concurrent`](@ref), or your own), [`TimeStepping`](@ref) |
| Process | [`Process`](@ref): [`simulate`](@ref), [`cyclic_steady_state`](@ref), [`performance`](@ref) |
| Design | [`add_variable!`](@ref), [`performance_gradient`](@ref), [`optimise!`](@ref) |
| Files | [`read_process`](@ref), [`write_process`](@ref), [`load_adsorbent`](@ref); see [Process files](@ref) |
| Results | [`ProcessResult`](@ref), [`write_results`](@ref), plots; see [Results and plots](@ref) |

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

A bed can be meshed across its radius as well as along it, in rings from the
axis to the wall, to resolve the temperature and composition across a wide
column (see [Beds meshed across their radius](@ref)). Give `radial_cells` to
[`add_bed!`](@ref), and the column the coefficients the 1D bed lumps
together: the effective radial conductivity of the packing, and the
coefficient of the film at the wall, which takes the place of
`inner_heat_transfer`:

```julia
column = Column(; length = 1.0, radius = 0.1445, outer_radius = 0.162, porosity = 0.37,
    radial_conductivity = 0.5, wall_heat_transfer = 30.0)
add_bed!(fs, :Bed, adsorbent, column; ncells = 30, radial_cells = 6)
```

Its ports are the rings at each end, so the devices, steps, performance and
results are those of any bed; the bed ends and the profiles of the results
are means over the cross-section, and `bed_fields.csv`,
[`plot_bed_section`](@ref) and [`plot_radial_profiles`](@ref) show the field
across it. The example [A column with a cooled wall, across its radius](@ref)
compares such a bed with a 1D one.

```@docs
Flowsheet
add_bed!
bed_dispersion
bed_radial_dispersion
compute_dispersion
compute_radial_dispersion
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

```@docs
Process
simulate(::Process)
initial_state
stages
cyclic_steady_state
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
```

## Index

```@index
Pages = ["user_interface.md", "results.md"]
```
