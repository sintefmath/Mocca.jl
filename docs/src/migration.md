```@meta
CurrentModule = Mocca
```

# Migrating from Mocca 0.1.0

Mocca 0.1.0 described a process by a constants struct, boundary conditions
and stage times. These were removed in favour of the objects in the [Reference](@ref):

| Removed | Use instead |
|---|---|
| `adsorptionConstants` with `FixedBed(constants)`, `setup_process_model(system, constants)`, `DualSiteLangmuir(constants)` | [`Gas`](@ref), [`Adsorbent`](@ref), [`Column`](@ref) and [`add_bed!`](@ref) |
| `setup_boundary_conditions`, `setup_forces(model, stage_times, bcs)` and the `AdsorptionBC`, `PressurisationBC`, `BlowdownBC`, `EvacuationBC` boundary conditions | a flowsheet with devices and a [`Cycle`](@ref), run as a [`Process`](@ref) |
| `setup_mocca_case` | [`read_process`](@ref) and [Process files](@ref) |
| `four_stage_vsa_flowsheet`, `two_bed_vsa_flowsheet`, `lpp_vsa_flowsheet`, `dual_adsorbent_vsa_flowsheet` | [`four_stage_vsa_process`](@ref), [`two_bed_vsa_process`](@ref), [`lpp_vsa_process`](@ref), [`dual_adsorbent_vsa_process`](@ref) |
| `wet_flue_gas_isotherm`, `wet_flue_gas_constants`, `wet_flue_gas_system`, `wet_flue_gas_initial_state` | [`load_adsorbent`](@ref), [`krishnamurthy_2014_gas`](@ref), [`krishnamurthy_2014_column`](@ref) |
| `Column()` for the unit entity | `Unit()` |

The constants structs remain for moving inputs over: `HaghpanahConstants`,
or `Mocca.parse_input` for a JSON file of Mocca 0.1.0, gives a struct that
[`convert_constants`](@ref) turns into the gas, feed, adsorbent and column of
a [`Flowsheet`](@ref).

```@docs
convert_constants
```
