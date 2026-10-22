```@meta
CurrentModule = Mocca
```

# Process files

A process can be read from and written to a JSON file with
[`read_process`](@ref) and [`write_process`](@ref), and an adsorbent with
[`load_adsorbent`](@ref) and [`save_adsorbent`](@ref). Mocca has a library of
both from published studies:

```julia
using Mocca
library_processes()     # "haghpanah_2013_four_stage_vsa", ...
library_adsorbents()    # "zeolite_13x_haghpanah_2013", ...
process = read_process("haghpanah_2013_four_stage_vsa")
write_process("my_vsa.json", process)
```

## Conventions

- Numbers are plain SI: Pa, K, m, s, mol, kg, J. A number may also be a
  string with a unit, such as `"0.1 bar"`, `"25 °C"`, `"14.45 cm"`,
  `"-36.6 kJ/mol"` or `"100 cm/s"`.
- Compositions and per-component values are objects keyed by component name,
  such as `{"CO2": 0.15, "N2": 0.85}`. Components left out of a composition
  get a trace (``10^{-10}``).
- Unknown keys are errors, so that a misspelt setting is not silently
  ignored. Errors say where they are, such as
  `MoccaInputError in step 3 (blowdown) of the cycle, to: ...`.

## A process file

```json
{
  "format": "mocca-process",
  "version": 1,
  "name": "four-step VSA",
  "gas": {
    "components": ["CO2", "N2"],
    "molar_masses": {"N2": 0.028},
    "heat_capacities": {"CO2": 697.5687, "N2": 1096.4},
    "viscosity": 1.72e-5,
    "density": 1.22638310956,
    "diffusivity": 1.6e-5
  },
  "feed": {"composition": {"CO2": 0.15, "N2": 0.85}, "temperature": "25 °C", "pressure": "1 bar"},
  "adsorbents": {"13X": "zeolite_13x_haghpanah_2013"},
  "columns": {
    "column": {"length": 1.0, "radius": 0.1445, "outer_radius": 0.162, "porosity": 0.37,
               "axial_conductivity": 0.0903,
               "wall": {"density": 7800, "heat_capacity": 502, "conductivity": 16},
               "inner_heat_transfer": 8.6, "outer_heat_transfer": 2.5, "ambient_temperature": 298.15}
  },
  "units": [
    {"name": "Bed", "type": "fixed_bed", "adsorbent": "13X", "column": "column", "ncells": 30,
     "thermal": "with_wall", "dispersion_velocity": "adsorption"}
  ],
  "devices": [
    {"name": "V_feed", "role": "feed", "at": ["Bed", "bottom"]},
    {"name": "V_product", "role": "product", "at": ["Bed", "top"]},
    {"name": "V_vacuum", "role": "extract", "at": ["Bed", "bottom"]}
  ],
  "cycle": {"steps": [
    {"type": "pressurise", "duration": 15},
    {"type": "adsorb", "duration": 15, "velocity": 1.0},
    {"type": "blowdown", "duration": 30, "to": "0.2 bar"},
    {"type": "evacuate", "duration": 40, "to": "0.1 bar"}
  ]},
  "time_stepping": {"max_dt": 1.0, "first_dt": 0.05, "keep_step_pattern": "auto"},
  "options": {"key_component": "CO2", "ambient_pressure": 1e5, "pump_efficiency": 0.72,
              "heat_capacity_ratio": 1.4, "ramp_rate": 0.5},
  "initial_state": {"Bed": {"pressure": "1 bar"}},
  "parameters": {"Bed": {"AdsorbentDensity": 1130}},
  "variables": [
    {"name": "P_L", "places": [{"step": "evacuation", "field": "to"}],
     "bounds": ["0.01 bar", "0.5 bar"], "description": "evacuation pressure"}
  ]
}
```

The sections:

| Section | Contents |
|---|---|
| `components` | Optional: components not known to Mocca, each with `name`, `molar_mass`, `heat_capacity` and `viscosity`, and optionally `formula`, `cas` and `shomate` (the five coefficients). |
| `gas` | The [`Gas`](@ref): its `components` in order (the first is the one to capture), and optional overrides. Heat capacities are constant, `heat_capacities` per component; or with `"thermodynamics": "shomate"` from the Shomate equation, with optional `shomate` coefficients `[A, B, C, D, E]` per component. |
| `feed` | The [`Feed`](@ref): `composition`, `temperature`, `pressure` and optionally `light_product`. |
| `adsorbents` | By name: an adsorbent object (as in an adsorbent file), the name of a library adsorbent, or the path of an adsorbent file relative to the process file. |
| `columns` | By name: [`Column`](@ref) objects. |
| `units` | Beds, in order, as in [`add_bed!`](@ref). `thermal` is `"with_wall"`, `"adiabatic"` or `{"type": "isothermal", "temperature": ...}`. |
| `devices` | Devices, in order, as in [`add_device!`](@ref): a `role`, and `at` (a `[unit, port]` pair) or, for a link, `from` and `to`. An optional `boundary` sets `pressure`, `temperature` and `composition`. |
| `cycle` | `steps`, a list of steps; or `sequences`, a list of `{"label": ..., "steps": [...]}` run side by side, and an optional `cycle_time`. |
| `time_stepping` | See [`TimeStepping`](@ref). |
| `options` | See [`Process`](@ref). Without `heat_capacity_ratio` the pumps take that of the feed. |
| `initial_state` | Per bed: any of `pressure`, `temperature`, `wall_temperature` and `composition`. |
| `parameters` | Per unit: Jutul parameters to override. |
| `variables` | Design variables (see [`add_variable!`](@ref)), each with `places`: `{"step": ..., "field": ..., "sequence": ...}`, `{"feed": "pressure"}` or `{"unit": ..., "parameter": ...}`. |

## Steps

Each step has a `type` and its fields as keys:

| `type` | Step | Fields |
|---|---|---|
| `pressurise` | [`Pressurise`](@ref) | `duration`, `to`, `from`, `source` (`"feed"` or `"light_product"`), `bed`, `device`, `ramp_rate`, `conductance`, `name` |
| `adsorb` | [`Adsorb`](@ref) | `duration`, one of `velocity`, `superficial_velocity` and `flow_rate`, `pressure`, `bed`, `feed_device`, `outlet`, `conductance`, `name` |
| `blowdown` | [`Blowdown`](@ref) | `duration`, `to`, `from`, `through`, `device`, `bed`, `ramp_rate`, `conductance`, `name` |
| `evacuate` | [`Evacuate`](@ref) | as `blowdown` |
| `purge` | [`Purge`](@ref) | `duration`, one of `velocity`, `superficial_velocity` and `flow_rate`, `pressure`, `bed`, `conductance`, `name` |
| `equalise` | [`Equalise`](@ref) | `duration`, `beds`, `device`, `conductance`, `name` |
| `idle` | [`Idle`](@ref) | `duration`, `bed`, `name` |
| `concurrent` | [`Concurrent`](@ref) | `steps`, `name` |

Step types defined with [`register_step!`](@ref) are read and written the
same way.

## An adsorbent file

```json
{
  "format": "mocca-adsorbent",
  "version": 1,
  "name": "zeolite 13X",
  "reference": "Haghpanah et al. (2013)",
  "density": 1130.0,
  "heat_capacity": 1070.0,
  "particle_diameter": 0.002,
  "particle_porosity": 0.35,
  "tortuosity": 3.0,
  "adsorbed_heat_capacity": {"CO2": 697.5687, "N2": 1096.4},
  "isotherm": {
    "type": "dual_site_langmuir",
    "reference_temperature": 298.15,
    "basis": "volume",
    "components": {
      "CO2": {"qsb": 3489.44, "b0": 8.65e-7, "dUb": -36641.21, "qsd": 2872.35, "d0": 2.63e-8, "dUd": -35690.66},
      "N2": {"qsb": 6613.551, "b0": 2.5e-6, "dUb": -15800.0}
    }
  }
}
```

`basis` is `"volume"` for capacities in mol per m³ of adsorbent or `"mass"`
for mol/kg. Isotherm types defined with [`register_isotherm!`](@ref) are read
the same way.

## The earlier format

`Mocca.parse_input` reads the JSON files of Mocca 0.1.0 into constants
structs, and [`convert_constants`](@ref) turns those into the objects above,
to move such inputs over. Mocca no longer simulates them directly.
