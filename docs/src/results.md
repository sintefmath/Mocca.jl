```@meta
CurrentModule = Mocca
```

# Results and plots

[`simulate`](@ref) and [`cyclic_steady_state`](@ref) return a
[`ProcessResult`](@ref). It holds the state of every bed and device after each
time step the solver took:

| Field | Contents |
|---|---|
| `result.states` | One state per time step; `result.states[k][:Bed]` is the state of unit `:Bed` |
| `result.timesteps` | Length of each time step [s]; `cumsum(result.timesteps)` gives the times |
| `result.cycles` | Number of cycles simulated |
| `result.converged`, `result.iterations`, `result.history` | For a cyclic steady state: whether it converged, in how many iterations, and the change per iteration |

A result of [`simulate`](@ref) covers every cycle simulated, from time 0. A
result of [`cyclic_steady_state`](@ref) covers only the steady cycle, with its
times measured from the start of that cycle.

## Performance

[`performance`](@ref) returns a NamedTuple; the formulas are on the
[Model equations](@ref) page.

| Field | Meaning | Unit |
|---|---|---|
| `purity` | Mole fraction of the key component in the extract | – |
| `recovery` | Fraction of the key component fed that leaves in the extract | – |
| `specific_energy` | Pump and compressor work per tonne of key component in the extract | kWh/t |
| `productivity` | Key component in the extract per m³ of adsorbent per second | mol/(m³·s) |
| `energy` | Pump and compressor work over the cycle | J |
| `cycle_time` | Duration of the cycle | s |
| `feed`, `extract`, `product` | Moles of each component through the feed, extract and light product devices over the cycle | mol |
| `adsorbent_volume`, `adsorbent_mass` | Adsorbent in all beds | m³, kg |

```julia
perf = performance(result)
perf.purity, perf.recovery
perf.extract[1] * 44.01e-3       # kg of CO₂ captured per cycle
```

## Variables of beds and devices

The state of a fixed bed holds these variables. Values per cell are ordered
from the bottom of the bed (cell 1, ``z = 0``) to the top (cell `ncells`), and
values per component are in the order of the [`Gas`](@ref).

| Variable | Size | Meaning | Unit |
|---|---|---|---|
| `:Pressure` | ncells | Gas pressure | Pa |
| `:Temperature` | ncells | Gas and sorbent temperature | K |
| `:WallTemperature` | ncells | Wall temperature (`WithWall` beds only) | K |
| `:y` | components × ncells | Gas mole fractions | – |
| `:AdsorbedConcentration` | components × ncells | Loading ``q_i`` | mol/m³ adsorbent |

To convert a loading to mol/kg, divide it by the adsorbent `density`. The cell
centres are at `fs[:Bed].data_domain[:cell_centroids][1, :]` [m], for a
flowsheet `fs` (`process.flowsheet`).

The state of a flow device holds the flow through it and the gas at each of
its two ports:

| Variable | Size | Meaning | Unit |
|---|---|---|---|
| `:MolarFlow` | 1 | Molar flow, positive from inlet to outlet | mol/s |
| `:InletPressure`, `:OutletPressure` | 1 | Pressure at each port | Pa |
| `:InletTemperature`, `:OutletTemperature` | 1 | Temperature at each port | K |
| `:InletComposition`, `:OutletComposition` | components × 1 | Mole fractions at each port | – |

For example, the CO₂ mole fraction leaving the top of the bed, and the flow
through the vacuum line, through time:

```julia
t = cumsum(result.timesteps)
y_CO2_top = [s[:Bed][:y][1, end] for s in result.states]
F_vacuum = [s[:V_vacuum][:MolarFlow][1] for s in result.states]
```

[`stream_totals`](@ref) adds up the moles of each component that passed
through a device:

```julia
n = stream_totals(result.states, result.timesteps, :V_vacuum)   # mol of CO₂, mol of N₂
```

## Export to a spreadsheet

[`write_results`](@ref) writes a result to a folder of CSV files, in bara, °C,
mol%, kg/h, Nm³/h and mol/kg, with the units in the column headers:

```julia
result = cyclic_steady_state(process)
write_results("results/base_case", result)
```

| File | Contents |
|---|---|
| `summary.csv` | Feed, stage durations, purity, recovery, productivity, specific energy, the mass balance of each component over the cycle, and the beds |
| `streams.csv` | Stream table: one row per device, with its composition, mean flow and amount per cycle |
| `stage_streams.csv` | The stream table for each stage |
| `stream_history.csv` | Flow, composition, temperature and pressures of every device at every time step |
| `bed_history.csv` | Pressure, temperatures, gas composition and loading at both ends of each bed at every time step |
| `bed_profiles.csv` | The same along each bed, at the start of the cycle and the end of each stage |
| `process.json` | The process, to read back with [`read_process`](@ref) |

The time series have one row per time step and device (or bed end), with the
cycle and stage it belongs to, so a pivot table or a filter picks out one
device or stage. For a spreadsheet set up with a decimal comma, use
`write_results(dir, result; delim = ';', decimal = ',')`. The
[Getting started](@ref) example writes these files and reads some of them
back.

For other variables, write the CSV file yourself:

```julia
t = cumsum(result.timesteps)
open("bed_top.csv", "w") do io
    println(io, "time_s,pressure_Pa,temperature_K,y_CO2,q_CO2_mol_per_m3")
    for (tk, s) in zip(t, result.states)
        b = s[:Bed]
        println(io, join((tk, b[:Pressure][end], b[:Temperature][end], b[:y][1, end],
            b[:AdsorbedConcentration][1, end]), ","))
    end
end
```

## Plots

Mocca draws its plots with [Makie](https://docs.makie.org). Each plot
function returns a figure, which shows in VS Code's plot pane or a notebook.
From a terminal, save it to a file instead: add CairoMakie to your
environment (see [Installation](@ref)) and use `save`:

```julia
using CairoMakie
bed = process.flowsheet[:Bed]
bed_states = [s[:Bed] for s in result.states]
f = plot_cell(bed_states, bed, result.timesteps, 30)   # history of cell 30, the top of a 30-cell bed
save("bed_top.png", f)
```

| Function | Shows |
|---|---|
| [`plot_cell`](@ref) | Every variable of one cell through time, such as the outlet |
| [`plot_state`](@ref) | Every variable along the bed at one time |
| [`plot_profiles`](@ref) | Mole fractions and loadings along the bed at several times, such as the end of each step |
| [`plot_flowsheet_streams`](@ref) | Flow and composition through each device, with the steps shaded |
| [`plot_cell_comparison`](@ref), [`plot_state_comparison`](@ref) | Several simulations of the same bed overlaid |

```@docs
ProcessResult
write_results
plot_cell
plot_state
plot_profiles
plot_flowsheet_streams
plot_cell_comparison
plot_state_comparison
stream_totals
vacuum_pump_energy
specific_energy_kwh_per_tonne
```
