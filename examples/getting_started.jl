#=
# Getting started

This example goes through a whole study of a process: read it from Mocca's
library, change its operating settings, run it to cyclic steady state,
measure its performance, look at the results, export them to a spreadsheet,
save the process to a file to edit and share, and check that the cycle is
steady.

The process is the four-step vacuum swing adsorption (VSA) cycle of
[Haghpanah et al. (2013)](https://doi.org/10.1021/ie302658y), which captures
CO₂ from a flue gas of 15% CO₂ in N₂ on zeolite 13X:

* **Pressurisation:** the top of the bed is closed and feed enters at the
  bottom, bringing the bed up to the feed pressure.
* **Adsorption:** feed enters the bottom; CO₂ is adsorbed and nitrogen leaves
  the top.
* **Blowdown:** the bottom is closed and the bed is taken down to an
  intermediate pressure through the top, which removes most of the nitrogen.
* **Evacuation:** the top is closed and the bed is evacuated through the
  bottom, which gives the CO₂ product.
=#

import Mocca
using Mocca.Units

# # 1. Read a process from the library
# The library holds the processes and adsorbents of published studies:
println(Mocca.library_processes())
println(Mocca.library_adsorbents())

#-
# A process holds everything needed to run it: the gas and feed, the adsorbent
# and column, the flowsheet of beds and devices, the cycle of steps, the time
# stepping and the design variables.
process = Mocca.read_process("haghpanah_2013_four_stage_vsa")

# # 2. Change the operating settings
# The design variables of the process are named settings. Setting one changes
# every place it acts: the evacuation pressure `P_L` is also where the next
# pressurisation starts. Values are in SI units; `Mocca.Units` has
# conversions such as `bar`.
Mocca.variable_names(process)

#-
process[:P_L] = 0.05bar
process[:t_ads] = 30.0
process[:P_L], process[:t_ads]

# The cycle is compiled into stages, each a set of device settings:
for st in Mocca.stages(process)
    println(rpad(st.name, 16), st.duration, " s, opens ", join(keys(st.settings), ", "))
end

# # 3. Run it to cyclic steady state
# `cyclic_steady_state` repeats the cycle until the state of the bed at the end
# of a cycle is the state it started from, accelerating the approach with
# Newton's method. The result holds the steady cycle.
result = Mocca.cyclic_steady_state(process)
println("converged: $(result.converged), cycles: $(result.iterations)")

# # 4. Measure its performance
# Purity and recovery of CO₂ in the extract (the vacuum line), the pump and
# compressor energy per tonne of CO₂, and the productivity per m³ of adsorbent:
perf = Mocca.performance(result)
println("purity $(round(100perf.purity, digits = 1))%, recovery $(round(100perf.recovery, digits = 1))%, ",
    "$(round(perf.specific_energy, digits = 1)) kWh/t, productivity $(round(perf.productivity, sigdigits = 3)) mol/(m³·s)")

# # 5. Look at the results
# The flow and composition through each device over the cycle, with the stages
# shaded. Flow is positive from a device's inlet to its outlet, so the vacuum
# line is positive while it draws gas out of the bed.
stages = Mocca.stages(process)
f_streams = Mocca.plot_flowsheet_streams(result.states, process.model, result.timesteps; stages = stages)

# Mole fractions and loadings along the bed at the end of each stage. Most of
# the CO₂ is held near the feed end, at the bottom. Blowdown removes most of
# the nitrogen, and at the end of evacuation the gas in the bed is nearly pure
# CO₂.
bed = process.flowsheet[:Bed]
f_profiles = Mocca.plot_profiles(result.states, result.timesteps, bed, cumsum([st.duration for st in stages]);
    labels = [st.name for st in stages], unit = :Bed)

# Each state in `result.states` holds the state of every bed and device after a
# time step, in SI units, for any other quantity. For example, the range of the
# bed temperature over the cycle:
T_min = minimum(minimum(s[:Bed][:Temperature]) for s in result.states)
T_max = maximum(maximum(s[:Bed][:Temperature]) for s in result.states)
println("bed temperature from $(round(T_min - 273.15, digits = 1)) to $(round(T_max - 273.15, digits = 1)) °C")

# # 6. Export to a spreadsheet
# `write_results` writes the result to a folder of CSV files that open in a
# spreadsheet, in bara, °C, mol%, kg/h, Nm³/h and mol/kg, with the units in
# the column headers. It also writes the process, to rerun it later.
dir = joinpath(mktempdir(), "base_case")
files = Mocca.write_results(dir, result)
foreach(println, basename.(files))

# `summary.csv` holds the feed, the cycle, the performance and the mass balance
# of each component over the cycle:
print(read(joinpath(dir, "summary.csv"), String))

# `streams.csv` is a stream table, one row per device, with its composition,
# flow, amount per cycle, temperature and pressures. Some of its columns:
rows = split.(readlines(joinpath(dir, "streams.csv")), ',')
columns = indexin(["Device", "CO2 [mol%]", "N2 [mol%]", "Mean flow [kg/h]", "Total [mol/cycle]"], rows[1])
for row in rows
    println(join(rpad.(row[columns], 20)))
end

# `stage_streams.csv` breaks it down by stage. `stream_history.csv` and
# `bed_history.csv` hold the devices and both ends of each bed at every time
# step, and `bed_profiles.csv` the state along each bed at the end of each
# stage, with the cycle and stage of each row so that a filter or a pivot
# table picks out one device or stage:
foreach(println, first(eachline(joinpath(dir, "bed_history.csv")), 4))

# For a spreadsheet set up with a decimal comma, pass `delim = ';', decimal = ','`.

# # 7. Save the process, edit it and read it back
# A process file holds every layer of the process and reads back as the same
# process:
file = Mocca.write_process(joinpath(mktempdir(), "vsa.json"), process)
same = Mocca.read_process(file)
same[:P_L] / bar

# Process files can be written by hand. Numbers are SI, or strings with a unit,
# such as `"0.1 bar"` or `"25 °C"`. Compositions and per-component values are
# keyed by component name, steps by their `type`, and an adsorbent may be the
# name of a library adsorbent. The [Process files](@ref) page describes the
# format.
text = """
{
  "name": "a short cycle",
  "gas": {"components": ["CO2", "N2"], "molar_masses": {"N2": 0.028}, "heat_capacities": {"CO2": 697.5687, "N2": 1096.4}},
  "feed": {"composition": {"CO2": 0.15, "N2": 0.85}, "pressure": "1 bar", "temperature": "25 °C"},
  "adsorbents": {"13X": "zeolite_13x_haghpanah_2013"},
  "columns": {"column": {"length": "1 m", "radius": "14.45 cm", "outer_radius": "16.2 cm", "porosity": 0.37}},
  "units": [{"name": "Bed", "adsorbent": "13X", "column": "column", "ncells": 30}],
  "devices": [
    {"name": "V_feed", "role": "feed", "at": ["Bed", "bottom"]},
    {"name": "V_product", "role": "product", "at": ["Bed", "top"]},
    {"name": "V_vacuum", "role": "extract", "at": ["Bed", "bottom"]}
  ],
  "cycle": {"steps": [
    {"type": "pressurise", "duration": "10 s"},
    {"type": "adsorb", "duration": "10 s", "velocity": "1 m/s"},
    {"type": "blowdown", "duration": "20 s", "to": "0.2 bar"},
    {"type": "evacuate", "duration": "30 s", "to": "0.1 bar"}
  ]},
  "time_stepping": {"max_dt": 1.0, "first_dt": 0.05}
}
"""
file = joinpath(mktempdir(), "short.json")
write(file, text)
short = Mocca.read_process(file)

# Mistakes are reported with where they are. A blowdown that goes up in
# pressure:
try
    write(file, replace(text, "\"to\": \"0.2 bar\"" => "\"to\": \"2 bar\""))
    Mocca.read_process(file)
catch e
    println(sprint(showerror, e))
end

# # 8. Check the steady state
# `simulate` runs a given number of cycles from a given state, rather than to
# steady state. Started from the start of the steady cycle, `result.state0`,
# the bed comes back to the same state at the end of every cycle. The largest
# relative change of any variable over three cycles:
result_cycles = Mocca.simulate(process; cycles = 3, state0 = result.state0)
Mocca.cycle_change(process.model, result.state0, result_cycles.states[end])

# The top of the bed, where the nitrogen product leaves, over the three cycles:
bed_states = [s[:Bed] for s in result_cycles.states]
f_outlet = Mocca.plot_cell(bed_states, bed, result_cycles.timesteps, 30)
