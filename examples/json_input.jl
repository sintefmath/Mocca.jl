# # Quick start: processes from files
# A whole process — gas, feed, adsorbents, columns, beds, devices, cycle,
# time stepping and design variables — can be read from a JSON file. Mocca has
# a library of processes and adsorbents from published studies.

import Mocca
using Mocca.Units

println(Mocca.library_processes())
println(Mocca.library_adsorbents())

# # Read a process from the library
# The four-step VSA cycle of Haghpanah et al. (2013), with the design
# variables of their optimisation study:
process = Mocca.read_process("haghpanah_2013_four_stage_vsa")

# # Change it
# Design variables set every place they act: the evacuation pressure is also
# where the next pressurisation starts.
process[:P_L] = 0.05bar
process[:t_ads] = 30.0
Mocca.stages(process)

# # Run it to cyclic steady state and measure it
result = Mocca.cyclic_steady_state(process)
perf = Mocca.performance(result)
println("purity $(round(100perf.purity, digits = 1))%, recovery $(round(100perf.recovery, digits = 1))%, ",
    "$(round(perf.specific_energy, digits = 1)) kWh/t, productivity $(round(perf.productivity, sigdigits = 3)) mol/(m³·s)")

# # Write it, edit it by hand, read it back
# Numbers are SI, or strings with a unit, such as `"0.1 bar"` or `"25 °C"`.
# Compositions and per-component values are keyed by component name, steps
# by their `type`, and an adsorbent may be the name of a library adsorbent.
file = Mocca.write_process(joinpath(mktempdir(), "vsa.json"), process)

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

# Mistakes are reported with where they are:
try
    write(file, replace(text, "\"to\": \"0.2 bar\"" => "\"to\": \"2 bar\""))
    Mocca.read_process(file)
catch e
    println(sprint(showerror, e))
end

# # Simulate a few cycles and plot the outlet of the bed
result = Mocca.simulate(short; cycles = 3)
bed_states = [s[:Bed] for s in result.states]
f = Mocca.plot_cell(bed_states, short.flowsheet[:Bed], result.timesteps, 30)
