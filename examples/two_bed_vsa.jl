#=
# Two-bed VSA with pressure equalisation

This example builds a process of two beds that work in turn: while one bed is
pressurised and fed, the other is evacuated. Between the two halves of the
cycle, a pressure equalisation joins the tops of the beds, so that the bed
that has just been fed passes some of its gas to the bed that has just been
evacuated, instead of sending it to the vacuum pump.

    feed ──V_feed_A──┐               ┌──V_product_A── raffinate
                     bottom  A  top ─┤
    CO2 ◄─V_vacuum_A─┘               │
                                    V_eq
    feed ──V_feed_B──┐               │
                     bottom  B  top ─┤
    CO2 ◄─V_vacuum_B─┘               └──V_product_B── raffinate

The cycle has six stages:

| Stage | Bed A | Bed B |
|---|---|---|
| 1 | pressurisation | evacuation |
| 2 | adsorption | evacuation |
| 3 | equalisation, A to B | equalisation, A to B |
| 4 | evacuation | pressurisation |
| 5 | evacuation | adsorption |
| 6 | equalisation, B to A | equalisation, B to A |

The beds and gas are those of [Haghpanah et al.
(2013)](https://doi.org/10.1021/ie302658y), as in [Getting
started](getting_started.md). The same process is ready-made as
`Mocca.two_bed_vsa_process()`.
=#

import Mocca
using Mocca.Units

# # 1. The flowsheet
# Two beds, each with a feed, a product and an extract device, and a link
# between their tops. A link is a valve with a bed at each end; the
# equalisation steps find it by the beds it joins.
gas = Mocca.haghpanah_2013_gas()
feed = Mocca.Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), temperature = 298.15, pressure = 1bar)
adsorbent = Mocca.load_adsorbent("zeolite_13x_haghpanah_2013")
column = Mocca.haghpanah_2013_column()

fs = Mocca.Flowsheet(gas, feed)
for b in (:A, :B)
    Mocca.add_bed!(fs, b, adsorbent, column; ncells = 20)
    Mocca.add_device!(fs, Symbol(:V_feed_, b), :feed; at = (b, :bottom))
    Mocca.add_device!(fs, Symbol(:V_product_, b), :product; at = (b, :top))
    Mocca.add_device!(fs, Symbol(:V_vacuum_, b), :extract; at = (b, :bottom))
end
Mocca.add_device!(fs, :V_eq, :link; from = (:A, :top), to = (:B, :top))
fs

# # 2. The cycle
# With several beds, each step names its bed, and `Concurrent` runs steps on
# different beds at the same time, as one stage. `Equalise` opens the link
# between two beds, and gas flows from the one at the higher pressure. The
# size of the link valve, its `conductance` [mol/(s·Pa)], sets how far the
# pressures come together in the time given.
#
# Pressures carry over from step to step within each bed. After an
# equalisation both beds are taken to be at the mean of their pressures,
# which is where the next pressurisation and evacuation ramps start. The
# evacuation of each bed runs through two stages, the second carrying on along
# the ramp of the first.
t_press, t_ads, t_eq = 15.0, 15.0, 10.0
P_eq = 0.55bar
evacuate(b, t, from) = Mocca.Evacuate(t; bed = b, from = from, to = 0.1bar)
P_mid = 0.1bar + (P_eq - 0.1bar) * exp(-0.5 * t_press)   # the ramp at 0.5 s⁻¹ after the first stage
cycle = Mocca.Cycle([
    Mocca.Concurrent(Mocca.Pressurise(t_press; bed = :A, from = P_eq), evacuate(:B, t_press, P_eq)),
    Mocca.Concurrent(Mocca.Adsorb(t_ads; bed = :A, velocity = 1.0), evacuate(:B, t_ads, P_mid)),
    Mocca.Equalise(t_eq; beds = (:A, :B), conductance = 1e-6),
    Mocca.Concurrent(Mocca.Pressurise(t_press; bed = :B, from = P_eq), evacuate(:A, t_press, P_eq)),
    Mocca.Concurrent(Mocca.Adsorb(t_ads; bed = :B, velocity = 1.0), evacuate(:A, t_ads, P_mid)),
    Mocca.Equalise(t_eq; beds = (:B, :A), conductance = 1e-6),
])

# The process lists each stage and the devices it opens:
process = Mocca.Process(fs, cycle; name = "two-bed VSA")

# # 3. Cyclic steady state and performance
# Two beds approach their steady state more slowly than one, so the first
# steady state takes a few minutes.
result = Mocca.cyclic_steady_state(process)
perf = Mocca.performance(result)
println("converged: $(result.converged) after $(result.iterations) cycles")
println("purity $(round(100perf.purity, digits = 1))%, recovery $(round(100perf.recovery, digits = 1))%, ",
    "$(round(perf.specific_energy, digits = 1)) kWh/t")

# The streams through the devices over the cycle. The link carries gas from A
# to B in the first equalisation (positive flow, from its inlet at A) and back
# in the second:
f_streams = Mocca.plot_flowsheet_streams(result.states, process.model, result.timesteps;
    stages = Mocca.stages(process), devices = [:V_feed_A, :V_vacuum_A, :V_product_A, :V_eq])

# The top of each bed through the cycle. The beds run half a cycle apart, and
# their pressures meet at about half a bar in each equalisation:
f_beds = Mocca.plot_cell_comparison(["A" => ([s[:A] for s in result.states], result.timesteps),
    "B" => ([s[:B] for s in result.states], result.timesteps)], fs[:A], 20)

# # 4. The stream table of each stage
# `write_results` writes the stream table of each stage to `stage_streams.csv`.
# It shows how much gas moves between the beds during equalisation, and what
# it is:
dir = joinpath(mktempdir(), "two_bed")
Mocca.write_results(dir, result)
lines = readlines(joinpath(dir, "stage_streams.csv"))
header = split(lines[1], ',')
columns = indexin(["Stage", "Device", "From", "To", "CO2 [mol%]", "Total [mol/stage]"], header)
for line in filter(contains("equalisation"), lines)
    println(join(split(line, ',')[columns], ", "))
end

# About 1.3 mol of gas, with 10% CO₂, moves from A to B in the first
# equalisation, and the same back in the second (a negative amount: against
# the direction of the link).
#
# Stage names with a comma, such as "A adsorption, B evacuation", are quoted
# in the file, so only the equalisation rows are split this simply here. A
# spreadsheet or a CSV reader such as CSV.jl reads them all.
#
# The mass balance over the cycle in `summary.csv` closes over both beds:
foreach(println, filter(startswith("Mass balance"), readlines(joinpath(dir, "summary.csv"))))

# # 5. What the equalisation is worth
# A design variable that sets the conductance of the link in both
# equalisation steps. Nearly closing it leaves the gas of the fed bed to the
# vacuum pump. The steady state for the new setting starts from the one
# already found.
Mocca.add_variable!(process, :C_eq, ("equalisation A to B", :conductance), ("equalisation B to A", :conductance))
process[:C_eq] = 1e-9
result_no_eq = Mocca.cyclic_steady_state(process)
println("converged: $(result_no_eq.converged) after $(result_no_eq.iterations) cycles")
no_eq = Mocca.performance(result_no_eq)
for (label, p) in ("with equalisation" => perf, "without" => no_eq)
    println(rpad(label, 20), "purity $(round(100p.purity, digits = 1))%, recovery $(round(100p.recovery, digits = 1))%, ",
        "$(round(p.specific_energy, digits = 1)) kWh/t")
end

# Equalisation sends the nitrogen-rich gas left in the fed bed to the other
# bed, rather than into the CO₂ product through the vacuum line, which raises
# the purity by about 16 percentage points and saves pump work.
