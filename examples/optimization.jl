# # Optimization of a cyclic process

# We optimise the four-stage vacuum swing adsorption (VSA) cycle of
# Haghpanah et al. (2013) at cyclic steady state: maximise the CO₂ recovery
# while keeping the purity of the CO₂ product above 90%, by tuning
# * the feed velocity,
# * the low pressure, and
# * the intermediate pressure.
#
# Each evaluation finds the cyclic steady state, starting from the last one,
# and its gradient with respect to the design variables from adjoints (see the
# [gradients](flowsheet_vsa_gradients.md) example). A gradient-based optimiser
# (BFGS within the bounds) then needs far fewer evaluations than a population
# of designs would. For more details about the VSA modelling, see the
# [process from its parts](process_from_parts.md) example.

import Mocca
using Mocca.Units

# # The process and its design variables
# A coarse grid keeps the example quick; the bounds are those of Table 5 of
# the paper.
process = Mocca.four_stage_vsa_process(ncells = 20)
Mocca.add_variable!(process, :v_feed, ("adsorption", :velocity); bounds = (0.1, 2.0), description = "interstitial feed velocity [m/s]")
Mocca.add_variable!(process, :P_L, ("evacuation", :to); bounds = (0.05bar, 0.5bar), description = "low pressure")
Mocca.add_variable!(process, :P_I, ("blowdown", :to); bounds = (0.05bar, 0.5bar), description = "intermediate pressure")

start = Mocca.performance(process)
println("start: purity $(round(100start.purity, digits = 1))%, recovery $(round(100start.recovery, digits = 1))%")

# # Run the optimisation
# The purity target is a penalty, so the optimum may miss it slightly. The
# blowdown must end above the evacuation pressure: P_L − P_I ≤ −0.01 bar.
report(values, r) = println("  ", join(["$k = $(round(v, sigdigits = 4))" for (k, v) in pairs(values)], ", "),
    ": purity $(round(100r.metrics.purity, digits = 2))%, recovery $(round(100r.metrics.recovery, digits = 2))%")
best = Mocca.optimise!(process; maximise = :recovery, at_least = (purity = 0.9,),
    linear_constraints = [(Dict(:P_L => 1.0, :P_I => -1.0), -0.01bar)], max_it = 10, log = report)

# The optimised design and its performance:
println("evaluations: $(best.evaluations)")
println(best.values)
println("purity $(round(100best.metrics.purity, digits = 2))%, recovery $(round(100best.metrics.recovery, digits = 2))%, ",
    "$(round(best.metrics.specific_energy, digits = 1)) kWh/t")

# The process is left at the best design. To share it, export its steady
# state to a spreadsheet; `write_results` also saves the process, as
# `process.json`, to rerun it with `read_process`:
dir = joinpath(mktempdir(), "optimised")
Mocca.write_results(dir, Mocca.cyclic_steady_state(process))
foreach(println, filter(startswith("Performance"), readlines(joinpath(dir, "summary.csv"))))
