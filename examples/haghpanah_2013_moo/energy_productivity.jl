# Energy-productivity Pareto front of the four-step cycle at a fixed
# evacuation pressure (Haghpanah et al. 2013, section 7.2, Figure 12), by the
# ε-constraint method: minimise energy with purity and recovery held at 90% or
# more and productivity above a target, stepping the target up along the
# front. With `pvsa`, the feed may be compressed (P_H from 1 to 10 bar).
#
#   julia --project examples/haghpanah_2013_moo/energy_productivity.jl P_L vsa|pvsa [targets...]
#
# with P_L in bar and productivity targets in mol CO2 per m³ adsorbent per s.
# Writes every evaluation and the front to results/haghpanah_2013_moo/.
include(joinpath(@__DIR__, "problem.jl"))

P_L = parse(Float64, get(ARGS, 1, "0.02")) * BAR
pvsa = get(ARGS, 2, "vsa") == "pvsa"
targets = length(ARGS) > 2 ? parse.(Float64, ARGS[3:end]) : [0.0, 0.3, 0.4, 0.5, 0.6, 0.7]
bounds = design_bounds(pvsa = pvsa, P_L = P_L)
tag = @sprintf("%s_PL_%.3f", pvsa ? "pvsa" : "vsa", P_L / BAR)

# Start: near the minimum-energy point of the paper's Table 6 for P_L = 0.02 bar
d_start = [94.89, 122.84, 189.46, max(0.18BAR, P_L + 0.02BAR), P_L, 0.58, 1.0BAR]
pvsa && (d_start[7] = 1.5BAR; d_start[4] = max(d_start[4], 0.2BAR))
u = scaled_design(clamp.(d_start, bounds.lo, bounds.hi), bounds)
warm = WarmStart()

evals = open(joinpath(RESULTS, "$(tag)_evaluations.csv"), "w")
write_row(evals, ["target", "objective", DESIGN_NAMES..., "purity", "recovery", "productivity", "energy", "seconds", "iterations", "jacobians"])
front = open(joinpath(RESULTS, "$(tag)_front.csv"), "w")
write_row(front, ["target", DESIGN_NAMES..., "purity", "recovery", "productivity", "energy", "evaluations", "seconds"])

# Energy is scaled to O(1) for the optimiser
scale = Dict(:energy => 100.0, :productivity => 0.5, :purity => 1.0, :recovery => 1.0)
for target in targets
    n = Ref(0)
    log(d, e, v, s) = (n[] += 1; m = e.metrics; write_row(evals, [target, v, d..., m.purity, m.recovery, m.productivity, m.energy, s, e.css.cycles, e.css.jacobian_builds]))
    t = @elapsed best = optimise(u, bounds, warm; maximise = :energy, sign = -1, energy = true, scale = scale,
        at_least = (purity = 0.9, recovery = 0.9, productivity = target), max_it = 20, log = log)
    global u = best.u
    m = best.metrics
    write_row(front, [target, best.d..., m.purity, m.recovery, m.productivity, m.energy, n[], t])
    @printf("productivity ≥ %.2f: energy %.1f kWh/t, productivity %.3f mol/m³/s, purity %.2f%%, recovery %.2f%% (%d evaluations, %.0f min)\n  %s\n",
        target, m.energy, m.productivity, 100m.purity, 100m.recovery, n[], t / 60, design_string(best.d))
    flush(stdout)
end
close(evals)
close(front)
