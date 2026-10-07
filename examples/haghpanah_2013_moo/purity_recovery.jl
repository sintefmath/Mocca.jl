# Purity-recovery Pareto front of the four-step VSA cycle (Haghpanah et al.
# 2013, section 7.1.1, Figure 10a), by the ε-constraint method: maximise
# recovery with purity held above a target, stepping the target up along the
# front. Each point starts from the previous optimum and its steady state.
#
#   julia --project examples/haghpanah_2013_moo/purity_recovery.jl P_L_min [vsa|pvsa] [targets...]
#
# with P_L_min, the lower bound on the evacuation pressure, in bar (the paper
# uses 0.05, 0.03, 0.02 and 0.01 for VSA, and 0.13, 0.07 and 0.05 for PVSA,
# section 7.1.2 and Figure 10b, where the feed pressure P_H is also free), and
# purity targets as fractions. Writes every evaluation and the front to
# results/haghpanah_2013_moo/.
include(joinpath(@__DIR__, "problem.jl"))

P_L_min = parse(Float64, get(ARGS, 1, "0.02")) * BAR
pvsa = get(ARGS, 2, "vsa") == "pvsa"
first_target = length(ARGS) > 1 && ARGS[2] in ("vsa", "pvsa") ? 3 : 2
targets = length(ARGS) >= first_target ? parse.(Float64, ARGS[first_target:end]) : [0.80, 0.85, 0.88, 0.90, 0.92, 0.94, 0.96]
bounds = design_bounds(pvsa = pvsa, P_L_min = P_L_min)
tag = @sprintf("%s_PLmin_%.3f", pvsa ? "pvsa" : "vsa", P_L_min / BAR)

# Start: a long cycle at a low evacuation pressure and moderate feed velocity;
# for PVSA, a moderately compressed feed
d_start = [80.0, 120.0, 180.0, max(0.18BAR, P_L_min + 0.05BAR), P_L_min + 0.002BAR, 0.6, 1.0BAR]
pvsa && (d_start[4] = max(0.5BAR, P_L_min + 0.1BAR); d_start[7] = 3.0BAR)
u = scaled_design(clamp.(d_start, bounds.lo, bounds.hi), bounds)
warm = WarmStart()

evals = open(joinpath(RESULTS, "$(tag)_evaluations.csv"), "w")
write_row(evals, ["target", "objective", DESIGN_NAMES..., "purity", "recovery", "productivity", "seconds", "iterations", "jacobians"])
front = open(joinpath(RESULTS, "$(tag)_front.csv"), "w")
write_row(front, ["target", DESIGN_NAMES..., "purity", "recovery", "productivity", "evaluations", "seconds"])

for target in targets
    n = Ref(0)
    log(d, e, v, s) = (n[] += 1; m = e.metrics; write_row(evals, [target, v, d..., m.purity, m.recovery, m.productivity, s, e.css.cycles, e.css.jacobian_builds]))
    t = @elapsed best = optimise(u, bounds, warm; maximise = :recovery, at_least = (purity = target,), max_it = 20, log = log)
    global u = best.u
    m = best.metrics
    write_row(front, [target, best.d..., m.purity, m.recovery, m.productivity, n[], t])
    @printf("purity ≥ %.2f: purity %.2f%%, recovery %.2f%% (%d evaluations, %.0f min)\n  %s\n",
        target, 100m.purity, 100m.recovery, n[], t / 60, design_string(best.d))
    flush(stdout)
end
close(evals)
close(front)
