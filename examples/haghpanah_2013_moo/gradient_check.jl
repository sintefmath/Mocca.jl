# Check the design gradients of problem.jl against finite differences, at the
# minimum-energy VSA point of Haghpanah et al. (2013), Table 6 (P_L = 0.02 bar),
# and compare its performance with the paper's.
include(joinpath(@__DIR__, "problem.jl"))

d0 = [94.89, 122.84, 189.46, 0.18BAR, 0.02BAR, 0.58, 1.0BAR]
# Stretch the steps of this schedule when the durations change, as the
# duration gradients do
REFERENCE_DURATIONS[] = [T_PRESS, d0[1], d0[2], d0[3]]
warm = WarmStart()
t = @elapsed res0 = evaluate(d0, warm; energy = true)
println("Table 6, P_L = 0.02 bar: ", design_string(d0))
@printf("  Mocca: purity %.2f%%, recovery %.2f%%, energy %.1f kWh/t, productivity %.3f mol/m³/s (%.0f s, %d Newton iterations, %d Jacobians)\n",
    100 * res0.metrics.purity, 100 * res0.metrics.recovery, res0.metrics.energy, res0.metrics.productivity, t, res0.css.cycles, res0.css.jacobian_builds)
println("  paper: purity ≥ 90%, recovery ≥ 90%, energy 148.96 kWh/t, productivity 0.49 mol/m³/s")

t = @elapsed res1 = evaluate(d0, warm; energy = true)
@printf("Second evaluation from the warm start: %.0f s\n", t)

# The energy has a kink in P_H here: the feed compressor only works above its
# 1 bar suction pressure, which is P_H. The adjoint gives the derivative for
# P_H increasing, the side the bounds allow (a forward difference gives
# 1.155e-3), so the central difference for it disagrees.
steps = [1.0, 1.0, 1.0, 0.002BAR, 0.0005BAR, 0.005, 0.01BAR]
for i in 1:7
    h = 1e-2 * steps[i]
    dp = copy(d0); dp[i] += h
    dm = copy(d0); dm[i] -= h
    mp = evaluate(dp, warm; energy = true, gradient = false).metrics
    mm = evaluate(dm, warm; energy = true, gradient = false).metrics
    println(DESIGN_NAMES[i])
    for k in keys(res0.metrics)
        fd = (mp[k] - mm[k]) / (2h)
        @printf("  %-13s adjoint % .6e   finite difference % .6e   rel. diff %.1e\n", k, res0.gradients[k][i], fd,
            abs(res0.gradients[k][i] - fd) / max(abs(fd), 1e-300))
    end
    flush(stdout)
end
