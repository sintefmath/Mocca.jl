#=
# Gradients at cyclic steady state

This example computes how the performance of the four-stage VSA cycle of
Haghpanah et al. (2013) at cyclic steady state changes with its operating
settings: the feed velocity, the blowdown and evacuation pressures, and the
durations of the four stages. These are the gradients an optimiser needs, and here they
come from adjoints rather than from re-running the cycle to steady state once
per setting.

The cycle maps its initial state `x` to its final state `Φ(x, p)`, where `p`
holds the settings. Cyclic steady state is the fixed point `x* = Φ(x*, p)`,
which moves when `p` changes. For an objective `J` over one cycle, the implicit
function theorem gives

    dJ/dp = ∂J/∂p + μᵀ ∂Φ/∂p,   where   (I − ∂Φ/∂x)ᵀ μ = ∂J/∂x.

Every term comes from adjoint solves over a single cycle: GMRES finds `μ`
from products with `(∂Φ/∂x)ᵀ`, and a last adjoint gives the derivatives with
respect to all parameters and stage settings together.
=#

import Jutul
import Mocca
using CairoMakie

# # 1. The flowsheet
#
# The bed, devices and stages of the [flowsheet VSA](flowsheet_vsa.md)
# example, with 30 cells. Each stage starts with short steps (`first_dt`).
# Without them, the 1 s step at the abrupt start of a stage can converge to a
# slightly different solution for a slightly different start, and the cycle
# map jumps. A gradient of a map that jumps is meaningless.
#
# The stage durations are design variables too. With `reference_durations`,
# a stage of a different duration keeps the steps of the reference schedule,
# stretched in proportion, which is how the duration gradients treat it.
constants = Mocca.HaghpanahConstants{Float64}()
reference = [15.0, 15.0, 30.0, 40.0]
function vsa_case(constants; durations = reference)
    fs, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = 30, stage_durations = durations)
    model = Mocca.setup_flowsheet_model(fs)
    parameters = Mocca.setup_flowsheet_parameters(fs)
    forces, timesteps = Mocca.setup_schedule(fs, model, stages; num_cycles = 1, max_dt = 1.0, first_dt = 0.05,
        reference_durations = reference)
    bed0 = Mocca.setup_process_state(fs[:Bed]; Pressure = constants.p_low, Temperature = 298.15,
        WallTemperature = constants.T_a, y = [1e-10, 1.0 - 1e-10])
    state0 = Mocca.setup_flowsheet_state(fs; Bed = bed0)
    return (; fs, stages, model, parameters, forces, timesteps, state0)
end
case = vsa_case(constants)
(; model, parameters, forces, timesteps) = case;

# # 2. Cyclic steady state
#
# Anderson acceleration gets close in about a hundred cycles. Newton's method
# then converges it tightly in a few iterations, using the Jacobian of the
# cycle map from one adjoint per bed variable. The gradient is only as
# accurate as the steady state, since the slowest mode of this cycle decays by
# only about 1% per cycle.
loose = Mocca.simulate_to_cyclic_steady_state(model, case.state0, parameters, forces, timesteps; tol = 1e-4, max_cycles = 300)
css = Mocca.newton_cyclic_steady_state(model, loose.state0, parameters, forces, timesteps; tol = 1e-10)
println("Anderson: $(loose.cycles) cycles; Newton: $(css.cycles) iterations, cycle change ", css.history)

# # 3. Objectives
#
# Adjoint objectives are sums over the steps of the cycle. The metrics of the
# process are ratios of four such sums: the CO₂ and N₂ in the extract, the CO₂
# fed, and the vacuum pump work for blowdown (through `V_product`) and
# evacuation (through `V_vacuum`).
pump_product = Mocca.VacuumPumpObjective(:V_product)
pump_vacuum = Mocca.VacuumPumpObjective(:V_vacuum)
objectives = (
    co2_extract = Mocca.StreamObjective(:V_vacuum; component = 1),
    n2_extract = Mocca.StreamObjective(:V_vacuum; component = 2),
    co2_feed = Mocca.StreamObjective(:V_feed; component = 1),
    energy = (m, s, dt, info, f) -> pump_product(m, s, dt, info, f) + pump_vacuum(m, s, dt, info, f),
)
g = Mocca.cyclic_steady_state_gradient(model, css.state0, parameters, forces, timesteps, objectives; with_forces = true)
println("GMRES iterations: ", map(r -> r.iterations, g))

# Productivity is the CO₂ produced per kilogram of sorbent per second, so it
# also depends on the cycle time.
M_CO2 = constants.molecularMassOfCO2
m_s = Mocca.sorbent_mass(case.fs[:Bed], parameters[:Bed])
T_cycle = sum(reference)
a, b, f, E = (g[k].objective for k in keys(objectives))
metrics = (purity = a / (a + b), recovery = a / f, energy_kWh_per_t = E / 3.6e6 / (a * M_CO2 / 1000),
    productivity = a / (m_s * T_cycle))
for (k, v) in pairs(metrics)
    println(rpad(k, 18), round(v, sigdigits = 4))
end

# # 4. Gradients with respect to the operating settings
#
# `g.co2_extract.forces[i]` has the layout of the forces of stage `i`, with
# each number replaced by the derivative of the objective with respect to it,
# and `g.co2_extract.durations[i]` is its derivative with respect to the
# duration of stage `i`. The stages are pressurisation, adsorption, blowdown
# and evacuation. Some
# settings appear in two stages: the evacuation pressure PL ends the
# evacuation ramp and starts the pressurisation ramp, and the intermediate
# pressure PI ends the blowdown ramp and starts the evacuation ramp. The
# derivative with respect to such a setting is the sum over the stages that
# use it.
A = π * constants.r_in^2
settings = (
    v_feed = (r -> r.forces[2][:V_feed].law.rate * A, constants.v_feed),
    p_low = (r -> r.forces[4][:V_vacuum].outlet.condition.pressure.stop + r.forces[1][:V_feed].inlet.condition.pressure.start,
        constants.p_low),
    p_intermediate = (r -> r.forces[3][:V_product].outlet.condition.pressure.stop + r.forces[4][:V_vacuum].outlet.condition.pressure.start,
        constants.p_intermediate),
    t_press = (r -> r.durations[1], reference[1]),
    t_ads = (r -> r.durations[2], reference[2]),
    t_blow = (r -> r.durations[3], reference[3]),
    t_evac = (r -> r.durations[4], reference[4]),
)

# The gradients of the ratios follow by the quotient rule, from the gradients
# of the four sums. A stage duration also changes the cycle time, by the same
# amount, which only productivity sees:
function metric_gradients(d, dT)
    da, db, df, dE = d.co2_extract, d.n2_extract, d.co2_feed, d.energy
    return (purity = (da * b - a * db) / (a + b)^2,
        recovery = (da * f - a * df) / f^2,
        energy_kWh_per_t = (dE * a - E * da) / a^2 / 3.6e6 / (M_CO2 / 1000),
        productivity = (da * T_cycle - a * dT) / (m_s * T_cycle^2))
end
is_duration(k) = startswith(String(k), "t_")
sensitivity = NamedTuple{keys(settings)}(Tuple(
    metric_gradients(map(get, g), is_duration(k) ? 1.0 : 0.0) for (k, (get, value)) in pairs(settings)))
# Relative sensitivities, `(x/J) dJ/dx`: the percentage change in each metric
# for a 1% change in each setting.
println(rpad("", 18), join([rpad(k, 18) for k in keys(metrics)]))
for (k, s) in pairs(sensitivity)
    x = settings[k][2]
    println(rpad(k, 18), join([rpad(round(x * s[m] / metrics[m], sigdigits = 3), 18) for m in keys(metrics)]))
end

# # 5. Check against finite differences
#
# Change the feed velocity and the adsorption time by ±0.01%, converge each
# case to steady state with Newton's method from the steady state above, and
# difference the metrics.
function metrics_at(constants; durations = reference)
    c = vsa_case(constants; durations = durations)
    r = Mocca.newton_cyclic_steady_state(c.model, css.state0, c.parameters, c.forces, c.timesteps; tol = 1e-10)
    s, ts = r.states, r.timesteps
    n_ext = Mocca.stream_totals(s, ts, :V_vacuum)
    n_feed = Mocca.stream_totals(s, ts, :V_feed)
    W = Mocca.vacuum_pump_energy(s, ts, :V_product) + Mocca.vacuum_pump_energy(s, ts, :V_vacuum)
    return (purity = Mocca.purity(n_ext), recovery = Mocca.recovery(n_ext, n_feed),
        energy_kWh_per_t = Mocca.specific_energy_kwh_per_tonne(W, n_ext, M_CO2),
        productivity = Mocca.productivity(n_ext, sum(durations), m_s))
end
h = 1e-4 * constants.v_feed
mp = metrics_at(Mocca.HaghpanahConstants{Float64}(v_feed = constants.v_feed + h))
mm = metrics_at(Mocca.HaghpanahConstants{Float64}(v_feed = constants.v_feed - h))
println("feed velocity")
for k in keys(metrics)
    fd = (mp[k] - mm[k]) / (2h)
    println("  ", rpad(k, 18), "adjoint ", round(sensitivity.v_feed[k], sigdigits = 8), "   finite difference ", round(fd, sigdigits = 8))
end
h = 1e-4 * reference[2]
mp = metrics_at(constants; durations = reference .+ [0, h, 0, 0])
mm = metrics_at(constants; durations = reference .- [0, h, 0, 0])
println("adsorption time")
for k in keys(metrics)
    fd = (mp[k] - mm[k]) / (2h)
    println("  ", rpad(k, 18), "adjoint ", round(sensitivity.t_ads[k], sigdigits = 8), "   finite difference ", round(fd, sigdigits = 8))
end

# # 6. Relative sensitivities
fig = Figure(size = (1100, 340))
for (j, m) in enumerate(keys(metrics))
    ax = Axis(fig[1, j], title = String(m), xticks = (1:length(settings), [String(k) for k in keys(settings)]),
        xticklabelrotation = π / 6, ylabel = j == 1 ? "% change per 1% change" : "")
    vals = [settings[k][2] * sensitivity[k][m] / metrics[m] for k in keys(settings)]
    barplot!(ax, 1:length(vals), vals)
    hlines!(ax, [0.0], color = :black, linewidth = 0.5)
end
fig
