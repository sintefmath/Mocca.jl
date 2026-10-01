#=
# Gradients at cyclic steady state

This example computes how the performance of the four-stage VSA cycle of
Haghpanah et al. (2013) at cyclic steady state changes with its operating
settings: the feed velocity, the blowdown and evacuation pressures, and the
valve conductances. These are the gradients an optimiser needs, and here they
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
constants = Mocca.HaghpanahConstants{Float64}()
function vsa_case(constants)
    fs, stages = Mocca.four_stage_vsa_flowsheet(constants; ncells = 30)
    model = Mocca.setup_flowsheet_model(fs)
    parameters = Mocca.setup_flowsheet_parameters(fs)
    forces, timesteps = Mocca.setup_schedule(fs, model, stages; num_cycles = 1, max_dt = 1.0, first_dt = 0.05)
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

M_CO2 = constants.molecularMassOfCO2
a, b, f, E = (g[k].objective for k in keys(objectives))
metrics = (purity = a / (a + b), recovery = a / f, energy_kWh_per_t = E / 3.6e6 / (a * M_CO2 / 1000))
for (k, v) in pairs(metrics)
    println(rpad(k, 18), round(v, sigdigits = 4))
end

# # 4. Gradients with respect to the operating settings
#
# `g.co2_extract.forces[i]` has the layout of the forces of stage `i`, with
# each number replaced by the derivative of the objective with respect to it.
# The stages are pressurisation, adsorption, blowdown and evacuation. Some
# settings appear in two stages: the evacuation pressure PL ends the
# evacuation ramp and starts the pressurisation ramp, and the intermediate
# pressure PI ends the blowdown ramp and starts the evacuation ramp. The
# derivative with respect to such a setting is the sum over the stages that
# use it.
A = π * constants.r_in^2
settings = (
    v_feed = (gf -> gf[2][:V_feed].law.rate * A, constants.v_feed, "feed velocity [m/s]"),
    p_low = (gf -> gf[4][:V_vacuum].outlet.condition.pressure.stop + gf[1][:V_feed].inlet.condition.pressure.start,
        constants.p_low, "evacuation pressure PL [Pa]"),
    p_intermediate = (gf -> gf[3][:V_product].outlet.condition.pressure.stop + gf[4][:V_vacuum].outlet.condition.pressure.start,
        constants.p_intermediate, "blowdown pressure PI [Pa]"),
    C_vacuum = (gf -> gf[4][:V_vacuum].law.conductance, forces[end][:V_vacuum].law.conductance, "vacuum valve conductance"),
)

# The gradients of the ratios follow by the quotient rule, from the gradients
# of the four sums:
function metric_gradients(d)
    da, db, df, dE = d.co2_extract, d.n2_extract, d.co2_feed, d.energy
    return (purity = (da * b - a * db) / (a + b)^2,
        recovery = (da * f - a * df) / f^2,
        energy_kWh_per_t = (dE * a - E * da) / a^2 / 3.6e6 / (M_CO2 / 1000))
end
sensitivity = map(settings) do (get, value, label)
    metric_gradients(map(r -> get(r.forces), g))
end
# Relative sensitivities, `(x/J) dJ/dx`: the percentage change in each metric
# for a 1% change in each setting.
println(rpad("", 18), join([rpad(k, 18) for k in keys(metrics)]))
for (k, s) in pairs(sensitivity)
    x = settings[k][2]
    println(rpad(k, 18), join([rpad(round(x * s[m] / metrics[m], sigdigits = 3), 18) for m in keys(metrics)]))
end

# The vacuum valve has the half-cell conductance of the bed end, so it holds
# the bed bottom at the vacuum line's pressure whatever its exact value, and
# the metrics hardly depend on it. The pressures set where the cycle operates.

# # 5. Check against finite differences
#
# Change the feed velocity by ±0.01%, converge each case to steady state with
# Newton's method from the steady state above, and difference the metrics.
function metrics_at(constants)
    c = vsa_case(constants)
    r = Mocca.newton_cyclic_steady_state(c.model, css.state0, c.parameters, c.forces, c.timesteps; tol = 1e-10)
    s, ts = r.states, r.timesteps
    n_ext = Mocca.stream_totals(s, ts, :V_vacuum)
    n_feed = Mocca.stream_totals(s, ts, :V_feed)
    W = Mocca.vacuum_pump_energy(s, ts, :V_product) + Mocca.vacuum_pump_energy(s, ts, :V_vacuum)
    return (purity = Mocca.purity(n_ext), recovery = Mocca.recovery(n_ext, n_feed),
        energy_kWh_per_t = Mocca.specific_energy_kwh_per_tonne(W, n_ext, M_CO2))
end
h = 1e-4 * constants.v_feed
mp = metrics_at(Mocca.HaghpanahConstants{Float64}(v_feed = constants.v_feed + h))
mm = metrics_at(Mocca.HaghpanahConstants{Float64}(v_feed = constants.v_feed - h))
for k in keys(metrics)
    fd = (mp[k] - mm[k]) / (2h)
    println(rpad(k, 18), "adjoint ", round(sensitivity.v_feed[k], sigdigits = 8), "   finite difference ", round(fd, sigdigits = 8))
end

# # 6. Relative sensitivities
fig = Figure(size = (900, 320))
for (j, m) in enumerate(keys(metrics))
    ax = Axis(fig[1, j], title = String(m), xticks = (1:length(settings), [String(k) for k in keys(settings)]),
        xticklabelrotation = π / 6, ylabel = j == 1 ? "% change per 1% change" : "")
    vals = [settings[k][2] * sensitivity[k][m] / metrics[m] for k in keys(settings)]
    barplot!(ax, 1:length(vals), vals)
    hlines!(ax, [0.0], color = :black, linewidth = 0.5)
end
fig
