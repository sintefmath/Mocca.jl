# The four-step VSA/PVSA optimisation problems of Haghpanah et al. (2013),
# "Multiobjective optimization of a four-step adsorption process for
# postcombustion CO2 capture via finite volume simulation", Ind. Eng. Chem.
# Res. 52, 4249-4265, section 7, solved with gradients at cyclic steady state
# instead of NSGA-II.
#
# Design variables (Table 5), with the pressurisation step fixed at 20 s:
#
#   t_ads, t_bd, t_evac [s], P_I, P_L [Pa], v_feed [m/s, interstitial], and P_H [Pa]
#   for the PVSA cycle (fixed at 1 bar for VSA).
#
# Included by the drivers in this folder.

using Mocca, Jutul, Printf, LinearAlgebra

# Several optimisations run side by side; one BLAS thread each avoids
# oversubscribing the cores
BLAS.set_num_threads(1)

const BAR = 1e5
const T_PRESS = 20.0
const NCELLS = 30
const MAX_DT = 1.0
const FIRST_DT = 0.05
const M_CO2 = 44.01e-3
const RESULTS = joinpath(@__DIR__, "..", "..", "results", "haghpanah_2013_moo")
mkpath(RESULTS)

const DESIGN_NAMES = (:t_ads, :t_bd, :t_evac, :P_I, :P_L, :v_feed, :P_H)

"""
Bounds of Table 5. `P_L_min` is the lower bound on the evacuation pressure,
which the paper varies between runs; `P_L` fixes it instead. `pvsa` frees P_H.
"""
function design_bounds(; pvsa = false, P_L_min = 0.05BAR, P_L = nothing)
    lo = [20.0, 30.0, 30.0, (pvsa ? 0.13 : 0.05)BAR, P_L_min, 0.1, 1.0BAR]
    hi = [100.0, 200.0, 200.0, (pvsa ? 3.0 : 0.5)BAR, 0.5BAR, 2.0, (pvsa ? 10.0 : 1.0)BAR]
    if !isnothing(P_L)
        lo[5] = hi[5] = P_L
    end
    free = findall(lo .< hi)
    return (lo = lo, hi = hi, free = free)
end

# Design vector from the free variables scaled to [0, 1], and back
full_design(u, b) = (d = copy(b.lo); d[b.free] .= b.lo[b.free] .+ u .* (b.hi[b.free] .- b.lo[b.free]); d)
scaled_design(d, b) = (d[b.free] .- b.lo[b.free]) ./ (b.hi[b.free] .- b.lo[b.free])

"Linear constraint P_I ≥ P_L + 0.01 bar (eq 65) in the scaled variables, as A u ≤ b."
function blowdown_constraint(b)
    a = zeros(length(b.free))
    i_PI, i_PL = findfirst(==(4), b.free), findfirst(==(5), b.free)
    rhs = b.lo[4] - b.lo[5] - 0.01BAR
    a[i_PI] = -(b.hi[4] - b.lo[4])
    isnothing(i_PL) ? (rhs = b.lo[4] - b.lo[5] - 0.01BAR) : (a[i_PL] = b.hi[5] - b.lo[5])
    return (A = reshape(a, 1, :), b = [rhs])
end

# ------------------------------------------------------------------------------
# The cycle for a design
# ------------------------------------------------------------------------------

# Stage durations whose steps every schedule stretches (see
# `reference_durations` in setup_schedule). The duration gradients are those
# of the schedule with each stage's steps stretched in proportion, so an
# optimisation keeps one reference: then the model the optimiser sees is the
# one the gradients describe. `nothing` uses the usual schedule.
const REFERENCE_DURATIONS = Ref{Any}(nothing)

# Rebuild the cycle Jacobian when an iteration gains less than this factor
const REFRESH_RATIO = Ref(0.5)

function vsa_case(d)
    t_ads, t_bd, t_evac, P_I, P_L, v_feed, P_H = d
    base = Mocca.HaghpanahConstants{Float64}()
    constants = Mocca.HaghpanahConstants{Float64}(p_high = P_H, p_intermediate = P_I, p_low = P_L, v_feed = base.Φ * v_feed)
    # The bed of the constants: heat of sorption of CO2 only, and the axial
    # dispersion at their reference velocity rather than v_feed
    k = Mocca.convert_constants(constants)
    fs = Mocca.Flowsheet(k.gas, k.feed)
    Mocca.add_bed!(fs, :Bed, k.adsorbent, k.column; ncells = NCELLS, dispersion_velocity = k.dispersion_velocity)
    Mocca.add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
    Mocca.add_device!(fs, :V_product, :product; at = (:Bed, :top))
    Mocca.add_device!(fs, :V_vacuum, :extract; at = (:Bed, :bottom))
    cycle = Mocca.Cycle([Mocca.Pressurise(T_PRESS), Mocca.Adsorb(t_ads; velocity = v_feed),
        Mocca.Blowdown(t_bd; to = P_I), Mocca.Evacuate(t_evac; to = P_L)])
    # As in the paper, the column starts full of nitrogen at P_L
    process = Mocca.Process(fs, cycle; ramp_rate = k.λ,
        initial_state = Dict(:Bed => (pressure = P_L, temperature = constants.T0, wall_temperature = constants.T_a,
            composition = [1e-10, 1.0 - 1e-10])))
    (; model, parameters) = process
    stages = Mocca.stages(process)
    forces, timesteps = Mocca.setup_schedule(process.flowsheet, model, stages; num_cycles = 1, max_dt = MAX_DT,
        first_dt = FIRST_DT, reference_durations = REFERENCE_DURATIONS[])
    state0 = Mocca.initial_state(process)
    A = π * constants.r_in^2
    V_ads = sum(parameters[:Bed][:SolidVolume])
    return (; process, stages, model, parameters, forces, timesteps, state0, A, V_ads, ε = base.Φ)
end

# Warm start: the last steady state and Jacobian, carried between evaluations
mutable struct WarmStart
    state::Any
    jacobian::Any
end
WarmStart() = WarmStart(nothing, nothing)

# Cyclic steady state for a design, by Newton's method from the warm start
# (reusing its Jacobian as long as it converges well), falling back to
# Anderson acceleration from scratch
function steady_state(case, warm::WarmStart; tol = 1e-9)
    (; model, parameters, forces, timesteps) = case
    if !isnothing(warm.state)
        try
            r = Mocca.newton_cyclic_steady_state(model, Mocca.restart_state(model, warm.state), parameters, forces, timesteps;
                tol = tol, max_iterations = 40, jacobian = warm.jacobian, refresh_ratio = REFRESH_RATIO[])
            if r.converged
                warm.state, warm.jacobian = r.state0, r.jacobian
                return r
            end
        catch e
            e isa InterruptException && rethrow()
        end
    end
    loose = Mocca.simulate_to_cyclic_steady_state(model, case.state0, parameters, forces, timesteps; tol = 1e-3, max_cycles = 500)
    r = Mocca.newton_cyclic_steady_state(model, loose.state0, parameters, forces, timesteps; tol = tol, max_iterations = 8)
    r.converged || @warn "Cyclic steady state not converged" r.history
    warm.state, warm.jacobian = r.state0, r.jacobian
    return r
end

# ------------------------------------------------------------------------------
# Objectives and metrics
# ------------------------------------------------------------------------------

# Sums over one cycle, as adjoint objectives. Energy (eqs 66-69): the feed
# compressor takes flue gas at 1 bar up to the bed inlet pressure in
# pressurisation and adsorption, and vacuum pumps discharge at 1 bar in
# blowdown and evacuation; η = 0.72, γ = 1.4.
const OBJECTIVES = (
    co2_extract = Mocca.StreamObjective(:V_vacuum; component = 1),
    n2_extract = Mocca.StreamObjective(:V_vacuum; component = 2),
    co2_feed = Mocca.StreamObjective(:V_feed; component = 1),
)
const FEED_COMPRESSOR = Mocca.CompressorObjective(:V_feed; suction_pressure = 1BAR)
const BLOWDOWN_PUMP = Mocca.VacuumPumpObjective(:V_product; discharge_pressure = 1BAR)
const EVACUATION_PUMP = Mocca.VacuumPumpObjective(:V_vacuum; discharge_pressure = 1BAR)
energy_objective(m, s, dt, info, f) = FEED_COMPRESSOR(m, s, dt, info, f) + BLOWDOWN_PUMP(m, s, dt, info, f) + EVACUATION_PUMP(m, s, dt, info, f)
const ALL_OBJECTIVES = merge(OBJECTIVES, (energy = energy_objective,))

# Derivative of a cycle sum with respect to the design variables, from its
# gradient with respect to each stage's settings and duration. Stages:
# 1 pressurisation, 2 adsorption, 3 blowdown, 4 evacuation.
function design_gradient(r, case)
    f = r.forces
    ramp(x) = x.outlet.condition.pressure
    g = zeros(7)
    g[1:3] .= r.durations[2:4]
    # P_I ends the blowdown ramp and starts the evacuation ramp
    g[4] = ramp(f[3][:V_product]).stop + ramp(f[4][:V_vacuum]).start
    # P_L ends the evacuation ramp and starts the pressurisation ramp
    g[5] = ramp(f[4][:V_vacuum]).stop + f[1][:V_feed].inlet.condition.pressure.start
    g[6] = f[2][:V_feed].law.rate * case.ε * case.A
    # P_H ends the pressurisation ramp, holds the product end in adsorption
    # and starts the blowdown ramp
    g[7] = f[1][:V_feed].inlet.condition.pressure.stop + f[2][:V_product].outlet.condition.pressure +
        ramp(f[3][:V_product]).start + f[2][:V_feed].inlet.condition.pressure
    return g
end

"""
    evaluate(d, warm; energy = false, gradient = true)

Process performance at cyclic steady state for design `d`, and its gradient
with respect to `d` if `gradient`: purity and recovery (eqs 60-61), energy
[kWh per tonne CO2] (eq 70) and productivity [mol CO2 per m³ adsorbent per s]
(eq 71).
"""
function evaluate(d, warm::WarmStart; energy = false, gradient = true)
    case = vsa_case(d)
    css = steady_state(case, warm)
    objectives = energy ? ALL_OBJECTIVES : OBJECTIVES
    T_cycle = T_PRESS + d[1] + d[2] + d[3]
    if gradient
        g = Mocca.cyclic_steady_state_gradient(case.model, css.state0, case.parameters, case.forces, case.timesteps, objectives;
            with_forces = true, jacobian = css.jacobian, targets = Mocca._no_targets(case.model), rtol = 1e-7)
        sums = map(r -> r.objective, g)
        dsums = map(r -> design_gradient(r, case), g)
    else
        states, steps, step_forces = Mocca._simulate_steps(case.model, css.state0, case.parameters, case.forces, case.timesteps)
        sums = map(G -> Mocca._evaluate_objective(G, case.model, states, steps, step_forces), objectives)
        dsums = map(_ -> zeros(7), objectives)
    end
    a, b, f = sums.co2_extract, sums.n2_extract, sums.co2_feed
    da, db, df = dsums.co2_extract, dsums.n2_extract, dsums.co2_feed
    dT = [1.0, 1.0, 1.0, 0, 0, 0, 0]
    m = (purity = a / (a + b), recovery = a / f, productivity = a / (case.V_ads * T_cycle))
    dm = (purity = (da .* b .- a .* db) ./ (a + b)^2, recovery = (da .* f .- a .* df) ./ f^2,
        productivity = (da .* T_cycle .- a .* dT) ./ (case.V_ads * T_cycle^2))
    if energy
        E, dE = sums.energy, dsums.energy
        to_kwh_per_t = 1 / 3.6e6 / (M_CO2 / 1000)
        m = merge(m, (energy = E / a * to_kwh_per_t,))
        dm = merge(dm, (energy = (dE .* a .- E .* da) ./ a^2 .* to_kwh_per_t,))
    end
    return (metrics = m, gradients = dm, css = css)
end

# ------------------------------------------------------------------------------
# One point of a Pareto front: maximise one metric, with others held above
# targets by quadratic penalties (the paper also penalises its constraints)
# ------------------------------------------------------------------------------

"""
    optimise(u0, bounds, warm; maximise, sign = 1, at_least = (;), weight = 1e3,
        energy = false, max_it = 25, log = nothing)

Maximise `sign` × metric `maximise` (use `sign = -1` to minimise) subject to
each metric `k` in `at_least` being at least `at_least[k]`, over the scaled
design variables, starting from `u0`. The constraints are quadratic penalties
with `weight`, so the result may violate them slightly; it is reported as is.
"""
function optimise(u0, bounds, warm; maximise, sign = 1, at_least = NamedTuple(), scale = Dict{Symbol, Float64}(),
        weight = 1e3, energy = false, max_it = 25, log = nothing)
    best = Ref{Any}(nothing)
    d0 = full_design(u0, bounds)
    REFERENCE_DURATIONS[] = [T_PRESS, d0[1], d0[2], d0[3]]
    function f(u)
        d = full_design(u, bounds)
        t_eval = @elapsed e = evaluate(d, warm; energy = energy)
        s_obj = get(scale, maximise, 1.0)
        v = sign * e.metrics[maximise] / s_obj
        g = sign .* e.gradients[maximise] ./ s_obj
        for (k, target) in pairs(at_least)
            s = get(scale, k, 1.0)
            short = (target - e.metrics[k]) / s
            if short > 0
                v -= weight * short^2
                g = g .+ 2weight * short .* e.gradients[k] ./ s
            end
        end
        gu = g[bounds.free] .* (bounds.hi[bounds.free] .- bounds.lo[bounds.free])
        if isnothing(best[]) || v > best[].v
            best[] = (v = v, u = copy(u), d = d, metrics = e.metrics)
        end
        isnothing(log) || log(d, e, v, t_eval)
        return (v, gu)
    end
    Jutul.unit_box_bfgs(u0, f; maximize = true, lin_ineq = blowdown_constraint(bounds), max_it = max_it,
        grad_tol = 1e-4, obj_change_tol = 1e-5, max_initial_update = 0.05, stepIncreaseTol = 3.0, print = 0)
    return best[]
end

# ------------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------------

function write_row(io, values)
    println(io, join([x isa AbstractFloat ? @sprintf("%.8g", x) : string(x) for x in values], ","))
    flush(io)
end

function design_string(d)
    return @sprintf("t_ads %.1f s, t_bd %.1f s, t_evac %.1f s, P_I %.3f bar, P_L %.3f bar, v %.3f m/s, P_H %.2f bar",
        d[1], d[2], d[3], d[4] / BAR, d[5] / BAR, d[6], d[7] / BAR)
end
