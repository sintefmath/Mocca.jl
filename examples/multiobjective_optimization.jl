# # Multi-objective Optimization of a VSA Process
#
# This example demonstrates multi-objective optimization of a cyclic Vacuum Swing
# Adsorption (VSA) simulation model, maximizing both CO₂ **purity** and **recovery**
# as described in [Haghpanah et al. 2013](http://dx.doi.org/10.1021/ie302658y).
#
# We optimize 7 parameters using a hybrid algorithm:
# * **Gradient-based inner loop**: `p_high`, `p_intermediate`, `p_low`, `v_feed` — these
#   parameters appear in the boundary conditions and can be differentiated through using
#   Jutul's adjoint-based optimization.
# * **Derivative-free outer loop**: `t_ads`, `t_blow`, `t_evac` — stage durations affect
#   the simulation structure (timesteps and forces) and cannot be differentiated through,
#   so we use a Nelder–Mead simplex method for global search.
#
# The multi-objective problem is handled via weighted-sum scalarization:
# ```math
# J = \alpha \times \text{purity} + (1 - \alpha) \times \text{recovery}
# ```
# where `α ∈ [0, 1]` controls the trade-off. Different values of `α` trace out
# points on the Pareto front.

# # Import modules

import Mocca
import Jutul
import Jutul: si_unit
import Jutul.DictOptimization: optimize, DictParameters, free_optimization_parameter!

# # Performance metrics
#
# **Recovery** is the fraction of the input CO₂ that is captured during evacuation:
# ```math
# \text{recovery} = \frac{CO_{2,\text{out}}^{\text{evac}}}{CO_{2,\text{in}}^{\text{total}}}
# ```
#
# **Purity** is the mole fraction of CO₂ in the evacuated product stream:
# ```math
# \text{purity} = \frac{CO_{2,\text{out}}^{\text{evac}}}{CO_{2,\text{out}}^{\text{evac}} + N_{2,\text{out}}^{\text{evac}}}
# ```
#
# Both metrics are computed from mass fluxes during the last cycle
# (once the process has approximately reached cyclic steady state).

function compute_purity_recovery(model, states, step_infos, forces, t_stage, num_cycles)
    total_co2_flux_in = 0.0
    total_co2_flux_out = 0.0
    total_flux_out = 0.0

    start_time_last_cycle = sum(t_stage) * (num_cycles - 1)

    for (step_info, state, force_outer) in zip(step_infos, states, forces)
        dt = step_info[:dt]
        time = step_info[:time]

        if time >= start_time_last_cycle
            force = force_outer.bc

            if force isa Mocca.PressurisationBC
                mass_flux = Mocca.mass_flux_left(state, model, time, force)
                total_co2_flux_in -= mass_flux[Mocca.CO2INDEX] * dt
            end

            if force isa Mocca.AdsorptionBC
                mass_flux = Mocca.mass_flux_left(state, model, time, force)
                total_co2_flux_in -= mass_flux[Mocca.CO2INDEX] * dt
            end

            if force isa Mocca.EvacuationBC
                mass_flux = Mocca.mass_flux_left(state, model, time, force)
                total_co2_flux_out -= mass_flux[Mocca.CO2INDEX] * dt
                for i in eachindex(mass_flux)
                    total_flux_out -= mass_flux[i] * dt
                end
            end
        end
    end

    recovery = total_co2_flux_in > 0 ? total_co2_flux_out / total_co2_flux_in : 0.0
    purity   = total_flux_out > 0    ? total_co2_flux_out / total_flux_out    : 0.0
    return (purity = purity, recovery = recovery)
end;

# # Case setup factory
#
# We create a factory function that returns a `setup_case` closure for a given
# set of stage durations. The returned function is called by Jutul's optimizer
# to build a fresh simulation case from the current parameter iterate.

function make_setup_case(t_stage, num_cycles)
    function setup_case(prm, step_info = missing)
        param_dict_symb = Dict(Symbol(k) => v for (k, v) in prm)
        RealT = valtype(param_dict_symb)

        input = Mocca.haghpanah_cyclic_input()
        input["processSpecification"]["stage_durations"] = [15.0, t_stage...]
        input["processSpecification"]["num_cycles"] = num_cycles

        constants, info = Mocca.parse_input(input; typeT = RealT)
        for (k, v) in param_dict_symb
            setproperty!(constants, Symbol(k), v)
        end

        case, = Mocca.setup_mocca_case(constants, info)
        return case
    end
    return setup_case
end;

# # Inner gradient-based optimization
#
# For a fixed set of stage durations, this function optimizes the four
# continuous parameters `p_high`, `p_intermediate`, `p_low`, and `v_feed`
# using Jutul's adjoint-based gradient optimizer.

function inner_gradient_optimization(t_stage;
    num_cycles = 3,
    α = 0.5,
    max_it = 5
)
    bar = si_unit(:bar)

    constants_ref, = Mocca.parse_input(Mocca.haghpanah_cyclic_input(); typeT = Float64)

    prm_guess = Dict(
        "p_high"         => constants_ref.p_high,
        "v_feed"         => constants_ref.v_feed,
        "p_intermediate" => constants_ref.p_intermediate,
        "p_low"          => constants_ref.p_low
    )

    dprm = DictParameters(prm_guess)
    free_optimization_parameter!(dprm, "p_high";         abs_min = 0.5bar, abs_max = 3.0bar)
    free_optimization_parameter!(dprm, "v_feed";         abs_min = 0.1,    abs_max = 2.0)
    free_optimization_parameter!(dprm, "p_intermediate"; abs_min = 0.01bar, abs_max = 0.5bar)
    free_optimization_parameter!(dprm, "p_low";          abs_min = 0.01bar, abs_max = 0.5bar)

    full_t_stage = [15.0, t_stage...]

    function objective_func(model, state0, states, step_infos, forces, input_data)
        metrics = compute_purity_recovery(model, states, step_infos, forces,
                                          full_t_stage, num_cycles)
        return α * metrics.purity + (1 - α) * metrics.recovery
    end

    wrapped_objective = Jutul.WrappedGlobalObjective(objective_func)
    setup_fn = make_setup_case(t_stage, num_cycles)

    try
        prm_opt = optimize(dprm, wrapped_objective, setup_fn;
            max_it = max_it, maximize = true, info_level = -1)
        final_obj = dprm.history.objectives[end]
        return (objective = final_obj, parameters = prm_opt, history = dprm)
    catch e
        @warn "Inner optimization failed for t_stage=$t_stage" exception = e
        return (objective = 0.0, parameters = prm_guess, history = dprm)
    end
end;

# # Nelder–Mead simplex method
#
# A derivative-free optimizer for the outer loop over stage durations.
# The method maintains a simplex of `n+1` vertices in `n`-dimensional space
# and iteratively improves the worst vertex through reflection, expansion,
# contraction, or shrinkage operations.
#
# This is suitable here because:
# * Stage durations are not differentiable through (they change simulation structure).
# * We have only 3 outer parameters, so the simplex has just 4 vertices.
# * Each vertex evaluation runs a full inner gradient optimization.

function nelder_mead_maximize(f, x0, lower, upper;
    max_iter = 20,
    tol = 1e-4,
    nm_α = 1.0,   # reflection coefficient
    nm_γ = 2.0,   # expansion coefficient
    nm_ρ = 0.5,   # contraction coefficient
    nm_σ = 0.5    # shrink coefficient
)
    n = length(x0)
    clamp_to_bounds(x) = clamp.(x, lower, upper)

    # Initialize simplex with n+1 vertices
    simplex = Vector{Vector{Float64}}(undef, n + 1)
    simplex[1] = clamp_to_bounds(collect(Float64, x0))
    for i in 1:n
        xi = copy(simplex[1])
        xi[i] += 0.1 * (upper[i] - lower[i])
        simplex[i + 1] = clamp_to_bounds(xi)
    end

    # Evaluate all vertices
    fvals = [f(x) for x in simplex]
    best_history = Float64[]

    for iter in 1:max_iter
        # Sort: best (highest) first for maximization
        order = sortperm(fvals, rev = true)
        simplex = simplex[order]
        fvals = fvals[order]

        push!(best_history, fvals[1])
        println("Outer iteration $iter: best objective = $(round(fvals[1]; digits=6))")

        # Centroid of the n best vertices (exclude worst)
        x_c = sum(simplex[1:n]) ./ n

        # Reflection
        x_r = clamp_to_bounds(x_c .+ nm_α .* (x_c .- simplex[end]))
        f_r = f(x_r)

        if fvals[1] >= f_r > fvals[n]
            simplex[end] = x_r
            fvals[end] = f_r
        elseif f_r > fvals[1]
            # Expansion
            x_e = clamp_to_bounds(x_c .+ nm_γ .* (x_r .- x_c))
            f_e = f(x_e)
            if f_e > f_r
                simplex[end] = x_e
                fvals[end] = f_e
            else
                simplex[end] = x_r
                fvals[end] = f_r
            end
        else
            # Contraction
            x_co = clamp_to_bounds(x_c .+ nm_ρ .* (simplex[end] .- x_c))
            f_co = f(x_co)
            if f_co > fvals[end]
                simplex[end] = x_co
                fvals[end] = f_co
            else
                # Shrink towards best vertex
                for i in 2:(n + 1)
                    simplex[i] = clamp_to_bounds(simplex[1] .+ nm_σ .* (simplex[i] .- simplex[1]))
                    fvals[i] = f(simplex[i])
                end
            end
        end

        # Convergence check
        frange = maximum(fvals) - minimum(fvals)
        if frange < tol
            println("Converged after $iter iterations (objective range = $frange)")
            break
        end
    end

    best_idx = argmax(fvals)
    return (x = simplex[best_idx], objective = fvals[best_idx], history = best_history)
end;

# # Running the hybrid optimization
#
# We now combine the inner gradient-based optimization with the outer
# Nelder–Mead search. The weight `α` controls the purity–recovery trade-off:
# * `α = 1.0`: maximize purity only
# * `α = 0.0`: maximize recovery only
# * `α = 0.5`: equal weight to both objectives

α_weight = 0.5
num_cycles = 3
max_inner_it = 5
max_outer_it = 10;

# Initial guess and bounds for stage durations `[t_ads, t_blow, t_evac]` in seconds.
# `t_pressurisation = 15 s` is kept fixed.
t_stage_init  = [15.0, 30.0, 40.0]
t_stage_lower = [ 5.0, 10.0, 10.0]
t_stage_upper = [60.0, 120.0, 120.0];

# Define the outer objective: for each candidate `[t_ads, t_blow, t_evac]`,
# run the inner gradient optimization over `[p_high, p_intermediate, p_low, v_feed]`
# and return the best weighted objective.
function outer_objective(t_stage)
    println("  Stage durations: t_ads=$(round(t_stage[1];digits=1))s, " *
            "t_blow=$(round(t_stage[2];digits=1))s, t_evac=$(round(t_stage[3];digits=1))s")
    result = inner_gradient_optimization(t_stage;
        num_cycles = num_cycles, α = α_weight, max_it = max_inner_it)
    println("  → weighted objective = $(round(result.objective; digits=6))")
    return result.objective
end;

# Run the Nelder–Mead outer optimization
result = nelder_mead_maximize(outer_objective, t_stage_init, t_stage_lower, t_stage_upper;
    max_iter = max_outer_it)

# # Retrieve final results
#
# Run one last inner optimization at the best stage durations to obtain the
# optimal pressure and velocity parameters together with purity and recovery.

final = inner_gradient_optimization(result.x;
    num_cycles = num_cycles, α = α_weight, max_it = max_inner_it);

# Display the 7 optimized parameters
println("\n===== Multi-objective Optimization Results =====")
println("Weight α = $α_weight  (objective = α·purity + (1-α)·recovery)\n")
println("Stage durations (global search):")
println("  t_pressurisation = 15.0 s  (fixed)")
println("  t_adsorption     = $(round(result.x[1]; digits=2)) s")
println("  t_blowdown       = $(round(result.x[2]; digits=2)) s")
println("  t_evacuation     = $(round(result.x[3]; digits=2)) s")
println("\nPressure / velocity parameters (gradient-based):")
for (k, v) in final.parameters
    println("  $k = $(round(v; digits=6))")
end
println("\nBest weighted objective = $(round(result.objective; digits=6))")

# We can also plot the outer optimization history
Mocca.plot_optimization_history(final.history;
    yscale = identity, ylabel = "α·Purity + (1-α)·Recovery")
