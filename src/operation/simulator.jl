# Mocca's convergence criteria are relative errors over a step (see
# balances/convergence.jl). Over the many steps and cycles of a cyclic process they
# add up, so they are held tighter than Jutul's default of 1e-3.
const DEFAULT_NONLINEAR_TOLERANCE = 1e-6

# Jutul gives up on a step after halving it 5 times. Opening a valve onto a bed
# at another pressure can take more: the start of a cycle of the dual-adsorbent
# example needs 7. Allowing more cuts changes nothing for a step that converges
# with fewer. Newton's method can also wander so far in such a step that the
# Jacobian is singular; that too cuts the step (`failure_cuts_timestep`) rather
# than ending the simulation.
const DEFAULT_MAX_TIMESTEP_CUTS = 10

function setup_process_simulator(model, state0, parameters;
        timestep_selector_cfg = nothing,
        initial_dt = 1.0,
        nonlinear_tolerance = DEFAULT_NONLINEAR_TOLERANCE,
        max_timestep_cuts = DEFAULT_MAX_TIMESTEP_CUTS,
        failure_cuts_timestep = true,
        kwargs...
    )

    # Set up simulator
    sim = Jutul.Simulator(model; state0 = state0, parameters = parameters)

    # Set up timestep selectors
    t_base = Jutul.TimestepSelector(initial_absolute = initial_dt)
    timesteppers = Vector{Any}()
    push!(timesteppers, t_base)
    if !isnothing(timestep_selector_cfg)
        for (k, v) in pairs(timestep_selector_cfg)
            if model isa Jutul.MultiModel && haskey(model.models, k)
                # Per-unit form, e.g. (Bed = (Pressure = 1e4, y = 0.05),)
                for (k2, v2) in pairs(v)
                    push!(timesteppers, Jutul.VariableChangeTimestepSelector(k2, v2; model = k, relative = false))
                end
            else
                push!(timesteppers, Jutul.VariableChangeTimestepSelector(k, v, relative = false))
            end
        end
    end

    # Set up config
    cfg = Jutul.simulator_config(sim;
        timestep_selectors = timesteppers,
        nonlinear_tolerance = nonlinear_tolerance,
        max_timestep_cuts = max_timestep_cuts,
        failure_cuts_timestep = failure_cuts_timestep,
        kwargs...
    )

    return (sim, cfg)
end

function simulate_process(state0, model, dt, parameters, forces; kwargs...)
    case = MoccaCase(model, dt, forces, state0 = state0, parameters = parameters)
    simulate_process(case; kwargs...)
end

function simulate_process(case::MoccaCase;
    simulator = missing,
    config = missing,
    kwargs...
)
    (; model, forces, state0, parameters, dt) = case

    if ismissing(simulator)
        (sim, cfg_new) = setup_process_simulator(model, state0, parameters; kwargs...)
        config = cfg_new
        extra_arg = NamedTuple()
    else
        sim = simulator
        @assert !ismissing(config) "If simulator is provided, config must also be provided"
        # May have been passed kwarg that should be accounted for
        if length(kwargs) > 0
            config = copy(config)
            for (k, v) in kwargs
                config[k] = v
            end
        end
        extra_arg = (state0 = case.state0, parameters = case.parameters)
    end

    result = Jutul.simulate!(sim, dt;
        config = config,
        forces = forces,
        extra_arg...
    )

    states, timesteps = Jutul.expand_to_ministeps(result)
    return states, timesteps
end

# Simulate `timesteps` from `x`, with every step the solver took, and whether
# it got to the end: Jutul gives up on a step after cutting it, and returns
# the states up to there.
function _simulate_schedule(model, x, parameters, forces, timesteps; kwarg...)
    case = MoccaCase(model, timesteps, forces; state0 = x, parameters = parameters)
    states, ts = simulate_process(case; output_substates = true, kwarg...)
    completed = !isempty(ts) && isapprox(sum(ts), sum(timesteps))
    return (states = states, timesteps = ts, completed = completed)
end
