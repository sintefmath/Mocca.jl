function setup_process_simulator(model, state0, parameters;
        timestep_selector_cfg = nothing,
        initial_dt = 1.0,
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
        sim = Jutul.Simulator(model; state0 = state0, parameters = parameters)
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
