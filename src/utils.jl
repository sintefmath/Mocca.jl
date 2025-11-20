export mocca_domain
"""
    mocca_domain(mesh::Jutul.CartesianMesh, system::AdsorptionSystem)

Set up a data domain for carbon capture simulation.
"""
function mocca_domain(mesh, system::AdsorptionSystem; kwarg...)
    domain = JutulDarcy.reservoir_domain(mesh, porosity = system.p.Φ, permeability = system.permeability)
    domain[:diffusion_coefficient] = system.dispersion
    domain[:thermal_conductivity] = system.p.K_z
    nc = Jutul.number_of_cells(mesh)

    dx = map(i -> first(Jutul.cell_dims(mesh, i)), 1:nc)
    domain[:dx] = dx

    for (k, v) in kwarg
        domain[k] = v
    end
    return domain
end

# TODO: This causes problem for adjoint simulation. Is it needed?
function Jutul.select_linear_solver(model::AdsorptionModel; kwarg...)
    #return Jutul.LUSolver(; kwarg...)
    return nothing
end

function setup_adsorption_simulator(model, state0, parameters;
        var_tstep_cfg = nothing,
        initial_dt = 1.0,
        kwargs...
    )

    # Set up simulator
    sim = Jutul.Simulator(model; state0 = state0, parameters = parameters)

    # Set up timestep selectors
    t_base = Jutul.TimestepSelector(initial_absolute = initial_dt)
    timesteppers = Vector{Any}()
    push!(timesteppers, t_base)
    if !isnothing(var_tstep_cfg)
        for (k, v) in pairs(var_tstep_cfg)
            t_i = Jutul.VariableChangeTimestepSelector(k, v, relative=false)
            push!(timesteppers, t_i)
        end
    end

    # Set up config
    cfg = Jutul.simulator_config(sim;
        timestep_selectors = timesteppers,
        kwargs...
    )

    return (sim, cfg)
end

function simulate_adsorption(state0, model, dt, parameters, forces; kwargs...)
    case = JutulCase(model, dt, forces, state0 = state0, parameters = parameters)
    simulate_adsorption(case; kwargs...)
end

function simulate_adsorption(case::JutulCase;
    simulator = missing,
    config = missing,
    kwargs...
)
    (; model, forces, state0, parameters, dt) = case

    if ismissing(simulator)
        sim = Jutul.Simulator(model; state0 = state0, parameters = parameters)
        (sim, cfg_new) = setup_adsorption_simulator(model, state0, parameters; kwargs...)
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
