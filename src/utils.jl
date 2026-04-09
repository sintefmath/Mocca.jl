"""
    mocca_domain(mesh; porosity, permeability, r_in, r_out, ...)

Set up a DataDomain for an adsorption column simulation.
Cell-level data (porosity, permeability) goes on Cells.
Column-level data (geometry, material, transport, heat transfer) goes on a
single Column entity.
"""
function mocca_domain(mesh;
    porosity, permeability,
    r_in, r_out,
    diffusion_coefficient, thermal_conductivity,
    adsorbent_density, adsorbent_heat_capacity,
    wall_density, wall_heat_capacity, wall_conductivity,
    fluid_viscosity, fluid_density,
    inner_htc, outer_htc, ambient_temperature,
    kwarg...
)
    domain = JutulDarcy.reservoir_domain(mesh, porosity = porosity, permeability = permeability)

    # Register Column entity (count = 1)
    domain.entities[Column()] = 1

    # Column-level geometry
    domain[:r_in, Column()] = r_in
    domain[:r_out, Column()] = r_out

    # Column-level transport coefficients
    domain[:diffusion_coefficient, Column()] = diffusion_coefficient
    domain[:thermal_conductivity, Column()] = thermal_conductivity

    # Column-level material properties
    domain[:adsorbent_density, Column()] = adsorbent_density
    domain[:adsorbent_heat_capacity, Column()] = adsorbent_heat_capacity
    domain[:wall_density, Column()] = wall_density
    domain[:wall_heat_capacity, Column()] = wall_heat_capacity
    domain[:wall_conductivity, Column()] = wall_conductivity
    domain[:fluid_viscosity, Column()] = fluid_viscosity
    domain[:fluid_density, Column()] = fluid_density
    domain[:inner_htc, Column()] = inner_htc
    domain[:outer_htc, Column()] = outer_htc
    domain[:ambient_temperature, Column()] = ambient_temperature

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
