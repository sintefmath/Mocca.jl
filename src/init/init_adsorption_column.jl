function initial_adsorbed_concentration(model; kwargs...)

    # Extract values from kwargs
    if haskey(kwargs, :Pressure) && haskey(kwargs, :y) && haskey(kwargs, :Temperature)
        p_init = kwargs[:Pressure]
        y_init = kwargs[:y]
        temperature_init = kwargs[:Temperature]
        R = model.system.p.R
        ncells = Jutul.number_of_cells(model.domain)

        # Ensure vectors have correct length for cell-wise operations
        p_vec = isa(p_init, Number) ? fill(p_init, ncells) : p_init
        temp_vec = isa(temperature_init, Number) ? fill(temperature_init, ncells) : temperature_init

        cTot = p_vec ./ (R * temp_vec)
        c = y_init' .* cTot

        q_init = map(1:ncells) do i
            qstar = compute_equilibrium(model.system, c[i,:], temp_vec[i])
        end
        q_init = stack(q_init) # Convert Vector of SVectors to Matrix

    else
        error("Need to set Pressure, y and Temperature in order to determine initial adsorbed concentration in the column")
    end
    return q_init
end

function setup_adsorption_model(system::AdsorptionSystem;
    ncells = 100
)
    constants = system.p
    dx = sqrt(pi*constants.r_in^2)
    mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))
    domain = mocca_domain(mesh, system)

    model = Jutul.SimulationModel(domain, system, general_ad = true)
    return model
end

function setup_adsorption_state(model; kwargs...)
    system = model.system
    g = JutulDarcy.physical_representation(model.data_domain)
    ncells = prod(g.dims)

    q_init = initial_adsorbed_concentration(model; kwargs...)

    state0 = Jutul.setup_state(model;
        AdsorbedConcentration = q_init,
        kwargs...)

    return state0
end

function setup_adsorption_parameters(model; kwargs...)
    system = model.system
    volumes = model.data_domain[:volumes]
    solid_volume = volumes * (1 - system.p.Φ)
    fluid_vol = volumes * system.p.Φ

    parameters = Jutul.setup_parameters(model;
        SolidVolume=solid_volume,
        FluidVolume=fluid_vol,
        kwargs...
    )

    return parameters
end

function setup_dcb_forces(model)
    constants = model.system.p

    ncells = Jutul.number_of_cells(model.domain)
    bc = Mocca.AdsorptionBC(y_feed = constants.y_feed, PH = constants.p_high, v_feed = constants.v_feed,
                                T_feed = constants.T_feed, cell_left = 1, cell_right = ncells);
    dcb_forces = Jutul.setup_forces(model, bc=bc);

    return dcb_forces
end

function setup_cyclic_forces(model, stage_times;
        num_cycles = 3,
        max_dt = 1.0
    )
    cycle_time = sum(stage_times)
    step_end = cumsum(stage_times)

    constants = model.system.p
    ncells = Jutul.number_of_cells(model.domain)

    d_press = Mocca.PressurisationBC(y_feed=constants.y_feed, PH=constants.p_high, PL=constants.p_low,
        λ=constants.λ, T_feed=constants.T_feed, cell_left=1, cell_right=ncells,
        cycle_time=cycle_time, previous_step_end=0)

    d_ads = Mocca.AdsorptionBC(y_feed=constants.y_feed, PH=constants.p_high, v_feed=constants.v_feed,
        T_feed=constants.T_feed, cell_left=1, cell_right=ncells)

    d_blow = Mocca.BlowdownBC(PH=constants.p_high, PI=constants.p_intermediate,
        λ=constants.λ, cell_left=1, cell_right=ncells,
        cycle_time=cycle_time, previous_step_end=step_end[2])


    d_evac = Mocca.EvacuationBC(PL=constants.p_low, PI=constants.p_intermediate,
        λ=constants.λ, cell_left=1, cell_right=ncells,
        cycle_time=cycle_time, previous_step_end=step_end[3])

    bcs = [d_press, d_ads, d_blow, d_evac]

    timesteps = Float64[]
    sim_forces = []

    for j = 1:num_cycles
        for i in eachindex(stage_times)
            numsteps = stage_times[i] / max_dt
            append!(timesteps, repeat([max_dt], Int(floor(numsteps))))
            append!(sim_forces, repeat([Jutul.setup_forces(model, bc=bcs[i])], Int(floor(numsteps))))
        end
    end

    return (sim_forces, timesteps)
end
