function initial_adsorbed_concentration(model, t_init, p_init, y_init)
    # Check that necessary state variables are provided
    !ismissing(t_init) || error("The key Temperature must be present to initialize AdsorbedConcentration")
    !ismissing(p_init) || error("The key Pressure must be present to initialize AdsorbedConcentration")
    !ismissing(y_init) || error("The key y must be present to initialize AdsorbedConcentration")
    R = model.system.p.R
    ncells = Jutul.number_of_cells(model.domain)

    # Ensure vectors have correct length for cell-wise operations
    p_vec = isa(p_init, Number) ? fill(p_init, ncells) : p_init
    temp_vec = isa(t_init, Number) ? fill(t_init, ncells) : t_init

    cTot = p_vec ./ (R * temp_vec)
    c = y_init' .* cTot

    q_init = map(1:ncells) do i
        qstar = compute_equilibrium(model.system, c[i, :], temp_vec[i])
    end
    q_init = stack(q_init) # Convert Vector of SVectors to Matrix
    return q_init
end

function setup_process_model(system::AdsorptionSystem;
    ncells = 100
)
    constants = system.p
    dx = sqrt(pi*constants.r_in^2)
    mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))
    domain = mocca_domain(mesh, system)

    model = Jutul.SimulationModel(domain, system, general_ad = true)
    return model
end

function setup_process_state(model; kwargs...)
    # If not provided, the initial adsorbed concentration is calculated from the other state variables
    if haskey(kwargs, :AdsorbedConcentration)
        q_init = kwargs[:AdsorbedConcentration]
    else
        T0 = get(kwargs, :Temperature, missing)
        p_init = get(kwargs, :Pressure, missing)
        y_init = get(kwargs, :y, missing)
        q_init = initial_adsorbed_concentration(model, T0, p_init, y_init)
    end

    state0 = Jutul.setup_state(model;
        AdsorbedConcentration = q_init,
        kwargs...)

    return state0
end

function setup_process_parameters(model; kwargs...)
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

function setup_forces(model,stage_times,stage_names;
    num_cycles=1,max_dt = 1.0)

    constants = model.system.p;
    ncells = Jutul.number_of_cells(model.domain);
    
    bcs = setup_stage_bcs(constants,stage_times,stage_names,ncells);

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

function setup_stage_bcs(constants,stage_times,stage_names,ncells)


	cycle_time = sum(stage_times);
	step_end = cumsum(stage_times);

	bcs = [];

	for i in 1:length(stage_times)
		if stage_names[i] == "pressurisation"
			bc = Mocca.PressurisationBC(y_feed = constants.y_feed, PH = constants.p_high, PL = constants.p_low,
				λ = constants.λ, T_feed = constants.T_feed, cell_left = 1, cell_right = ncells,
				cycle_time = cycle_time, previous_step_end = i == 1 ? 0 : step_end[i-1])
		elseif stage_names[i] == "adsorption"
			bc = Mocca.AdsorptionBC(y_feed = constants.y_feed, PH = constants.p_high, v_feed = constants.v_feed,
				T_feed = constants.T_feed, cell_left = 1, cell_right = ncells)
		elseif stage_names[i] == "blowdown"
			bc = Mocca.BlowdownBC(PH = constants.p_high, PI = constants.p_intermediate,
				λ = constants.λ, cell_left = 1, cell_right = ncells,
				cycle_time = cycle_time, previous_step_end = i == 1 ? 0 : step_end[i-1])
		elseif stage_names[i] == "evacuation"
			bc = Mocca.EvacuationBC(PL = constants.p_low, PI = constants.p_intermediate,
				λ = constants.λ, cell_left = 1, cell_right = ncells,
				cycle_time = cycle_time, previous_step_end = i == 1 ? 0 : step_end[i-1])
		else
			error("Boundary condition type $(stage_names[i]) not recognized")
		end
		push!(bcs, bc)
	end

	return bcs
end