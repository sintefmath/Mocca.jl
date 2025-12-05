using Mocca

function setup_mocca_case(inputStruct::ConstantsStruct)

	# We set up a two component adsorption system. This system type is associated
	# with the appropriate equations and primary and secondary variables.
	if inputStruct.system_type == "TwoComponentAdsorptionSystem"
		system = Mocca.TwoComponentAdsorptionSystem(inputStruct);
	else
		error("System type $(inputStruct.system_type) not recognized")
	end

	# Define the model
	model = Mocca.setup_process_model(system; ncells = inputStruct.ncells);

	# # Setup the initial state

	# The final thing required to create the simulator is the intial state of the
	# system.


	P_init = inputStruct.P_init;
	T_init = inputStruct.T0;
	Tw_init = inputStruct.T_a;

	# To avoid numerical errors we set the initial CO2 concentration to be very
	# small instead of 0.

	if sum(inputStruct.y_init) != 1.0 
		error("Initial concentration must sum to 1.0")
	end

	# Now we can initialise the state in the column
	state0 = Mocca.setup_process_state(model;
		Pressure = P_init,
		Temperature = T_init,
		WallTemperature = Tw_init,
		y = inputStruct.y_init
	)

	parameters = Mocca.setup_process_parameters(model);

	# # Setup the timestepping and boundary conditions

	stage_types = inputStruct.stage_types;
	stage_durations = inputStruct.stage_durations;
	num_cycles = inputStruct.num_cycles;
	maxdt = inputStruct.maxdt;

	# Define the full cyclic simulation by stacking subsequent stages in time
	# for a specified number of cycles
	

	sim_forces, timesteps = setup_forces(model,stage_durations,stage_types;
	num_cycles=num_cycles, max_dt = maxdt);


	# # Simulate
	# IF timestepping config is provided then setup timesteppers
	if hasproperty(inputStruct, :timestep_selectors)
		ts_select = inputStruct.timestep_selectors;
		timestep_selector_cfg = (y = ts_select["y"]["change"],
			Temperature = ts_select["Temperature"]["change"],
			Pressure = ts_select["Pressure"]["change"],
		)
	else
		timestep_selector_cfg = nothing
	end


	# We define the simulation setup with initial states and parameters, a linear solver
	# and other configurable options


	case = Mocca.MoccaCase(model, timesteps, sim_forces; 
	state0 = state0, parameters = parameters)

	return case, timestep_selector_cfg, inputStruct.info_level

end

