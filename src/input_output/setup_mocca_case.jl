using Mocca

function setup_mocca_case(constants::ConstantsStruct, info::InfoStruct)

	# We set up a two component adsorption system. This system type is associated
	# with the appropriate equations and primary and secondary variables.
	if info.system_type == "TwoComponentAdsorptionSystem"
		system = Mocca.TwoComponentAdsorptionSystem(constants);
	else
		error("System type $(info.system_type) not recognized")
	end

	# Define the model
	model = Mocca.setup_process_model(system; ncells = info.ncells);
	push!(model.output_variables, :CellDx)

	# # Setup the initial state

	# The final thing required to create the simulator is the intial state of the
	# system.


	P_init = constants.P_init;
	T_init = constants.T0;
	Tw_init = hasproperty(constants, :Tw_init) ? constants.Tw_init : constants.T_a;
	y_init = constants.y_init;

	# To avoid numerical errors we set the initial CO2 concentration to be very
	# small instead of 0.

	if sum(constants.y_init) != 1.0 
		error("Initial concentration must sum to 1.0")
	end

	# Now we can initialise the state in the column
	state0 = Mocca.setup_process_state(model;
		Pressure = P_init,
		Temperature = T_init,
		WallTemperature = Tw_init,
		y = y_init
	)

	parameters = Mocca.setup_process_parameters(model);

	# # Setup the timestepping and boundary conditions

	stage_types = info.stage_types;
	stage_durations = info.stage_durations;
	num_cycles = info.num_cycles;
	maxdt = info.maxdt;

	# Define the full cyclic simulation by stacking subsequent stages in time
	# for a specified number of cycles
	

	sim_forces, timesteps = Mocca.setup_forces(model,stage_durations,stage_types;
	num_cycles=num_cycles, max_dt = maxdt);


	# # Simulate
	# IF timestepping config is provided then setup timesteppers
	if ~isempty(info.timestep_selectors)
		ts_select = info.timestep_selectors;
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

	return case, timestep_selector_cfg, info.info_level

end

