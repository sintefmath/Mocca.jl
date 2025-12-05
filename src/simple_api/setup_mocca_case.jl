
using Jutul
using JutulDarcy
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

	bar = Jutul.si_unit(:bar);
	P_init = inputStruct.P_init*bar;
	T0 = inputStruct.T0;
	Tw_init = inputStruct.Tw_init;

	# To avoid numerical errors we set the initial CO2 concentration to be very
	# small instead of 0.

	if sum(inputStruct.y_init) != 1.0 
		error("Initial concentration must sum to 1.0")
	end

	# Now we can initialise the state in the column
	state0 = Mocca.setup_process_state(model;
		Pressure = P_init,
		Temperature = T0,
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
		ts_config = inputStruct.timestep_selectors;
		timesteppers = []
		push!(timesteppers, Jutul.TimestepSelector(initial_absolute = 1.0))
		for (k, v) in ts_config
			push!(timesteppers, Jutul.VariableChangeTimestepSelector(Symbol(k), v["change"], relative = v["relative"]))
		end
	end


	# We define the simulation setup with initial states and parameters, a linear solver
	# and other configurable options
	sim = Jutul.Simulator(model; state0 = state0, parameters = parameters)

	lsolve = Jutul.LUSolver()

	cfg = Jutul.simulator_config(sim;
		timestep_selectors = timesteppers,
		output_substates = true,
		linear_solver = lsolve,
		info_level = inputStruct.info_level,
	);

	case = Mocca.MoccaCase(model, timesteps, sim_forces; 
	state0 = state0, parameters = parameters)

	return case, cfg, sim

end

