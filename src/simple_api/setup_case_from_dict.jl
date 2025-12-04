function setup_case_from_dict(input_dict::Dict{Symbol, Any})

	using Jutul
	using JutulDarcy
	using Mocca

	constants = parse_constants_from_dict(input_dict)

	# We calculate the permability and dispersion which are used to specify the model
	permeability = Mocca.compute_permeability(constants);
	axial_dispersion = Mocca.calc_dispersion(constants);

	# We set up a two component adsorption system. This system type is associated
	# with the appropriate equations and primary and secondary variables.
	if input_dict[:system_type] == "TwoComponentAdsorptionSystem"
		system = Mocca.TwoComponentAdsorptionSystem(; permeability = permeability, dispersion = axial_dispersion, p = constants);
	else
		error("System type $(input_dict[:system_type]) not recognized")
	end

	# Jutul uses finite volume discretisation in space. To model a 1D cylindrical column
	# we setup a cartesian grid with dimensions ncells x 1 x 1.
	# To ensure we have the correct interface area between cells we set each dimension
	# to the square root of the inner column area. We can then define the simulation domain.

	ncells = input_dict[:ncells];
	dx = sqrt(pi*constants.r_in^2);
	mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx));
	domain = Mocca.mocca_domain(mesh, system);

	# The domain also contains the mass diffusion coefficient to calculate mass
	# transport between cells and the thermal conductivity to calculate heat
	# transfer.

	# # Create the model
	# Now we can assemble the model which contains the domain and the system of equations.
	model = Jutul.SimulationModel(domain, system, general_ad = true);

	# # Setup the initial state

	# The final thing required to create the simulator is the intial state of the
	# system.

	bar = Jutul.si_unit(:bar);
	P_init = input_dict[:P_init]*bar;
	T_init = input_dict[:T_init];
	Tw_init = input_dict[:Tw_init];

	# To avoid numerical errors we set the initial CO2 concentration to be very
	# small instead of 0.

	if sum(input_dict[:y_init]) != 1.0
		error("Initial concentration must sum to 1.0")
	end

	y_init = hcat(input_dict[:y_init]);

	# Now we can initialise the state in the column
	state0, prm = Mocca.initialise_state_AdsorptionColumn(P_init, T_init, Tw_init, y_init, model);

	# # Setup the timestepping and boundary conditions

	# For the DCB we are only running the adsorption stage of a VSA process.
	# We will use a total time of 5000 seconds with a single report step


	stage_names = input_dict[:stages];
	stage_times = input_dict[:t_stage];
	numcycles = input_dict[:numcycles];
	maxdt = input_dict[:maxdt];
	bcs = setup_forces(stage_times,stage_names);

	# Define the full cyclic simulation by stacking subsequent stages in time
	# for a specified number of cycles
	

	sim_forces, timesteps = setup_forces(model,stage_times,stage_names;
	num_cycles=numcycles, max_dt = maxdt);


	# # Simulate

	# Set up timesteppers based on target changes with an initial timestep of 1 day
	t_c = Jutul.VariableChangeTimestepSelector(:y, 0.01, relative = false)
	t_t = Jutul.VariableChangeTimestepSelector(:Temperature, 10.0, relative = false)
	t_p = Jutul.VariableChangeTimestepSelector(:Pressure, 10.0, relative = false);
	t_base = Jutul.TimestepSelector(initial_absolute = 1.0)
	timesteppers = [t_base, t_c, t_t, t_p];

	# We define the simulation setup with initial states and parameters, a linear solver
	# and other configurable options
	sim = Jutul.Simulator(model; state0 = state0, parameters = prm)

	lsolve = Jutul.LUSolver()

	cfg = Jutul.simulator_config(sim;
		timestep_selectors = timesteppers,
		output_substates = true,
		linear_solver = lsolve,
		info_level = input_dict[:info_level],
	);

	return (sim = sim, timesteps = timesteps, sim_forces = sim_forces, cfg = cfg)

end

