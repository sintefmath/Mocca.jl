function haghpanah_DCB(constants)

    
    # We calculate the permability and dispersion which are used to specify the model
    permeability = Mocca.compute_permeability(constants);
    axial_dispersion = Mocca.calc_dispersion(constants);

    # We set up a two component adsorption system. This system type is associated
    # with the appropriate equations and primary and secondary variables.

    system = Mocca.TwoComponentAdsorptionSystem(; permeability = permeability, dispersion = axial_dispersion, p = constants);

    # Jutul uses finite volume discretisation in space. To model a 1D cylindrical column
    # we setup a cartesian grid with dimensions ncells x 1 x 1.
    # To ensure we have the correct interface area between cells we set each dimension
    # to the square root of the inner column area. We can then define the simulation domain.

    ncells = 200;
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
    P_init = 1*bar;
    T_init = 298.15;
    Tw_init = constants.T_a;

    # To avoid numerical errors we set the initial CO2 concentration to be very
    # small instead of 0.

    yCO2 = fill(1e-10, ncells)
    y_init = hcat(yCO2, 1 .- yCO2);

    # Now we can initialise the state in the column
    state0, prm = Mocca.initialise_state_AdsorptionColumn(P_init, T_init, Tw_init, y_init, model);

    # # Setup the timestepping and boundary conditions

    # For the DCB we are only running the adsorption stage of a VSA process.
    # We will use a total time of 5000 seconds with a single report step

    t_ads = 5000;
    maxdt = 5000.0;
    numsteps = Int(floor(t_ads / maxdt));
    timesteps = fill(maxdt, numsteps);

    # We set up boundary conditions for an adsorption stage. AdsorptionBC sets a fixed
    # velocity, concentration and temperature at the inlet, and fixed pressure at
    # the outlet. By convention we assume the inlet bc is applied on the left hand
    # side and the outlet bc is applied on the right hand side.

    bc = Mocca.AdsorptionBC(y_feed = constants.y_feed, PH = constants.p_high, v_feed = constants.v_feed,
                                    T_feed = constants.T_feed, cell_left = 1, cell_right = ncells);

    sim_forces = Jutul.setup_forces(model, bc=bc);

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
        info_level = 0
    );

    return (sim=sim, timesteps=timesteps, sim_forces=sim_forces, cfg=cfg)
end