# # History matching of simulation models

# We demonstrate history matching of a Direct Column Breakthrough (DCB) simulation model in Mocca. 
# We leverage the powerful and flexible optimization functionality of Jutul to set up and perform the history matching.
# For more details about the DCB modelling, see the [Simulate DCB](simulate_DCB.md) example.

# Import necessary modules
import Jutul
import Mocca

# We create a function for setting up new simulation cases from the value of the parameter we wish to tune
function setup_case(prm, step_info=missing)
    RealT = typeof(prm["v_feed"])
    ncells = 200

    constants = Mocca.HaghpanahConstants{RealT}(h_in = 0.0, h_out = 0.0, v_feed = prm["v_feed"])
    system = Mocca.TwoComponentAdsorptionSystem(constants)
    model = Mocca.setup_adsorption_model(system; ncells = ncells);

    bar = Jutul.si_unit(:bar)
    P_init = 1 * bar
    T_init = 298.15
    Tw_init = constants.T_a

    yCO2_2 = 1e-10
    y_init = [yCO2_2, 1.0 - yCO2_2] # [CO2, N2]

    state0 = Mocca.setup_adsorption_state(model;
        Pressure=P_init,
        Temperature=T_init,
        WallTemperature=Tw_init,
        y=y_init
    )
    parameters = Mocca.setup_adsorption_parameters(model)

    t_ads = 5000.0
    maxdt = 5000.0
    numsteps = Int(floor(t_ads / maxdt))
    timesteps = fill(maxdt, numsteps)

    sim_forces = Mocca.setup_dcb_forces(model)

    case = Mocca.MoccaCase(model, timesteps, sim_forces; state0=state0, parameters=parameters)
    return case
end;

# # Create synthetic reference data
constants_ref = Mocca.HaghpanahConstants{Float64}(h_in=0.0, h_out=0.0)
prm_ref = Dict("v_feed" => constants_ref.v_feed)
case_ref = setup_case(prm_ref);

# Configure simulator which will be used in the history matching
timestep_selector_cfg = (y=0.01, Temperature=10.0, Pressure=10.0)
sim, cfg = Mocca.setup_adsorption_simulator(case_ref.model, case_ref.state0, case_ref.parameters;
    timestep_selector_cfg = timestep_selector_cfg,
    initial_dt = 1.0,
    output_substates = true,
    info_level = -1
);

# Run reference simulation to generate and generate "ground truth" data from the result
states, timesteps_out = Mocca.simulate_adsorption(case_ref;
    simulator = sim,
    config = cfg
);

times_ref = cumsum(timesteps_out)
total_time = times_ref[end]
last_cell_idx = Jutul.number_of_cells(case_ref.model.domain)
qCO2_ref = map(s -> getindex(s[:AdsorbedConcentration], 1, last_cell_idx), states)
qCO2_ref_by_time = Jutul.get_1d_interpolator(times_ref, qCO2_ref);

# # Setting up and solving the optimization problem
# We define a suitable objective function to quantify the match between our simulations and the reference solution.
# Here we choose deviation of adsorbed CO2 in the last grid cell.
function objective_function(model, state, dt, step_info, forces)
    current_time = step_info[:time]
    q_co2 = getindex(state[:AdsorbedConcentration], 1, last_cell_idx)
    q_co2_ref = qCO2_ref_by_time(current_time)
    v = dt/total_time*(q_co2 - q_co2_ref)^2
    return v
end;

# Perturb the known parameter ``v_{feed}`` to form our initial guess for the optimization
prm_guess = Dict("v_feed" => constants_ref.v_feed+0.2)

# Activate ``v_{feed}`` as a free parameter
dprm = Jutul.DictOptimization.DictParameters(prm_guess)
Jutul.DictOptimization.free_optimization_parameter!(dprm, "v_feed"; rel_min = 0.1, rel_max = 10.0)

# Run the optimization
prm_opt = Jutul.DictOptimization.optimize(dprm, objective_function, setup_case;
    config = cfg,
    max_it = 10,
    obj_change_tol = 1e-3
);

# We can see a clear reduction of the objective function value throughout the optimization iterations,
# indicating a close match between the reference solution and our simulation.
f = Mocca.plot_optimization_history(dprm)

# We can look at the optimization result:
dprm

# and see that the value matches the reference parameter value
constants_ref.v_feed
