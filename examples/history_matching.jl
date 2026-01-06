# # History matching of simulation models

# We demonstrate history matching of a Direct Column Breakthrough (DCB) simulation model in Mocca. 
# We leverage the powerful and flexible optimization functionality of Jutul to set up and perform the history matching.
# For more details about the DCB modelling, see the [Simulate DCB](simulate_DCB.md) example.

# Import necessary modules
import Jutul
import Mocca

# We create a function for setting up new simulation cases from the value of the parameter we wish to tune
function setup_case(prm, step_info=missing)

    param_dict_symb = Dict(Symbol(k) => v for (k, v) in prm)
    RealT = valtype(param_dict_symb)
    constants, info = Mocca.parse_input(Mocca.haghpanah_DCB_input(); typeT=RealT)

    for (k, v) in param_dict_symb
        print(k)
        print(v)
        setproperty!(constants, Symbol(k), v)
    end
    case,  = Mocca.setup_mocca_case(constants, info)

    return case
end;

# # Create synthetic reference data
constants_ref, = Mocca.parse_input(Mocca.haghpanah_DCB_input(); typeT=Float64)

prm_ref = Dict("v_feed" => constants_ref.v_feed)
case_ref = setup_case(prm_ref);

# Configure simulator which will be used in the history matching
timestep_selector_cfg = (y=0.01, Temperature=10.0, Pressure=10.0)
sim, cfg = Mocca.setup_process_simulator(case_ref.model, case_ref.state0, case_ref.parameters;
    timestep_selector_cfg = timestep_selector_cfg,
    initial_dt = 1.0,
    output_substates = true,
    info_level = -1
);

# Run reference simulation to generate and generate "ground truth" data from the result
states, timesteps_out = Mocca.simulate_process(case_ref;
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
