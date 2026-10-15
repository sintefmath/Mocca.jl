# # History matching of simulation models

# We demonstrate history matching of a Direct Column Breakthrough (DCB) simulation model in Mocca.
# We leverage the powerful and flexible optimization functionality of Jutul to set up and perform the history matching.
# For more details about the DCB modelling, see the [breakthrough curve](breakthrough.md) example.

# Import necessary modules
import Jutul
import Mocca
using Mocca.Units

# # The breakthrough
# CO₂/N₂ fed at a constant interstitial velocity into an adiabatic column of
# zeolite 13X full of nitrogen, as in the breakthrough example but on a
# coarser grid.
gas = Mocca.haghpanah_2013_gas()
feed = Mocca.Feed(gas; composition = (CO2 = 0.15, N2 = 0.85), temperature = 298.15, pressure = 1bar)
fs = Mocca.Flowsheet(gas, feed)
Mocca.add_bed!(fs, :Bed, Mocca.load_adsorbent("zeolite_13x_haghpanah_2013"), Mocca.haghpanah_2013_column();
    ncells = 50, thermal = Mocca.Adiabatic())
Mocca.add_device!(fs, :V_feed, :feed; at = (:Bed, :bottom))
Mocca.add_device!(fs, :V_product, :product; at = (:Bed, :top))

v_ref = 1.0 # interstitial velocity [m/s]
duration = 1000.0
process = Mocca.Process(fs, Mocca.Cycle([Mocca.Adsorb(duration; velocity = v_ref)]))

# We create a function for setting up new simulation cases from the value of
# the parameter we wish to tune, the feed velocity. The feed is the flow
# controller on `V_feed`, so a case for velocity `v` keeps the model, initial
# state and parameters of the process, and gives the controller the flow rate
# of `v`. The rate is the only thing that depends on `v`, and it accepts the
# automatic differentiation numbers the optimizer passes in.
bed = process.flowsheet[:Bed]
stage = only(Mocca.stages(process))
state0 = Mocca.initial_state(process)
function setup_case(prm, step_info = missing)
    settings = copy(stage.settings)
    rate = Mocca.interstitial_flow_rate(bed, prm["velocity"])
    settings[:V_feed] = merge(settings[:V_feed], (law = Mocca.VolumetricFlow(rate),))
    forces, dt = Mocca.setup_schedule(process.flowsheet, process.model, [Mocca.Stage(stage.name, stage.duration, settings)];
        max_dt = 20.0)
    return Mocca.MoccaCase(process.model, dt, forces; state0 = state0, parameters = process.parameters)
end;

# # Create synthetic reference data
case_ref = setup_case(Dict("velocity" => v_ref));

# Configure simulator which will be used in the history matching. The time
# step selector limits the change of the bed's variables per step.
sim, cfg = Mocca.setup_process_simulator(case_ref.model, case_ref.state0, case_ref.parameters;
    timestep_selector_cfg = (Bed = (y = 0.01, Temperature = 10.0, Pressure = 10.0),),
    initial_dt = 1.0,
    output_substates = true,
    info_level = -1
);

# Run reference simulation to generate "ground truth" data from the result
states, timesteps_out = Mocca.simulate_process(case_ref;
    simulator = sim,
    config = cfg
);

times_ref = cumsum(timesteps_out)
total_time = times_ref[end]
last_cell_idx = Jutul.number_of_cells(bed.domain)
qCO2_ref = map(s -> s[:Bed][:AdsorbedConcentration][1, last_cell_idx], states)
qCO2_ref_by_time = Jutul.get_1d_interpolator(times_ref, qCO2_ref);

# # Setting up and solving the optimization problem
# We define a suitable objective function to quantify the match between our simulations and the reference solution.
# Here we choose deviation of adsorbed CO2 in the last grid cell.
function objective_function(model, state, dt, step_info, forces)
    current_time = step_info[:time]
    q_co2 = state[:Bed][:AdsorbedConcentration][1, last_cell_idx]
    q_co2_ref = qCO2_ref_by_time(current_time)
    v = dt/total_time*(q_co2 - q_co2_ref)^2
    return v
end;

# Perturb the known velocity to form our initial guess for the optimization
prm_guess = Dict("velocity" => 1.5 * v_ref);

# Activate the velocity as a free parameter
dprm = Jutul.DictOptimization.DictParameters(prm_guess)
Jutul.DictOptimization.free_optimization_parameter!(dprm, "velocity"; rel_min = 0.001, rel_max = 100.0)

# Run the optimization
prm_opt = Jutul.DictOptimization.optimize(dprm, objective_function, setup_case;
    config = cfg,
    max_it = 10,
    obj_change_tol = 1e-6,
    solution_history = true
);

# We can see a clear reduction of the objective function value throughout the optimization iterations,
# indicating a close match between the reference solution and our simulation.
f = Mocca.plot_optimization_history(dprm)

# We can look at the optimization result:
dprm

# and see that the value matches the reference velocity
v_ref

# Check size of objective mismatch
dprm.history.objectives
