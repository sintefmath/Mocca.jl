#=
# Cyclic Vacuum Swing Adsorption simulation
This example shows how to set up and run a cyclic vacuum swing adsorption
simulation as described in [Haghpanah et al. 2013](http://dx.doi.org/10.1021/ie302658y)

This simulation inolves injection of a two component flue gas (CO2 and N2)
into a column of Zeolite 13X. The CO2 is preferentially adsorbed onto the zeolite.
Pressure in the column is then reduced to enable desorption of CO2 for collection.

This is a four stage process comprising:
* Pressurisation: Where the RHS of the column is closed and flue gas is injected from the LHS, at velocity ``v_{feed}``, to bring the column pressure up to ``P_H`` (Pressure High).
* Adsorption: Where both ends of the column are open and flue gas is injected from the RHS at velocity ``v_{feed}``. Pressure at the LHS is ``P_H``.
* Blowdown: Where the LHS of the column is closed and the column is evacuated at ``P_I`` (Pressure Intermediate).
* Evacuation: Where the RHS of the column is closed and the column is evacuated from the LHS at ``P_L`` (Pressure Low).

Each stage is modelled using the same governing equations but with different boundary
conditions. Adsorption onto Zeolite 13X is modelled with a dual-site Langmuir adsorption isotherm.
=#

# First we load the necessary modules
import Jutul: si_unit
import Mocca

# We define parameters, and set up the system and domain as in the [Simulate DCB](simulate_DCB.md) example.
constants = Mocca.HaghpanahConstants{Float64}()
system = Mocca.TwoComponentAdsorptionSystem(constants);

# # Create the model
# Now we can assemble the model which contains the domain and the system of equations.
ncells = 200
model = Mocca.setup_adsorption_model(system; ncells = ncells);

# # Setup the initial state and parameters
# Initial values for pressure and temperature of the system
bar = si_unit(:bar);
P_init = 1*bar;
T_init = 298.15;
Tw_init = constants.T_a;

# To avoid numerical errors we set the initial CO2 concentration to be very
# small and not exactly zero
yCO2_2 = 1e-10
y_init = [yCO2_2, 1.0 - yCO2_2] # [CO2, N2]

state0 = Mocca.setup_adsorption_state(model;
    Pressure = P_init,
    Temperature = T_init,
    WallTemperature = Tw_init,
    y = y_init
)
parameters = Mocca.setup_adsorption_parameters(model);

# # Set up the stage timings and boundary conditions
# Here we have 4 stages and we specify a duration in seconds that we will run each stage.
t_press = 15
t_ads = 15
t_blow = 30
t_evac= 40
stage_times = [t_press, t_ads, t_blow, t_evac];

# Set up cyclic boundary conditions and timesteps for the simulation
sim_forces, timesteps = Mocca.setup_cyclic_forces(model, stage_times;
    num_cycles = 3,
    max_dt = 1.0
);

# # Simulate
# Now we are ready to run the simulation
case = Mocca.MoccaCase(model, timesteps, sim_forces; state0 = state0, parameters = parameters)
states, timesteps_out = Mocca.simulate_adsorption(case;
    output_substates = true,
    info_level = -1
);

# # Visualisation
# We plot primary variables at the outlet through time
outlet_cell = ncells
f_outlet = Mocca.plot_cell(states, model, timesteps_out, outlet_cell)

# We also plot primary variables along the column at the end of the simulation
f_column = Mocca.plot_state(states[end], model)
