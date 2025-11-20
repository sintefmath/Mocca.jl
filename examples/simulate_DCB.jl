# # Direct Column Breakthrough simulation

# This example shows how to set up and run a direct column breakthrough (DCB) simulation
# as described in [Haghpanah et al. 2013](http://dx.doi.org/10.1021/ie302658y).
# This simulation inolves injection of a two component flue gas (CO2 and N2)
# into a column of Zeolite 13X initially filled with N2.

# Adsorption onto Zeolite 13X is modelled with a dual-site Langumuir adsorption isotherm.
# Injection flow rate is fixed at the right hand side of the column and the left hand side of the column
# is open.
# There is no heat transfer between the column and the column wall.
#-
# First we load the necessary modules
import Jutul
import Mocca

#=
Then we define parameters which we want. We have defined a structure containing
parameters from Haghpanah et al. 2013 which we load now.
As we are doing a DCB simulation we will set the heat transfer coefficient between
the column and the wall and the wall and the outside to zero.
=#
constants = Mocca.HaghpanahConstants{Float64}(h_in = 0.0, h_out = 0.0);

# We set up a two component adsorption system. This system type is associated
# with the appropriate equations and primary and secondary variables.
system = Mocca.TwoComponentAdsorptionSystem(constants)

# # Define the model
# Next we need to make the model. This model contains information about
# the domain (grid) over which we will solve the equations and information
# about the system of equations which we are solving.
#-
# Because JutulDarcy has its roots in reservoir simulation we formulate
# our velocity equation in the following form:
#
# ```math
# v = -\frac{k}{\mu}(\nabla P)
# ```
#
# where `k` is known as the permeability.
# In this instance we use the pressure drop equation for plug flow in a packed
# bed:
#
# ```math
# v=-\frac{4}{150}\frac{\epsilon^3}{(1-\epsilon)^2} r_{i n}^2 \frac{1}{\mu}(\nabla P)
# ```
#
# The permeability of the system is then given by:
#
# ```math
# k = \frac{4}{150}\frac{\epsilon^3}{(1-\epsilon)^2} r_{i n}^2
# ```
ncells = 200
model = Mocca.setup_adsorption_model(system; ncells = ncells)

# # Setup the initial state
bar = Jutul.si_unit(:bar);
P_init = 1*bar;
T_init = 298.15;
Tw_init = constants.T_a;

# To avoid numerical errors we set the initial CO2 concentration to be very
# small instead of 0.
yCO2_2 = 1e-10
y_init = [yCO2_2, 1.0 - yCO2_2] # [CO2, N2]

state0 = Mocca.setup_adsorption_state(model;
    Pressure = P_init,
    Temperature = T_init,
    WallTemperature = Tw_init,
    y = y_init
)
parameters = Mocca.setup_adsorption_parameters(model)

# # Setup the timestepping and boundary conditions

# For the DCB we are only running the adsorption stage of a VSA process.
# We will use a total time of 5000 seconds with a single report step
t_ads = 5000;
maxdt = 5000.0;
numsteps = Int(floor(t_ads / maxdt));
timesteps = fill(maxdt, numsteps);

# TODO: function to set up adsorption forces?

# We set up boundary conditions for an adsorption stage. AdsorptionBC sets a fixed
# velocity, concentration and temperature at the inlet, and fixed pressure at
# the outlet. By convention we assume the inlet bc is applied on the left hand
# side and the outlet bc is applied on the right hand side.
bc = Mocca.AdsorptionBC(y_feed = constants.y_feed, PH = constants.p_high, v_feed = constants.v_feed,
                                T_feed = constants.T_feed, cell_left = 1, cell_right = ncells);
sim_forces = Jutul.setup_forces(model, bc=bc);


var_tstep_cfg = (y = 0.01, Temperature = 10.0, Pressure = 10.0)

case = Jutul.JutulCase(model, timesteps, sim_forces; state0 = state0, parameters = parameters)
states, sub_timesteps = Mocca.simulate_adsorption(case;
    var_tstep_cfg = var_tstep_cfg,
    output_substates = true,
);

# We plot primary variables at the outlet through time
outlet_cell = ncells
f_outlet = Mocca.plot_cell(states, model, sub_timesteps, outlet_cell)

# We also plot primary variables along the column at the end of the simulation
f_column = Mocca.plot_state(states[end], model)
