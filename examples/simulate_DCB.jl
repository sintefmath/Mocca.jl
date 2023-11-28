# # Direct Column Breakthrough simulation
# This example shows how to setup and run a direct column breakthrough simulation
# as described in [Haghpanah et al. 2013](http://dx.doi.org/10.1021/ie302658y)
# This simulation inolves injection of a two component flue gas (CO2 and N2)
# into a column of Zeolite 13X initially filled with N2. 
# Adsorption onto Zeolite 13X is modelled with a dual-site Langumuir adsorption isotherm.
# Injection flow rate is fixed at the rhs of the column and the lhs of the column
# is open.
# There is no heat transfer between the column and the column wall
#-
# First we load the necessary modules
import Jutul
import JutulDarcy
import Mocca

#=
Then we define parameters which we want. We have defined a structure containing
parameters from Haghpanah et al. 2013 which we load now. 
As we are doing a DCB simulation we will set the heat transfer coefficient between 
the column and the wall and the wall and the outside to 0.
=#

constants = Mocca.HaghpanahConstants{Float64}(h_in=0.0,h_out=0.0)

# # Define the model
# Next we need to make the model. This model contains information about
# the domain (grid) over which we will solve the equations and information
# about the system of equations which we are solving.
#-
# Because JutulDarcy has its roots in reservoir simulation we formulate 
# our velocity equation in the following form:
# ``q = -\frac{k}{\mu}\frac{\partial{P}\partial{x}}``
# where `k` is known as the permeability. 
# In this instance we use the pressure drop equation for plug flow in a packed 
# bed:
#
# ``v=-\frac{4}{150}\left(\frac{\epsilon}{1-\epsilon}\right)^2 r_{i n}^2 \frac{1}{\mu}(\nabla P)``
#
# The permeability of the system is then given by:
# ```math
# k = \frac{4}{150}\left(\frac{\epsilon}{1-\epsilon}\right)^2 r_{i n}^2
# ```
# We add the parameters and the desired permeability as an input.
permeability = Mocca.compute_permeability(constants)
axial_dispersion = Mocca.calc_dispersion(constants)



# -
# We setup a two component adsorption system. This system type is associated 
# with the appropriate equations and primary and secondary variables. 

system = Mocca.TwoComponentAdsorptionSystem(; permeability = permeability, dispersion = axial_dispersion, p = constants)

# Jutul uses finite volume discretisation in space. To model a 1D cylindrical column
# we setup a cartesian grid with ncells x 1 x 1 dimensions.
# To ensure we have the correct interface area between cells we set each dimension 
# to the square root of the inner column area.

ncells = 200
dx = sqrt(pi*constants.r_in^2)
mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))


#-
# The domain also contains the mass diffusion coefficient to calculate mass 
# transport between cells and the thermal conductivity to calculate heat 
# transfer.

domain = Mocca.mocca_domain(mesh, system)

# # Create the model
# Now we can assemble the model which contains the domain and the system of equations.
model = Jutul.SimulationModel(domain, system, general_ad = true)

# # Setup the initial state

# The final thing required to create the simulator is the intial state of the 
# system. 

bar = 1e5
P_init = 1.0*bar
T_init = 298.15
Tw_init = constants.T_a

# To avoid numerical errors we set the initial CO2 concentration to be very 
# small instead of 0.

yCO2 = ones(ncells)*1e-10
y_init = hcat(yCO2, 1 .- yCO2)

state0, prm = Mocca.initialise_state_AdsorptionColumn(P_init, T_init, Tw_init, y_init, model)

# # Setup the timestepping and boundary conditions 

# For the DCB we are only running the adsorption stage of a VSA process. 
# We will use a total time of 5000 seconds with 1 second timesteps.

t_ads = 5000
maxdt = 1.0


numsteps = Int(floor(t_ads / maxdt))
timesteps = fill(maxdt, numsteps)



# We setup boundary conditions for an adsorption stage. AdsorptionBC sets a fixed 
# velocity and concentration and temperature at the inlet and fixed pressure at
# the outlet. By convention we assume the inlet bc is applied on the left hand 
# side and the outlet bc is applied on the right hand side. 

bc = Mocca.AdsorptionBC(y_feed = constants.y_feed, PH = constants.p_high, v_feed = constants.v_feed,
                                T_feed = constants.T_feed, cell_left = 1, cell_right = ncells) 


sim_forces = Jutul.setup_forces(model,bc=bc)

# # Simulate
# Now we are ready to run the simulation. For more verbose simulation output 
# info_level can be set to a number higher than 0.


states, report = Jutul.simulate(
    state0,
    model,
    timesteps,
    linear_solver = Jutul.LUSolver(),
    forces=sim_forces,
    parameters = prm,
    info_level = 0
)

# # Plot

# We plot primary variables at the outlet through time
outlet_cell = ncells
f_outlet = Mocca.plot_cell(states,model,timesteps,outlet_cell)
display(f_outlet)

# We plot primary variables along the column at the end of the simulation
f_column = Mocca.plot_state(states[end], model)
display(f_column)