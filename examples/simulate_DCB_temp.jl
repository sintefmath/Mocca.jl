# # Direct Column Breakthrough simulation
# This example shows how to setup and run a direct column breakthrough simulation
# as described in [Haghpanah et al. 2013](dx.doi.org/10.1021/ie302658y)
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

# Then we define parameters which we want. We have defined a structure containing
# parameters from Haghpanah et al. 2013 which we load now. 
# As we are doing a DCB simulation we will set the heat transfer coefficient between 
# the column and the wall and the wall and the outside to 0.

parameters = Mocca.HaghpanahParameters(h_in=0,h_out=0)



# Then we need to make the model. This model contains information about
# the domain (mesh)) which we will solve the equations over and a information
# about the system of equations which we are solving.

# In this instance we will use the same system as in Haghpanah, which is 
# a two component adsorption system. This system type is associated with 
# the appropriate equations and primary and secondary variables. We also 
# add the parameters and the desired velocity model as an input.
permeability = Mocca.compute_permeability(parameters)
axial_dispersion = Mocca.calc_dispersion(parameters)

system = Mocca.TwoComponentAdsorptionSystem(; permeability = permeability, dispersion = axial_dispersion)

# Jutul uses finite volume discretisation in space. To model a 1D cylindrical column
# we setup a cartesian grid with ncells x 1 x 1 dimensions.
# To ensure we have the correct interface area between cells we set 
# ``dx = \sqrt(\pi*r_{in})``
# where `r_{in}` is the inner radius of the column.

ncells = 200

dx = sqrt(pi*parameters.r_in^2)
mesh = Jutul.CartesianMesh((ncells, 1, 1), (parameters.L, dx, dx))

# Because JutulDarcy has its roots in reservoir simulation we need to formulate 
# our velocity equation in the following form:
# ``q = -\frac{k}{\mu}\frac{\partial{P}\partial{x}}``
# where `k` is known as the permeability. 
# In this instance we use the pressure drop equation for plug flow in a packed 
# bed. The velocity model specified when instantiating the system automatically
# calculates the permeability required for input to the domain.
#-
# The domain also contains the mass diffusion coefficient to calculate mass 
# transport between cells and the thermal conductivity to calculate heat 
# transfer

domain = JutulDarcy.reservoir_domain(mesh, porosity = parameters.Φ, permeability = system.permeability)
domain[:diffusion_coefficient] = system.dispersion
domain[:thermal_conductivity] = parameters.K_z  #TODO : do we need this here? And what about line above?

# # Create the model
# Now we can assemble the model which contains the domain and the system of equations.
model = Jutul.SimulationModel(domain, system, general_ad = true)

# # Setup the initial state

# The final thing required to create the simulator is the intial state of the system
# #WRITE
# #TODO: can this be put in functions
barsa = 1e5 #TODO: see if this exists
P_init = 0.4*barsa 
T_init = 298.15
Tw_init = parameters.T_a

yCO2 = ones(ncells)*1e-10
y_init = hcat(yCO2, 1 .- yCO2)

state0, prm = Mocca.initialise_state_AdsorptionColumn(P_init, T_init, Tw_init, y_init, model)

# # Make the simulator and setup the timestepping and boundary conditions 

# We combine everything into a Jutul.Simulator instance.

sim = Jutul.Simulator(model, state0=state0, parameters=prm)


# # Setup schedule

# For the DCB we are only running the adsorption stage of a VSA process. 
# We will use a total time of 5000 seconds with 1 second timesteps.

t_ads = 5000
timesteps = []
sim_forces = []
maxdt = 1.0

#WRITE : write

bc = Mocca.AdsorptionBC(y_feed = parameters.y_feed, PH = parameters.p_high, v_feed = parameters.v_feed,
                                T_feed = parameters.T_feed, cell_left = 1, cell_right = ncells) 


numsteps = t_ads / maxdt
append!(timesteps,repeat([maxdt],Int(floor(numsteps))))
append!(sim_forces,repeat([Jutul.setup_forces(sim.model,bc=bc)],Int(floor(numsteps))))


# # Simulate
#WRITE 
states, report = Jutul.simulate(
    sim,
    timesteps,
    forces=sim_forces,
)


# # Plot
#WRITE : 
Mocca.plot_outlet(model,states,timesteps)
