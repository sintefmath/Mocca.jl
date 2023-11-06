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

constants = Mocca.HaghpanahConstants(h_in=0,h_out=0)

# # Define the model
# Next we need to make the model. This model contains information about
# the domain (mesh) which we will solve the equations over and a information
# about the system of equations which we are solving.
#-
# Because JutulDarcy has its roots in reservoir simulation we need to store some paramters 



# we need to formulate 
# our velocity equation in the following form:
# ``q = -\frac{k}{\mu}\frac{\partial{P}\partial{x}}``
# where `k` is known as the permeability. 
# In this instance we use the pressure drop equation for plug flow in a packed 
# bed:
#
# ``v=-\frac{4}{150}\left(\frac{\epsilon}{1-\epsilon}\right)^2 r_{i n}^2 \frac{1}{\mu}(\nabla P)``
#
# In this case the permeability of the system is given by:
# ```math
# k = \frac{4}{150}\left(\frac{\epsilon}{1-\epsilon}\right)^2 r_{i n}^2
# ```
# We add the parameters and the desired permeability as an input.
permeability = Mocca.compute_permeability(constants)
axial_dispersion = Mocca.calc_dispersion(constants)



# -
# In this instance we will use the same system as in Haghpanah, which is 
# a two component adsorption system. This system type is associated with 
# the appropriate equations and primary and secondary variables. 

system = Mocca.TwoComponentAdsorptionSystem(; permeability = permeability, dispersion = axial_dispersion, p = constants)

# Jutul uses finite volume discretisation in space. To model a 1D cylindrical column
# we setup a cartesian grid with ncells x 1 x 1 dimensions.
# To ensure we have the correct interface area between cells we set 
# ``dx = \sqrt(\pi*r_{in})``
# where `r_{in}` is the inner radius of the column.

ncells = 10

dx = sqrt(pi*constants.r_in^2)
mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))


#-
# The domain also contains the mass diffusion coefficient to calculate mass 
# transport between cells and the thermal conductivity to calculate heat 
# transfer

domain = JutulDarcy.reservoir_domain(mesh, porosity = constants.Φ, permeability = system.permeability)
domain[:diffusion_coefficient] = system.dispersion
domain[:thermal_conductivity] = constants.K_z  #TODO : do we need this here? And what about line above?

# # Create the model
# Now we can assemble the model which contains the domain and the system of equations.
model = Jutul.SimulationModel(domain, system, general_ad = true)

# # Setup the initial state

# The final thing required to create the simulator is the intial state of the system
# #WRITE
# #TODO: can this be put in functions
barsa = 1e5 #TODO: see if this exists
P_init = 1.0*barsa 
T_init = 298.15
Tw_init = constants.T_a

yCO2 = ones(ncells)*1e-10
y_init = hcat(yCO2, 1 .- yCO2)

state0, prm = Mocca.initialise_state_AdsorptionColumn(P_init, T_init, Tw_init, y_init, model)

# # Setup the timestepping and boundary conditions 

# For the DCB we are only running the adsorption stage of a VSA process. 
# We will use a total time of 5000 seconds with 1 second timesteps.

t_ads = 50
timesteps = []
sim_forces = []
maxdt = 1.0

#WRITE : write

bc = Mocca.AdsorptionBC(y_feed = constants.y_feed, PH = constants.p_high, v_feed = constants.v_feed,
                                T_feed = constants.T_feed, cell_left = 1, cell_right = ncells) 


numsteps = Int(floor(t_ads / maxdt))
timesteps = fill(maxdt, numsteps)
# append!(sim_forces,repeat([],Int(floor(numsteps))))

sim_forces = Jutul.setup_forces(model,bc=bc)

# # Simulate
#WRITE 
states, report = Jutul.simulate(
    state0,
    model,
    timesteps,
    forces=sim_forces,
    parameters = prm,
    info_level = 0
)


# # Plot
#WRITE : 

outlet_cell = ncells
Mocca.plot_cell(states,model,timesteps,outlet_cell)


Mocca.plot_state(states[end], model)


# Mocca.plot_outlet(model,states,timesteps)
# pvars = [:Pressure]
# cell = ncells
# #Mocca.plot_cell(model,states,timesteps,cell,pvars)

# nc = size(states[end][:Pressure], 1)
# x = model.data_domain[:cell_centroids][1,:]
# t = cumsum(timesteps)

# # plot pressure
# using CairoMakie
# symbol = :Pressure
# cell = nc
# f = CairoMakie.Figure()
# ax = CairoMakie.Axis(f, title=String(symbol), xlabel=CairoMakie.L"t")
# CairoMakie.lines!(ax, t, Float64.([result[symbol][cell] for result in states]), color=:darkgray)
# display(f)