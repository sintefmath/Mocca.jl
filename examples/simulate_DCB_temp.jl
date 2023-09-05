# First we load the necessary modules 

# Then we define parameters which we want. We have defined a structure containing
# parameters from Haghpanah et al. 2013 which we can load now.

parameters = HaghpanahParameters()


# As we are doing a DCB simulation we will set the heat transfer coefficient between 
# the column and the wall and the wall and the outside to 0.

parameters.h_in = 0
parameters.h_out = 0


# Then we need to make the model. This model contains information about
# the domain (mesh)) which we will solve the equations over and a information
# about the system of equations which we are solving.

# In this instance we will use the same system as in Haghpanah, which is 
# a two component adsorption system. This system type is associated with 
# the appropriate equations and primary and secondary variables.

system = TwoComponentAdsorptionSystem()

# The column is modelled as a 1D simulation. This is represented in Jutul as a domain
# with one cell in the y and z domains. To get the right areas #WRITE 

dx = sqrt(pi*parameters.r_in^2)
mesh = Jutul.CartesianMesh((ncells, 1, 1), (parameters.L, dx, dx))

# The domain also includes information about transmissibility between cells
# To calculate this we need to know the #WRITE


domain = JutulDarcy.reservoir_domain(mesh, porosity=system.p.Φ, permeability=perm)
domain[:diffusion_coefficient] = axial_dispersion(system)
domain[:thermal_conductivity] = system.p.K_z 


# Now we can assemble the model
model = Jutul.SimulationModel(domain, system)

# The final thing required to create the simulator is the intial state of the system
# #WRITE
# #TODO: can this be put in functions


# initPressure = 0.4*barsa 
initPressure = 1*barsa   #DEBUG  
initT = 298.15

p_init = ones(ncells)*initPressure
temperature_init = ones(ncells)*initT

yCO2 = ones(ncells)*1e-10

y_init = hcat(yCO2, 1 .- yCO2)

cTot = p_init ./ (R * temperature_init)
c = y_init .* cTot
qN2 = ones(ncells)
for i in 1:ncells
    qstar = compute_equilibrium(system, c[i,:], temperature_init[i])
    qN2[i] = qN2[i]*qstar[2]
end

qCO2 = ones(ncells)*0
q_init = hcat(qCO2, qN2)

walltemperature_init = ones(ncells)*parameters.T_a

state0 = Jutul.setup_state(model,
    Pressure = p_init,
    y = y_init',
    AdsorbedConcentration = q_init',
    Temperature = temperature_init,
    WallTemperature = walltemperature_init)

# Make simulator
sim = Jutul.Simulator(model, state0=state0, parameters=parameters)


# Setup schedule

# Simulate