using Mocca
import Jutul
import JutulDarcy


## Setup parameters


parameters = HaghpanahParameters()

system = AdsorptionSystem(forcing_term_coefficient=forcing_term_coefficient, p=parameters)
perm = compute_permeability(system)

barsa = 1e5


## Setup grid

dx = sqrt(pi*parameters.r_in^2)

mesh = Jutul.CartesianMesh((ncells, 1, 1), (parameters.L, dx, dx))

domain = JutulDarcy.reservoir_domain(mesh, porosity=system.p.Φ, permeability=perm)
domain[:diffusion_coefficient] = axial_dispersion(system)
domain[:thermal_conductivity] = system.p.K_z # TODO: Check this and the one above


# TODO: Figure out a better way to compute the volumes
volumes = ones(ncells) * prod(mesh.deltas)
@assert volumes == domain[:volumes]
solid_volume = volumes * (1 - system.p.Φ)
fluid_volume = volumes * system.p.Φ

parameters = Jutul.setup_parameters(model,
    solidVolume=solid_volume,
    fluidVolume=fluid_volume
)





## Setup model



model = Jutul.SimulationModel(domain, system, general_ad=general_ad)


## Setup initial state



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

# Make schedule

# Simulate

# Plot










ncells = 200

## Intialise Haghpanah parameters


simulator, state0, parameters =
initialize_Haghpanah_model(forcing_term_coefficient=1.0, ncells = ncells)       


sim = Jutul.Simulator(model, state0=state0, parameters=parameters)

## Setup BCs
pars = simulator.model.system.p


# Set timesteps
t_ads = 300

t_stage = [t_ads]
numcycles = 1

timesteps = []
sim_forces = []
maxdt = 1.0



# d_ads: Optimize v_feed []
d_ads = Mocca.AdsorptionBC(y_feed = pars.y_feed, PH = pars.p_high, v_feed = pars.v_feed,
                                T_feed = pars.T_feed, cell_left = 1, cell_right = ncells) #TODO: Don't hardcode end cell!                               


bcs = [d_ads]

for j = 1:numcycles
    for i in eachindex(t_stage)
        numsteps = t_stage[i] / maxdt
        append!(timesteps,repeat([maxdt],Int(floor(numsteps))))
        append!(sim_forces,repeat([Jutul.setup_forces(simulator.model,bc=bcs[i])],Int(floor(numsteps))))
    end
end


states, report = Jutul.simulate(
    simulator,
    timesteps,
    info_level=0,
    forces=sim_forces,
    # max_nonlinear_iterations=0,
    # max_timestep_cuts = 0
)

Mocca.plot_outlet(model,states,timesteps)