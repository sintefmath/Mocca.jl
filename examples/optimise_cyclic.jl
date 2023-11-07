#=
# Cyclic Vacuum Swing Adsorption simulation
This example shows how to setup and run a cyclic vacuum swing adsorption 
simulation as described in [Haghpanah et al. 2013](dx.doi.org/10.1021/ie302658y)
-
This simulation inolves injection of a two component flue gas (CO2 and N2)
into a column of Zeolite 13X. The CO2 is preferentially adsorbed onto the zeolite.
Pressure in the column is then reduced to enable desorption of CO2 for collection.
-
This is a four stage process comprising:
* Pressurisation: where the rhs of the column is closed and flue gas is injected from
the lhs, at velocity `v_feed`, to bring the column pressure up to `PH` (Pressure High)
* Adsorption: Where both ends of the column are open and flue gas is injected from the 
rhs at velocity `v_feed`. Pressure at the lhs is `PH`
* Blowdown: where the lhs of the column is closed and the column is evacuated at
`PI` (Pressure Intermediate) 
* Evacuation: where the rhs of the column is closed and the column is evacuated from the 
lhs at `PL` (Pressure Low).

Each stage is modelled using the same governing equations but with different boundary
conditions.

Adsorption onto Zeolite 13X is modelled with a dual-site Langumuir adsorption isotherm.
=#

# First we load the necessary modules
import Jutul
import JutulDarcy

import Mocca

# Then we define parameters which we want. We have defined a structure containing
# parameters from Haghpanah et al. 2013 which we load now. 

constants = Mocca.HaghpanahConstants()


# Then we need to make the model. This model contains information about
# the domain (mesh)) which we will solve the equations over and a information
# about the system of equations which we are solving.

# In this instance we will use the same system as in Haghpanah, which is 
# a two component adsorption system. This system type is associated with 
# the appropriate equations and primary and secondary variables. We also 
# add the parameters and the desired velocity model as an input.
permeability = Mocca.compute_permeability(constants)
axial_dispersion = Mocca.calc_dispersion(constants)

system = Mocca.TwoComponentAdsorptionSystem(; permeability = permeability, dispersion = axial_dispersion, p = constants)

# Jutul uses finite volume discretisation in space. To model a 1D cylindrical column
# we setup a cartesian grid with ncells x 1 x 1 dimensions.
# To ensure we have the correct interface area between cells we set 
# ``dx = \sqrt(\pi*r_{in})``
# where `r_{in}` is the inner radius of the column.

ncells = 200

dx = sqrt(pi*constants.r_in^2)
mesh = Jutul.CartesianMesh((ncells, 1, 1), (constants.L, dx, dx))

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



# Set timesteps
t_press = 15
t_ads = 15
t_blow = 30
t_evac= 40

t_stage = [t_press, t_ads, t_blow, t_evac]

cycle_time = sum(t_stage)
step_end = cumsum(t_stage)


# cell_right = simulator.model.data_domain


# d_press: Optimize PL [0.05 - 0.5] bar but maybe go even lower (0.01?)
d_press = Mocca.PressurisationBC(y_feed = constants.y_feed, PH = constants.p_high, PL = constants.p_low,
                                λ = constants.λ, T_feed = constants.T_feed, cell_left = 1, cell_right = ncells,
                                cycle_time = cycle_time, previous_step_end = 0)
# d_ads: Optimize v_feed []
d_ads = Mocca.AdsorptionBC(y_feed = constants.y_feed, PH = constants.p_high, v_feed = constants.v_feed,
                                T_feed = constants.T_feed, cell_left = 1, cell_right = ncells) #TODO: Don't hardcode end cell!                               
# d_blow: Optimize PI [0.05 - 0.5] bar (should be 0.01 bar diff between d_press and d_blow)


d_blow = Mocca.BlowdownBC(PH = constants.p_high, PI = constants.p_intermediate,
                            λ = constants.λ, cell_left = 1, cell_right = ncells,
                            cycle_time = cycle_time, previous_step_end = step_end[2]) 
                            
                            
                            #TODO: Don't hardcode end cell!                               
# d_evac: Optimize PI [0.05 - 0.5] bar (should be 0.01 bar diff between d_press and d_blow)
d_evac = Mocca.EvacuationBC(PL = constants.p_low, PI = constants.p_intermediate,
                            λ = constants.λ, cell_left = 1, cell_right = ncells,
                            cycle_time = cycle_time, previous_step_end = step_end[3]) 
                            #TODO: Don't hardcode end cell!                               


# numstages = 4
bcs = [d_press, d_ads, d_blow, d_evac]






numcycles = 3

timesteps = []
sim_forces = []
maxdt = 1





for j = 1:numcycles
    for i in eachindex(t_stage)
        numsteps = t_stage[i] / maxdt
        append!(timesteps,repeat([maxdt],Int(floor(numsteps))))
        append!(sim_forces,repeat([Jutul.setup_forces(model,bc=bcs[i])],Int(floor(numsteps))))
    end
end



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

# (states_mat, times_mat) = Mocca.get_matlab_states("data/$datapath")

# sims_all = [(states, timesteps)]
# Mocca.plot_pvars_outlet(model, sims_all)


##
# using GLMakie, Jutul
# GLMakie.activate!(inline = false)
# plot_interactive(simulator.model, states)

##
function mocca_purity_objective(model, state, dt, step_no, forces)
    function local_purity(bc::Mocca.EvacuationBC)
        y = state[:y]
        p = state[:Pressure][1]
        
        time = 0.0
        for i in 1:step_no
            time += timesteps[i]
        end
        q = Mocca.pressure_left(bc, time)
        y_co2 = y[1, 1]
        y_n2 = y[2, 1]
        return q*(y_co2/(y_co2 + y_n2))
    end

    function local_purity(bc)
        0.0
    end

    return local_purity(forces.bc)
end

function mocca_recovery_objective(model, state, dt, step_no, forces)
    0
end
# function mocca_objective(model, state, dt, step_no, forces; purity = true, recovery = true)
#     obj = 0
#     if purity
#         obj += mocca_purity_objective(model, state, dt, step_no, forces)
#     end
#     if recovery
#         obj += mocca_recovery_objective(model, state, dt, step_no, forces)
#     end
#     return obj
# end

obj = 0.0
for (i, state) in enumerate(states)
    global obj += mocca_purity_objective(model, state, timesteps[i], i, sim_forces[i])
end
obj