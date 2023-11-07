using Mocca
import Jutul
import JutulDarcy
import MAT


ncells = 200

## Intialise Haghpanah parameters
simulator, state0, parameters =
initialize_Haghpanah_model(forcing_term_coefficient=1.0, ncells = ncells)       

## Setup BCs
pars = simulator.model.system.p


# Set timesteps
t_ads = 5000

t_stage = [t_ads]
numcycles = 1

timesteps = []
sim_forces = []
maxdt = 1



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



Mocca.plot_outlet(simulator.model,states,timesteps)