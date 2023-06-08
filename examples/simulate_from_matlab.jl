using Mocca
import Jutul
import JutulDarcy
import MAT

## Read in MATLAB packed problem
datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp.mat"

## Intialise parameters from MATLAB
simulator, state0, parameters =
    initialize_from_matlab("data/$datapath",
        forcing_term_coefficient=1.0)

## Setup BCs
pars = simulator.model.system.p


d = Mocca.PressurisationBC(y_feed = pars.y_feed, PH = pars.p_high, PL = pars.p_low,
                                λ = pars.λ, T_feed = pars.T_feed, cell_left = 1)

d = Mocca.AdsorptionBC(y_feed = pars.y_feed, PH = pars.p_high, v_feed = pars.v_feed,
                                T_feed = pars.T_feed, cell_left = 1, cell_right = 30) #TODO: Don't hardcode end cell!                               

forces = Jutul.setup_forces(simulator.model, bc=d)


times_matlab = collect(Iterators.flatten(MAT.matread("data/VSA_Comparison_HAG_n30_nc1_julia_comp.mat")["results"]["time"]))
times_matlab_zero = zeros(length(times_matlab) + 1)
times_matlab_zero[2:end] = times_matlab

timesteps = times_matlab - times_matlab_zero[1:end-1]
@show timesteps

sim_forces = repeat([forces], length(timesteps))


states, report = Jutul.simulate(
    simulator,
    timesteps,
    info_level=0,
    forces=sim_forces,
     max_nonlinear_iterations=0,
    max_timestep_cuts = 0
)#000)
##
display(Mocca.plot_states(states))
##
#display(Mocca.plot_outlet(cumsum(timesteps), states))
display(Mocca.plot_against_matlab_mat(states, 
    # "data/VSA_Comparison_HAG_n30_nc1_julia_comp.mat", 
    "data/$datapath",
    cumsum(timesteps)[end], 
    cumsum(timesteps)))
##
# plot_states(states)
