using Mocca
import Jutul
import JutulDarcy
import MAT

## Read in MATLAB packed problem
datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp.mat"
datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp_adsorption.mat"
datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp_blowdown.mat"
# datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp_adsorptionbc_sloping.mat"

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

d = Mocca.BlowdownBC(PH = pars.p_high, PI = pars.p_intermediate,
                            λ = pars.λ, cell_right = 30) #TODO: Don't hardcode end cell!                               

forces = Jutul.setup_forces(simulator.model, bc=d)


times_matlab = collect(Iterators.flatten(MAT.matread("data/$datapath")["results"]["time"]))
times_matlab_zero = zeros(length(times_matlab) + 1)
times_matlab_zero[2:end] = times_matlab

timesteps = times_matlab - times_matlab_zero[1:end-1]

@show timesteps

sim_forces = repeat([forces], length(timesteps))


states, report = Jutul.simulate(
    simulator,
    timesteps,
    info_level=4,
    forces=sim_forces,
   # max_nonlinear_iterations=0,
   # max_timestep_cuts = 0
)
display(Mocca.plot_states(states))

display(Mocca.plot_against_matlab_mat(states, 
    "data/$datapath",
    cumsum(timesteps)[end], 
    cumsum(timesteps)))
