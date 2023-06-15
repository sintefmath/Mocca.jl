using Mocca
import Jutul
import JutulDarcy
import MAT

## Read in MATLAB packed problem
datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp.mat"
datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp_adsorption.mat"
datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp_blowdown.mat"
datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp_evacuation.mat"
datapath = "VSA_Comparison_HAG_n30_nc1_julia_comp_all.mat"


# # Intialise parameters from MATLAB
# simulator, state0, parameters =
#     initialize_from_matlab("data/$datapath",
#         forcing_term_coefficient=1.0)

## Intialise Haghpanah parameters
simulator, state0, parameters =
initialize_Haghpanah_model(forcing_term_coefficient=1.0)       

## Setup BCs
pars = simulator.model.system.p

# d_press: Optimize PL [0.05 - 0.5] bar but maybe go even lower (0.01?)
d_press = Mocca.PressurisationBC(y_feed = pars.y_feed, PH = pars.p_high, PL = pars.p_low,
                                λ = pars.λ, T_feed = pars.T_feed, cell_left = 1)
# d_ads: Optimize v_feed []
d_ads = Mocca.AdsorptionBC(y_feed = pars.y_feed, PH = pars.p_high, v_feed = pars.v_feed,
                                T_feed = pars.T_feed, cell_left = 1, cell_right = 30) #TODO: Don't hardcode end cell!                               
# d_blow: Optimize PI [0.05 - 0.5] bar (should be 0.01 bar diff between d_press and d_blow)
d_blow = Mocca.BlowdownBC(PH = pars.p_high, PI = pars.p_intermediate,
                            λ = pars.λ, cell_right = 30) #TODO: Don't hardcode end cell!                               
# d_evac: Optimize PI [0.05 - 0.5] bar (should be 0.01 bar diff between d_press and d_blow)
d_evac = Mocca.EvacuationBC(PL = pars.p_low, PI = pars.p_intermediate,
                            λ = pars.λ, cell_left = 1, cell_right = 30) #TODO: Don't hardcode end cell!                               


# forces = Jutul.setup_forces(simulator.model, bc=d_evac)


# numstages = 4
bcs = [d_press, d_ads, d_blow, d_evac]

# Set timesteps
t_press = 15
t_ads = 15
t_blow = 30
t_evac= 40

t_stage = [t_press, t_ads, t_blow, t_evac]


timesteps = []
sim_forces = []
maxdt = 1.0
for i in eachindex(t_stage)
    numsteps = t_stage[i] / maxdt
    append!(timesteps,repeat([maxdt],Int(floor(numsteps))))
    append!(sim_forces,repeat([Jutul.setup_forces(simulator.model,bc=bcs[i])],Int(floor(numsteps))))
end

# for i in eachindex(t_stage)
#     numsteps = t_stage[i] / maxdt
#     forces = Jutul.setup_forces(simulator.model,bc=bcs[i])
#     for j in 1:numsteps
#         push!(timesteps, maxdt)
#         push!(sim_forces, forces)
#     end
# end


# times_matlab = collect(Iterators.flatten(MAT.matread("data/$datapath")["results"]["time"]))
# times_matlab_zero = zeros(length(times_matlab) + 1)
# times_matlab_zero[2:end] = times_matlab

# timesteps = times_matlab - times_matlab_zero[1:end-1]

# @show timesteps

# sim_forces = repeat([forces], length(timesteps))


states, report = Jutul.simulate(
    simulator,
    timesteps,
    info_level=4,
    forces=sim_forces,
   # max_nonlinear_iterations=0,
   # max_timestep_cuts = 0
)
display(Mocca.plot_states(states))

# display(Mocca.plot_against_matlab_mat(states, 
#     "data/$datapath",
#     cumsum(timesteps)[end], 
#     cumsum(timesteps)))
##
using GLMakie, Jutul
GLMakie.activate!(inline = false)
plot_interactive(simulator.model, states)

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

model = simulator.model
obj = 0.0
for (i, state) in enumerate(states)
    global obj += mocca_purity_objective(model, state, timesteps[i], i, sim_forces[i])
end
obj