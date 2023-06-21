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
maxdt = 0.5


times_matlab = collect(Iterators.flatten(MAT.matread("data/$datapath")["results"]["time"]))
ts_matlab = diff(vcat(0,times_matlab))
stepStart = [1,24,47,85]
stepEnd = [23,46,84,132]
numSteps = [23,23,38,48]



# cell_right = simulator.model.data_domain


# d_press: Optimize PL [0.05 - 0.5] bar but maybe go even lower (0.01?)
d_press = Mocca.PressurisationBC(y_feed = pars.y_feed, PH = pars.p_high, PL = pars.p_low,
                                λ = pars.λ, T_feed = pars.T_feed, cell_left = 1, cell_right = ncells,
                                t_stage = t_stage)
# d_ads: Optimize v_feed []
d_ads = Mocca.AdsorptionBC(y_feed = pars.y_feed, PH = pars.p_high, v_feed = pars.v_feed,
                                T_feed = pars.T_feed, cell_left = 1, cell_right = ncells) #TODO: Don't hardcode end cell!                               
# d_blow: Optimize PI [0.05 - 0.5] bar (should be 0.01 bar diff between d_press and d_blow)
d_blow = Mocca.BlowdownBC(PH = pars.p_high, PI = pars.p_intermediate,
                            λ = pars.λ, cell_left = 1, cell_right = ncells,
                            t_stage = t_stage) #TODO: Don't hardcode end cell!                               
# d_evac: Optimize PI [0.05 - 0.5] bar (should be 0.01 bar diff between d_press and d_blow)
d_evac = Mocca.EvacuationBC(PL = pars.p_low, PI = pars.p_intermediate,
                            λ = pars.λ, cell_left = 1, cell_right = ncells,
                            t_stage = t_stage) #TODO: Don't hardcode end cell!                               


# forces = Jutul.setup_forces(simulator.model, bc=d_evac)


# numstages = 4
bcs = [d_ads]




# for j = 1:numcycles
#         for i in eachindex(t_stage)
#             ns = numSteps[i]
#             append!(timesteps,ts_matlab[stepStart[i]:stepEnd[i]])
#             append!(sim_forces,repeat([Jutul.setup_forces(simulator.model,bc=bcs[i])],ns))
#         end
# end


for j = 1:numcycles
    for i in eachindex(t_stage)
        numsteps = t_stage[i] / maxdt
        append!(timesteps,repeat([maxdt],Int(floor(numsteps))))
        append!(sim_forces,repeat([Jutul.setup_forces(simulator.model,bc=bcs[i])],Int(floor(numsteps))))
    end
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
    info_level=0,
    forces=sim_forces,
   # max_nonlinear_iterations=0,
   # max_timestep_cuts = 0
)
#display(Mocca.plot_outlet(states))

# display(Mocca.plot_against_matlab_mat(states, 
#     "data/$datapath",
#     cumsum(timesteps)[end], 
#     cumsum(timesteps)))



# ## Plotting
# states_all = []
# push!(states_all,states)

# 

# plot_pvars_spatial(model, states_all)

# outputfile = "mocca_jl_nc200_DCB.mat"
# mocca_jl_nc200_DCB_dt1.mat
outputfile = "mocca_jl_nc200_DCB_dt0_5.mat"

states_jl = Dict()
states_jl["times"] = timesteps

for i in eachindex(states)
    state = states[i]
    state_new = Dict()  
    for (key, value) in state
        state_new[String(key)] = value
    end
    num = lpad(i,3,"0")
    statename = "s_$num"
    states_jl[statename] = state_new
end




MAT.matwrite("data/$outputfile",states_jl)

