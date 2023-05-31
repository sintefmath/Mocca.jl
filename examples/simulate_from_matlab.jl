using Mocca
import Jutul
import JutulDarcy
import MAT


simulator =
    initialize_from_matlab("data/only_pressurisation.mat",
        forcing_term_coefficient=1.0)

g = Jutul.physical_representation(simulator.model)
model = simulator.model

# TODO: Is the transmisibility for the bc computed correctly here? Seems to match...
# d = Mocca.PressurationBC(trans= 2 * Mocca.compute_permeability(model.system) / Mocca.compute_dx(model, 1) * (pi * model.system.p.r_in^2))
d = Mocca.AdsorptionBC(trans= 2 * Mocca.compute_permeability(model.system) / Mocca.compute_dx(model, 1) * (pi * model.system.p.r_in^2))

forces = Jutul.setup_forces(simulator.model, bc=d)

numberoftimesteps = 15000
dt = 15.0 / numberoftimesteps
times_matlab = collect(Iterators.flatten(MAT.matread("data/VSA_Comparison_HAG_n30_nc1_julia_comp.mat")["results"]["time"]))
times_matlab_zero = zeros(length(times_matlab) + 1)
times_matlab_zero[2:end] = times_matlab

timesteps = times_matlab - times_matlab_zero[1:end-1]
@show timesteps

# timesteps = repeat([dt], numberoftimesteps)[1:1000]


nc = size(simulator.storage.primary_variables.Pressure, 1)
@show nc
# d = JutulDarcy.FlowBoundaryCondition(
#     nc,
#     2 * simulator.storage.primary_variables.Pressure[1].value, # TODO: Do this nicer
#     trans_flow = g.trans[1],
#     fractional_flow = (0.5, 0.5),
# )
# forces = Jutul.setup_forces(simulator.model, sources = [], bc = d)
states, report = Jutul.simulate(
    simulator,
    timesteps,
    info_level=0,
    forces=forces,
    max_nonlinear_iterations=20,
)#000)
##
#display(Mocca.plot_states(states))
#display(Mocca.plot_outlet(cumsum(timesteps), states))
display(Mocca.plot_against_matlab_mat(states, 
    "data/VSA_Comparison_HAG_n30_nc1_julia_comp.mat", 
    cumsum(timesteps)[end], 
    cumsum(timesteps)))
##
# plot_states(states)