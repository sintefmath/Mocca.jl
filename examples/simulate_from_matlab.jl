using Mocca
import Jutul
import JutulDarcy 
simulator =
    initialize_from_matlab("data/VSA_Hag_Simplified.mat", 
        forcing_term_coefficient = 1.0)

g = Jutul.physical_representation(simulator.model)
model = simulator.model
d = Mocca.PressurationBC(trans = Mocca.compute_permeability(model.system)/Mocca.compute_dx(model, 1))
forces = Jutul.setup_forces(simulator.model, bc=d)

numberoftimesteps = 1_000
dt = 1.0 / numberoftimesteps
timesteps = repeat([dt], numberoftimesteps)

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
    info_level = 0,
    forces = forces,
    max_nonlinear_iterations = 10000,
)#000)

display(Mocca.plot_states(states))
