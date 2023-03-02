using Mocca
import Jutul
simulator = initialize_from_matlab("data/VSA_Hag_Simplified.mat", forcing_term_coefficient=0.0)
forces = Jutul.setup_forces(simulator.model)

numberoftimesteps = 100
dt = 1.0 / numberoftimesteps
timesteps = repeat([dt], numberoftimesteps)

states, report = Jutul.simulate(simulator, timesteps, info_level=3, forces=forces, max_nonlinear_iterations=10000)#000)

display(plot_states(states))