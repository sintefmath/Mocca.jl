using JutulDarcy, Jutul
function solve_bl(;nc = 100, time = 1.0, nstep = 100, general_ad = false)
    T = time
    tstep = repeat([T/nstep], nstep)
    G = get_1d_reservoir(nc, general_ad = general_ad)
    nc = number_of_cells(G)
    timesteps = tstep*3600*24 # Convert time-steps from days to seconds

    bar = 1e5
    p0 = 100*bar
    # Define system and realize on grid
    sys = ImmiscibleSystem((LiquidPhase(), VaporPhase()))
    model = SimulationModel(G, sys)
    kr = BrooksCoreyRelPerm(sys, [2.0, 2.0], [0.2, 0.2])
    replace_variables!(model, RelativePermeabilities = kr)
    tot_time = sum(timesteps)
    irate = 500*sum(G.grid.pore_volumes)/tot_time
    src  = [SourceTerm(1, irate, fractional_flow = [1.0, 0.0]), 
            SourceTerm(nc, -irate, fractional_flow = [1.0, 0.0])]
    forces = setup_forces(model, sources = src)

    parameters = setup_parameters(model, PhaseViscosities = [1e-3, 5e-3]) # 1 and 5 cP
    state0 = setup_state(model, Pressure = p0, Saturations = [0.0, 1.0])
    # Simulate and return
    sim = Simulator(model, state0 = state0, parameters = parameters)
    states, report = simulate(sim, timesteps, forces = forces)
    return states, model, report
end
## Perform test
n, n_f = 100, 1000
states, model, report = solve_bl(nc = n)
states_refined, = solve_bl(nc = n_f, nstep = 1000);
## Plot results
using GLMakie
x = range(0, stop = 1, length = n)
x_f = range(0, stop = 1, length = n_f)
f = Figure()
ax = Axis(f[1, 1], ylabel = "Saturation", title = "Buckley-Leverett")
for i in 1:6:length(states)
    lines!(ax, x, states[i][:Saturations][1, :], color = :darkgray)
end
# Plot refined reference
lines!(ax, x_f, states_refined[end][:Saturations][1, :], color = :red)
display(f)
