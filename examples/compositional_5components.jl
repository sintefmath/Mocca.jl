
using MultiComponentFlash
n2_ch4 = MolecularProperty(0.0161594, 4.58e6, 189.515, 9.9701e-05, 0.00854)
co2 = MolecularProperty(0.04401, 7.3866e6, 304.200, 9.2634e-05, 0.228)
c2_5 = MolecularProperty(0.0455725, 4.0955e6, 387.607, 2.1708e-04, 0.16733)
c6_13 = MolecularProperty(0.117740, 3.345e6, 597.497, 3.8116e-04, 0.38609)
c14_24 = MolecularProperty(0.248827, 1.768e6, 698.515, 7.2141e-04, 0.80784)

bic = [0.11883 0.00070981 0.00077754 0.01 0.011;
       0.00070981 0.15 0.15 0.15 0.15;
       0.00077754 0.15 0 0 0;
       0.01 0.15 0 0 0;
       0.011 0.15 0 0 0]

mixture = MultiComponentMixture([n2_ch4, co2, c2_5, c6_13, c14_24], A_ij = bic, names = ["N2-CH4", "CO2", "C2-5", "C6-13", "C14-24"])
eos = GenericCubicEOS(mixture, PengRobinson())

using Jutul, JutulDarcy, JutulViz
nx = ny =20
nz = 2
bar = 1e5
dims = (nx, ny, nz)
g = CartesianMesh(dims, (1000.0, 1000.0, 1.0))
nc = number_of_cells(g)
#Darcy = 9.869232667160130e-13
K = repeat([5.0e-14], 1, nc)
res = discretized_domain_tpfv_flow(tpfv_geometry(g), porosity = 0.25, permeability = K)
## Set up a vertical well in the first corner, perforated in all layers
prod = setup_vertical_well(g, K, nx, ny, name = :Producer)
## Set up an injector in the opposite corner, perforated in all layers
inj = setup_vertical_well(g, K, 1, 1, name = :Injector)

## Plot the permeability (scaled to Darcy) and the wells
fig, ax, p = plot_cell_data(g, K[:])
plot_well!(ax, g, inj, textscale = 0.1)
plot_well!(ax, g, prod, color = :darkblue, textscale = 0.1)


rhoLS, rhoVS = 1000.0, 100.0
rhoS = [rhoLS, rhoVS]
L, V = LiquidPhase(), VaporPhase()
# Define system and realize on grid
sys = MultiPhaseCompositionalSystemLV(eos, (L, V))
model, parameters = setup_reservoir_model(res, sys, wells = [inj, prod], reference_densities = rhoS, block_backend = true);
kr = BrooksCoreyRelPerm(sys, 2.0, 0.0, 1.0)
model = replace_variables!(model, RelativePermeabilities = kr)
T0 = repeat([387.45], 1, nc)
parameters[:Reservoir][:Temperature] = T0
state0 = setup_reservoir_state(model, Pressure = 225*bar, OverallMoleFractions = [0.463, 0.01640, 0.20520, 0.19108, 0.12432]);


day = 24*3600.0
dt = repeat([2.0]*day, 365)
rate_target = TotalRateTarget(0.0015)
I_ctrl = InjectorControl(rate_target, [0, 1, 0, 0, 0], density = rhoVS)
bhp_target = BottomHolePressureTarget(100*bar)
P_ctrl = ProducerControl(bhp_target)

controls = Dict()
controls[:Injector] = I_ctrl
controls[:Producer] = P_ctrl
forces = setup_reservoir_forces(model, control = controls)

sim, config = setup_reservoir_simulator(model, state0, parameters, info_level = -1);
states, reports = simulate!(sim, dt, forces = forces, config = config);


## Once the simulation is done, we can plot the states
f, = plot_interactive(g, map(x -> x[:Reservoir], states))
display(f)
## Plot the wells
# wd = full_well_outputs(sim.model, states, forces)
# time = report_times(reports)
# plot_well_results(wd, time)
