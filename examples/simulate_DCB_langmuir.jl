# # Direct Column Breakthrough simulation with Langmuir.jl backend
#
# This example shows how to run a DCB simulation using Langmuir.jl for isotherm
# calculations with IAST (Ideal Adsorbed Solution Theory) for competitive multicomponent
# adsorption.
#
# This simulation involves injection of a two component flue gas (CO2 and N2)
# into a column of Zeolite 13X initially filled with N2.
#
# Adsorption onto Zeolite 13X is modelled with a dual-site Langmuir adsorption isotherm
# using the Langmuir.jl package for isotherm calculations.
# Injection flow rate is fixed at the right hand side of the column and the left hand side
# is open. There is no heat transfer between the column and the column wall.

# First we load the necessary modules
import Jutul: si_unit
import Mocca
using Printf

#=
Define parameters from Haghpanah et al. 2013.
As we are doing a DCB simulation we set the heat transfer coefficient between
the column and the wall and the wall and the outside to zero.
=#
constants = Mocca.HaghpanahConstants{Float64}(h_in = 0.0, h_out = 0.0);

# We set up a flexible adsorption system with Langmuir.jl backend.
# This system uses Langmuir.jl for isotherm calculations with IAST for competitive adsorption.
# The API is similar to TwoComponentAdsorptionSystem.
system = Mocca.FlexibleAdsorptionSystemWithLangmuir(constants; use_iast = true);

# # Define the model
# Next we need to make the model. This model contains information about
# the domain (grid) over which we will solve the equations and information
# about the system of equations which we are solving.
ncells = 200
model = Mocca.setup_process_model(system; ncells = ncells);

# # Setup the initial state and parameters of the simulation
# Initial values for pressure and temperature of the system
bar = si_unit(:bar);
P_init = 1*bar;
T_init = 298.15;
Tw_init = constants.T_a;

# To avoid numerical errors we set the initial CO2 concentration to be very
# small and not exactly zero
yCO2_init = 1e-10
y_init = [yCO2_init, 1.0 - yCO2_init] # [CO2, N2]

state0 = Mocca.setup_process_state(model;
    Pressure = P_init,
    Temperature = T_init,
    WallTemperature = Tw_init,
    y = y_init
)
parameters = Mocca.setup_process_parameters(model);

# # Setup the timestepping and boundary conditions

# For the DCB we are only running the adsorption stage of a VSA process.
# We will use a total time of 5000 seconds with a single report step
t_ads = 5000.0;
maxdt = 5000.0;

# For the DCB we set up boundary conditions for just an adsorption stage. This sets a fixed
# velocity, concentration and temperature at the inlet, and fixed pressure at
# the outlet. By convention we assume the inlet bc is applied on the left hand
# side and the outlet bc is applied on the right hand side.
sim_forces, timesteps = Mocca.setup_forces(model, [t_ads], ["adsorption"];
    num_cycles = 1, max_dt = maxdt);

# Specify target change of the different state variables for dynamic timestepping
timestep_selector_cfg = (y = 0.01, Temperature = 10.0, Pressure = 10.0);

# # Run the simulation
# Collect everything into a fully specified simulation case and start the simulation
case = Mocca.MoccaCase(model, timesteps, sim_forces; state0 = state0, parameters = parameters)
states, timesteps_out = Mocca.simulate_process(case;
    timestep_selector_cfg = timestep_selector_cfg,
    output_substates = true,
);

# # Comparison with original Haghpanah implementation
# To verify the Langmuir.jl backend produces correct results, we compare
# with the original dual-site Langmuir implementation

println("\nRunning comparison with original implementation...")

# Create original system
system_original = Mocca.TwoComponentAdsorptionSystem(constants)
model_original = Mocca.setup_process_model(system_original; ncells = ncells)

state0_original = Mocca.setup_process_state(model_original;
    Pressure = P_init,
    Temperature = T_init,
    WallTemperature = Tw_init,
    y = y_init
)
parameters_original = Mocca.setup_process_parameters(model_original)

sim_forces_original, _ = Mocca.setup_forces(model_original, [t_ads], ["adsorption"];
    num_cycles = 1, max_dt = maxdt)

case_original = Mocca.MoccaCase(model_original, timesteps, sim_forces_original;
    state0 = state0_original, parameters = parameters_original)

states_original, timesteps_original = Mocca.simulate_process(case_original;
    timestep_selector_cfg = timestep_selector_cfg,
    output_substates = true,
)

# Compare final outlet compositions
outlet_cell = ncells
y_langmuir = states[end][:y][:, outlet_cell]
y_original = states_original[end][:y][:, outlet_cell]

println("\nOutlet composition comparison:")
println("  Component │  Langmuir.jl │  Original  │  Error")
println("  ──────────┼──────────────┼────────────┼─────────")
for i in 1:2
    comp_name = i == 1 ? "CO2" : "N2 "
    error_pct = abs(y_langmuir[i] - y_original[i]) / y_original[i] * 100
    @printf("  %s       │   %.6f   │  %.6f  │  %.3f%%\n", 
            comp_name, y_langmuir[i], y_original[i], error_pct)
end

# # Visualisation and Comparison
# Create comparison plot showing both implementations
println("\nGenerating comparison plot...")

import CairoMakie
CairoMakie.activate!()

# Get outlet data for both implementations
times_langmuir = cumsum(timesteps_out)
times_original = cumsum(timesteps_original)

y_CO2_langmuir = [states[i][:y][1, outlet_cell] for i in 1:length(states)]
y_CO2_original = [states_original[i][:y][1, outlet_cell] for i in 1:length(states_original)]

T_langmuir = [states[i][:Temperature][outlet_cell] for i in 1:length(states)]
T_original = [states_original[i][:Temperature][outlet_cell] for i in 1:length(states_original)]

# Create comparison figure with the two most important variables
fig = CairoMakie.Figure(size=(1200, 400))

# CO2 breakthrough comparison
ax1 = CairoMakie.Axis(fig[1, 1], 
    xlabel="Time [s]", 
    ylabel="CO₂ Mole Fraction [-]",
    title="Outlet CO₂ Breakthrough Comparison")
CairoMakie.lines!(ax1, times_langmuir, y_CO2_langmuir, 
    label="Langmuir.jl (IAST)", linewidth=2.5, color=:blue)
CairoMakie.lines!(ax1, times_original, y_CO2_original, 
    label="Original Haghpanah", linewidth=2, linestyle=:dash, color=:red)
CairoMakie.axislegend(ax1, position=:rb)

# Temperature comparison
ax2 = CairoMakie.Axis(fig[1, 2], 
    xlabel="Time [s]", 
    ylabel="Temperature [K]",
    title="Outlet Temperature Comparison")
CairoMakie.lines!(ax2, times_langmuir, T_langmuir, 
    label="Langmuir.jl (IAST)", linewidth=2.5, color=:blue)
CairoMakie.lines!(ax2, times_original, T_original, 
    label="Original Haghpanah", linewidth=2, linestyle=:dash, color=:red)
CairoMakie.axislegend(ax2, position=:rt)

display(fig)
