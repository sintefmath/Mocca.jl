#=
# Cyclic VSA simulation with explicit parameter setup

This example sets up a four-stage vacuum swing adsorption (VSA) process
by constructing each component explicitly.  It shows how to compose an
isotherm, mass transfer model, system, domain, and boundary conditions
from individual physical parameters. This should act as a minimal example
for users who want to set up a custom simulation.

The column separates a CO₂/N₂ mixture on Zeolite 13X.
=#

import Jutul
import Jutul: si_unit
import Mocca
using StaticArrays

# # 1. Define the physics

# Dual-site Langmuir isotherm for CO₂/N₂ on Zeolite 13X
isotherm = Mocca.DualSiteLangmuir(
    qsb = SVector(3489.44, 6613.551),       # saturation loading, site b  [mol/m³]
    b0  = SVector(8.65e-7, 2.5e-6),         # pre-exponential factor      [m³/mol]
    ΔUb = SVector(-36641.21, -15800.0),     # adsorption energy, site b   [J/mol]
    qsd = SVector(2872.35, 0.0),            # saturation loading, site d  [mol/m³]
    d0  = SVector(2.63e-8, 0.0),            # pre-exponential factor      [m³/mol]
    ΔUd = SVector(-35690.66, 0.0),          # adsorption energy, site d   [J/mol]
)

# Particle and transport properties
D_m = 1.6e-5    # molecular diffusivity [m²/s]
d_p = 2e-3      # particle diameter     [m]

mass_transfer = Mocca.LinearDrivingForce(
    D_m = D_m,
    τ   = 3.0,       # tortuosity        [-]
    ϵ_p = 0.35,      # particle porosity [-]
    d_p = d_p,
)

# Assemble the adsorption system
system = Mocca.TwoComponentAdsorptionSystem(
    isotherm            = isotherm,
    mass_transfer       = mass_transfer,
    molecular_masses    = SVector(44.01e-3, 28.0e-3),     # [kg/mol]
    heat_capacity_gas      = SVector(697.5687, 1096.4),   # C_pg [J/(kg·K)]
    heat_capacity_adsorbed = SVector(697.5687, 1096.4),   # C_pa [J/(kg·K)]
)

# # 2. Build the column domain

ncells   = 200
L        = 1.0       # column length     [m]
r_in     = 0.1445    # inner radius      [m]
r_out    = 0.162     # outer wall radius [m]
porosity = 0.37      # bed porosity      [-]

permeability = Mocca.compute_permeability(porosity, d_p)
dispersion   = Mocca.compute_dispersion(D_m, 1.0, d_p)

mesh = Mocca.column_mesh(ncells, L, r_in)

domain = Mocca.mocca_domain(mesh;
    porosity             = porosity,
    permeability         = permeability,
    r_in                 = r_in,
    r_out                = r_out,
    diffusion_coefficient = dispersion,
    thermal_conductivity = 0.0903,     # K_z  [W/(m·K)]
    adsorbent_density    = 1130.0,     # ρ_s  [kg/m³]
    adsorbent_heat_capacity = 1070.0,  # C_ps [J/(kg·K)]
    wall_density         = 7800.0,     # ρ_w  [kg/m³]
    wall_heat_capacity   = 502.0,      # C_pw [J/(kg·K)]
    wall_conductivity    = 16.0,       # K_w  [W/(m·K)]
    fluid_viscosity      = 1.72e-5,    # μ    [Pa·s]
    fluid_density        = 1.22638310956,  # ρ_g  [kg/m³]
    inner_htc            = 8.6,        # h_in [W/(m²·K)]
    outer_htc            = 2.5,        # h_out[W/(m²·K)]
    ambient_temperature  = 298.15,     # T_a  [K]
)

# # 3. Create the model, initial state and parameters

model = Mocca.setup_process_model(system, domain)

bar = si_unit(:bar)
state0 = Mocca.setup_process_state(model;
    Pressure        = 1 * bar,
    Temperature     = 298.15,
    WallTemperature = 298.15,
    y = [1e-10, 1.0 - 1e-10],
)
parameters = Mocca.setup_process_parameters(model)

# # 4. Boundary conditions for each stage

y_feed  = SVector(0.15, 0.85)
T_feed  = 298.15
p_high  = 1e5
p_low   = 0.1e5
p_inter = 0.2e5
λ       = 0.5

stage_times = [15.0, 15.0, 30.0, 40.0]

bcs = [
    Mocca.PressurisationBC(y_feed = y_feed, PH = p_high, PL = p_low, λ = λ, T_feed = T_feed),
    Mocca.AdsorptionBC(y_feed = y_feed, PH = p_high, v_feed = 0.37, T_feed = T_feed),
    Mocca.BlowdownBC(PH = p_high, PI = p_inter, λ = λ),
    Mocca.EvacuationBC(PL = p_low, PI = p_inter, λ = λ),
]

# # 5. Assemble forces and simulate

sim_forces, timesteps = Mocca.setup_forces(model, stage_times, bcs;
    num_cycles = 3, max_dt = 1.0)

case = Mocca.MoccaCase(model, timesteps, sim_forces;
    state0 = state0, parameters = parameters)

states, timesteps_out = Mocca.simulate_process(case;
    output_substates = true,
    info_level = 0,
)

# # 6. Visualise

f_outlet = Mocca.plot_cell(states, model, timesteps_out, ncells)
f_column = Mocca.plot_state(states[end], model)
