using Mocca

# Quick start example
# We import and run a DCB adsorption simulation from a JSON file.

# Setup filepath to JSON input
json_dir = joinpath(dirname(pathof(Mocca)), "../models/json/")

# Load input parameters from JSON 
filepath = joinpath(json_dir, "haghpanah_DCB_input_simple.json")
input_pars = Mocca.parse_input(filepath)

# Option 2: Load input from the detailed JSON format
# filepath = joinpath(json_dir, "haghpanah_DCB_input.json")
# input_pars = Mocca.parse_input(filepath)

# Option 3: Get input pars from a predefined Julia dictionary
(constants, info )= Mocca.parse_input(haghpanah_DCB_input())
    
# Setup simulation case and timestep configuration
case, ts_config = Mocca.setup_mocca_case(constants, info)

# Run simulation
states, timesteps = Mocca.simulate_process(case; timestep_selector_cfg = ts_config,
    output_substates = true, info_level = 0)

# Export results to CSV
Mocca.export_cell_results(joinpath(Mocca.moccaResultsDir, "haghpanah_DCB_results.csv"),
    case, states, timesteps; format="csv")

# Plot results at the outlet
f = Mocca.plot_outlet(case, states, timesteps)
display(f)
