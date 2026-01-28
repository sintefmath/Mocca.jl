```@meta
CurrentModule = Mocca
```

# Quick start example
Run the following code to quickly setup and run a
Direct Column breakthrough adsorption simulation using predefined input parameters.

The example uses some utility functions which simplify the simulation setup.
To see the steps used in more detail, please refer to the
[Direct Column Breakthrough simulation](simulate_DCB.md) example.

```julia
using Mocca

# Import and load input parameters
json_dir = joinpath(dirname(pathof(Mocca)), "../models/json/")
filepath = joinpath(json_dir, "haghpanah_DCB_input_simple.json")
(constants, info ) = Mocca.parse_input(filepath)

# Setup and run simulation
case, ts_config = Mocca.setup_mocca_case(constants, info)
states, timesteps = Mocca.simulate_process(case; timestep_selector_cfg = ts_config,
    output_substates = true, info_level = 0)

# Save results to CSV and plot
Mocca.export_cell_results(joinpath(Mocca.moccaResultsDir, "haghpanah_DCB_results.csv"),
    case, states, timesteps; format="csv")
f = Mocca.plot_outlet(case, states, timesteps)
display(f)
```
