using Mocca
#=
# Mocca Quick Start Example
This example shows how to set up and run a cyclic vacuum swing adsorption
simulation using a JSON input file to define the simulation parameters.

The input we use is the same as in the [DCB simulation](simulate_DCB.md)
i.e. the direct column breakthrough experiment from [Haghpanah et al. 2013](http://dx.doi.org/10.1021/ie302658y)
=#

## Load input parameters from simple JSON input file

# Mocca.parse_input can read input parameters from a JSON file or from a Julia dictionary and output a 
# Mocca.adsorptionInput type structure.

# Here we read from a JSON file with a simplified format [haghpanah_DCB_input_simple.json](../models/json/haghpanah_DCB_input_simple.json)
filepath = joinpath(@__DIR__, "../models/json/haghpanah_DCB_input_simple.json")
input_pars = Mocca.parse_input(filepath)

## Load input parameters from detailed JSON  file 
# We could also read from the detailed JSON file using the same function  
# See [haghpanah_DCB_input.json](../models/json/haghpanah_DCB_input.json)
# Uncomment the lines below to use this file instead.   

# filepath = joinpath(@__DIR__, "../models/json/haghpanah_DCB_input.json")
# input_pars = Mocca.parse_input(filepath)

# Load input parameters from Julia dictionary (optional)

# Alternatively we could read from a Julia dictionary directly. 
# See [models/julia/model_library.jl](../models/julia/model_library.jl) for an example dictionary.
# Uncomment the lines below to use this method.

# input_pars = Mocca.parse_input(haghpanah_DCB_input())

## Setup simulation case

# Once loaded we setup a Mocca simulation case conatining the model, initial state,
# parameters, boundary conditions and timestepping configuration.
case, ts_config = Mocca.setup_mocca_case(input_pars)

## Run simulation

states, timesteps = Mocca.simulate_process(case; timestep_selector_cfg = ts_config,
    output_substates = true, info_level = 0);

## Export output to CSV and plot results

Mocca.export_cell_results(joinpath(Mocca.moccaResultsDir, "haghpanah_DCB_results.csv"), case, states, timesteps; format="csv")

f = Mocca.plot_outlet(case, states, timesteps);

display(f)