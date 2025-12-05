# Plan for Mocca IO

## File structure

Mocca
- json
-- haghpanah_input.json
-- mocca_adsorption_schema.json

- input_output
-- parse_input : Reads JSON or Dict to make adsorptionInput Struct 
--- (Later split into separate structs in an inputDict => Subfunctions: parseConstants, parseSimulation, parseBoundaryConditions, parseSolver, setupInputDict)

-- setup_mocca_case : Takes adsorptionInput and makes a case
-- export_results : exports simulation output to csv (for now..)
-- plot_outlet : plots values at outlet


## Input file
haghpanah_input.json

with schema:
mocca_adsorption_schema.json

## quick_start.jl

fpath =  'input.json'

adsorptionInput = parse_input(fpath)

case, config = setup_mocca_case(adsorptionInput)

results, timesteps = Mocca.simulate_process(case; config = config, output_substates = true)

f = Mocca.plot_outlet(case,states)

Mocca.export_cell_results("testout.csv", case, result; format="csv")


## TODO 4/12/25

Change parse_PSA_constants=> parse_input so that it parses a dict and creates an adsorption case