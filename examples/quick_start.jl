import Jutul
import JutulDarcy
import Mocca

filepath = joinpath(@__DIR__, "../json/haghpanah_constants.json")
# This creates a constants struct from the JSON file
input_pars = Mocca.parse_input(filepath)
input_pars.h_in = 0.0
input_pars.h_out = 0.0

case, ts_config = Mocca.setup_mocca_case(input_pars)

states, timesteps = Mocca.simulate_process(case; timestep_selector_cfg = ts_config,
    output_substates = true, info_level = 0);


Mocca.export_cell_results("testout_2.csv", case, states, timesteps; format="csv")

f = Mocca.plot_outlet(case, states, timesteps);

display(f)