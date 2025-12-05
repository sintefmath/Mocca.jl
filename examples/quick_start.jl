import Jutul
import JutulDarcy
import Mocca

filepath = joinpath(@__DIR__, "../json/haghpanah_constants.json")
# This creates a constants struct from the JSON file
input_pars = Mocca.parse_input(filepath)
input_pars.h_in = 0.0
input_pars.h_out = 0.0

case, cfg, sim = Mocca.setup_mocca_case(input_pars)

states, timesteps_out = Mocca.simulate_process(case; simulator = sim, config = cfg)

f = Mocca.plot_outlet(case,states)

Mocca.export_cell_results("testout.csv", case, states; format="csv")