import Jutul
import JutulDarcy
import Mocca

filepath = joinpath(@__DIR__, "../json/haghpanah_constants.json")
constants = Mocca.parse_PSA_constants(filepath)
constants.h_in = 0.0
constants.h_out = 0.0

case = Mocca.haghpanah_DCB(constants);
result = Mocca.simulate_case(case);

f = Mocca.plot_outlet(case,result)

