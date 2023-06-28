
module MoccaMaker
using Test
import Jutul
import JutulDarcy
using Mocca
using MAT
using Logging
data_from_mrst = matread("data/VSA_Hag_Simplified.mat")
@info "data is " data_from_mrst
field(name) = data_from_mrst["problem"]["SimulatorSetup"]["state0"][name]
nc = size(field("qCO2"), 1)

initial_temperature = field("T")
# timesteps = tstep*3600*24 # Convert time-steps from days to seconds
sys = AdsorptionFlowSystem(forcing_term_coefficient = 0.0)

bar = 1e5
p0 = 0.4 * bar

yCO2 = field("yCO2")

@show yCO2
# yCO2 = 1e-10
initY = hcat(yCO2, 1.0 .- yCO2)'
@show size(initY)

@show initY
p0 = field("pressure")
@show p0
ctot = p0 / sys.R ./ initial_temperature
c = zeros(size(initY))

c[1, :] = ctot .* initY[1,:]
c[2, :] = ctot .* initY[2,:]
equilinit = [compute_equilibrium(sys, c[:, i], initial_temperature[i]) for i in 1:nc]

@info "from mocca" equilinit
@info "from mrst" field("qstar")

@test first.(equilinit) ≈ field("qstar")[1]
@test [eq[2] for  eq in equilinit] ≈ field("qstar")[2]

k = [compute_ki(sys, c[:, i], equilinit[i]) for i in 1:nc]

@show k

second(arr) = arr[2]
@test all(isapprox.(first.(k), 0.0036,  atol=0.0001))
@test all(isapprox.(second.(k), 2.9573,  atol=0.0001))
end
