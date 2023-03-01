module MoccaMaker
using Mocca
import Jutul
import JutulDarcy
using StaticArrays
using Parameters
using Tullio
using CairoMakie
using MakiePublication
using LinearAlgebra

time = 1.0
nc = 10
# nc = 2
nstep = 80
# nstep = 100000
general_ad = true
T = time

initial_temperature = 298
tstep = repeat([T / nstep], nstep)
# timesteps = tstep*3600*24 # Convert time-steps from days to seconds
timesteps = tstep
sys = AdsorptionFlowSystem(forcing_term_coefficient = 1.0)
ϵ = sys.Φ
r_in = sys.d_p / 2.0 #0.289 / 2.0
perm = 4 / 150 * ((ϵ / (1 - ϵ))^2) * r_in^2
# perm = 1e-14
# perm = 5.0625e-09

@info "Using perm " perm

G = JutulDarcy.get_1d_reservoir(nc, general_ad=general_ad, poro=ϵ, perm=perm)

ctx = Jutul.DefaultContext()
model = Jutul.SimulationModel(G, sys, context=ctx)

g = Jutul.physical_representation(model.domain)
# pv = vol * poro
pv = g.pore_volumes
volumes = pv / ϵ
solid_volume = volumes * (1 - ϵ)

# axial Dispersion
dp = sys.d_p
Dm = sys.D_m
V0_inter = 0.03653         # Interstitial inlet velocity [m/s]
V0 = V0_inter * ϵ         # Inlet velocity [m/s]

DL = 0.7 * Dm + 0.5 * V0 * dp
bar = 1e5
p0 = 1 * bar
parameters = Jutul.setup_parameters(model, Temperature=initial_temperature,
    solidVolume=solid_volume,
    axialDispersion=DL,
    fluidViscosity=1.72e-5,
    PhaseMassDensities=1.0)

# TODO: Find a nicer way to specify trans on the boundary
d = JutulDarcy.FlowBoundaryCondition(nc, 2 * p0, trans_flow=g.trans[1], fractional_flow=(0.5, 0.5))
# forces = Jutul.setup_forces(model, sources = [], bc = d)
forces = Jutul.setup_forces(model)
irate = 500 * sum(g.pore_volumes) / time
# src  = [JutulDarcy.SourceTerm(1, irate, fractional_flow = [1.0, 0.0]), 
#     JutulDarcy.SourceTerm(nc, -irate, fractional_flow = [1.0, 0.0])]
# forces = JutulDarcy.setup_forces(model, sources = src)
@info "parameter set" parameters g.trans

yCO2 = 1e-16
# yCO2 = 1e-10
initY = [yCO2, 1 - yCO2]
ctot = p0 / sys.R / initial_temperature
c = ctot .* initY
equilinit = compute_equilibrium(sys, c, initial_temperature) # TODO: Should this still be zero for CO2?
# equilinit = [0, equilinit[2]]
# equilinit = [0, 0.0]
@show equilinit


p_init = p0
p_init = collect(LinRange(p0, 0.899999 * p0, nc))
# p_init = [100*bar, 99.9999*bar]
state0 = Jutul.setup_state(model,
    Pressure=p_init,
    y=initY,
    adsorptionRates=equilinit)
# Simulate and return
sim = Jutul.Simulator(model, state0=state0, parameters=parameters)
states, report = Jutul.simulate(sim, timesteps, info_level=3, forces=forces, max_nonlinear_iterations=10)#000)
# states, report = Jutul.simulate(sim, timesteps, info_level = 3, forces=forces, max_timestep_cuts = 0)
# states, report = Jutul.simulate(sim, timesteps, info_level = 4, forces=forces, max_timestep_cuts = 0, max_nonlinear_iterations = 0)

# with_theme(theme_web()) do 
#     f = CairoMakie.Figure()
#     ax = CairoMakie.Axis(f[1, 1], ylabel = "y", title = "Adsorption")
#     x = range(0, stop = 1, length = nc)
#     for i in 1:6:length(states)
#         CairoMakie.lines!(ax, x, states[i][:y][2, :], color = :darkgray)
#     end
#     display(f)
# end

with_theme(theme_web()) do
    f = CairoMakie.Figure()
    x = collect(LinRange(0.0, 1.0, nc))
    key_to_label = Dict(
        :y => "y",
        :Pressure => "p",
        :adsorptionRates => "q"
    )
    for (nsymb, symbol) in enumerate([:y, :Pressure, :adsorptionRates])
        @show symbol
        if size(states[end][symbol], 2) == 1
            ax = CairoMakie.Axis(f[nsymb, 1], title=String(symbol), xlabel=L"x", ylabel=L"%$(key_to_label[symbol])")
            for j in 1:length(states)
                CairoMakie.lines!(ax, x, states[j][symbol][:], color=:darkgray)
            end
            CairoMakie.lines!(ax, x, states[end][symbol][:], color=:red)

        else
            for i in 1:size(states[end][symbol], 1)
                ax = CairoMakie.Axis(f[nsymb, i], title=String(symbol), xlabel=L"x", ylabel=L"%$(key_to_label[symbol])_%$i")
                for j in 1:length(states)
                    CairoMakie.lines!(ax, x, states[j][symbol][i, :], color=:darkgray)
                end
                CairoMakie.lines!(ax, x, states[end][symbol][i, :], color=:red)

            end
        end
    end
    CairoMakie.resize!(f.scene, (2*400, 3*400))
    display(f)
end
end
MoccaMaker.model