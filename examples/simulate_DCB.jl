using Mocca
import Jutul
import JutulDarcy
import MAT


ncells = 200

## Intialise Haghpanah parameters
simulator, state0, parameters =
initialize_Haghpanah_model(forcing_term_coefficient=1.0, ncells = ncells)       

## Setup BCs
pars = simulator.model.system.p


# Set timesteps
t_ads = 300

t_stage = [t_ads]
numcycles = 1

timesteps = []
sim_forces = []
maxdt = 1



# d_ads: Optimize v_feed []
d_ads = Mocca.AdsorptionBC(y_feed = pars.y_feed, PH = pars.p_high, v_feed = pars.v_feed,
                                T_feed = pars.T_feed, cell_left = 1, cell_right = ncells) #TODO: Don't hardcode end cell!                               


bcs = [d_ads]

for j = 1:numcycles
    for i in eachindex(t_stage)
        numsteps = t_stage[i] / maxdt
        append!(timesteps,repeat([maxdt],Int(floor(numsteps))))
        append!(sim_forces,repeat([Jutul.setup_forces(simulator.model,bc=bcs[i])],Int(floor(numsteps))))
    end
end


states, report = Jutul.simulate(
    simulator,
    timesteps,
    info_level=0,
    forces=sim_forces,
    # max_nonlinear_iterations=0,
    # max_timestep_cuts = 0
)

using CairoMakie

key_to_label = Dict(
    :y => "y",
    :Pressure => "p",
    :AdsorbedConcentration => "q",
    :Temperature => "T",
    :WallTemperature => "T_{wall}"
)



f = Figure()
ax = Axis(f[1, 1])
t = Float64.(cumsum(timesteps))



for (nsymb, symbol) in enumerate([:y, :Pressure, :AdsorbedConcentration, :Temperature, :WallTemperature])
    @show symbol

    if size(states[end][symbol], 2) == 1
        y = Float64.([result[symbol][end] for result in states])
        ax = Axis(f[nsymb, 1], title=String(symbol), xlabel=CairoMakie.L"t", ylabel=CairoMakie.L"%$(key_to_label[symbol])")
        lines!(ax, t, y, color=:darkgray)
    else
        for i in 1:size(states[end][symbol], 1)
            y = Float64.([val[:y][i, end] for val in states])
            ax = Axis(f[nsymb, i], title=String(symbol), xlabel=CairoMakie.L"t", ylabel=CairoMakie.L"%$(key_to_label[symbol])_%$i")
            lines!(ax, t, y, color=:darkgray)
        end
    end
end


resize!(f.scene, (2 * 400, 3 * 400))
f