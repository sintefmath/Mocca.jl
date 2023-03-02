import CairoMakie
import MakiePublication
function plot_states(states)
    return MakiePublication.with_theme(MakiePublication.theme_web()) do
        f = CairoMakie.Figure()
        nc = size(states[end][:Pressure], 1)
        x = collect(LinRange(0.0, 1.0, nc))
        key_to_label = Dict(
            :y => "y",
            :Pressure => "p",
            :adsorptionRates => "q"
        )

        for (nsymb, symbol) in enumerate([:y, :Pressure, :adsorptionRates])
            @show symbol
            # Truncating to float16 seems to be needed due to some weird cairomakie bug:
            # https://discourse.julialang.org/t/range-step-cannot-be-zero/66948/10
            # TODO: Fix the above
            if size(states[end][symbol], 2) == 1
                ax = CairoMakie.Axis(f[nsymb, 1], title=String(symbol), xlabel=CairoMakie.L"x", ylabel=CairoMakie.L"%$(key_to_label[symbol])")
                for j in eachindex(states)
                    CairoMakie.lines!(ax, x, Float16.(states[j][symbol][:]), color=:darkgray)
                end
                CairoMakie.lines!(ax, x, Float16.(states[end][symbol][:]), color=:red)

            else
                for i in 1:size(states[end][symbol], 1)
                    ax = CairoMakie.Axis(f[nsymb, i], title=String(symbol), xlabel=CairoMakie.L"x", ylabel=CairoMakie.L"%$(key_to_label[symbol])_%$i")
                    for j in eachindex(states)
                        CairoMakie.lines!(ax, x, Float16.(states[j][symbol][i, :]), color=:darkgray)
                    end

                    CairoMakie.lines!(ax, x, Float16.(states[end][symbol][i, :]), color=:red)

                end
            end
        end
        CairoMakie.resize!(f.scene, (2*400, 3*400))
        return f
    end
end