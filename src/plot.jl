using CairoMakie

function units_dict()
    unit_label = Dict(
        :y => "[-]",
        :Pressure => "[Pa]",
        :AdsorbedConcentration => "[mol/m^{3}]",
        :Temperature => "[K]",
        :WallTemperature => "[K]"
    )

    return unit_label
end


function plot_cell(states, model, timesteps, cell)

    pvars = model.primary_variables.keys
    comp_names = model.system.component_names


    units = units_dict()
    t = cumsum(timesteps)

    f = Figure(resolution = (900, 600))
    ga = f[2,3] = GridLayout()
    r = 1

    for (i, symbol) in enumerate(pvars)

        c = i
        if i > 3
            c = i - 3
            r = 2
        end

        ax = Axis(f[r,c],
            title=String(symbol),
            xlabel=L"t",
            ylabel=L"%$(units[symbol])")

        if size(states[end][symbol], 2) == 1
           lines!(ax, t, Float64.([result[symbol][cell] for result in states]))
        else
            for k in 1:size(states[end][symbol], 1)
                lines!(ax, t, Float64.([result[symbol][k, cell] for result in states]), label = comp_names[k])
            end
            leg = Legend(f[2,3], ax, tellwidth=false)
            
        end

        
    end

    return f
end

