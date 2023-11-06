using CairoMakie

function units_dict()
    unit_label = Dict(
        :y => "[-]",
        :Pressure => "[Pa]",
        :AdsorbedConcentration => "[mol m^{-3}]",
        :Temperature => "[K]",
        :WallTemperature => "[K]"
    )

    return unit_label
end


function plot_cell(states, model, timesteps, cell)

    pvars = model.primary_variables.keys


    units = units_dict()
    t = cumsum(timesteps)

    j = range(;  stop = ceil(size(pvars)[1]/2))

    f = Figure()
    for (i, symbol) in enumerate(pvars)
        ax = Axis(f[i, 1],
            title=String(symbol),
            xlabel=L"t",
            ylabel=L"%$(units[symbol])")

        if size(states[end][symbol], 2) == 1
           lines!(ax, t, Float64.([result[symbol][cell] for result in states]))
        else
            for k in 1:size(states[end][symbol], 1)
                lines!(ax, t, Float64.([result[symbol][k, cell] for result in states]))
            end
        end
    end


    return f
end

