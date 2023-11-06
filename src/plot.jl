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
    t = Float64.(cumsum(timesteps))

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
            xlabel=L"t\; [s]",
            ylabel=L"%$(units[symbol])")

        if size(states[end][symbol], 2) == 1
           lines!(ax, t, [result[symbol][cell] for result in states])
        else
            for k in 1:size(states[end][symbol], 1)
                lines!(ax, t, [result[symbol][k, cell] for result in states], label = comp_names[k])
            end
            leg = Legend(f[2,3], ax, tellwidth=false)
            
        end

        
    end

    return f
end


function plot_state(state, model)

    pvars = model.primary_variables.keys
    comp_names = model.system.component_names

    units = units_dict()
    x = model.data_domain[:cell_centroids][1,:]

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
            xlabel=L"x\; [m]",
            ylabel=L"%$(units[symbol])")

        if size(state[symbol], 2) == 1
           lines!(ax, x, state[symbol])
        else
            for k in 1:size(state[symbol], 1)
                lines!(ax, x, state[symbol][k,:], label = comp_names[k])
            end
            leg = Legend(f[2,3], ax, tellwidth=false)
            
        end

        
    end

    return f
end


