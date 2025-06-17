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

function prettyVarNames(vars::Vector{Symbol})
    pretty_names = Dict()

    for (i, symb) in enumerate(vars)
        var = String(symb)
        words = split(var, r"(?=[A-Z])")
        pretty_names[symb] = join(words, " ")
    end

    return pretty_names
end


function plot_cell(states, model, timesteps, cell)

    pvars = model.primary_variables.keys
    comp_names = model.system.component_names

    pretty_names = prettyVarNames(pvars)
    units = units_dict()
    t = Float64.(cumsum(timesteps))

    f = Figure(size = (900, 600))
    ga = f[2,3] = GridLayout()
    r = 1

    for (i, symb) in enumerate(pvars)
        c = i
        if i > 3
            c = i - 3
            r = 2
        end
        ax = Axis(f[r,c],
            title=pretty_names[symb],
            xlabel=L"t\; [s]",
            ylabel=L"%$(units[symb])")

        if size(states[end][symb], 2) == 1
           lines!(ax, t, [result[symb][cell] for result in states])
        else
            for k in 1:size(states[end][symb], 1)
                lines!(ax, t, [result[symb][k, cell] for result in states], label = comp_names[k])
            end
            leg = Legend(f[2,3], ax, tellwidth=false)
            
        end

        
    end

    return f
end


function plot_state(state, model)

    pvars = model.primary_variables.keys
    comp_names = model.system.component_names
    pretty_names = prettyVarNames(pvars)

    units = units_dict()
    x = model.data_domain[:cell_centroids][1,:]

    f = Figure(size = (900, 600))
    ga = f[2,3] = GridLayout()
    r = 1

    for (i, symb) in enumerate(pvars)
        c = i
        if i > 3
            c = i - 3
            r = 2
        end
        ax = Axis(f[r,c],
            title = pretty_names[symb],
            xlabel=L"x\; [m]",
            ylabel=L"%$(units[symb])")

        if size(state[symb], 2) == 1
           lines!(ax, x, state[symb])
        else
            for k in 1:size(state[symb], 1)
                lines!(ax, x, state[symb][k,:], label = comp_names[k])
            end
            leg = Legend(f[2,3], ax, tellwidth=false)
            
        end

        
    end

    return f
end


function plot_optimization_history(dict_parameters::Jutul.DictParameters)
    vals = dict_parameters.history.val

    f = Figure(size= (900, 600))
    ax = Axis(f[1,1]; xlabel = "Iteration #", ylabel = "Objective error")
    ax.xticks = 1:length(vals)
    lim_min, lim_max = extrema(vals)
    lim_min *= 0.5
    lim_max *= 2.0
    ylims!(ax, lim_min, lim_max)
    ax.yscale = log10

    scatter!(ax, vals)
    return f
end