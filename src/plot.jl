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

function plot_outlet(case, states, timesteps_out)

    # Convenience function to plot outlet cell variables over time

    outlet_cell = size(states[1][:y][1,:],1)
    f_outlet = Mocca.plot_cell(states, case.model, timesteps_out, outlet_cell);

    return f_outlet
end


function plot_optimization_history(dict_parameters::Jutul.DictParameters;
    yscale = log10,
    ylabel = "Objective error"
)
    vals = dict_parameters.history.objectives

    f = Figure(size= (900, 600))
    ax = Axis(f[1,1]; xlabel = "Iteration #", ylabel = ylabel)
    ax.xticks = 1:length(vals)
    lim_min, lim_max = extrema(vals)
    lim_min *= 0.5
    lim_max *= 2.0
    if yscale == identity
        lim_min = 0.0
    end
    ylims!(ax, lim_min, lim_max)
    ax.yscale = yscale

    scatter!(ax, vals)
    return f
end

"""
    plot_cell_comparison(runs, model, cell; component = 1)

Overlay the history of `cell` in several simulations of the same unit, one
panel per primary variable. `runs` is a vector of `label => (states,
timesteps)`, where `states` belong to that unit alone (for a flowsheet, pass
`[s[:Bed] for s in states]`). For variables with one value per component only
`component` is shown.
"""
function plot_cell_comparison(runs, model, cell; component = 1)
    curves = map(runs) do (label, (states, timesteps))
        t = Float64.(cumsum(timesteps))
        label => (t, symb -> [_component_value(s[symb], component, cell) for s in states])
    end
    return _comparison_figure(model, curves, L"t\; [s]", component)
end

"""
    plot_state_comparison(states, model; component = 1)

Overlay profiles along the column from several states of the same unit.
`states` is a vector of `label => state`.
"""
function plot_state_comparison(states, model; component = 1)
    x = model.data_domain[:cell_centroids][1, :]
    curves = map(states) do (label, state)
        label => (x, symb -> [_component_value(state[symb], component, c) for c in eachindex(x)])
    end
    return _comparison_figure(model, curves, L"x\; [m]", component)
end

_component_value(v::AbstractVector, component, cell) = v[cell]
_component_value(v::AbstractMatrix, component, cell) = v[component, cell]

function _comparison_figure(model, curves, xlabel, component)
    pvars = model.primary_variables.keys
    comp_name = model.system.component_names[component]
    pretty_names = prettyVarNames(pvars)
    units = units_dict()
    colors = Makie.wong_colors()
    styles = [:solid, :dash, :dot, :dashdot]
    f = Figure(size = (900, 600))
    ax = nothing
    for (i, symb) in enumerate(pvars)
        r, c = fldmod1(i, 3)
        title = pretty_names[symb]
        if model.primary_variables[symb] isa Jutul.VectorVariables
            title = "$title ($comp_name)"
        end
        ax = Axis(f[r, c], title = title, xlabel = xlabel, ylabel = L"%$(units[symb])")
        for (j, (label, (x, values))) in enumerate(curves)
            lines!(ax, x, values(symb),
                color = colors[mod1(j, length(colors))],
                linestyle = styles[mod1(j, length(styles))],
                label = label)
        end
    end
    Legend(f[2, 3], ax, tellwidth = false)
    return f
end

"""
    plot_flowsheet_streams(states, model, timesteps; stages = nothing, component = 1)

Molar flow through each flow device of a flowsheet (positive from inlet to
outlet), and the mole fraction of `component` in each stream while it flows.
With `stages` (as passed to `setup_schedule`), the stages of each cycle are
shaded.
"""
function plot_flowsheet_streams(states, model::Jutul.MultiModel, timesteps;
        stages = nothing,
        component = 1,
        devices = [k for (k, m) in pairs(model.models) if m isa FlowDeviceModel],
        flow_tol = 1e-6
    )
    t = Float64.(cumsum(timesteps))
    comp_name = model.models[first(devices)].system.component_names[component]
    colors = Makie.wong_colors()
    f = Figure(size = (900, 600))
    ax_F = Axis(f[1, 1], title = "Molar flow through devices", ylabel = L"[mol/s]")
    ax_y = Axis(f[2, 1], title = "$comp_name mole fraction in flowing streams", xlabel = L"t\; [s]", ylabel = L"[-]")
    linkxaxes!(ax_F, ax_y)
    if !isnothing(stages)
        for ax in (ax_F, ax_y)
            _shade_stages!(ax, stages, t[end], colors)
        end
        Legend(f[2, 2], [PolyElement(color = (colors[i], 0.25)) for i in eachindex(stages)],
            [st.name for st in stages], "Stage", tellheight = false)
    end
    F_max = maximum(abs(s[d][:MolarFlow][1]) for s in states for d in devices)
    for (j, d) in enumerate(devices)
        F = [s[d][:MolarFlow][1] for s in states]
        y = map(states) do s
            Fd = s[d][:MolarFlow][1]
            abs(Fd) > flow_tol * F_max || return NaN
            y_up = Fd >= 0 ? s[d][:InletComposition] : s[d][:OutletComposition]
            y_up[component, 1]
        end
        color = colors[mod1(length(colors) - j + 1, length(colors))]
        lines!(ax_F, t, F, color = color, label = String(d))
        lines!(ax_y, t, y, color = color, label = String(d))
    end
    Legend(f[1, 2], ax_F, "Device", tellheight = false)
    return f
end

function _shade_stages!(ax, stages, t_end, colors)
    t = 0.0
    while t < t_end
        for (i, st) in enumerate(stages)
            t >= t_end && break
            vspan!(ax, t, min(t + st.duration, t_end), color = (colors[mod1(i, length(colors))], 0.12))
            t += st.duration
        end
    end
    return ax
end
