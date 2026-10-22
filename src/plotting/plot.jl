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

# Colours shared with the documentation: one per kind of step and one per
# component, so that a reader learns them once. Anything else falls back to
# Makie's default palette.
const STEP_COLOURS = [
    "lpp" => "#c79ad9", "light product" => "#c79ad9", "purge" => "#c79ad9",
    "pressur" => "#8a63b8", "adsor" => "#3d8c47", "blowdown" => "#c48a1f",
    "evacua" => "#3f63c9", "equali" => "#8d8a94", "idle" => "#b9b6bf",
]
const COMPONENT_COLOURS = Dict("CO2" => "#c4362e", "N2" => "#4a5a75", "H2O" => "#2b8cbe")
const SINGLE_COLOUR = "#5b2f73"

# The colour of the kind of step that a stage or label name starts with, or
# `fallback`
function _step_colour(name, fallback)
    n = lowercase(String(name))
    hits = [(first(findfirst(k, n)), c) for (k, c) in STEP_COLOURS if occursin(k, n)]
    isempty(hits) && return fallback
    return last(hits[argmin(first.(hits))])
end

_component_colour(name, k) = get(COMPONENT_COLOURS, String(name), Makie.wong_colors()[mod1(k, 7)])

function prettyVarNames(vars::Vector{Symbol})
    pretty_names = Dict()

    for (i, symb) in enumerate(vars)
        var = String(symb)
        words = split(var, r"(?=[A-Z])")
        pretty_names[symb] = join(words, " ")
    end

    return pretty_names
end


"""
    plot_cell(states, model, timesteps, cell)

History of every variable of `cell` of a bed through time, one panel per
variable. `cell` is a port, such as `:top` for the end the light product
leaves by or `:bottom` for the feed end, whose values are the means over the
cross-section (see [`cross_section_mean`](@ref)); a cell `(i, j)`, `i`-th
along the column and in ring `j` from the axis (1 for a bed of one radial
cell); or the index of a cell, which for a bed of one radial cell is its
position along the column. `states` belong to the bed alone: for a process,
pass `[s[:Bed] for s in result.states]`, with `model` the bed
(`process.flowsheet[:Bed]`) and the `timesteps` of the result.
"""
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

        values = [_at_cell(model, result[symb], cell) for result in states]
        if states[end][symb] isa AbstractVector
           lines!(ax, t, values, color = SINGLE_COLOUR)
        else
            for k in 1:size(states[end][symb], 1)
                lines!(ax, t, [v[k] for v in values], label = comp_names[k],
                    color = _component_colour(comp_names[k], k))
            end
            leg = Legend(f[2,3], ax, tellwidth=false)
            
        end

        
    end

    return f
end


"""
    plot_state(state, model)

Every variable along a bed at one time, one panel per variable, averaged over
the cross-section (see [`axial_profile`](@ref)). `state` is the state of the
bed alone, such as `result.states[end][:Bed]` at the end of a simulation, and
`model` the bed (`process.flowsheet[:Bed]`).
"""
function plot_state(state, model)

    pvars = model.primary_variables.keys
    comp_names = model.system.component_names
    pretty_names = prettyVarNames(pvars)

    units = units_dict()
    x = axial_positions(model)

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

        v = axial_profile(model, state[symb])
        if v isa AbstractVector
           lines!(ax, x, v, color = SINGLE_COLOUR)
        else
            for k in 1:size(v, 1)
                lines!(ax, x, v[k,:], label = comp_names[k], color = _component_colour(comp_names[k], k))
            end
            leg = Legend(f[2,3], ax, tellwidth=false)
            
        end

        
    end

    return f
end

function plot_outlet(case, states, timesteps_out)

    # Convenience function to plot outlet cell variables over time

    outlet_cell = size(states[1][:y][1,:],1)
    f_outlet = plot_cell(states, case.model, timesteps_out, outlet_cell);

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
panel per primary variable. `cell` is a port, a cell `(i, j)` or the index of
a cell, as for [`plot_cell`](@ref). `runs` is a vector of `label => (states,
timesteps)`, where `states` belong to that unit alone (for a flowsheet, pass
`[s[:Bed] for s in states]`). For variables with one value per component only
`component` is shown.
"""
function plot_cell_comparison(runs, model, cell; component = 1)
    curves = map(runs) do (label, (states, timesteps))
        t = Float64.(cumsum(timesteps))
        label => (t, symb -> [_component_value(_at_cell(model, s[symb], cell), component) for s in states])
    end
    return _comparison_figure(model, curves, L"t\; [s]", component)
end

"""
    plot_state_comparison(states, model; component = 1)

Overlay profiles along the column, averaged over the cross-section (see
[`axial_profile`](@ref)), from several states of the same unit. `states` is a
vector of `label => state`.
"""
function plot_state_comparison(states, model; component = 1)
    x = axial_positions(model)
    curves = map(states) do (label, state)
        label => (x, symb -> (v = axial_profile(model, state[symb]); [_component_value(v, component, i) for i in eachindex(x)]))
    end
    return _comparison_figure(model, curves, L"x\; [m]", component)
end

_component_value(v::AbstractVector, component, cell) = v[cell]
_component_value(v::AbstractMatrix, component, cell) = v[component, cell]
_component_value(v::Number, component) = v
_component_value(v::AbstractVector, component) = v[component]

# The value of a cell variable `v` of `model` at `cell`, as plot_cell takes
# it: a number for a variable with one value per cell, and a vector for one
# with a column per cell
_at_cell(model, v, cell::Integer) = v isa AbstractVector ? v[cell] : v[:, cell]
_at_cell(model, v, port::Symbol) = cross_section_mean(model, v, port_cells(model, port))
function _at_cell(model, v, (i, j)::Tuple{Integer, Integer})
    layers = axial_layers(model)
    nz, nr = length(layers), length(first(layers))
    (1 <= i <= nz && 1 <= j <= nr) || throw(ArgumentError("cell ($i, $j) is outside the bed of $nz × $nr cells"))
    return _at_cell(model, v, i + (j - 1) * nz)
end

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
        Legend(f[2, 2], [PolyElement(color = (_step_colour(st.name, colors[mod1(i, length(colors))]), 0.25)) for (i, st) in enumerate(stages)],
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
            vspan!(ax, t, min(t + st.duration, t_end), color = (_step_colour(st.name, colors[mod1(i, length(colors))]), 0.12))
            t += st.duration
        end
    end
    return ax
end

"""
    plot_profiles(states, timesteps, model, times; labels, unit = nothing, components = all)

Profiles along a bed at the given `times` [s], measured from the start of
`states`, averaged over the cross-section (see [`axial_profile`](@ref)): the
gas mole fraction of each of `components` (top row) and its loading (bottom
row). For a dual-site Langmuir isotherm the loading is shown relative to the
component's saturation capacity `qsb + qsd`, as in the bed profile figures of
Haghpanah et al. and Krishnamurthy et al. For flowsheet states, `unit` names
the bed and `model` is the bed's model. Each time uses the last state at or
before it.
"""
function plot_profiles(states, timesteps, model, times;
        labels = ["t = $(round(t, digits = 1)) s" for t in times],
        unit = nothing,
        components = eachindex(model.system.component_names))
    t = cumsum(timesteps)
    unit_states = isnothing(unit) ? states : [s[unit] for s in states]
    names = model.system.component_names
    x = axial_positions(model)
    x = x ./ sum(model.data_domain[:dx][first(l)] for l in axial_layers(model))
    q_scale = saturation_capacity(model.system.isotherm)
    colors = Makie.wong_colors()
    f = Figure(size = (300 * length(components), 600))
    ax = nothing
    for (col, i) in enumerate(components)
        ax_y = Axis(f[1, col], title = "$(names[i]) in gas", ylabel = L"y\; [-]")
        q_label = isnothing(q_scale) ? L"q\; [mol/m^{3}]" : L"q/q_s\; [-]"
        ax = Axis(f[2, col], title = "$(names[i]) adsorbed", xlabel = "Dimensionless bed length [-]", ylabel = q_label)
        for (j, (tj, label)) in enumerate(zip(times, labels))
            k = findlast(<=(tj + 1e-6), t)
            isnothing(k) && error("No state at or before t = $tj s")
            s = unit_states[k]
            q = axial_profile(model, s[:AdsorbedConcentration])[i, :]
            isnothing(q_scale) || (q = q ./ q_scale[i])
            color = _step_colour(label, colors[mod1(j, length(colors))])
            lines!(ax_y, x, axial_profile(model, s[:y])[i, :], color = color, label = label)
            lines!(ax, x, q, color = color, label = label)
        end
    end
    Legend(f[1:2, length(components) + 1], ax)
    return f
end


# Edges [m] of the rings of a bed from the axis to the wall: those of a bed
# meshed across its radius, or the axis and the wall for one radial cell
function _ring_edges(model)
    d = model.data_domain
    _is_axisymmetric(d) || return [0.0, first(d[:r_in, Unit()])]
    return [0.0; cumsum(collect(Jutul.physical_representation(d).deltas[2]))]
end

# Cells along the column and across the radius of a bed
function _bed_shape(model)
    layers = axial_layers(model)
    return (length(layers), length(first(layers)))
end

# Values of `variable` of a bed state on its cells, one component of a
# variable with one value per component, as a matrix of cells along the
# column by rings from the axis
function _bed_field(model, s, variable, component)
    v = s[variable]
    return reshape(v isa AbstractVector ? v : v[component, :], _bed_shape(model)...)
end

_field_title(model, variable, component, v) =
    (name = prettyVarNames([variable])[variable]; v isa AbstractVector ? name : "$name, $(model.system.component_names[component])")

"""
    plot_bed_section(state, model; variable = :Temperature, component = 1, unit = nothing,
        colormap = :viridis, colorrange = nothing, radial_scale = 1)

A cut through a bed along its axis at one time: `variable` over the radius
and the length, the bed upright with its bottom port at the foot and mirrored
about the axis, so that it spans the diameter, with the wall from the inner to
the outer radius. For `variable = :Temperature` the wall shows its own
temperature, and otherwise it is grey. For a variable with one value per
component, `component` is shown. A bed of one radial cell is uniform across.

`state` is the state of the bed alone, or, with `unit` naming the bed, a
state of the flowsheet. `radial_scale` widens the plot, for a column much
longer than it is wide; the axis keeps the true radius.
"""
function plot_bed_section(state, model; variable = :Temperature, component = 1, unit = nothing,
        colormap = :viridis, colorrange = nothing, radial_scale = 1.0)
    s = isnothing(unit) ? state : state[unit]
    nz, nr = _bed_shape(model)
    F = _bed_field(model, s, variable, component)
    layers = axial_layers(model)
    z_edges = [0.0; cumsum(model.data_domain[:dx][first.(layers)])]
    r_edges = _ring_edges(model)
    r_in, r_out = last(r_edges), first(model.data_domain[:r_out, Unit()])
    L = last(z_edges)
    x_edges = [-reverse(r_edges[2:end]); r_edges]
    M = permutedims([reverse(F, dims = 2) F])
    wall = variable == :Temperature && has_wall(model.system) && haskey(s, :WallTemperature) && r_out > r_in
    Tw = wall ? s[:WallTemperature][[i + (nr - 1) * nz for i in 1:nz]] : Float64[]
    lo, hi = extrema([vec(F); Tw])
    cr = something(colorrange, lo == hi ? (lo - 0.5, hi + 0.5) : (lo, hi))

    f = Figure(size = (450, 700))
    ax = Axis(f[1, 1], title = _field_title(model, variable, component, s[variable]),
        xlabel = L"r\; [m]", ylabel = L"z\; [m]", xticks = LinearTicks(3))
    hm = heatmap!(ax, x_edges, z_edges, M; colormap = colormap, colorrange = cr)
    if r_out > r_in
        for xs in ([r_in, r_out], [-r_out, -r_in])
            if wall
                heatmap!(ax, xs, z_edges, reshape(Tw, 1, :); colormap = colormap, colorrange = cr)
            else
                poly!(ax, Rect(xs[1], 0.0, xs[2] - xs[1], L), color = (:grey, 0.6))
            end
        end
    end
    limits!(ax, -r_out, r_out, 0.0, L)
    Colorbar(f[1, 2], hm; label = L"%$(units_dict()[variable])")
    # The width of the bed from its height, so that the colour bar is as tall
    # as the bed, and the figure fitted to it
    colsize!(f.layout, 1, Aspect(1, 2 * r_out * radial_scale / L))
    resize_to_layout!(f)
    return f
end

"""
    plot_radial_profiles(states, timesteps, model, times; positions = [0.25, 0.5, 0.75],
        variable = :Temperature, component = 1, labels, unit = nothing)

Profiles of `variable` across the radius of a bed, from the axis to the wall,
one panel for each relative position along the column in `positions` (at the
nearest cells), and one line for each of `times` [s], measured from the start
of `states`, each the last state at or before it. For `variable =
:Temperature` the temperature of the wall around the outer ring is marked
with a square at the middle of the wall. For a variable with one value per
component, `component` is shown. For flowsheet states, `unit` names the bed
and `model` is the bed's model.
"""
function plot_radial_profiles(states, timesteps, model, times; positions = [0.25, 0.5, 0.75],
        variable = :Temperature, component = 1, labels = ["t = $(round(t, digits = 1)) s" for t in times], unit = nothing)
    t = cumsum(timesteps)
    unit_states = isnothing(unit) ? states : [s[unit] for s in states]
    nz, nr = _bed_shape(model)
    layers = axial_layers(model)
    z = axial_positions(model) ./ sum(model.data_domain[:dx][first.(layers)])
    r_edges = _ring_edges(model)
    rc = (r_edges[1:(end - 1)] .+ r_edges[2:end]) ./ 2
    r_in, r_out = last(r_edges), first(model.data_domain[:r_out, Unit()])
    wall = variable == :Temperature && has_wall(model.system)
    colors = Makie.wong_colors()
    f = Figure(size = (300 * length(positions) + 200, 420))
    Label(f[0, 1:length(positions)], _field_title(model, variable, component, first(unit_states)[variable]), font = :bold, tellwidth = false)
    ax = nothing
    entries = Any[]
    for (col, pos) in enumerate(positions)
        i = argmin(abs.(z .- pos))
        ax = Axis(f[1, col], title = "z/L = $(round(z[i], digits = 2))", xlabel = L"r\; [m]", xticks = LinearTicks(4),
            ylabel = col == 1 ? L"%$(units_dict()[variable])" : "")
        for (k, (tk, label)) in enumerate(zip(times, labels))
            n = findlast(<=(tk + 1e-6), t)
            isnothing(n) && error("No state at or before t = $tk s")
            s = unit_states[n]
            F = _bed_field(model, s, variable, component)
            color = _step_colour(label, colors[mod1(k, length(colors))])
            scatterlines!(ax, rc, F[i, :], color = color)
            col == 1 && push!(entries, [LineElement(color = color), MarkerElement(color = color, marker = :circle)] => label)
            if wall && haskey(s, :WallTemperature)
                scatter!(ax, [(r_in + r_out) / 2], [s[:WallTemperature][i + (nr - 1) * nz]], color = color, marker = :rect)
            end
        end
        wall && r_out > r_in && vlines!(ax, [r_in], color = :grey, linestyle = :dash)
    end
    wall && push!(entries, MarkerElement(color = :grey, marker = :rect) => "wall")
    Legend(f[1, length(positions) + 1], first.(entries), last.(entries))
    return f
end
