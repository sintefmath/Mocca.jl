using GLMakie
using JSON

const DEFAULT_CYCLIC_JSON = joinpath(@__DIR__, "..", "models", "json", "haghpanah_cyclic_input_simple.json")

# Section display names and their JSON keys with parameter ordering
const PARAM_SECTIONS = [
    ("Physical Constants",    "physicalConstants",   ["molecularMassOfCO2", "molecularMassOfN2", "R"]),
    ("DSL Isotherm",          "dslPars",             ["b0", "d0", "ΔUbi", "ΔUdi", "qsbi", "qsdi"]),
    ("Adsorbent Properties",  "adsorbentProps",      ["ϵ_p", "D_m", "τ", "d_p", "V0_inter", "ρ_s", "C_pa", "C_ps"]),
    ("Column Properties",     "columnProps",         ["Φ", "K_z", "K_w", "r_in", "r_out", "h_in", "h_out", "ρ_w", "C_pw", "L"]),
    ("Feed Properties",       "feedProps",           ["fluid_viscosity", "ρ_g", "C_pg", "T_feed", "v_feed", "y_feed"]),
    ("Boundary Conditions",   "boundaryConditions",  ["T_a", "p_high", "p_intermediate", "p_low", "λ"]),
    ("Initial Conditions",    "initialConditions",   ["P_init", "T0", "Tw_init", "y_init"]),
    ("Process Specification", "processSpecification",["stage_types", "stage_durations", "num_cycles"]),
    ("Simulation",            "simulation",          ["system_type", "ncells", "maxdt", "timestep_selectors"]),
    ("Solver",                "solver",              ["linear_solver", "info_level"]),
]

"""
    load_default_parameters()

Load the default haghpanah cyclic simulation parameters from the simple JSON format.
Returns a `Dict{String, Any}`.
"""
function load_default_parameters()
    return JSON.parse(read(DEFAULT_CYCLIC_JSON, String); dicttype=Dict{String, Any})
end

"""
    load_parameters_from_file(filepath::String)

Load parameters from a JSON file. Supports both simple and detailed formats.
Returns a `Dict{String, Any}`.
"""
function load_parameters_from_file(filepath::String)
    return JSON.parse(read(filepath, String); dicttype=Dict{String, Any})
end

"""
    export_parameters_to_file(filepath::String, params::Dict)

Export the current parameter dictionary to a JSON file in simple format.
"""
function export_parameters_to_file(filepath::String, params::Dict)
    open(filepath, "w") do io
        JSON.print(io, params, 2)
    end
end

"""
    value_to_string(val)

Convert a parameter value to a string representation for display in the GUI.
"""
function value_to_string(val)
    if isa(val, AbstractVector)
        return join(string.(val), ", ")
    elseif isa(val, Dict)
        return "{}"
    else
        return string(val)
    end
end

"""
    string_to_value(str::String, original_val)

Parse a string back into the appropriate type matching the original value.
"""
function string_to_value(str::String, original_val)
    if isa(original_val, AbstractVector)
        parts = split(str, ",")
        if isa(original_val, AbstractVector{<:AbstractString}) || (!isempty(original_val) && isa(first(original_val), AbstractString))
            return String.(strip.(parts))
        else
            return parse.(Float64, strip.(parts))
        end
    elseif isa(original_val, Integer)
        return parse(Int, strip(str))
    elseif isa(original_val, AbstractFloat)
        return parse(Float64, strip(str))
    elseif isa(original_val, AbstractString)
        return strip(str)
    elseif isa(original_val, Dict)
        return original_val
    else
        return str
    end
end

"""
    read_params_from_textboxes(textboxes::Dict, params::Dict)

Read current values from all GUI Textbox widgets and update the parameter dictionary.
Returns a new Dict with updated values.
"""
function read_params_from_textboxes(textboxes::Dict, params::Dict)
    updated = deepcopy(params)
    for (section_label, section_key, param_keys) in PARAM_SECTIONS
        if !haskey(updated, section_key)
            continue
        end
        for pkey in param_keys
            entry_key = (section_key, pkey)
            if haskey(textboxes, entry_key) && haskey(updated[section_key], pkey)
                original_val = params[section_key][pkey]
                text = textboxes[entry_key].stored_string[]
                if !isnothing(text)
                    try
                        updated[section_key][pkey] = string_to_value(text, original_val)
                    catch e
                        @warn "Could not parse value for $section_key.$pkey: $text" exception=e
                    end
                end
            end
        end
    end
    return updated
end

"""
    set_textbox_text!(tb::Textbox, text::String)

Set the displayed and stored text of a GLMakie Textbox widget.
"""
function set_textbox_text!(tb, text::String)
    tb.displayed_string[] = text
    tb.stored_string[] = text
end

"""
    populate_textboxes!(textboxes::Dict, params::Dict)

Update all GUI Textbox widgets with values from the parameter dictionary.
"""
function populate_textboxes!(textboxes::Dict, params::Dict)
    for (section_label, section_key, param_keys) in PARAM_SECTIONS
        if !haskey(params, section_key)
            continue
        end
        for pkey in param_keys
            entry_key = (section_key, pkey)
            if haskey(textboxes, entry_key) && haskey(params[section_key], pkey)
                val = params[section_key][pkey]
                set_textbox_text!(textboxes[entry_key], value_to_string(val))
            end
        end
    end
end

"""
    run_simulation_from_params(params::Dict)

Run a Mocca simulation using the given parameter dictionary.
Returns `(case, states, timesteps)`.
"""
function run_simulation_from_params(params::Dict)
    constants, info = Mocca.parse_input(params)
    case, ts_config, info_level = Mocca.setup_mocca_case(constants, info)
    states, timesteps = Mocca.simulate_process(case;
        timestep_selector_cfg = ts_config,
        output_substates = true,
        info_level = info_level
    )
    return case, states, timesteps
end

"""
    plot_results!(fig, case, states, timesteps)

Plot simulation results (outlet cell variables over time) into the given Figure layout.
Clears existing plot axes in the right-side grid and draws new ones.
"""
function plot_results!(plot_grid, case, states, timesteps)
    # Clear previous content
    empty!(plot_grid)

    model = case.model
    pvars = model.primary_variables.keys
    comp_names = model.system.component_names
    unit_label = units_dict()
    pretty_names = prettyVarNames(pvars)
    outlet_cell = size(states[1][:y][1,:], 1)
    t = Float64.(cumsum(timesteps))

    r = 1
    for (i, symb) in enumerate(pvars)
        c = i
        if i > 3
            c = i - 3
            r = 2
        end
        ax = Axis(plot_grid[r, c],
            title = pretty_names[symb],
            xlabel = "t [s]",
            ylabel = string(unit_label[symb]))

        if size(states[end][symb], 2) == 1
            lines!(ax, t, [result[symb][outlet_cell] for result in states])
        else
            for k in 1:size(states[end][symb], 1)
                lines!(ax, t, [result[symb][k, outlet_cell] for result in states], label = comp_names[k])
            end
            Legend(plot_grid[2, 3], ax, tellwidth=false)
        end
    end
end

# Web-inspired breakpoints (in pixels) for responsive layout
const BREAKPOINT_SM = 768
const BREAKPOINT_MD = 992
const BREAKPOINT_LG = 1200
const BREAKPOINT_XL = 1400

"""
    responsive_colsizes!(fig, window_width)

Adjust main figure column sizes based on window width using web-inspired breakpoints.
At smaller sizes the parameter panel gets proportionally more space so fields remain usable.
"""
function responsive_colsizes!(fig, window_width)
    if window_width < BREAKPOINT_SM
        colsize!(fig.layout, 1, Relative(0.55))
        colsize!(fig.layout, 2, Relative(0.45))
    elseif window_width < BREAKPOINT_MD
        colsize!(fig.layout, 1, Relative(0.50))
        colsize!(fig.layout, 2, Relative(0.50))
    elseif window_width < BREAKPOINT_LG
        colsize!(fig.layout, 1, Relative(0.45))
        colsize!(fig.layout, 2, Relative(0.55))
    elseif window_width < BREAKPOINT_XL
        colsize!(fig.layout, 1, Relative(0.40))
        colsize!(fig.layout, 2, Relative(0.60))
    else
        colsize!(fig.layout, 1, Relative(0.35))
        colsize!(fig.layout, 2, Relative(0.65))
    end
end

"""
    launch_gui()

Launch the Mocca parameter editor and simulation GUI using GLMakie.
Loads the default haghpanah cyclic simulation parameters on startup.

The layout is responsive to window size:
- Parameters are arranged in two columns so all fields are visible without scrolling.
- Column ratios adjust automatically when the window is resized, following
  web-inspired breakpoints (sm 768, md 992, lg 1200, xl 1400 px).
"""
function launch_gui()
    # Load default parameters
    params = Ref(load_default_parameters())

    fig = Figure(size = (1400, 1100), title = "Mocca Simulator")

    # === Top bar: toolbar + file I/O + status (spans full width) ===
    top_bar = fig[1, 1:2] = GridLayout()
    rowsize!(fig.layout, 1, Auto())

    btn_run = Button(top_bar[1, 1]; label="Run Simulation", fontsize=12)
    btn_reset = Button(top_bar[1, 2]; label="Reset Defaults", fontsize=12)
    Label(top_bar[1, 3], text="JSON path:", fontsize=11, halign=:right)
    filepath_tb = Textbox(top_bar[1, 4]; stored_string=DEFAULT_CYCLIC_JSON, fontsize=10)
    btn_import = Button(top_bar[1, 5]; label="Import", fontsize=11)
    btn_export = Button(top_bar[1, 6]; label="Export", fontsize=11)

    status_label = Label(top_bar[2, 1:6],
        text="Ready — Default haghpanah cyclic parameters loaded.",
        fontsize=12, halign=:left)

    # === Content: parameter editor (left) + plot area (right) ===
    param_layout = fig[2, 1] = GridLayout()
    right_layout = fig[2, 2] = GridLayout()

    # Set initial responsive column sizes
    responsive_colsizes!(fig, 1400)

    # === Parameter editing: split into 2 columns for better vertical fit ===
    # Col 1: Physical Constants(4), DSL Isotherm(7), Adsorbent Props(9), Column Props(11) ≈ 31 rows
    # Col 2: Feed Props(7), Boundary Cond(6), Initial Cond(5), Process Spec(4), Simulation(5), Solver(3) ≈ 30 rows
    param_col1 = param_layout[1, 1] = GridLayout()
    param_col2 = param_layout[1, 2] = GridLayout()

    split_at = 4  # first 4 sections in column 1, remaining in column 2

    textboxes = Dict{Tuple{String,String}, Textbox}()

    for (col_idx, sections) in enumerate([
            PARAM_SECTIONS[1:split_at],
            PARAM_SECTIONS[split_at+1:end]
        ])
        target = col_idx == 1 ? param_col1 : param_col2
        row = 1

        for (section_label, section_key, param_keys) in sections
            if !haskey(params[], section_key)
                continue
            end

            # Section header
            Label(target[row, 1:2], text=section_label,
                fontsize=13, halign=:left, color=:steelblue)
            row += 1

            for pkey in param_keys
                if !haskey(params[][section_key], pkey)
                    continue
                end
                val = params[][section_key][pkey]

                Label(target[row, 1], text=pkey,
                    fontsize=11, halign=:right)
                tb = Textbox(target[row, 2];
                    stored_string=value_to_string(val),
                    fontsize=11)
                textboxes[(section_key, pkey)] = tb
                row += 1
            end
        end
    end

    # === Plot area (right side) ===
    plot_grid = right_layout[1, 1] = GridLayout()
    Label(plot_grid[1, 1], text="Simulation results will appear here after running.",
        fontsize=14, color=:gray)

    # === Responsive layout: adjust column ratios on window resize ===
    on(fig.scene.events.window_area) do area
        responsive_colsizes!(fig, area.widths[1])
    end

    # === Button callbacks ===

    # Import JSON from file path in textbox
    on(btn_import.clicks) do _
        try
            filepath = filepath_tb.stored_string[]
            if !isnothing(filepath) && filepath != "" && isfile(filepath)
                params[] = load_parameters_from_file(filepath)
                populate_textboxes!(textboxes, params[])
                status_label.text[] = "Loaded parameters from: $(basename(filepath))"
            else
                status_label.text[] = "File not found: $(something(filepath, ""))"
            end
        catch e
            status_label.text[] = "Error importing: $e"
        end
    end

    # Export JSON to file path in textbox
    on(btn_export.clicks) do _
        try
            filepath = filepath_tb.stored_string[]
            if !isnothing(filepath) && filepath != ""
                current_params = read_params_from_textboxes(textboxes, params[])
                if !endswith(filepath, ".json")
                    filepath = filepath * ".json"
                end
                export_parameters_to_file(filepath, current_params)
                params[] = current_params
                status_label.text[] = "Parameters exported to: $(basename(filepath))"
            else
                status_label.text[] = "Enter a file path to export."
            end
        catch e
            status_label.text[] = "Error exporting: $e"
        end
    end

    # Run Simulation
    on(btn_run.clicks) do _
        try
            status_label.text[] = "Running simulation... please wait."

            current_params = read_params_from_textboxes(textboxes, params[])
            params[] = current_params

            case, states, timesteps = run_simulation_from_params(current_params)

            plot_results!(plot_grid, case, states, timesteps)

            status_label.text[] = "Simulation complete. Results plotted."
        catch e
            status_label.text[] = "Simulation error: $e"
            @error "Simulation failed" exception=(e, catch_backtrace())
        end
    end

    # Reset Defaults
    on(btn_reset.clicks) do _
        params[] = load_default_parameters()
        populate_textboxes!(textboxes, params[])
        status_label.text[] = "Reset to default haghpanah cyclic parameters."
    end

    display(fig)
    return fig
end
