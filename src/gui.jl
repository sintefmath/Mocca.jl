using Gtk4
using JSON
using CairoMakie

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
    read_params_from_entries(entries::Dict, params::Dict)

Read current values from all GUI entry widgets and update the parameter dictionary.
Returns a new Dict with updated values.
"""
function read_params_from_entries(entries::Dict, params::Dict)
    updated = deepcopy(params)
    for (section_label, section_key, param_keys) in PARAM_SECTIONS
        if !haskey(updated, section_key)
            continue
        end
        for pkey in param_keys
            entry_key = (section_key, pkey)
            if haskey(entries, entry_key) && haskey(updated[section_key], pkey)
                original_val = params[section_key][pkey]
                text = entries[entry_key][]
                try
                    updated[section_key][pkey] = string_to_value(text, original_val)
                catch e
                    @warn "Could not parse value for $section_key.$pkey: $text" exception=e
                end
            end
        end
    end
    return updated
end

"""
    populate_entries!(entries::Dict, params::Dict)

Update all GUI entry widgets with values from the parameter dictionary.
"""
function populate_entries!(entries::Dict, params::Dict)
    for (section_label, section_key, param_keys) in PARAM_SECTIONS
        if !haskey(params, section_key)
            continue
        end
        for pkey in param_keys
            entry_key = (section_key, pkey)
            if haskey(entries, entry_key) && haskey(params[section_key], pkey)
                val = params[section_key][pkey]
                buf = Gtk4.G_.get_buffer(entries[entry_key])
                buf[] = value_to_string(val)
            end
        end
    end
end

"""
    create_parameter_panel(params::Dict)

Create a scrollable panel with organized parameter editing fields.
Returns `(scroll_window, entries_dict)` where `entries_dict` maps
`(section_key, param_key)` to `GtkEntry` widgets.
"""
function create_parameter_panel(params::Dict)
    entries = Dict{Tuple{String,String}, GtkEntry}()
    main_box = GtkBox(:v, 4)

    for (section_label, section_key, param_keys) in PARAM_SECTIONS
        if !haskey(params, section_key)
            continue
        end

        # Create grid for this section's parameters
        grid = GtkGrid()
        Gtk4.G_.set_column_spacing(grid, 8)
        Gtk4.G_.set_row_spacing(grid, 4)

        row = 0
        for pkey in param_keys
            if !haskey(params[section_key], pkey)
                continue
            end
            val = params[section_key][pkey]

            label = GtkLabel(pkey)
            Gtk4.G_.set_halign(label, Gtk4.Align_START)
            Gtk4.G_.set_hexpand(label, false)

            entry = GtkEntry()
            Gtk4.G_.set_hexpand(entry, true)
            buf = Gtk4.G_.get_buffer(entry)
            buf[] = value_to_string(val)

            grid[1, row] = label
            grid[2, row] = entry
            entries[(section_key, pkey)] = entry
            row += 1
        end

        # Wrap in expander for collapsible sections
        expander = Gtk4.G_.Expander_new(section_label)
        Gtk4.G_.set_expanded(expander, true)
        Gtk4.G_.set_child(expander, grid)
        push!(main_box, expander)
    end

    # Wrap in scrolled window
    scroll = GtkScrolledWindow()
    Gtk4.G_.set_child(scroll, main_box)
    Gtk4.G_.set_hexpand(scroll, false)
    Gtk4.G_.set_vexpand(scroll, true)
    Gtk4.G_.set_min_content_width(scroll, 380)

    return scroll, entries
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
    render_plot_to_file(case, states, timesteps)

Render the outlet plot to a temporary PNG file using CairoMakie.
Returns the file path.
"""
function render_plot_to_file(case, states, timesteps)
    f = Mocca.plot_outlet(case, states, timesteps)
    tmpfile = tempname() * ".png"
    CairoMakie.save(tmpfile, f; px_per_unit=2)
    return tmpfile
end

"""
    launch_gui()

Launch the Mocca parameter editor and simulation GUI.
Loads the default haghpanah cyclic simulation parameters on startup.
"""
function launch_gui()
    # Load default parameters
    params = Ref(load_default_parameters())

    app = GtkApplication("org.mocca.simulator", 0)

    signal_connect(app, "activate") do app
        # Create main window
        win = GtkApplicationWindow(app, "Mocca Simulator")
        Gtk4.G_.set_default_size(win, 1200, 800)

        # Main vertical layout
        main_vbox = GtkBox(:v, 4)

        # === Toolbar ===
        toolbar = GtkBox(:h, 8)

        btn_import = GtkButton("Import JSON")
        btn_export = GtkButton("Export JSON")
        btn_run = GtkButton("Run Simulation")
        btn_reset = GtkButton("Reset Defaults")

        # Status label
        status_label = GtkLabel("Ready — Default haghpanah cyclic parameters loaded.")
        Gtk4.G_.set_hexpand(status_label, true)
        Gtk4.G_.set_halign(status_label, Gtk4.Align_START)

        push!(toolbar, btn_import)
        push!(toolbar, btn_export)
        push!(toolbar, btn_run)
        push!(toolbar, btn_reset)
        push!(toolbar, status_label)

        push!(main_vbox, toolbar)

        # === Content area: parameter editor (left) + plot (right) ===
        hpane = GtkPaned(:h)

        # Left: parameter editor
        param_scroll, entries = create_parameter_panel(params[])
        hpane[1] = param_scroll

        # Right: plot display area
        plot_box = GtkBox(:v, 4)
        Gtk4.G_.set_hexpand(plot_box, true)
        Gtk4.G_.set_vexpand(plot_box, true)

        plot_label = GtkLabel("Simulation results will appear here after running.")
        Gtk4.G_.set_vexpand(plot_label, true)
        Gtk4.G_.set_hexpand(plot_label, true)
        push!(plot_box, plot_label)

        # Scrolled window for plot
        plot_scroll = GtkScrolledWindow()
        Gtk4.G_.set_child(plot_scroll, plot_box)
        Gtk4.G_.set_hexpand(plot_scroll, true)
        Gtk4.G_.set_vexpand(plot_scroll, true)

        hpane[2] = plot_scroll
        Gtk4.G_.set_position(hpane, 400)

        push!(main_vbox, hpane)
        Gtk4.G_.set_vexpand(hpane, true)

        # Set main_vbox as window child
        Gtk4.G_.set_child(win, main_vbox)

        # === Signal handlers ===

        # Import JSON
        signal_connect(btn_import, "clicked") do _
            open_dialog("Import JSON Parameter File", win, ["*.json"]) do filepath
                if filepath != "" && isfile(filepath)
                    try
                        params[] = load_parameters_from_file(filepath)
                        populate_entries!(entries, params[])
                        Gtk4.G_.set_label(status_label, "Loaded parameters from: $(basename(filepath))")
                    catch e
                        Gtk4.G_.set_label(status_label, "Error loading file: $e")
                    end
                end
            end
        end

        # Export JSON
        signal_connect(btn_export, "clicked") do _
            save_dialog("Export JSON Parameter File", win, ["*.json"]) do filepath
                if filepath != ""
                    try
                        current_params = read_params_from_entries(entries, params[])
                        if !endswith(filepath, ".json")
                            filepath = filepath * ".json"
                        end
                        export_parameters_to_file(filepath, current_params)
                        params[] = current_params
                        Gtk4.G_.set_label(status_label, "Parameters exported to: $(basename(filepath))")
                    catch e
                        Gtk4.G_.set_label(status_label, "Error exporting: $e")
                    end
                end
            end
        end

        # Run Simulation
        signal_connect(btn_run, "clicked") do _
            try
                Gtk4.G_.set_label(status_label, "Running simulation... please wait.")

                # Read current parameters from GUI
                current_params = read_params_from_entries(entries, params[])
                params[] = current_params

                # Run simulation
                case, states, timesteps = run_simulation_from_params(current_params)

                # Render plot
                plot_file = render_plot_to_file(case, states, timesteps)

                # Update plot display
                # Remove old children from plot_box
                for child in collect(plot_box)
                    delete!(plot_box, child)
                end

                pic = GtkPicture(plot_file)
                Gtk4.G_.set_can_shrink(pic, true)
                Gtk4.G_.set_hexpand(pic, true)
                Gtk4.G_.set_vexpand(pic, true)
                push!(plot_box, pic)

                Gtk4.G_.set_label(status_label, "Simulation complete. Results plotted.")
            catch e
                Gtk4.G_.set_label(status_label, "Simulation error: $e")
                @error "Simulation failed" exception=(e, catch_backtrace())
            end
        end

        # Reset Defaults
        signal_connect(btn_reset, "clicked") do _
            params[] = load_default_parameters()
            populate_entries!(entries, params[])
            Gtk4.G_.set_label(status_label, "Reset to default haghpanah cyclic parameters.")
        end

        show(win)
    end

    # Run the application
    Gtk4.run(app)
end
