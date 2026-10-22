# Results written to CSV files for a spreadsheet: a summary with the
# performance and the mass balance, a stream table, the streams stage by stage
# and step by step, and the beds through time and along their length. Values
# are in process engineering units (bara, °C, mol%, kg/h, Nm³/h, mol/kg)
# rather than SI.

# Molar volume [m³/mol] of an ideal gas at normal conditions, 0 °C and 1 atm
const NORMAL_MOLAR_VOLUME = GAS_CONSTANT * 273.15 / 101325.0

# Flows below this fraction of the largest flow through any device count as
# no flow: a closed device, or a stream whose composition is undefined
const NO_FLOW_FRACTION = 1e-6

const RESULT_FILES = ("summary.csv", "streams.csv", "stage_streams.csv", "stream_history.csv",
    "bed_history.csv", "bed_profiles.csv")

"""
    write_results(dir, r::ProcessResult; sigdigits = 8, delim = ',', decimal = '.') → Vector{String}

Write result `r` to CSV files in directory `dir` (created if needed, files
overwritten), and the process to `process.json` (see [`write_process`](@ref))
to rerun it. Returns the paths written.

The files hold, with the units in the column headers:

| File | One row per | Contents |
|---|---|---|
| `summary.csv` | quantity | Feed, cycle, purity, recovery, productivity, energy, mass balance and beds, as `Section, Quantity, Value, Unit` |
| `streams.csv` | device | Stream table over the last cycle: composition, mean flow, amount per cycle, temperature and pressures |
| `stage_streams.csv` | stage and device | The same, for each stage of the last cycle, for devices with flow |
| `stream_history.csv` | time step and device | Flow, composition, temperature and pressures through time |
| `bed_history.csv` | time step, bed and end | Pressure, temperatures, gas composition and loadings at the bottom and top of each bed |
| `bed_profiles.csv` | snapshot, bed and cell | The same along each bed, at the start of the last cycle and the end of each of its stages |

Pressures are absolute, in bar (bara); temperatures in °C; compositions in
mol%; loadings in mol per kg of adsorbent; normal volumes (Nm³) at 0 °C and
1 atm. The composition and temperature of a stream are those of the gas
upstream of the device, and are left empty while nothing flows. Flows are
positive from the device's inlet (`From`) to its outlet (`To`). Amounts of a
stream that flows both ways are net amounts, with the gas that flowed back
under `Reverse flow`; the composition of net amounts is left empty if some
components went each way.

The time-series files hold one row per time step of the solver, with the cycle
and the stage it belongs to, so that they can be filtered or pivoted. For
spreadsheets set up with a decimal comma, pass `delim = ';', decimal = ','`.
"""
function write_results(dir::AbstractString, r::ProcessResult; sigdigits = 8, delim = ',', decimal = '.')
    decimal in ('.', ',') || input_error("write_results", "decimal", "use '.' or ',', got '$decimal'")
    delim == decimal && input_error("write_results", "delim", "the delimiter and the decimal mark are both '$delim'; use delim = ';' with decimal = ','")
    mkpath(dir)
    tables = (_summary_table(r), _stream_table(r), _stage_stream_table(r), _stream_history_table(r),
        _bed_history_table(r), _bed_profile_table(r))
    paths = [_write_csv(joinpath(dir, f), t; sigdigits = sigdigits, delim = delim, decimal = decimal)
             for (f, t) in zip(RESULT_FILES, tables)]
    json = joinpath(dir, "process.json")
    try
        push!(paths, write_process(json, r.process))
    catch e
        e isa MoccaInputError || rethrow()
        @warn "Could not write the process to $json: $(sprint(showerror, e))"
    end
    return paths
end

# ------------------------------------------------------------------------------
# Writing
# ------------------------------------------------------------------------------

function _write_csv(path, table; sigdigits, delim, decimal)
    open(path, "w") do io
        println(io, join((_csv_field(h, delim) for h in table.header), delim))
        for row in table.rows
            println(io, join((_csv_field(_csv_value(x, sigdigits, decimal), delim) for x in row), delim))
        end
    end
    return path
end

_csv_value(x::AbstractFloat, sigdigits, decimal) =
    isfinite(x) ? replace(string(round(Float64(x); sigdigits = sigdigits)), '.' => decimal) : ""
_csv_value(x::Integer, sigdigits, decimal) = string(x)
_csv_value(x::Bool, sigdigits, decimal) = x ? "yes" : "no"
_csv_value(::Nothing, sigdigits, decimal) = ""
_csv_value(x, sigdigits, decimal) = string(x)

function _csv_field(s::AbstractString, delim)
    needs_quotes = occursin(delim, s) || occursin('"', s) || occursin('\n', s) || s != strip(s)
    return needs_quotes ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end

_bara(P) = P / 1e5
_degC(T) = T - 273.15
_percent(x) = 100 * x
# Composition [mol%] of net amounts `n`, undefined if they flow both ways
function _fractions_percent(n)
    (iszero(sum(n)) || (any(<(0), n) && any(>(0), n))) && return fill(NaN, length(n))
    return _percent.(n ./ sum(n))
end

# ------------------------------------------------------------------------------
# Times, stages and streams
# ------------------------------------------------------------------------------

# Cycle, times and stage of each time step of `r`, the stage found from the
# middle of the step, since the stages are split into whole steps
function _step_times(r::ProcessResult)
    sts = stages(r.process)
    T = sum(st -> st.duration, sts)
    stage_end = cumsum([st.duration for st in sts])
    t = cumsum(r.timesteps)
    return map(eachindex(r.timesteps)) do k
        mid = t[k] - r.timesteps[k] / 2
        c = fld(mid, T)
        j = something(findfirst(>(mid - c * T), stage_end), length(sts))
        (time = t[k], cycle = Int(c) + 1, time_in_cycle = t[k] - c * T, stage_number = j, stage = sts[j].name)
    end
end

_start_of_last_cycle(r::ProcessResult) = (r.cycles - 1) * cycle_time(r.process)

_devices(fs) = [d for d in unit_names(fs) if is_device(fs, d)]

_role_name(fs, d) = (role = device_role(fs, d); isnothing(role) ? "" : string(role))

# What a device's port is connected to: a unit port, or a boundary
function _port_label(fs, d, port)
    for c in attachments(fs, d)
        c.device.port == port && return "$(c.unit.name) $(c.unit.port)"
    end
    return "boundary"
end

_flow_tol(r::ProcessResult) =
    NO_FLOW_FRACTION * maximum((abs(s[d][:MolarFlow][1]) for s in r.states for d in _devices(r.process.flowsheet)); init = 0.0)

# Flow, and the composition and temperature of the gas upstream
function _upstream(s)
    F = s[:MolarFlow][1]
    forward = F >= 0
    y = forward ? s[:InletComposition][:, 1] : s[:OutletComposition][:, 1]
    T = forward ? s[:InletTemperature][1] : s[:OutletTemperature][1]
    return F, y, T
end

# Totals of the stream through device `d` over the steps in `window`; the
# amounts match stream_totals
function _stream_summary(r::ProcessResult, d, window, tol)
    N = number_of_components(r.process.flowsheet.gas)
    n = zeros(N)
    backflow = 0.0
    flowing_time = 0.0
    weight = 0.0
    FT = 0.0
    P_in = [Inf, -Inf]
    P_out = [Inf, -Inf]
    for k in window
        s = r.states[k][d]
        dt = r.timesteps[k]
        F, y, T = _upstream(s)
        n .+= F .* y .* dt
        backflow += max(-F, 0.0) * dt
        abs(F) > tol || continue
        flowing_time += dt
        weight += abs(F) * dt
        FT += abs(F) * dt * T
        P_in .= (min(P_in[1], s[:InletPressure][1]), max(P_in[2], s[:InletPressure][1]))
        P_out .= (min(P_out[1], s[:OutletPressure][1]), max(P_out[2], s[:OutletPressure][1]))
    end
    nan_if_none(x) = flowing_time > 0 ? x : NaN
    return (moles = n, backflow = backflow, duration = sum(r.timesteps[window]), flowing_time = flowing_time,
        temperature = weight > 0 ? FT / weight : NaN,
        inlet_pressure = nan_if_none.(P_in), outlet_pressure = nan_if_none.(P_out))
end

# Columns shared by the stream table and the stage stream table, with the
# amounts per `per` (cycle or stage)
function _stream_header(gas, per)
    names = gas.names
    return vcat(["Device", "Role", "From", "To"],
        ["$c [mol%]" for c in names],
        ["Mean flow [kg/h]", "Mean flow [Nm3/h]", "Mean flow [mol/s]"],
        ["Total [mol/$per]", "Total [kg/$per]", "Reverse flow [mol/$per]"],
        ["$c [mol/$per]" for c in names],
        ["$c [kg/$per]" for c in names],
        ["Flowing time [s]", "Mean flow while flowing [Nm3/h]", "Mean temperature [degC]",
         "Min inlet pressure [bara]", "Max inlet pressure [bara]",
         "Min outlet pressure [bara]", "Max outlet pressure [bara]"])
end

function _stream_row(fs, d, st)
    M = fs.gas.molar_masses
    n = st.moles
    mass = n .* M
    F = sum(n) / st.duration
    return vcat(Any[string(d), _role_name(fs, d), _port_label(fs, d, :inlet), _port_label(fs, d, :outlet)],
        _fractions_percent(n),
        [sum(mass) / st.duration * 3600, F * NORMAL_MOLAR_VOLUME * 3600, F],
        [sum(n), sum(mass), st.backflow], n, mass,
        [st.flowing_time, st.flowing_time > 0 ? sum(n) / st.flowing_time * NORMAL_MOLAR_VOLUME * 3600 : NaN,
         _degC(st.temperature)],
        _bara.(st.inlet_pressure), _bara.(st.outlet_pressure))
end

# ------------------------------------------------------------------------------
# Tables
# ------------------------------------------------------------------------------

function _stream_table(r::ProcessResult)
    fs = r.process.flowsheet
    w = _last_cycle(r)
    tol = _flow_tol(r)
    rows = [_stream_row(fs, d, _stream_summary(r, d, w, tol)) for d in _devices(fs)]
    return (header = _stream_header(fs.gas, "cycle"), rows = rows)
end

function _stage_stream_table(r::ProcessResult)
    fs = r.process.flowsheet
    w = _last_cycle(r)
    info = _step_times(r)
    tol = _flow_tol(r)
    rows = Vector{Any}[]
    start = 0.0
    for (j, st) in enumerate(stages(r.process))
        steps = [k for k in w if info[k].stage_number == j]
        for d in _devices(fs)
            isempty(steps) && break
            s = _stream_summary(r, d, steps, tol)
            s.flowing_time > 0 || continue
            push!(rows, vcat(Any[j, st.name, start, st.duration], _stream_row(fs, d, s)))
        end
        start += st.duration
    end
    header = vcat(["Stage number", "Stage", "Stage start [s]", "Stage duration [s]"], _stream_header(fs.gas, "stage"))
    return (header = header, rows = rows)
end

function _stream_history_table(r::ProcessResult)
    fs = r.process.flowsheet
    gas = fs.gas
    info = _step_times(r)
    tol = _flow_tol(r)
    rows = Vector{Any}[]
    for (k, s) in enumerate(r.states), d in _devices(fs)
        F, y, T = _upstream(s[d])
        flowing = abs(F) > tol
        Fmass = F * sum(y .* gas.molar_masses)
        push!(rows, vcat(Any[info[k].cycle, info[k].time, info[k].time_in_cycle, info[k].stage, string(d), _role_name(fs, d),
                F, Fmass * 3600, F * NORMAL_MOLAR_VOLUME * 3600, flowing ? _degC(T) : NaN],
            flowing ? _percent.(y) : fill(NaN, length(y)),
            [_bara(s[d][:InletPressure][1]), _bara(s[d][:OutletPressure][1])]))
    end
    header = vcat(["Cycle", "Time [s]", "Time in cycle [s]", "Stage", "Device", "Role",
            "Flow [mol/s]", "Flow [kg/h]", "Flow [Nm3/h]", "Temperature [degC]"],
        ["$c [mol%]" for c in gas.names], ["Inlet pressure [bara]", "Outlet pressure [bara]"])
    return (header = header, rows = rows)
end

function _bed_header(gas)
    return vcat(["Pressure [bara]", "Temperature [degC]", "Wall temperature [degC]"],
        ["y $c [mol%]" for c in gas.names], ["q $c [mol/kg]" for c in gas.names])
end

# Values of a bed state at `cells`, averaged over them (see
# cross_section_mean), with the loading per kg of adsorbent. The wall
# temperature is that of the cells the wall touches.
function _bed_values(model, s, cells, density)
    m(v) = cross_section_mean(model, v, cells)
    Tw = haskey(s, :WallTemperature) ? _degC(cross_section_mean(model, s[:WallTemperature], intersect(cells, wall_cells(model)))) : NaN
    return vcat(Any[_bara(m(s[:Pressure])), _degC(m(s[:Temperature])), Tw],
        _percent.(m(s[:y])), m(s[:AdsorbedConcentration]) ./ density)
end

_density(p::Process, b) = first(p.parameters[b][:AdsorbentDensity])

function _bed_history_table(r::ProcessResult)
    p = r.process
    fs = p.flowsheet
    info = _step_times(r)
    rows = Vector{Any}[]
    for (k, s) in enumerate(r.states), b in beds(fs), port in keys(ports(fs[b]))
        push!(rows, vcat(Any[info[k].cycle, info[k].time, info[k].time_in_cycle, info[k].stage, string(b), string(port)],
            _bed_values(fs[b], s[b], port_cells(fs[b], port), _density(p, b))))
    end
    header = vcat(["Cycle", "Time [s]", "Time in cycle [s]", "Stage", "Bed", "End"], _bed_header(fs.gas))
    return (header = header, rows = rows)
end

function _bed_profile_table(r::ProcessResult)
    p = r.process
    fs = p.flowsheet
    w = _last_cycle(r)
    t = cumsum(r.timesteps)
    t0 = _start_of_last_cycle(r)
    snapshots = Any[("start of cycle", 0.0, r.state0)]
    stage_end = 0.0
    for st in stages(p)
        stage_end += st.duration
        k = findlast(k -> t[k] <= t0 + stage_end + 1e-6, collect(w))
        isnothing(k) || push!(snapshots, ("end of $(st.name)", stage_end, r.states[w[k]]))
    end
    rows = Vector{Any}[]
    for (label, tc, s) in snapshots, b in beds(fs)
        layers = axial_layers(fs[b])
        x = axial_positions(fs[b])
        L = sum(fs[b].data_domain[:dx][first(l)] for l in layers)
        for (i, cells) in enumerate(layers)
            push!(rows, vcat(Any[label, tc, string(b), i, x[i], x[i] / L], _bed_values(fs[b], s[b], cells, _density(p, b))))
        end
    end
    header = vcat(["Snapshot", "Time in cycle [s]", "Bed", "Cell", "Position [m]", "Relative position [-]"], _bed_header(fs.gas))
    return (header = header, rows = rows)
end

function _summary_table(r::ProcessResult)
    p = r.process
    fs = p.flowsheet
    gas = fs.gas
    names = gas.names
    k = p.options.key_component
    perf = performance(r)
    T = perf.cycle_time
    rows = Vector{Any}[]
    add(section, quantity, value, unit = "") = push!(rows, Any[section, quantity, value, unit])

    add("Process", "Name", p.name)
    add("Process", "Components", join(names, " / "))
    add("Process", "Key component", names[k])
    if isnothing(r.converged)
        add("Process", "Result", "simulation")
        add("Process", "Cycles simulated", r.cycles)
    else
        add("Process", "Result", "cyclic steady state")
        add("Process", "Converged", r.converged)
        add("Process", "Cycles to steady state", r.iterations)
        isempty(r.history) || add("Process", "Final cycle-to-cycle change", last(r.history), "-")
    end
    add("Process", "Mocca version", string(pkgversion(Interface)))
    add("Process", "Written", Libc.strftime("%Y-%m-%d %H:%M", time()))

    add("Feed", "Pressure", _bara(fs.feed.pressure), "bara")
    add("Feed", "Temperature", _degC(fs.feed.temperature), "degC")
    for (c, y) in zip(names, fs.feed.composition)
        add("Feed", "$c", _percent(y), "mol%")
    end
    F_feed = sum(perf.feed) / T
    add("Feed", "Mean feed flow", sum(perf.feed .* gas.molar_masses) / T * 3600, "kg/h")
    add("Feed", "Mean feed flow", F_feed * NORMAL_MOLAR_VOLUME * 3600, "Nm3/h")
    add("Feed", "Mean feed flow", F_feed, "mol/s")

    add("Cycle", "Cycle time", T, "s")
    for (j, st) in enumerate(stages(p))
        add("Cycle", "Stage $j: $(st.name)", st.duration, "s")
    end

    M_k = gas.molar_masses[k]
    add("Performance", "$(names[k]) purity in extract", _percent(perf.purity), "mol%")
    add("Performance", "$(names[k]) recovery", _percent(perf.recovery), "%")
    add("Performance", "$(names[k]) in extract", perf.extract[k] * M_k / T * 3600, "kg/h")
    add("Performance", "Productivity", perf.productivity, "mol/(m3 adsorbent.s)")
    add("Performance", "Productivity", perf.extract[k] / (perf.adsorbent_mass * T), "mol/(kg adsorbent.s)")
    add("Performance", "Specific energy", perf.specific_energy, "kWh/t $(names[k])")
    add("Performance", "Specific energy", perf.specific_energy * 3.6e-3, "GJ/t $(names[k])")
    add("Performance", "Pump and compressor work", perf.energy / 1e3, "kJ/cycle")
    add("Performance", "Mean pump and compressor power", perf.energy / T / 1e3, "kW")

    # Mass balance over the last cycle: what enters and leaves through the
    # boundary devices, and what the beds gain
    w = _last_cycle(r)
    totals(roles) = sum((stream_totals(r.states, r.timesteps, d; window = w) for role in roles for d in devices_with_role(fs, role));
        init = zeros(length(names)))
    n_in = totals(INLET_ROLES)
    n_out = totals(OUTLET_ROLES)
    inventory(s) = sum((component_inventory(fs[b], s[b], p.parameters[b]) for b in beds(fs)); init = zeros(length(names)))
    n_acc = inventory(r.states[last(w)]) .- inventory(r.state0)
    for (c, a, b, acc) in zip(vcat(names, "Total"), vcat(n_in, sum(n_in)), vcat(n_out, sum(n_out)), vcat(n_acc, sum(n_acc)))
        add("Mass balance", "$c in", a, "mol/cycle")
        add("Mass balance", "$c out", b, "mol/cycle")
        add("Mass balance", "$c accumulated in beds", acc, "mol/cycle")
        add("Mass balance", "$c imbalance", iszero(a) ? NaN : _percent((a - b - acc) / a), "% of in")
    end

    add("Adsorbent", "Volume", perf.adsorbent_volume, "m3")
    add("Adsorbent", "Mass", perf.adsorbent_mass, "kg")
    for b in beds(fs)
        spec = get(fs.specs, b, nothing)
        spec isa BedSpec || continue
        isempty(spec.adsorbent.name) || add("Beds", "$b adsorbent", spec.adsorbent.name)
        add("Beds", "$b length", spec.column.length, "m")
        add("Beds", "$b inner diameter", 2 * spec.column.radius, "m")
        add("Beds", "$b bed porosity", spec.column.porosity, "-")
        add("Beds", "$b particle diameter", spec.adsorbent.particle_diameter * 1e3, "mm")
        add("Beds", "$b adsorbent mass", sorbent_mass(fs[b], p.parameters[b]), "kg")
        add("Beds", "$b cells", spec.ncells)
    end
    return (header = ["Section", "Quantity", "Value", "Unit"], rows = rows)
end
