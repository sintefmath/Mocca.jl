# The Process: one object holding a flowsheet, its cycle, the initial state and
# solver options, with simulation, cyclic steady state and performance.

"""
Options of a [`Process`](@ref) that set how its performance is measured.
"""
struct ProcessOptions
    key_component::Int
    ambient_pressure::Float64
    pump_efficiency::Float64
    heat_capacity_ratio::Float64
    ramp_rate::Float64
end

"Last cyclic steady state and cycle Jacobian, to start the next solve from."
mutable struct WarmStart
    state::Any
    jacobian::Any
end
WarmStart() = WarmStart(nothing, nothing)

"""
    Process(fs, cycle; name = "", timestepping = TimeStepping(), key_component = 1,
        ambient_pressure = 101325.0, pump_efficiency = 0.72, heat_capacity_ratio = 1.4,
        ramp_rate = 0.5, initial_state = Dict(), parameters = Dict())

A cyclic adsorption process: the flowsheet `fs` (built with
[`Flowsheet(gas, feed)`](@ref Flowsheet), [`add_bed!`](@ref) and
[`add_device!`](@ref)) operated by `cycle` (see [`Cycle`](@ref)).

- `timestepping`: see [`TimeStepping`](@ref).
- `key_component`: the component to capture, by index or name, for
  [`performance`](@ref).
- `ambient_pressure` [Pa]: vacuum pumps discharge to it and feed compressors
  take gas from it; `pump_efficiency` and `heat_capacity_ratio` are those of
  both.
- `ramp_rate` [1/s]: of the pressure ramps of steps that do not set their own.
- `initial_state`: per bed, a NamedTuple with any of `pressure`,
  `temperature`, `wall_temperature` and `composition`. By default each bed
  starts at the lowest pressure of its cycle, at the feed temperature, full of
  light product.
- `parameters`: per unit, values of its Jutul parameters to override, such as
  `Dict(:Bed => Dict(:AdsorbentDensity => 1200.0))`.

Simulate it with [`simulate`](@ref), find its cyclic steady state with
[`cyclic_steady_state`](@ref), and measure it with [`performance`](@ref).
"""
mutable struct Process
    name::String
    flowsheet::Flowsheet
    cycle::Cycle
    timestepping::TimeStepping
    options::ProcessOptions
    initial_conditions::Dict{Symbol, Any}
    parameter_settings::Dict{Symbol, Dict{Symbol, Any}}
    variables::Vector{Any}
    model::Jutul.MultiModel
    parameters::Dict{Symbol, Any}
    compiled::Union{Nothing, CompiledCycle}
    reference_durations::Union{Nothing, Vector{Float64}}
    warm::WarmStart
end

function Process(fs::Flowsheet, cycle::Cycle; name = "", timestepping = TimeStepping(), key_component = 1,
        ambient_pressure = 101325.0, pump_efficiency = 0.72, heat_capacity_ratio = 1.4, ramp_rate = 0.5,
        initial_state = Dict{Symbol, Any}(), parameters = Dict{Symbol, Any}())
    _require_materials(fs, "Process")
    k = key_component isa Integer ? Int(key_component) : component_index(fs.gas, key_component)
    1 <= k <= number_of_components(fs.gas) || input_error("process", "key_component", "the gas has $(number_of_components(fs.gas)) components, got $k")
    options = ProcessOptions(k, ambient_pressure, pump_efficiency, heat_capacity_ratio, ramp_rate)
    compiled = compile_cycle(cycle, fs; ramp_rate = ramp_rate)
    _set_dispersion_velocities!(fs, compiled)
    _refresh_default_boundaries!(fs, compiled)
    for b in beds(fs)
        isempty(vcat([devices_at(fs, b, p) for p in keys(ports(fs.units[b]))]...)) && @warn "Bed $b has no devices attached"
    end
    model = setup_flowsheet_model(fs)
    ics = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in pairs(initial_state))
    for b in keys(ics)
        b in beds(fs) || input_error("initial_state", String(b), "$b is not a bed of the flowsheet")
    end
    settings = Dict{Symbol, Dict{Symbol, Any}}(Symbol(u) => Dict{Symbol, Any}(Symbol(n) => v for (n, v) in pairs(d)) for (u, d) in pairs(parameters))
    p = Process(String(name), fs, cycle, timestepping, options, ics, settings, Any[], model, Dict{Symbol, Any}(), compiled, nothing, WarmStart())
    p.parameters = _process_parameters(p)
    return p
end

function _process_parameters(p::Process)
    prm = setup_flowsheet_parameters(p.flowsheet)
    for (u, d) in p.parameter_settings
        haskey(prm, u) || input_error("parameters", String(u), "$u is not a unit of the flowsheet")
        for (n, v) in d
            haskey(prm[u], n) || input_error("parameters of $u", String(n), "unknown parameter; $u has $(join(sort(collect(keys(prm[u]))), ", "))")
            prm[u][n] .= v
        end
    end
    return prm
end

# Beds whose dispersion follows the adsorption velocity: rebuild them for it
function _set_dispersion_velocities!(fs, compiled)
    all_v = vcat(values(compiled.velocities)...)
    for b in beds(fs)
        spec = get(fs.specs, b, nothing)
        spec isa BedSpec && spec.dispersion_velocity == :adsorption && isnothing(spec.dispersion) || continue
        v = haskey(compiled.velocities, b) ? maximum(compiled.velocities[b]) : (isempty(all_v) ? nothing : maximum(all_v))
        if isnothing(v)
            @warn "No adsorption step sets a feed velocity, so the axial dispersion of bed $b uses 1 m/s; pass dispersion_velocity to add_bed! to set it"
            continue
        end
        _set_dispersion_velocity!(fs, b, v)
    end
end

# Default boundaries follow the cycle: extract and waste lines sit at the lowest
# pressure of their bed, the others at the feed pressure
function _refresh_default_boundaries!(fs, compiled)
    for (d, port) in fs.default_boundaries
        role = device_role(fs, d)
        c = fs.boundaries[(d, port)]
        P = fs.feed.pressure
        if role in (:extract, :waste)
            b = first(attachments(fs, d))[1]
            haskey(compiled.pressures, b) && (P = minimum(compiled.pressures[b]))
        end
        fs.boundaries[(d, port)] = PortCondition(pressure = P, temperature = fs.feed.temperature, composition = role_composition(fs, role))
    end
end

"Recompile the cycle after a change to its steps or the feed."
function _refresh!(p::Process)
    compiled = compile_cycle(p.cycle, p.flowsheet; ramp_rate = p.options.ramp_rate)
    _refresh_default_boundaries!(p.flowsheet, compiled)
    p.compiled = compiled
    return p
end

_compiled(p::Process) = isnothing(p.compiled) ? _refresh!(p).compiled : p.compiled

"""
    stages(p::Process) → Vector{Stage}

The stages of the process's cycle, with the device settings of each.
"""
stages(p::Process) = _compiled(p).stages

"Duration [s] of one cycle."
cycle_time(p::Process) = sum(st -> st.duration, stages(p))

"""
    reset_step_pattern!(p)

Take the current stage durations as the pattern of time steps that changed
durations stretch (see [`TimeStepping`](@ref)).
"""
reset_step_pattern!(p::Process) = (p.reference_durations = nothing; p)

function _keeps_step_pattern(p::Process)
    k = p.timestepping.keep_step_pattern
    k isa Bool && return k
    return any(v -> any(pl -> pl isa StepPlace && pl.field == :duration, v.places), p.variables)
end

function _reference_durations(p::Process)
    _keeps_step_pattern(p) || return nothing
    d = [st.duration for st in stages(p)]
    if isnothing(p.reference_durations) || length(p.reference_durations) != length(d)
        p.reference_durations = d
    end
    return p.reference_durations
end

"""
    cycle_schedule(p; cycles = 1) → (forces, timesteps)

Forces and time steps for `cycles` cycles of the process; see
[`setup_schedule`](@ref).
"""
function cycle_schedule(p::Process; cycles = 1)
    ts = p.timestepping
    return setup_schedule(p.flowsheet, p.model, stages(p); num_cycles = cycles, max_dt = ts.max_dt,
        first_dt = ts.first_dt, reference_durations = _reference_durations(p))
end

"""
    initial_state(p) → Dict

The initial state of every unit; see `initial_state` in [`Process`](@ref).
"""
function initial_state(p::Process)
    fs = p.flowsheet
    compiled = _compiled(p)
    bed_states = Dict{Symbol, Any}()
    for b in beds(fs)
        ic = get(p.initial_conditions, b, NamedTuple())
        model = fs.units[b]
        P = haskey(ic, :pressure) ? ic.pressure : (haskey(compiled.pressures, b) ? minimum(compiled.pressures[b]) : fs.feed.pressure)
        T = get(ic, :temperature, fs.feed.temperature)
        spec = get(fs.specs, b, nothing)
        Tw = get(ic, :wall_temperature, spec isa BedSpec ? spec.column.ambient_temperature : T)
        y = haskey(ic, :composition) ? composition(fs.gas, ic.composition; location = "initial state of $b") : fs.feed.light_product
        kw = Dict{Symbol, Any}(:Pressure => P, :Temperature => T, :y => y)
        haskey(model.primary_variables, :WallTemperature) && (kw[:WallTemperature] = Tw)
        bed_states[b] = setup_process_state(model; kw...)
    end
    return setup_flowsheet_state(fs; bed_states...)
end

# ------------------------------------------------------------------------------
# Simulation
# ------------------------------------------------------------------------------

"""
Result of simulating a [`Process`](@ref): the `states` and `timesteps` of
every step the solver took, the number of `cycles` they span, the initial
state `state0` of the last of them, and for a cyclic steady state whether it
`converged`, the change per iteration or cycle in `history` and the cycle
Jacobian if Newton's method found it.
"""
struct ProcessResult
    process::Process
    states::Vector{Any}
    timesteps::Vector{Float64}
    cycles::Int
    state0::Any
    converged::Union{Nothing, Bool}
    history::Vector{Float64}
    iterations::Int
    jacobian::Any
end

"""
    simulate(p::Process; cycles = 1, state0 = initial_state(p), info_level = -1, kwarg...) → ProcessResult

Simulate `cycles` cycles of the process. Extra keywords go to
`simulate_process`.
"""
function Jutul.simulate(p::Process; cycles = 1, state0 = nothing, info_level = -1, kwarg...)
    forces, dt = cycle_schedule(p; cycles = cycles)
    x0 = isnothing(state0) ? initial_state(p) : state0
    case = MoccaCase(p.model, dt, forces; state0 = x0, parameters = p.parameters)
    states, ts = simulate_process(case; output_substates = true, info_level = info_level, kwarg...)
    isapprox(sum(ts), sum(dt)) || error("The simulation stopped at t = $(sum(ts)) s of $(sum(dt)) s")
    n = cycles
    if cycles > 1
        w = cycle_window(ts, cycle_time(p), cycles)
        last0 = restart_state(p.model, states[first(w) - 1])
    else
        last0 = x0
    end
    return ProcessResult(p, collect(Any, states), ts, n, last0, nothing, Float64[], 0, nothing)
end

"""
    cyclic_steady_state(p; method = :newton, tol = nothing, max_cycles = 500,
        warm_start = true, info_level = -1, kwarg...) → ProcessResult

The cyclic steady state of the process, with the states of its last cycle:

- `:newton` (default): Anderson acceleration to a loose steady state, then
  Newton's method on the cycle map (see [`newton_cyclic_steady_state`](@ref))
  to `tol` (default `1e-9`). With `warm_start`, it starts from the last
  steady state found for this process, and its Jacobian, which is much faster
  when the settings have changed a little.
- `:anderson`: [`simulate_to_cyclic_steady_state`](@ref) to `tol` (default
  `1e-6`).
- `:mass_balance`: [`simulate_until_mass_balance`](@ref) on the key component,
  between the feed and light inlet devices and the outlet devices.
"""
function cyclic_steady_state(p::Process; method = :newton, tol = nothing, max_cycles = 500, warm_start = true, info_level = -1, kwarg...)
    forces, dt = cycle_schedule(p)
    model, prm = p.model, p.parameters
    has_warm = warm_start && !isnothing(p.warm.state)
    x0 = has_warm ? restart_state(model, p.warm.state) : initial_state(p)
    if method == :newton
        tol = something(tol, 1e-9)
        p.timestepping.first_dt < p.timestepping.max_dt ||
            @warn "Newton's method needs a smooth cycle map; start each stage with short steps (TimeStepping(first_dt = 0.05))" maxlog = 1
        if has_warm
            r = try
                newton_cyclic_steady_state(model, x0, prm, forces, dt; tol = tol, max_iterations = 40,
                    jacobian = p.warm.jacobian, info_level = info_level, kwarg...)
            catch e
                e isa InterruptException && rethrow()
                nothing
            end
            if !isnothing(r) && r.converged
                p.warm.state, p.warm.jacobian = r.state0, r.jacobian
                return ProcessResult(p, collect(Any, r.states), r.timesteps, 1, r.state0, true, r.history, r.cycles, r.jacobian)
            end
            x0 = initial_state(p)
        end
        loose = simulate_to_cyclic_steady_state(model, x0, prm, forces, dt; tol = 1e-3, max_cycles = max_cycles, info_level = info_level, kwarg...)
        r = newton_cyclic_steady_state(model, loose.state0, prm, forces, dt; tol = tol, max_iterations = 10, info_level = info_level, kwarg...)
        r.converged || @warn "Cyclic steady state not converged by Newton's method" r.history
        p.warm.state, p.warm.jacobian = r.state0, r.jacobian
        return ProcessResult(p, collect(Any, r.states), r.timesteps, 1, r.state0, r.converged, vcat(loose.history, r.history),
            loose.cycles + r.cycles, r.jacobian)
    elseif method == :anderson
        r = simulate_to_cyclic_steady_state(model, x0, prm, forces, dt; tol = something(tol, 1e-6), max_cycles = max_cycles, info_level = info_level, kwarg...)
        r.converged && (p.warm.state = r.state0)
        return ProcessResult(p, collect(Any, r.states), r.timesteps, 1, r.state0, r.converged, r.history, r.cycles, nothing)
    elseif method == :mass_balance
        fs = p.flowsheet
        inflow = Tuple(vcat(devices_with_role(fs, :feed), devices_with_role(fs, :light_inlet)))
        outflow = Tuple(vcat([devices_with_role(fs, r) for r in OUTLET_ROLES]...))
        kw = isnothing(tol) ? NamedTuple() : (tol = tol,)
        r = simulate_until_mass_balance(model, x0, prm, forces, dt; inflow = inflow, outflow = outflow,
            component = p.options.key_component, max_cycles = max_cycles, info_level = info_level, kw..., kwarg...)
        return ProcessResult(p, collect(Any, r.states), r.timesteps, 1, r.state0, r.converged, r.history, r.cycles, nothing)
    else
        input_error("cyclic_steady_state", "method", "use :newton, :anderson or :mass_balance, got $method")
    end
end

# ------------------------------------------------------------------------------
# Performance
# ------------------------------------------------------------------------------

_last_cycle(r::ProcessResult) = r.cycles == 1 ? eachindex(r.timesteps) : cycle_window(r.timesteps, cycle_time(r.process), r.cycles)

adsorbent_volume(p::Process) = sum(b -> sum(p.parameters[b][:SolidVolume]), beds(p.flowsheet))
adsorbent_mass(p::Process) = sum(b -> sorbent_mass(p.flowsheet.units[b], p.parameters[b]), beds(p.flowsheet))

"""
    performance(r::ProcessResult) → NamedTuple

Performance of the process over the last cycle of `r`, from the streams
through its devices (by their roles; see [`add_device!`](@ref)):

- `purity`: mole fraction of the key component in the extract;
- `recovery`: fraction of the key component fed that leaves in the extract;
- `specific_energy` [kWh per tonne of key component in the extract]: vacuum
  pumps on every outlet device whose gas leaves below ambient pressure, and
  compressors on every inlet device whose gas enters a bed above it, as in
  Haghpanah et al. (2013);
- `productivity` [mol of key component per m³ of adsorbent per s];
- `energy` [J per cycle], `cycle_time` [s], the component totals [mol per
  cycle] of the `feed`, `extract` and light `product`, and the
  `adsorbent_volume` [m³] and `adsorbent_mass` [kg].
"""
function performance(r::ProcessResult)
    p = r.process
    fs, o = p.flowsheet, p.options
    w = _last_cycle(r)
    states, ts = r.states, r.timesteps
    total(role) = (ds = devices_with_role(fs, role); N = number_of_components(fs.gas);
        isempty(ds) ? zeros(N) : sum(d -> stream_totals(states, ts, d; window = w), ds))
    n_feed, n_extract, n_product = total(:feed), total(:extract), total(:product)
    E = 0.0
    for role in OUTLET_ROLES, d in devices_with_role(fs, role)
        E += vacuum_pump_energy(states, ts, d; discharge_pressure = o.ambient_pressure, efficiency = o.pump_efficiency, γ = o.heat_capacity_ratio, window = w)
    end
    for role in INLET_ROLES, d in devices_with_role(fs, role)
        E += compressor_energy(states, ts, d; suction_pressure = o.ambient_pressure, efficiency = o.pump_efficiency, γ = o.heat_capacity_ratio, window = w)
    end
    k = o.key_component
    T = sum(ts[w])
    V = adsorbent_volume(p)
    return (purity = purity(n_extract; component = k), recovery = recovery(n_extract, n_feed; component = k),
        specific_energy = specific_energy_kwh_per_tonne(E, n_extract, fs.gas.molar_masses[k]; component = k),
        productivity = n_extract[k] / (V * T), energy = E, cycle_time = T,
        feed = n_feed, extract = n_extract, product = n_product, adsorbent_volume = V, adsorbent_mass = adsorbent_mass(p))
end

"""
    performance(p::Process; kwarg...) → NamedTuple

Performance at cyclic steady state; keywords go to
[`cyclic_steady_state`](@ref).
"""
performance(p::Process; kwarg...) = performance(cyclic_steady_state(p; kwarg...))

# ------------------------------------------------------------------------------
# Inspection
# ------------------------------------------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", p::Process)
    fs = p.flowsheet
    println(io, "Process", isempty(p.name) ? "" : " $(p.name)", ": ", join(fs.gas.names, "/"), " feed at ",
        _bar(fs.feed.pressure), ", ", fs.feed.temperature, " K, capturing ", fs.gas.names[p.options.key_component])
    print_flowsheet_diagram(io, fs)
    sts = try
        stages(p)
    catch e
        println(io, "  cycle: ", sprint(showerror, e))
        return
    end
    println(io, "  cycle of $(length(sts)) stages, $(round(sum(s -> s.duration, sts), sigdigits = 5)) s:")
    for st in sts
        open_ = [string(d) for d in fs.order if haskey(st.settings, d)]
        println(io, "    ", rpad(st.name, 34), lpad(string(round(st.duration, sigdigits = 5)), 8), " s   ", isempty(open_) ? "all closed" : join(open_, ", "))
    end
    if !isempty(p.variables)
        print(io, "  variables: ", join([string(v.name) for v in p.variables], ", "))
    end
end
