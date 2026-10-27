# Processes that run at steady state, such as moving bed temperature swing
# adsorption: the devices keep their settings, and the process is run until
# nothing changes any more, rather than through a cycle of steps.

"""
    SteadyProcess(fs, settings; name = "", duration = 3600.0,
        timestepping = TimeStepping(max_dt = 60.0, first_dt = 0.05), key_component = 1,
        ambient_pressure = 101325.0, pump_efficiency = 0.72, heat_capacity_ratio = nothing,
        ambient_temperature = 298.15, initial_state = Dict(), parameters = Dict())

A process on flowsheet `fs` whose devices keep the same settings, run until it
reaches steady state, such as a moving bed process (see
[`add_moving_bed!`](@ref)). `settings` gives each open device its flow law, or
a `NamedTuple` with the `law` and the [`PortCondition`](@ref) of a boundary
port (`inlet` or `outlet`); devices not given are closed:

```julia
SteadyProcess(fs, Dict(
    :V_feed => FixedMolarFlow(1878.0),
    :V_product => LinearValve(default_conductance(fs, :V_product)),
    :V_recycle => BackPressureController(101325.0, 1e-2)))
```

[`steady_state`](@ref) runs it, starting from `initial_state` and taking time
steps that grow from `timestepping.first_dt` to `timestepping.max_dt`, in
spans of `duration` until the state stops changing; [`simulate`](@ref) runs
one span. [`performance`](@ref) measures the steady state: devices with the
`:feed` role give the feed and those with the `:extract` role the product of
`key_component`.

`initial_state` gives units their starting `pressure`, `temperature` and
`composition` (by default those of the feed, with the light product), as a
`NamedTuple` by unit name; the sorbent starts at equilibrium with the gas, and
tubes and their fluid at the temperature of the unit. `parameters` overrides
unit parameters, as for a [`Process`](@ref).
"""
mutable struct SteadyProcess <: Process
    name::String
    flowsheet::Flowsheet
    settings::Dict{Symbol, Any}
    duration::Float64
    timestepping::TimeStepping
    options::ProcessOptions
    initial_conditions::Dict{Symbol, Any}
    parameter_settings::Dict{Symbol, Dict{Symbol, Any}}
    model::Jutul.MultiModel
    parameters::Dict{Symbol, Any}
    warm::WarmStart
    variables::Vector{Any}
    integrals::Vector{Pair{Symbol, Any}}
    measures::Vector{Pair{Symbol, Any}}
end

function SteadyProcess(fs::Flowsheet, settings; name = "", duration = 3600.0, timestepping = TimeStepping(max_dt = 60.0, first_dt = 0.05),
        key_component = 1, ambient_pressure = 101325.0, pump_efficiency = 0.72, heat_capacity_ratio = nothing,
        ambient_temperature = 298.15, initial_state = Dict{Symbol, Any}(), parameters = Dict{Symbol, Any}())
    fs = deepcopy(fs)
    duration > 0 || input_error("steady process", "duration", "must be positive, got $duration")
    k = key_component isa Integer ? Int(key_component) : component_index(fs.gas, key_component)
    1 <= k <= number_of_components(fs.gas) || input_error("steady process", "key_component", "the gas has $(number_of_components(fs.gas)) components, got $k")
    fs.feed.light_product_given || (fs.feed.light_product = light_product_composition(fs.feed.composition; component = k))
    isnothing(heat_capacity_ratio) || heat_capacity_ratio > 1 ||
        input_error("steady process", "heat_capacity_ratio", "must be greater than 1, got $heat_capacity_ratio")
    options = ProcessOptions(k, ambient_pressure, pump_efficiency, _maybe_float(heat_capacity_ratio), 0.5, ambient_temperature)
    s = Dict{Symbol, Any}()
    for (d, v) in pairs(settings)
        d = Symbol(d)
        haskey(fs, d) && is_device(fs, d) || input_error("settings", String(d), "$d is not a device of the flowsheet")
        s[d] = _device_setting(fs, d, v)
    end
    for u in holdup_units(fs)
        isempty(vcat([devices_at(fs, u, p) for p in keys(ports(fs[u]))]...)) && @warn "Unit $u has no devices attached"
    end
    model = setup_flowsheet_model(fs)
    ics = Dict{Symbol, Any}(Symbol(u) => v for (u, v) in pairs(initial_state))
    for u in keys(ics)
        u in holdup_units(fs) || input_error("initial_state", String(u), "$u is not a unit with holdup of the flowsheet")
    end
    prm_settings = Dict{Symbol, Dict{Symbol, Any}}(Symbol(u) => Dict{Symbol, Any}(Symbol(n) => v for (n, v) in pairs(d)) for (u, d) in pairs(parameters))
    p = SteadyProcess(String(name), fs, s, Float64(duration), timestepping, options, ics, prm_settings, model, Dict{Symbol, Any}(),
        WarmStart(), Any[], Pair{Symbol, Any}[], Pair{Symbol, Any}[])
    p.parameters = _steady_parameters(p)
    return p
end

# A device's settings as a NamedTuple with its law and any boundary conditions
function _device_setting(fs, d, v)
    loc = "settings of $d"
    v isa AbstractDeviceLaw && return (law = v,)
    v isa NamedTuple || input_error(loc, "", "expected a device law or a NamedTuple with law, inlet and outlet; got $(typeof(v))")
    for k in keys(v)
        k in (:law, :inlet, :outlet) || input_error(loc, String(k), "unknown setting; use law, inlet and outlet")
    end
    haskey(v, :law) && !(v.law isa AbstractDeviceLaw) && input_error(loc, "law", "expected a device law, got $(typeof(v.law))")
    for p in (:inlet, :outlet)
        haskey(v, p) || continue
        is_connected(fs, d, p) && input_error(loc, String(p), "the $p of $d is connected to a unit, so it has no boundary condition")
        v[p] isa PortCondition || input_error(loc, String(p), "expected a PortCondition")
    end
    return v
end

# The model and parameters of the process again, after a unit of its
# flowsheet was rebuilt
function _rebuild!(p::SteadyProcess)
    p.model = setup_flowsheet_model(p.flowsheet)
    p.parameters = _steady_parameters(p)
    return p
end

function _steady_parameters(p::SteadyProcess)
    fs = p.flowsheet
    prm = setup_flowsheet_parameters(fs)
    for b in beds(fs)
        get(fs.specs, b, nothing) isa BedSpec && !_dispersion_follows_state(fs, b) && (prm[b][:AxialDispersion] .= bed_dispersion(fs, b))
    end
    for (u, d) in p.parameter_settings
        haskey(prm, u) || input_error("parameters", String(u), "$u is not a unit of the flowsheet")
        for (n, v) in d
            haskey(prm[u], n) || input_error("parameters of $u", String(n), "unknown parameter; $u has $(join(sort(collect(keys(prm[u]))), ", "))")
            prm[u][n] .= v
        end
    end
    return prm
end

"The device settings of steady process `p`, as one [`Stage`](@ref) lasting `duration`."
operating_stage(p::SteadyProcess; duration = p.duration) = Stage(isempty(p.name) ? "steady operation" : p.name, duration, p.settings...)

"""
    stages(p::SteadyProcess) → Vector{Stage}

The one stage of a steady process, its [`operating_stage`](@ref).
"""
stages(p::SteadyProcess) = [operating_stage(p)]

# A steady process has no steps whose time steps it keeps
reset_step_pattern!(p::SteadyProcess) = p

# The settings of a device, and what a setting place of it holds: a field of
# its flow law, or the condition at a boundary port, its own or the
# flowsheet's
function _device_setting_of(p::SteadyProcess, x::SettingPlace)
    s = get(p.settings, x.device, nothing)
    isnothing(s) && input_error("variable", "places", "$(x.device) is not an open device of the process")
    if isnothing(x.port)
        haskey(s, :law) || input_error("variable", "places", "$(x.device) has no flow law set")
        hasfield(typeof(s.law), x.field) ||
            input_error("variable", "places", "a $(nameof(typeof(s.law))) has no field $(x.field); it has $(join(fieldnames(typeof(s.law)), ", "))")
        return (settings = s, value = getfield(s.law, x.field))
    end
    x.port in (:inlet, :outlet) || input_error("variable", "places", "a device port is :inlet or :outlet, got $(x.port)")
    is_connected(p.flowsheet, x.device, x.port) && input_error("variable", "places", "the $(x.port) of $(x.device) is connected to a unit, so it has no boundary condition")
    x.field in (:pressure, :temperature) || input_error("variable", "places", "a port condition sets :pressure or :temperature, got $(x.field)")
    c = haskey(s, x.port) ? s[x.port] : boundary_condition(p.flowsheet, x.device, x.port)
    v = getfield(c, x.field)
    v isa Real || input_error("variable", "places", "the $(x.field) at the $(x.port) of $(x.device) is $(repr(v)), not a number")
    return (settings = s, value = v, condition = c)
end

_get(p::SteadyProcess, x::SettingPlace) = _device_setting_of(p, x).value

function _set!(p::SteadyProcess, x::SettingPlace, v)
    d = _device_setting_of(p, x)
    s = d.settings
    if isnothing(x.port)
        law = s.law
        values = [f == x.field ? v : getfield(law, f) for f in fieldnames(typeof(law))]
        p.settings[x.device] = merge(s, (law = Base.typename(typeof(law)).wrapper(promote(values...)...),))
    else
        c = d.condition
        c = PortCondition(pressure = x.field == :pressure ? Float64(v) : c.pressure,
            temperature = x.field == :temperature ? Float64(v) : c.temperature, composition = c.composition)
        p.settings[x.device] = merge(s, NamedTuple{(x.port,)}((c,)))
    end
end

function initial_state(p::SteadyProcess)
    fs = p.flowsheet
    unit_states = Dict{Symbol, Any}()
    for u in holdup_units(fs)
        ic = get(p.initial_conditions, u, NamedTuple())
        P = get(ic, :pressure, fs.feed.pressure)
        T = get(ic, :temperature, fs.feed.temperature)
        y = haskey(ic, :composition) ? composition(fs.gas, ic.composition; location = "initial state of $u") : fs.feed.light_product
        unit_states[u] = setup_unit_state(fs[u]; pressure = P, temperature = T, composition = y,
            wall_temperature = get(ic, :wall_temperature, T))
    end
    return setup_flowsheet_state(fs; unit_states...)
end

# A steady state is the fixed point of one long backward-Euler step: a state
# x that the step from x returns satisfies (M(x) − M(x))/Δt + F(x) = 0, so
# F(x) = 0 whatever Δt. Its gradients are those of a cycle of that one step
# (see cyclic_steady_state_gradient), exact for any Δt. The step is not
# solved: it returns the steady state, so the adjoints linearise it there. (On
# fine grids the solver cannot take so long a step from the steady state at
# its tolerance.) A long step makes the step map ∂Φ/∂x nearly zero, so that
# the adjoint solves for the steady state converge in an iteration or two.
# The integrals over the step are Δt times the rates.
const STEADY_GRADIENT_STEP = 1e6

function _fixed_point(p::SteadyProcess, r)
    Δt = STEADY_GRADIENT_STEP
    forces, dt = setup_schedule(p.flowsheet, p.model, [operating_stage(p; duration = Δt)]; max_dt = Δt)
    x = restart_state(p.model, r.states[end])
    return (forces = forces, timesteps = dt, state0 = x, jacobian = nothing, result = r,
        steps = ([r.states[end]], dt, forces), measured = (states = r.states[end:end], timesteps = [1.0]), scale = 1 / Δt)
end

"""
Result of running a [`SteadyProcess`](@ref): the `states` and `timesteps` of
every step the solver took, and for [`steady_state`](@ref) whether it
`converged` and the change of the state over each span in `history`.
"""
struct SteadyResult <: ProcessResult
    process::SteadyProcess
    states::Vector{Any}
    timesteps::Vector{Float64}
    converged::Union{Nothing, Bool}
    history::Vector{Float64}
end

# Steps of a span: growing from first_dt on a cold start, otherwise of max_dt
function _span_schedule(p::SteadyProcess, duration; cold = true)
    ts = p.timestepping
    first_dt = cold ? ts.first_dt : ts.max_dt
    return setup_schedule(p.flowsheet, p.model, [operating_stage(p; duration = duration)]; max_dt = ts.max_dt, first_dt = first_dt)
end

"""
    simulate(p::SteadyProcess; duration = p.duration, state0 = initial_state(p), info_level = -1, kwarg...) → SteadyResult

Run the process for `duration` [s]. Extra keywords go to `simulate_process`.
"""
function Jutul.simulate(p::SteadyProcess; duration = p.duration, state0 = nothing, info_level = -1, kwarg...)
    x0 = isnothing(state0) ? initial_state(p) : state0
    forces, dt = _span_schedule(p, duration; cold = isnothing(state0))
    states, ts, ok = _simulate_schedule(p.model, x0, p.parameters, forces, dt; info_level = info_level, kwarg...)
    ok || throw(SimulationFailure("the simulation stopped at t = $(sum(ts)) s of $(sum(dt)) s"))
    return SteadyResult(p, collect(Any, states), ts, nothing, Float64[])
end

"""
    steady_state(p::SteadyProcess; tol = 1e-6, max_spans = 10, state0 = nothing,
        warm_start = true, progress = false, info_level = -1, kwarg...) → SteadyResult

Run the process in spans of `p.duration` until the largest relative change of
any unit's state over a span (see [`cycle_change`](@ref)) is below `tol`, for
at most `max_spans` spans. It starts from `state0` if given, otherwise from
the last steady state the process found if `warm_start`, which suits a
process whose settings change a little at a time, as in an optimisation, and
otherwise from [`initial_state`](@ref). Far from the last design, or near a
change of regime, a warm start may settle into another regime than a cold
one would. The result holds the states of every span; check its `converged`.
With `progress`, each span is printed with its change.
"""
function steady_state(p::SteadyProcess; tol = 1e-6, max_spans = 10, state0 = nothing, warm_start = true, progress = false,
        info_level = -1, kwarg...)
    isnothing(state0) && warm_start && !isnothing(p.warm.state) && (state0 = p.warm.state)
    x = isnothing(state0) ? initial_state(p) : state0
    cold = isnothing(state0)
    states, ts, history = Any[], Float64[], Float64[]
    converged = false
    for _ in 1:max_spans
        forces, dt = _span_schedule(p, p.duration; cold = cold)
        span, span_ts, ok = _simulate_schedule(p.model, x, p.parameters, forces, dt; info_level = info_level, kwarg...)
        ok || throw(SimulationFailure("the simulation stopped at t = $(sum(ts) + sum(span_ts)) s"))
        append!(states, span)
        append!(ts, span_ts)
        x_next = restart_state(p.model, span[end])
        push!(history, cycle_change(p.model, x, x_next))
        progress === false || println(progress === true ? stdout : progress, "  span $(length(history)): change $(_sci(history[end])) → $(_sci(tol))")
        x = x_next
        cold = false
        history[end] < tol && (converged = true; break)
    end
    converged && (p.warm.state = x)
    return SteadyResult(p, states, ts, converged, history)
end

"""
    performance(r::SteadyResult) → NamedTuple

Performance at the end of `r`, which should be steady, as rates. Its measures
(see [`performance_measures`](@ref)) are functions of integrals over one
second (see [`performance_integrals`](@ref)):

- `purity` and `recovery` of the key component in the streams of the
  `:extract` devices, against those of the `:feed` devices;
- `specific_energy` [kWh/t of key component], from the `power`;
- `productivity` [mol of key component per m³ of adsorbent per second];
- with heat exchangers, `specific_heat_duty` [J/kg of key component], the
  heating per key component captured, and `specific_equivalent_work` [kWh/t],
  the power and the `equivalent_work` of the heating together;
- any measure added with [`add_measure!`](@ref);
- `feed`, `extract` and `product` (the `:product` devices): molar flow of each
  component [mol/s];
- `power` [W]: vacuum pumps on the outlets and compressors on the inlets, as
  for a [`Process`](@ref), and the fans and blowers of links that move gas up
  in pressure, `blower_power`, such as recycle blowers;
- `heat_duty` [W] of each unit with a heat exchanger, heat given to the
  sorbent; `heating` and `cooling` [W], the heat given and taken by fluid
  coming from outside the process, hotter than the ambient temperature or not
  (units whose fluid comes from another unit only pass heat around); and
  `equivalent_work` [W], the work equivalent to the heating, weighted by the
  Carnot factor of the temperatures it is given at (see
  [`EquivalentWorkObjective`](@ref));
- `sorbent_inventory` [kg] and `circulation_time` [s] of the moving sorbent,
  and `adsorbent_volume` [m³] of every bed.
"""
function performance(r::SteadyResult)
    p = r.process
    fs = p.flowsheet
    state = r.states[end]
    prm = p.parameters
    N = number_of_components(fs.gas)
    rate(role) = (ds = devices_with_role(fs, role); isempty(ds) ? zeros(N) : sum(d -> stream_rates(state, d), ds))
    s = performance_sums(p, [state], [1.0])
    duties = Dict{Symbol, Float64}(u => heat_duty(fs, state, prm, u) for u in holdup_units(fs) if has_heat_exchanger(fs[u].system))
    moving = !isempty(moving_sorbent_units(fs))
    return merge(_measure_values(p, s), (feed = rate(:feed), extract = rate(:extract), product = rate(:product),
        power = s.energy, blower_power = s.blower_energy, heat_duty = duties, heating = s.heating, cooling = s.cooling,
        equivalent_work = s.equivalent_work, sorbent_inventory = s.sorbent_inventory,
        circulation_time = moving ? sorbent_circulation_time(fs, prm) : NaN, adsorbent_volume = s.adsorbent_volume))
end


function Base.show(io::IO, ::MIME"text/plain", p::SteadyProcess)
    fs = p.flowsheet
    println(io, "Steady process", isempty(p.name) ? "" : " $(p.name)", ": ", join(fs.gas.names, "/"), " feed at ",
        _bar(fs.feed.pressure), ", ", fs.feed.temperature, " K, capturing ", fs.gas.names[p.options.key_component])
    print_flowsheet_diagram(io, fs)
    open_ = [string(d) for d in unit_names(fs) if haskey(p.settings, d)]
    print(io, "  open devices: ", isempty(open_) ? "none" : join(open_, ", "))
end

function Base.show(io::IO, ::MIME"text/plain", r::SteadyResult)
    print(io, "SteadyResult: $(length(r.timesteps)) steps over $(round(sum(r.timesteps), sigdigits = 5)) s")
    isnothing(r.converged) || print(io, r.converged ? ", converged" : ", not converged",
        isempty(r.history) ? "" : " (change $(round(r.history[end], sigdigits = 3)))")
end
