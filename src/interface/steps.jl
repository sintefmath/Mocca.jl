# Cycles in process terms: step templates, an open set that users extend, and
# their compilation into the device settings of each Stage.
#
# A step acts on one bed (or, like Equalise, on several) and opens the devices
# it needs, found by their roles; every other device stays closed. Pressure
# carries over from step to step: a ramp starts from the pressure the bed was
# left at by the step before, around the cycle, unless `from` is given.

"""
    AbstractStep

A step of a cycle, written in process terms. The built-in steps are
[`Pressurise`](@ref), [`Adsorb`](@ref), [`Blowdown`](@ref),
[`Evacuate`](@ref), [`Purge`](@ref), [`Equalise`](@ref) and [`Idle`](@ref).

To add a kind of step, define a subtype with a `duration` field (and usually
`bed` and `name` fields) and a method of [`step_settings`](@ref). Optionally
add methods of [`end_pressure`](@ref) (if the step leaves its bed at a known
pressure), [`step_beds`](@ref) (if it acts on several beds) and
[`default_step_name`](@ref). Register it with [`register_step!`](@ref) to
read and write it in files; files give its fields as keywords. For example, a
rinse of the bed with heavy product at the feed pressure:

```julia
Base.@kwdef mutable struct Rinse <: Mocca.AbstractStep
    duration::Float64
    bed::Union{Nothing, Symbol} = nothing
    name::String = "rinse"
end

function Mocca.step_settings(s::Rinse, ctx::Mocca.StepContext)
    d = Mocca.find_device(ctx, (:feed,))
    return [d => (law = Mocca.valve(ctx, d), inlet = Mocca.boundary(ctx, d, ctx.fs.feed.pressure))]
end
Mocca.end_pressure(s::Rinse, ctx) = ctx.fs.feed.pressure
Mocca.register_step!("rinse", Rinse)
```
"""
abstract type AbstractStep end

"""
    StepContext

What a step sees when its settings are built: the flowsheet `fs`, its `bed`
(resolved: a step that names no bed acts on the flowsheet's only bed), the
bed's `pressure` at the start of the step (or `nothing` if no earlier step
set it), the start pressures of all beds (`pressures`), the default
`ramp_rate` [1/s] and a `location` for error messages. Helpers for step
authors: [`find_device`](@ref), [`valve`](@ref), [`boundary`](@ref),
[`start_pressure`](@ref) and [`ramp`](@ref).
"""
struct StepContext
    fs::Flowsheet
    bed::Union{Nothing, Symbol}
    pressure::Union{Nothing, Float64}
    pressures::Dict{Symbol, Float64}
    ramp_rate::Float64
    location::String
end

"""
    step_settings(step, ctx::StepContext) → Vector{Pair{Symbol, NamedTuple}}

The devices `step` opens and their settings, as `device => (law = ...,
inlet = ..., outlet = ...)`; see [`Stage`](@ref). Devices not listed are
closed.
"""
function step_settings end

"""
    end_pressure(step, ctx) → pressure or nothing

The pressure [Pa] `step` leaves its bed at, for the ramp of the next step, or
`nothing` if it leaves it as it was.
"""
end_pressure(step::AbstractStep, ctx) = nothing

"Pressures [Pa] each bed of `step` is left at, by bed."
function end_pressures(step::AbstractStep, ctx)
    p = end_pressure(step, ctx)
    return (isnothing(p) || isnothing(ctx.bed)) ? Dict{Symbol, Float64}() : Dict(ctx.bed => Float64(p))
end

"""
    step_beds(step) → Vector{Symbol}

The beds `step` acts on; by default its `bed` field, if set.
"""
step_beds(step::AbstractStep) = (hasfield(typeof(step), :bed) && !isnothing(step.bed)) ? [step.bed] : Symbol[]

step_duration(step::AbstractStep) = step.duration

"Name of a step of this kind, when it is not given one."
default_step_name(step::AbstractStep) = lowercase(string(nameof(typeof(step))))

step_name(step::AbstractStep) = (hasfield(typeof(step), :name) && !isempty(step.name)) ? step.name : default_step_name(step)

"""
    find_device(ctx, roles; port = nothing, device = nothing) → Symbol

The device of the step's bed with the first of `roles` that is present,
attached at `port` (`:bottom` or `:top`) if given. With `device`, check that
it is attached to the bed and return it.
"""
function find_device(ctx::StepContext, roles; port = nothing, device = nothing)
    fs, bed = ctx.fs, ctx.bed
    isnothing(bed) && input_error(ctx.location, "bed", "the flowsheet has several beds ($(join(beds(fs), ", "))); say which one the step acts on with bed = ...")
    if !isnothing(device)
        haskey(fs, device) && is_device(fs, device) || input_error(ctx.location, "device", "$device is not a device of the flowsheet")
        any(a -> a.unit.name == bed, attachments(fs, device)) || input_error(ctx.location, "device", "$device is not attached to $bed")
        return device
    end
    ports_ = isnothing(port) ? collect(keys(ports(fs[bed]))) : [port]
    isnothing(port) || haskey(ports(fs[bed]), port) || input_error(ctx.location, "through", "$bed has no port $port; its ports are $(join(keys(ports(fs[bed])), ", "))")
    for role in roles, p in ports_, d in devices_at(fs, bed, p)
        device_role(fs, d) == role && return d
    end
    where_ = isnothing(port) ? "" : " at the $port"
    input_error(ctx.location, "", "$bed has no $(join(string.(roles), " or ")) device$where_; add one with add_device!")
end

"""
    valve(ctx, device; conductance = nothing) → LinearValve

A fully open valve for `device`, with its [`default_conductance`](@ref)
unless `conductance` [m³/(s·Pa)] is given.
"""
valve(ctx::StepContext, d::Symbol; conductance = nothing) = LinearValve(isnothing(conductance) ? default_conductance(ctx.fs, d) : conductance)

"""
    boundary(ctx, device, pressure) → PortCondition

The boundary of `device` with `pressure` (a number or a profile such as
[`ExponentialRamp`](@ref)), and the temperature and composition of its
default boundary.
"""
function boundary(ctx::StepContext, d::Symbol, pressure)
    p = _boundary_port(ctx.fs, d)
    isnothing(p) && input_error(ctx.location, "", "device $d has no boundary port")
    c = boundary_condition(ctx.fs, d, p)
    return PortCondition(pressure = pressure, temperature = c.temperature, composition = c.composition)
end

function _boundary_port(fs::Flowsheet, d::Symbol)
    for p in DEVICE_PORTS
        is_connected(fs, d, p) || return p
    end
    return nothing
end

"""
    start_pressure(ctx, from) → Float64

`from` if given, else the pressure the bed was left at by the step before.
"""
function start_pressure(ctx::StepContext, from)
    isnothing(from) || return Float64(from)
    isnothing(ctx.pressure) && input_error(ctx.location, "from", "the pressure of $(ctx.bed) at the start of this step is unknown; give from = ..., or a step earlier in the cycle that sets it")
    return ctx.pressure
end

"`ExponentialRamp` from `p0` to `p1` at `rate` [1/s], by default the cycle's."
ramp(ctx::StepContext, p0, p1; rate = nothing) = ExponentialRamp(p0, p1, isnothing(rate) ? ctx.ramp_rate : rate)

function _check_rising(ctx, p0, p1)
    p1 >= p0 * (1 - 1e-12) || input_error(ctx.location, "to", "the step should raise the pressure of $(ctx.bed), but goes from $(_bar(p0)) to $(_bar(p1))")
end
function _check_falling(ctx, p0, p1)
    p1 <= p0 * (1 + 1e-12) || input_error(ctx.location, "to", "the step should lower the pressure of $(ctx.bed), but goes from $(_bar(p0)) to $(_bar(p1))")
end
_bar(p) = "$(round(p / 1e5, sigdigits = 4)) bar"

# The feed rate [m³/s] of a step given by one of an interstitial velocity, a
# superficial velocity or a volumetric flow rate
function _flow_rate(ctx, bed_model, velocity, superficial_velocity, flow_rate)
    given = count(!isnothing, (velocity, superficial_velocity, flow_rate))
    given == 1 || input_error(ctx.location, "velocity", "give one of velocity (interstitial), superficial_velocity or flow_rate")
    !isnothing(velocity) && return interstitial_flow_rate(bed_model, velocity)
    !isnothing(superficial_velocity) && return superficial_flow_rate(bed_model, superficial_velocity)
    return flow_rate
end

function _interstitial_velocity(bed_model, rate)
    A = compute_column_face_area(bed_model, NamedTuple(Jutul.setup_parameters(bed_model)))
    return rate / (first(bed_model.data_domain[:porosity]) * A)
end

_other_port(p::Symbol) = p == :bottom ? :top : (p == :top ? :bottom : error("unknown port $p"))

# ------------------------------------------------------------------------------
# Built-in steps
# ------------------------------------------------------------------------------

"""
    Pressurise(duration; to = feed pressure, from = nothing, source = :feed, bed = nothing,
        device = nothing, ramp_rate = nothing, conductance = nothing, name = "")

Raise the bed from `from` (by default where the step before left it) to `to`
[Pa], through the bed's feed device (`source = :feed`) or its light inlet
(`source = :light_product`, light product pressurisation). The pressure on
the far side of the valve follows an exponential ramp at `ramp_rate` [1/s].
"""
Base.@kwdef mutable struct Pressurise <: AbstractStep
    duration::Float64
    to::Union{Nothing, Float64} = nothing
    from::Union{Nothing, Float64} = nothing
    source::Symbol = :feed
    bed::Union{Nothing, Symbol} = nothing
    device::Union{Nothing, Symbol} = nothing
    ramp_rate::Union{Nothing, Float64} = nothing
    conductance::Union{Nothing, Float64} = nothing
    name::String = ""
end
Pressurise(duration; kwarg...) = Pressurise(; duration = duration, kwarg...)
default_step_name(s::Pressurise) = s.source == :light_product ? "LPP" : "pressurisation"

function step_settings(s::Pressurise, ctx::StepContext)
    s.source in (:feed, :light_product) || input_error(ctx.location, "source", "use :feed or :light_product, got $(s.source)")
    d = find_device(ctx, s.source == :feed ? (:feed,) : (:light_inlet,); device = s.device)
    p0 = start_pressure(ctx, s.from)
    p1 = something(s.to, ctx.fs.feed.pressure)
    _check_rising(ctx, p0, p1)
    return [d => (law = valve(ctx, d; conductance = s.conductance), inlet = boundary(ctx, d, ramp(ctx, p0, p1; rate = s.ramp_rate)))]
end
end_pressure(s::Pressurise, ctx) = something(s.to, ctx.fs.feed.pressure)

"""
    Adsorb(duration; velocity = nothing, superficial_velocity = nothing, flow_rate = nothing,
        pressure = feed pressure, bed = nothing, feed_device = nothing, outlet = nothing,
        conductance = nothing, name = "")

Feed the bed at constant flow, given as one of an interstitial `velocity`
[m/s], a `superficial_velocity` [m/s] or a `flow_rate` [m³/s at the bed
inlet], with light product leaving the other end at `pressure` [Pa] through
its product device, or through a link into another bed. A bed without a feed
device of its own is fed through a link opened by the step of the bed before
it, and takes no flow setting.
"""
Base.@kwdef mutable struct Adsorb <: AbstractStep
    duration::Float64
    velocity::Union{Nothing, Float64} = nothing
    superficial_velocity::Union{Nothing, Float64} = nothing
    flow_rate::Union{Nothing, Float64} = nothing
    pressure::Union{Nothing, Float64} = nothing
    bed::Union{Nothing, Symbol} = nothing
    feed_device::Union{Nothing, Symbol} = nothing
    outlet::Union{Nothing, Symbol} = nothing
    conductance::Union{Nothing, Float64} = nothing
    name::String = ""
end
Adsorb(duration; kwarg...) = Adsorb(; duration = duration, kwarg...)
default_step_name(::Adsorb) = "adsorption"

_has_role_at(fs, bed, role) = [(d, p) for p in keys(ports(fs[bed])) for d in devices_at(fs, bed, p) if device_role(fs, d) == role]

function step_settings(s::Adsorb, ctx::StepContext)
    fs, bed = ctx.fs, ctx.bed
    isnothing(bed) && find_device(ctx, (:feed,))  # reports the missing bed
    P = something(s.pressure, fs.feed.pressure)
    out = Pair{Symbol, Any}[]
    feeds = isnothing(s.feed_device) ? _has_role_at(fs, bed, :feed) : [(find_device(ctx, (:feed,); device = s.feed_device), nothing)]
    if !isempty(feeds)
        d = first(feeds[1])
        rate = _flow_rate(ctx, fs[bed], s.velocity, s.superficial_velocity, s.flow_rate)
        push!(out, d => (law = VolumetricFlow(rate), inlet = boundary(ctx, d, P)))
        in_port = only(a.unit.port for a in attachments(fs, d) if a.unit.name == bed)
    else
        any(!isnothing, (s.velocity, s.superficial_velocity, s.flow_rate)) &&
            input_error(ctx.location, "velocity", "$bed has no feed device, so it is fed through a link and takes no flow setting")
        links_in = [a.unit.port for d in devices_with_role(fs, :link) for a in attachments(fs, d) if a.unit.name == bed && a.device.port == :outlet]
        isempty(links_in) && input_error(ctx.location, "", "$bed has neither a feed device nor a link feeding it")
        in_port = links_in[1]
    end
    if isnothing(s.outlet)
        out_port = _other_port(in_port)
        cands = [d for d in devices_at(fs, bed, out_port) if device_role(fs, d) == :product]
        if isempty(cands)
            cands = [d for d in devices_at(fs, bed, out_port) if device_role(fs, d) == :link && any(a -> a.unit.name == bed && a.device.port == :inlet, attachments(fs, d))]
        end
        isempty(cands) && input_error(ctx.location, "outlet", "$bed has no product device or link at its $out_port")
        o = cands[1]
    else
        o = find_device(ctx, (:product, :link); device = s.outlet)
    end
    if device_role(fs, o) == :link
        push!(out, o => (law = valve(ctx, o; conductance = s.conductance),))
    else
        push!(out, o => (law = valve(ctx, o; conductance = s.conductance), outlet = boundary(ctx, o, P)))
    end
    return out
end
end_pressure(s::Adsorb, ctx) = something(s.pressure, ctx.fs.feed.pressure)

"""
    Blowdown(duration; to, from = nothing, through = :top, device = nothing, bed = nothing,
        ramp_rate = nothing, conductance = nothing, name = "")

Lower the bed to `to` [Pa] through an outlet device at its `through` end:
the product device at the top (co-current), or a waste or extract device at
the bottom (counter-current).
"""
Base.@kwdef mutable struct Blowdown <: AbstractStep
    duration::Float64
    to::Float64
    from::Union{Nothing, Float64} = nothing
    through::Symbol = :top
    device::Union{Nothing, Symbol} = nothing
    bed::Union{Nothing, Symbol} = nothing
    ramp_rate::Union{Nothing, Float64} = nothing
    conductance::Union{Nothing, Float64} = nothing
    name::String = ""
end
Blowdown(duration; kwarg...) = Blowdown(; duration = duration, kwarg...)

"""
    Evacuate(duration; to, from = nothing, through = :bottom, device = nothing, bed = nothing,
        ramp_rate = nothing, conductance = nothing, name = "")

Lower the bed to `to` [Pa] through an outlet device at its `through` end,
by default the extract (vacuum) device at the bottom, which takes the heavy
product.
"""
Base.@kwdef mutable struct Evacuate <: AbstractStep
    duration::Float64
    to::Float64
    from::Union{Nothing, Float64} = nothing
    through::Symbol = :bottom
    device::Union{Nothing, Symbol} = nothing
    bed::Union{Nothing, Symbol} = nothing
    ramp_rate::Union{Nothing, Float64} = nothing
    conductance::Union{Nothing, Float64} = nothing
    name::String = ""
end
Evacuate(duration; kwarg...) = Evacuate(; duration = duration, kwarg...)
default_step_name(::Evacuate) = "evacuation"

_depressurise_roles(::Blowdown) = (:product, :waste, :extract)
_depressurise_roles(::Evacuate) = (:extract, :waste, :product)

function step_settings(s::Union{Blowdown, Evacuate}, ctx::StepContext)
    d = find_device(ctx, _depressurise_roles(s); port = s.through, device = s.device)
    p0 = start_pressure(ctx, s.from)
    _check_falling(ctx, p0, s.to)
    return [d => (law = valve(ctx, d; conductance = s.conductance), outlet = boundary(ctx, d, ramp(ctx, p0, s.to; rate = s.ramp_rate)))]
end
end_pressure(s::Union{Blowdown, Evacuate}, ctx) = s.to

"""
    Purge(duration; velocity = nothing, superficial_velocity = nothing, flow_rate = nothing,
        pressure = nothing, bed = nothing, conductance = nothing, name = "")

Light product enters through the bed's light inlet at a constant flow (one of
an interstitial `velocity`, a `superficial_velocity` or a `flow_rate`) and
leaves through the extract device at the other end, held at `pressure` [Pa],
by default where the step before left the bed.
"""
Base.@kwdef mutable struct Purge <: AbstractStep
    duration::Float64
    velocity::Union{Nothing, Float64} = nothing
    superficial_velocity::Union{Nothing, Float64} = nothing
    flow_rate::Union{Nothing, Float64} = nothing
    pressure::Union{Nothing, Float64} = nothing
    bed::Union{Nothing, Symbol} = nothing
    conductance::Union{Nothing, Float64} = nothing
    name::String = ""
end
Purge(duration; kwarg...) = Purge(; duration = duration, kwarg...)

function step_settings(s::Purge, ctx::StepContext)
    fs, bed = ctx.fs, ctx.bed
    d_in = find_device(ctx, (:light_inlet,))
    in_port = only(a.unit.port for a in attachments(fs, d_in) if a.unit.name == bed)
    d_out = find_device(ctx, (:extract, :waste); port = _other_port(in_port))
    P = start_pressure(ctx, s.pressure)
    rate = _flow_rate(ctx, fs[bed], s.velocity, s.superficial_velocity, s.flow_rate)
    in_c = boundary_condition(fs, d_in, :inlet)
    return [d_in => (law = VolumetricFlow(rate), inlet = boundary(ctx, d_in, evaluate_profile(in_c.pressure, 0.0))),
        d_out => (law = valve(ctx, d_out; conductance = s.conductance), outlet = boundary(ctx, d_out, P))]
end
end_pressure(s::Purge, ctx) = isnothing(s.pressure) ? ctx.pressure : s.pressure

"""
    Equalise(duration; beds, device = nothing, conductance = nothing, name = "")

Open the link between the two `beds`, such as `(:A, :B)`, so that their
pressures equalise; gas flows from the higher to the lower. Both beds are then
taken to be at the mean of their start pressures, the start of the next ramp.
"""
Base.@kwdef mutable struct Equalise <: AbstractStep
    duration::Float64
    beds::Vector{Symbol}
    device::Union{Nothing, Symbol} = nothing
    conductance::Union{Nothing, Float64} = nothing
    name::String = ""
end
Equalise(duration; beds, kwarg...) = Equalise(; duration = duration, beds = collect(Symbol, beds), kwarg...)
step_beds(s::Equalise) = s.beds
default_step_name(s::Equalise) = "equalisation $(s.beds[1]) to $(s.beds[2])"

function _equalisation_link(s::Equalise, ctx)
    length(s.beds) == 2 || input_error(ctx.location, "beds", "equalisation joins two beds, got $(s.beds)")
    if !isnothing(s.device)
        return s.device
    end
    for d in devices_with_role(ctx.fs, :link)
        us = Set(a.unit.name for a in attachments(ctx.fs, d))
        us == Set(s.beds) && return d
    end
    input_error(ctx.location, "", "no link joins $(s.beds[1]) and $(s.beds[2]); add one with add_device!(fs, name, :link; from, to)")
end

function step_settings(s::Equalise, ctx::StepContext)
    d = _equalisation_link(s, ctx)
    return [d => (law = valve(ctx, d; conductance = s.conductance),)]
end

function end_pressures(s::Equalise, ctx)
    p = [get(ctx.pressures, b, NaN) for b in s.beds]
    any(isnan, p) && return Dict{Symbol, Float64}()
    m = (p[1] + p[2]) / 2
    return Dict(b => m for b in s.beds)
end

"""
    Idle(duration; bed = nothing, name = "")

Keep the bed closed.
"""
Base.@kwdef mutable struct Idle <: AbstractStep
    duration::Float64
    bed::Union{Nothing, Symbol} = nothing
    name::String = ""
end
Idle(duration; kwarg...) = Idle(; duration = duration, kwarg...)
step_settings(::Idle, ctx) = Pair{Symbol, Any}[]

"""
    Concurrent(steps...; name = "")

Steps on different beds that run at the same time, as one stage of the cycle.
They must have the same duration and open different devices. The stage is
named after its steps, such as `"A adsorption, B evacuation"`.
"""
struct Concurrent
    steps::Vector{Any}
    name::String
    Concurrent(steps::AbstractVector, name::AbstractString) = new(collect(Any, steps), String(name))
end
Concurrent(steps...; name = "") = Concurrent(collect(Any, steps), String(name))

# ------------------------------------------------------------------------------
# Registry, for files
# ------------------------------------------------------------------------------

const STEP_TYPES = Dict{String, Type}()

"""
    register_step!(name, T)

Make the step type `T` readable from and writable to process files under
`"type": name`. Its fields are given as keywords.
"""
register_step!(name::AbstractString, T::Type) = (STEP_TYPES[String(name)] = T; T)

for (n, T) in ("pressurise" => Pressurise, "adsorb" => Adsorb, "blowdown" => Blowdown, "evacuate" => Evacuate,
        "purge" => Purge, "equalise" => Equalise, "idle" => Idle)
    register_step!(n, T)
end

step_type_name(T::Type) = (k = findfirst(==(T), STEP_TYPES); isnothing(k) ? nothing : k)

# ------------------------------------------------------------------------------
# Cycles
# ------------------------------------------------------------------------------

"""
    Cycle(steps)
    Cycle(label => steps, label => steps, ...; cycle_time = nothing)

A cycle of steps (see [`AbstractStep`](@ref)), each run after the one
before. Elements may also be [`Concurrent`](@ref) groups and plain
[`Stage`](@ref)s.

With several labelled sequences, such as one per bed of a two-bed process
whose steps take different times, the sequences run side by side and are
combined with [`parallel_stages`](@ref): a shorter one is padded to
`cycle_time` (by default the longest) with its devices closed. Each sequence
tracks the pressures of its own beds.
"""
struct Cycle
    sequences::Vector{Pair{String, Vector{Any}}}
    cycle_time::Union{Nothing, Float64}
end
Cycle(steps::AbstractVector) = Cycle(["" => collect(Any, steps)], nothing)
Cycle(sequences::Pair{<:AbstractString, <:AbstractVector}...; cycle_time = nothing) =
    Cycle([String(first(s)) => collect(Any, last(s)) for s in sequences], isnothing(cycle_time) ? nothing : Float64(cycle_time))

"All steps of the cycle, in order, with the label of their sequence."
function cycle_steps(c::Cycle)
    out = Tuple{String, Any}[]
    for (label, elems) in c.sequences, e in elems
        if e isa Concurrent
            for s in e.steps
                push!(out, (label, s))
            end
        else
            push!(out, (label, e))
        end
    end
    return out
end

"""
    CompiledCycle

The [`Stage`](@ref)s of a [`Cycle`](@ref) on a flowsheet, with the
pressures each bed passes through and the adsorption velocities of each bed.
"""
struct CompiledCycle
    stages::Vector{Stage}
    pressures::Dict{Symbol, Vector{Float64}}
    velocities::Dict{Symbol, Vector{Float64}}
end

"""
    compile_cycle(cycle, fs; ramp_rate = 0.5) → CompiledCycle

The stages of `cycle` on the flowsheet `fs`. Pressure ramps run at
`ramp_rate` [1/s] unless a step sets its own.
"""
function compile_cycle(cycle::Cycle, fs::Flowsheet; ramp_rate = 0.5)
    pressures = Dict{Symbol, Vector{Float64}}()
    velocities = Dict{Symbol, Vector{Float64}}()
    seqs = Pair{String, Vector{Stage}}[]
    for (label, elems) in cycle.sequences
        isempty(elems) && input_error(isempty(label) ? "cycle" : "$label sequence", "", "has no steps")
        push!(seqs, label => _compile_sequence(elems, fs, Float64(ramp_rate), label, pressures, velocities))
    end
    if length(seqs) == 1
        stages = last(seqs[1])
    else
        kw = isnothing(cycle.cycle_time) ? NamedTuple() : (cycle_time = cycle.cycle_time,)
        stages = _at(() -> parallel_stages(seqs...; kw...), "cycle")
    end
    return CompiledCycle(stages, pressures, velocities)
end

_location(label, i, name) = "step $i ($name)" * (isempty(label) ? " of the cycle" : " of the $label sequence")

function _resolve_bed(fs, step, loc)
    bs = step_beds(step)
    if isempty(bs)
        hasfield(typeof(step), :bed) || return nothing
        all_beds = beds(fs)
        return length(all_beds) == 1 ? all_beds[1] : nothing
    end
    for b in bs
        b in holdup_units(fs) || input_error(loc, "bed", "$b is not a unit with holdup of the flowsheet; beds are $(join(beds(fs), ", "))")
    end
    return bs[1]
end

function _element_steps(e)
    e isa Concurrent && return e.steps
    return Any[e]
end

function _compile_sequence(elems, fs, ramp_rate, label, pressures_seen, velocities)
    # Pressures carry around the cycle, so the first steps start from where
    # the last ones leave the beds: a first pass finds those
    p = Dict{Symbol, Float64}()
    for _ in 1:2, (i, e) in enumerate(elems)
        start = copy(p)
        for s in _element_steps(e)
            s isa AbstractStep || continue
            loc = _location(label, i, step_name(s))
            bed = _resolve_bed(fs, s, loc)
            ctx = StepContext(fs, bed, isnothing(bed) ? nothing : get(start, bed, nothing), start, ramp_rate, loc)
            # A step may fail here because the pressure it starts from is not
            # known yet. The second pass calls it again with the pressures
            # known, so any error that remains is reported there.
            ends = try
                end_pressures(s, ctx)
            catch err
                err isa InterruptException && rethrow()
                Dict{Symbol, Float64}()
            end
            merge!(p, ends)
        end
    end
    stages = Stage[]
    for (i, e) in enumerate(elems)
        if e isa Stage
            push!(stages, e)
            continue
        end
        steps = _element_steps(e)
        e isa Concurrent || e isa AbstractStep || input_error(_location(label, i, "?"), "", "expected a step, Concurrent or Stage, got $(typeof(e))")
        start = copy(p)
        settings = Dict{Symbol, Any}()
        names = String[]
        durations = Float64[]
        multi = e isa Concurrent && length(unique(vcat([step_beds(s) for s in steps if s isa AbstractStep]...))) > 1
        for s in steps
            s isa AbstractStep || input_error(_location(label, i, "?"), "", "a Concurrent group holds steps, got $(typeof(s))")
            name = step_name(s)
            loc = _location(label, i, name)
            step_duration(s) > 0 || input_error(loc, "duration", "must be positive, got $(step_duration(s))")
            bed = _resolve_bed(fs, s, loc)
            ctx = StepContext(fs, bed, isnothing(bed) ? nothing : get(start, bed, nothing), start, ramp_rate, loc)
            for (d, v) in step_settings(s, ctx)
                haskey(settings, d) && input_error(loc, "", "device $d is already set by another step of the same stage")
                settings[d] = v
            end
            ends = end_pressures(s, ctx)
            merge!(p, ends)
            for b in step_beds(s) ∪ (isnothing(bed) ? Symbol[] : [bed])
                haskey(start, b) && push!(get!(pressures_seen, b, Float64[]), start[b])
                haskey(ends, b) && push!(get!(pressures_seen, b, Float64[]), ends[b])
            end
            _record_velocity!(velocities, fs, s, bed, settings)
            push!(names, multi && !isnothing(bed) ? "$bed $name" : name)
            push!(durations, step_duration(s))
        end
        all(d -> abs(d - durations[1]) <= 1e-9 * max(1, durations[1]), durations) ||
            input_error(_location(label, i, join(names, ", ")), "duration", "concurrent steps must have the same duration, got $(join(durations, ", ")) s")
        stage_name = e isa Concurrent && !isempty(e.name) ? e.name : join(names, ", ")
        push!(stages, Stage(stage_name, durations[1], settings))
    end
    return stages
end

_record_velocity!(velocities, fs, s, bed, settings) = nothing
function _record_velocity!(velocities, fs, s::Adsorb, bed, settings)
    isnothing(bed) && return
    for (d, v) in settings
        device_role(fs, d) == :feed && v.law isa VolumetricFlow || continue
        u = isnothing(s.velocity) ? _interstitial_velocity(fs[bed], v.law.rate) : s.velocity
        push!(get!(velocities, bed, Float64[]), u)
    end
end

# ------------------------------------------------------------------------------
# Time stepping
# ------------------------------------------------------------------------------

"""
    TimeStepping(; max_dt = 1.0, first_dt = 0.05, keep_step_pattern = :auto)

How each stage is split into time steps: steps of at most `max_dt` [s],
starting with `first_dt` and doubling, which carries the solver smoothly
through the abrupt start of a stage (and keeps the cycle map smooth for
Newton's method and gradients; see [`setup_schedule`](@ref)).

`keep_step_pattern` keeps the number of steps of each stage fixed when its
duration changes, stretching them instead, which is what the gradients with
respect to durations describe. `:auto` does this when the process has a
duration among its design variables; the pattern is then that of the
durations when the first such variable was added, or when
[`reset_step_pattern!`](@ref) was last called.
"""
struct TimeStepping
    max_dt::Float64
    first_dt::Float64
    keep_step_pattern::Union{Bool, Symbol}
end
function TimeStepping(; max_dt = 1.0, first_dt = min(0.05, max_dt), keep_step_pattern = :auto)
    max_dt > 0 || input_error("time stepping", "max_dt", "must be positive, got $max_dt")
    0 < first_dt <= max_dt || input_error("time stepping", "first_dt", "must be positive and at most max_dt = $max_dt, got $first_dt")
    keep_step_pattern isa Bool || keep_step_pattern == :auto || input_error("time stepping", "keep_step_pattern", "use true, false or :auto")
    return TimeStepping(Float64(max_dt), Float64(first_dt), keep_step_pattern)
end

function Base.show(io::IO, ::MIME"text/plain", c::Cycle)
    steps = cycle_steps(c)
    println(io, "Cycle of $(length(steps)) step$(length(steps) == 1 ? "" : "s")", length(c.sequences) > 1 ? " in $(length(c.sequences)) sequences" : "")
    for (label, s) in steps
        s isa AbstractStep || (println(io, "  ", s isa Stage ? "stage $(s.name), $(s.duration) s" : string(s)); continue)
        print(io, "  ", isempty(label) ? "" : "[$label] ", rpad(step_name(s), 16), lpad(string(round(step_duration(s), sigdigits = 5)), 8), " s")
        bs = step_beds(s)
        isempty(bs) || print(io, "  ", join(bs, ", "))
        println(io)
    end
end
