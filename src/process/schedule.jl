"""
    Stage(name, duration; settings...)

One step of an operating schedule. Each keyword names a device and gives its
settings for the stage as a `NamedTuple` with any of `law`, `inlet` and
`outlet` (see `Jutul.setup_forces` for a `FlowDeviceModel`). Devices not named
are closed. Time-dependent profiles in port conditions use time since the start
of the stage.

```julia
Stage("adsorption", 15.0;
    V_feed = (law = VolumetricFlow(v_feed * A),),
    V_product = (law = LinearValve(C),))
```

Settings can also be passed as `device => settings` pairs, which is convenient
when they are built programmatically: `Stage("adsorption", 15.0, :V_feed => (law = Closed(),))`.
"""
struct Stage
    name::String
    duration::Float64
    settings::Dict{Symbol, Any}
end
Stage(name, duration; settings...) = Stage(String(name), Float64(duration), Dict{Symbol, Any}(settings))
Stage(name, duration, settings::Pair{Symbol}...) = Stage(String(name), Float64(duration), Dict{Symbol, Any}(settings...))

"""
    setup_schedule(fs, model, stages; num_cycles = 1, max_dt = 1.0, first_dt = max_dt) → (forces, timesteps)

Forces and time steps for `num_cycles` repetitions of `stages` on the flowsheet
`fs` and its `MultiModel`. Each stage is split into steps no longer than
`max_dt`. With `first_dt < max_dt`, each stage starts with a step of
`first_dt`, doubling up to `max_dt`, which helps the solver through the abrupt
change when a valve opens or a flow controller starts. Boundary ports use the stage's port condition if given, otherwise the
flowsheet default from [`set_boundary!`](@ref).
"""
function setup_schedule(fs::Flowsheet, model::Jutul.MultiModel, stages::AbstractVector{Stage};
        num_cycles = 1, max_dt = 1.0, first_dt = max_dt)
    for st in stages, k in keys(st.settings)
        haskey(fs.units, k) && _is_device(fs.units[k]) || error("Stage $(st.name) sets $k, which is not a device in the flowsheet")
    end
    cycle_time = sum(st -> st.duration, stages)
    forces = Any[]
    timesteps = Float64[]
    t = 0.0
    for cycle in 1:num_cycles
        for st in stages
            f = _stage_forces(fs, model, st, t)
            dt = _stage_steps(st.duration, max_dt, first_dt)
            append!(timesteps, dt)
            append!(forces, fill(f, length(dt)))
            t += st.duration
        end
    end
    @assert isapprox(t, num_cycles * cycle_time)
    return (forces, timesteps)
end

# Steps for one stage: from `first_dt`, doubling up to `max_dt`, then equal
# steps no longer than `max_dt` to the end of the stage
function _stage_steps(duration, max_dt, first_dt)
    steps = Float64[]
    t = 0.0
    h = min(first_dt, max_dt)
    while h < max_dt && t + h < duration - 1e-9
        push!(steps, h)
        t += h
        h = min(2h, max_dt)
    end
    rest = duration - t
    n = max(1, ceil(Int, rest / max_dt - 1e-9))
    append!(steps, fill(rest / n, n))
    return steps
end

function _stage_forces(fs::Flowsheet, model::Jutul.MultiModel, st::Stage, t_start)
    per_unit = Dict{Symbol, Any}()
    for name in fs.order
        unit = fs.units[name]
        _is_device(unit) || continue
        s = get(st.settings, name, NamedTuple())
        law = get(s, :law, Closed())
        conds = map(DEVICE_PORTS) do p
            if is_connected(fs, name, p)
                haskey(s, p) && error("Stage $(st.name) sets a condition on port $p of $name, which is connected to a unit")
                nothing
            else
                shift_profile(get(s, p, fs.boundaries[(name, p)]), t_start)
            end
        end
        per_unit[name] = Jutul.setup_forces(unit; law = law, inlet = conds[1], outlet = conds[2])
    end
    return Jutul.setup_forces(model; per_unit...)
end

"""
    parallel_stages(label => stages, ...; cycle_time = longest sequence) → Vector{Stage}

Combine stage sequences that run at the same time on separate devices into one
sequence of stages. An example is two beds in series whose blowdown and
evacuation times differ:

```julia
parallel_stages("silica gel" => sg_stages, "13X" => zeolite_stages)
```

Each sequence is padded to `cycle_time` with an `idle` stage in which its
devices are closed. The combined sequence has a stage boundary wherever any
sequence has one, and each of its stages holds the settings of every
sequence's current stage. It is named after them, prefixed by their labels,
such as `"silica gel blowdown, 13X blowdown"`. Where a stage is split, its
pressure profiles keep running from the start of the original stage.

Each device may appear in only one sequence.
"""
function parallel_stages(sequences::Pair{<:AbstractString, <:AbstractVector{Stage}}...;
        cycle_time = maximum(seq -> sum(st -> st.duration, last(seq)), sequences))
    tol = 1e-9
    owners = Dict{Symbol, String}()
    for (label, seq) in sequences, st in seq, d in keys(st.settings)
        get!(owners, d, label) == label || error("Device $d appears in the sequences for both $(owners[d]) and $label")
    end
    padded = map(sequences) do (label, seq)
        rest = cycle_time - sum(st -> st.duration, seq)
        rest > -tol || error("The $label sequence lasts longer than the cycle time $cycle_time s")
        label => (rest > tol ? vcat(seq, Stage("idle", rest)) : collect(seq))
    end
    starts = [cumsum(vcat(0.0, [st.duration for st in seq])) for (_, seq) in padded]
    breaks = sort(vcat(starts...))
    breaks = breaks[vcat(true, diff(breaks) .> tol)]
    stages = Stage[]
    for (t0, t1) in zip(breaks[1:end-1], breaks[2:end])
        settings = Dict{Symbol, Any}()
        names = String[]
        for ((label, seq), s) in zip(padded, starts)
            k = findlast(<=(t0 + tol), s[1:end-1])
            st = seq[k]
            elapsed = t0 - s[k]
            for (d, v) in st.settings
                settings[d] = _shift_setting(v, -elapsed)
            end
            push!(names, "$label $(st.name)")
        end
        push!(stages, Stage(join(names, ", "), t1 - t0, settings))
    end
    return stages
end

_shift_setting(s::NamedTuple, dt) = map(v -> v isa PortCondition ? shift_profile(v, dt) : v, s)
