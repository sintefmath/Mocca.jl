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
"""
struct Stage
    name::String
    duration::Float64
    settings::Dict{Symbol, Any}
end
Stage(name, duration; settings...) = Stage(String(name), Float64(duration), Dict{Symbol, Any}(settings))

"""
    setup_schedule(fs, model, stages; num_cycles = 1, max_dt = 1.0) → (forces, timesteps)

Forces and time steps for `num_cycles` repetitions of `stages` on the flowsheet
`fs` and its `MultiModel`. Each stage is split into steps no longer than
`max_dt`. Boundary ports use the stage's port condition if given, otherwise the
flowsheet default from [`set_boundary!`](@ref).
"""
function setup_schedule(fs::Flowsheet, model::Jutul.MultiModel, stages::AbstractVector{Stage}; num_cycles = 1, max_dt = 1.0)
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
            n = max(1, ceil(Int, st.duration / max_dt - 1e-9))
            dt = st.duration / n
            append!(timesteps, fill(dt, n))
            append!(forces, fill(f, n))
            t += st.duration
        end
    end
    @assert isapprox(t, num_cycles * cycle_time)
    return (forces, timesteps)
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
