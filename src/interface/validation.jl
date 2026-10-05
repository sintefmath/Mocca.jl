"""
    validate(p::Process) → Vector{MoccaInputError}

Check the description of a process and return the problems found, or an
empty vector. Constructors already throw the first error they find; this
also reports problems that do not stop a simulation but likely give wrong
results:

- a cycle that does not compile;
- a unit with no devices, or a flowsheet without feed or extract devices;
- a bed whose axial dispersion was set for a velocity far from that of its
  adsorption step;
- design variables whose places cannot be read.
"""
function validate(p::Process)
    issues = MoccaInputError[]
    fs = p.flowsheet
    note(loc, field, msg) = push!(issues, MoccaInputError(String(loc), String(field), String(msg)))
    compiled = try
        compile_cycle(p.cycle, fs; ramp_rate = p.options.ramp_rate)
    catch e
        e isa MoccaInputError ? push!(issues, e) : note("cycle", "", sprint(showerror, e))
        nothing
    end
    for u in holdup_units(fs)
        isempty(vcat([devices_at(fs, u, q) for q in keys(ports(fs[u]))]...)) && note("unit $u", "", "no devices are attached")
    end
    isempty(devices_with_role(fs, :feed)) && note("flowsheet", "", "there is no feed device")
    isempty(devices_with_role(fs, :extract)) && note("flowsheet", "", "there is no extract device, so purity, recovery and energy cannot be measured")
    if !isnothing(compiled)
        for b in beds(fs)
            spec = get(fs.specs, b, nothing)
            spec isa BedSpec && isnothing(spec.dispersion) && haskey(compiled.velocities, b) || continue
            v = maximum(compiled.velocities[b])
            u = spec.velocity_used
            (v > 2u || v < u / 2) && note("bed $b", "dispersion_velocity",
                "its axial dispersion is set for $(round(u, sigdigits = 3)) m/s, but its adsorption step feeds at $(round(v, sigdigits = 3)) m/s")
        end
    end
    for v in p.variables, x in v.places
        try
            _get(p, x)
        catch e
            e isa MoccaInputError ? push!(issues, e) : note("variable $(v.name)", "", sprint(showerror, e))
        end
    end
    return issues
end
