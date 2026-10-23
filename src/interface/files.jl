# Process and adsorbent files in JSON.
#
# Numbers are plain SI. A number may also be written as a string with a unit,
# such as "0.1 bar" (see Mocca.Units). Compositions and per-component values
# are objects keyed by component name. Step and isotherm types are looked up
# by their "type" in registries that users extend.

const LIBRARY_DIR = normpath(joinpath(@__DIR__, "..", "..", "models", "library"))
const ADSORBENT_LIBRARY = joinpath(LIBRARY_DIR, "adsorbents")
const PROCESS_LIBRARY = joinpath(LIBRARY_DIR, "processes")
const BUILTIN_COMPONENTS = Set(keys(COMPONENTS))

const PROCESS_FORMAT = "mocca-process"
const ADSORBENT_FORMAT = "mocca-adsorbent"
const FILE_VERSION = 1

_obj() = JSON.Object{String, Any}()
function _obj(pairs_...)
    o = _obj()
    for (k, v) in pairs_
        isnothing(v) || (o[String(k)] = v)
    end
    return o
end

_num(x::Real, loc, key) = Float64(x)
_num(x::AbstractString, loc, key) = _at(() -> Units.parse_quantity(x), "$loc, $key")
_num(x, loc, key) = input_error(loc, String(key), "expected a number, got $(repr(x))")
_optnum(d, key, loc) = haskey(d, key) && !isnothing(d[key]) ? _num(d[key], loc, key) : nothing

function _check_keys(d, allowed, loc)
    d isa AbstractDict || input_error(loc, "", "expected an object, got $(repr(d))")
    for k in keys(d)
        k in allowed || input_error(loc, k, "unknown key; expected one of $(join(allowed, ", "))")
    end
end
_require(d, key, loc) = haskey(d, key) ? d[key] : input_error(loc, key, "is required")

_by_name(names, v) = (o = _obj(); for (n, x) in zip(names, v); o[n] = x; end; o)
_json_float(x) = isfinite(x) ? x : nothing

# ------------------------------------------------------------------------------
# Isotherms
# ------------------------------------------------------------------------------

const ISOTHERM_TYPES = Dict{String, Any}()

"""
    register_isotherm!(name, T; read)

Make the isotherm data type `T` (an [`AbstractIsothermData`](@ref)) readable
from files under `"type": name`: `read(dict)` builds it. Write it by adding a
method of `Mocca.isotherm_to_dict(data)` that returns a `Dict` with `"type" =>
name`.
"""
register_isotherm!(name::AbstractString, T::Type; read) = (ISOTHERM_TYPES[String(name)] = (T, read); T)

isotherm_to_dict(x) = input_error("adsorbent", "isotherm", "isotherms of type $(typeof(x)) cannot be written to a file; give the isotherm by component, such as with dual_site_langmuir")

function isotherm_to_dict(data::DualSiteLangmuirData)
    comps = _obj()
    for (i, n) in enumerate(data.components)
        e = _obj("qsb" => data.qsb[i], "b0" => data.b0[i], "dUb" => data.ΔUb[i])
        if data.qsd[i] != 0 || data.d0[i] != 0 || data.ΔUd[i] != 0
            e["qsd"], e["d0"], e["dUd"] = data.qsd[i], data.d0[i], data.ΔUd[i]
        end
        comps[n] = e
    end
    return _obj("type" => "dual_site_langmuir", "reference_temperature" => data.reference_temperature,
        "basis" => String(data.basis), "components" => comps)
end

function _read_dsl(d)
    loc = "dual_site_langmuir isotherm"
    _check_keys(d, ("type", "reference_temperature", "basis", "components"), loc)
    entries = [String(n) => NamedTuple{Tuple(Symbol.(collect(keys(e))))}(Tuple(_num(v, "$loc for $n", k) for (k, v) in e)) for (n, e) in _require(d, "components", loc)]
    return dual_site_langmuir(entries...; reference_temperature = something(_optnum(d, "reference_temperature", loc), 298.15),
        basis = Symbol(get(d, "basis", "volume")))
end
register_isotherm!("dual_site_langmuir", DualSiteLangmuirData; read = _read_dsl)

function isotherm_from_dict(d, loc)
    d isa AbstractDict || input_error(loc, "isotherm", "expected an object with a \"type\"")
    t = String(_require(d, "type", "$loc isotherm"))
    haskey(ISOTHERM_TYPES, t) || input_error(loc, "isotherm", "unknown isotherm type $t; known types are $(join(sort(collect(keys(ISOTHERM_TYPES))), ", "))")
    return _at(() -> last(ISOTHERM_TYPES[t])(d), loc)
end

# ------------------------------------------------------------------------------
# Adsorbents
# ------------------------------------------------------------------------------

_mass_transfer_to_dict(::Nothing) = nothing
function _mass_transfer_to_dict(mt::LinearDrivingForce)
    D = mt.D_m isa Number ? mt.D_m : collect(mt.D_m)
    return _obj("type" => "linear_driving_force", "molecular_diffusivity" => D, "tortuosity" => mt.τ,
        "particle_porosity" => mt.ϵ_p, "particle_diameter" => mt.d_p)
end
_mass_transfer_to_dict(mt) = input_error("adsorbent", "mass_transfer", "mass transfer of type $(typeof(mt)) cannot be written to a file")

function _mass_transfer_from_dict(d, loc)
    isnothing(d) && return nothing
    _check_keys(d, ("type", "molecular_diffusivity", "tortuosity", "particle_porosity", "particle_diameter"), "$loc mass transfer")
    get(d, "type", "") == "linear_driving_force" || input_error(loc, "mass_transfer", "only \"linear_driving_force\" is known")
    D = d["molecular_diffusivity"]
    D_m = D isa AbstractVector ? SVector{length(D)}(Float64.(D)) : _num(D, loc, "molecular_diffusivity")
    return LinearDrivingForce(D_m, _num(d["tortuosity"], loc, "tortuosity"), _num(d["particle_porosity"], loc, "particle_porosity"),
        _num(d["particle_diameter"], loc, "particle_diameter"))
end

"""
    adsorbent_to_dict(a; gas = nothing) → Dict

The JSON description of adsorbent `a`. An adsorbed-phase heat capacity given
as a vector needs the `gas` it is ordered by.
"""
function adsorbent_to_dict(a::Adsorbent; gas = nothing)
    cp = a.adsorbed_heat_capacity
    if cp isa AbstractVector && !(eltype(cp) <: Pair)
        isnothing(gas) && input_error("adsorbent $(a.name)", "adsorbed_heat_capacity", "a vector of heat capacities needs the gas to write it by component")
        cp = _by_name(gas.names, cp)
    elseif !isnothing(cp)
        cp = _by_name([String(k) for (k, _) in pairs_of(cp)], [v for (_, v) in pairs_of(cp)])
    end
    return _obj("format" => ADSORBENT_FORMAT, "version" => FILE_VERSION, "name" => a.name, "reference" => a.reference,
        "density" => a.density, "heat_capacity" => a.heat_capacity, "particle_diameter" => a.particle_diameter,
        "particle_porosity" => a.particle_porosity, "tortuosity" => a.tortuosity, "adsorbed_heat_capacity" => cp,
        "isotherm" => isotherm_to_dict(a.isotherm),
        "mass_transfer" => _mass_transfer_to_dict(a.mass_transfer))
end

function adsorbent_from_dict(d; name = nothing)
    loc = "adsorbent" * (isnothing(name) ? "" : " $name")
    _check_keys(d, ("format", "version", "name", "reference", "density", "heat_capacity", "particle_diameter", "particle_porosity",
        "tortuosity", "adsorbed_heat_capacity", "isotherm", "mass_transfer"), loc)
    get(d, "format", ADSORBENT_FORMAT) == ADSORBENT_FORMAT || input_error(loc, "format", "expected \"$ADSORBENT_FORMAT\"")
    cp = get(d, "adsorbed_heat_capacity", nothing)
    cp = isnothing(cp) ? nothing : Dict(String(k) => _num(v, loc, "adsorbed_heat_capacity") for (k, v) in cp)
    return Adsorbent(String(get(d, "name", something(name, "")));
        isotherm = isotherm_from_dict(_require(d, "isotherm", loc), loc),
        density = _num(_require(d, "density", loc), loc, "density"),
        heat_capacity = _num(_require(d, "heat_capacity", loc), loc, "heat_capacity"),
        particle_diameter = _num(_require(d, "particle_diameter", loc), loc, "particle_diameter"),
        particle_porosity = _num(_require(d, "particle_porosity", loc), loc, "particle_porosity"),
        tortuosity = something(_optnum(d, "tortuosity", loc), 3.0),
        adsorbed_heat_capacity = cp,
        mass_transfer = _mass_transfer_from_dict(get(d, "mass_transfer", nothing), loc),
        reference = String(get(d, "reference", "")))
end

"""
    load_adsorbent(name) → Adsorbent

An adsorbent from the library by name (see [`library_adsorbents`](@ref)), or
from a JSON file at the path `name`.
"""
function load_adsorbent(name::AbstractString)
    path = isfile(name) ? name : joinpath(ADSORBENT_LIBRARY, name * ".json")
    isfile(path) || input_error("adsorbent $name", "", "no such file, and not in the library; library adsorbents are $(join(library_adsorbents(), ", "))")
    return adsorbent_from_dict(JSON.parsefile(path); name = splitext(basename(path))[1])
end

"Names of the adsorbents in Mocca's library."
library_adsorbents() = sort([splitext(f)[1] for f in readdir(ADSORBENT_LIBRARY) if endswith(f, ".json")])

"""
    save_adsorbent(path, a; gas = nothing)

Write adsorbent `a` to a JSON file.
"""
save_adsorbent(path::AbstractString, a::Adsorbent; gas = nothing) = (JSON.json(path, adsorbent_to_dict(a; gas = gas); pretty = true); path)

# ------------------------------------------------------------------------------
# Steps
# ------------------------------------------------------------------------------

function step_to_dict(s::AbstractStep)
    name = step_type_name(typeof(s))
    isnothing(name) && input_error("cycle", "", "steps of type $(typeof(s)) are not registered; call register_step!(\"name\", $(typeof(s)))")
    o = _obj("type" => name)
    for f in fieldnames(typeof(s))
        v = getfield(s, f)
        (isnothing(v) || (f == :name && isempty(v))) && continue
        o[String(f)] = v isa Symbol ? String(v) : (v isa AbstractVector{Symbol} ? String.(v) : v)
    end
    return o
end

function _from_json(FT, v, loc, key)
    isnothing(v) && return nothing
    Symbol <: FT && v isa AbstractString && return Symbol(v)
    Vector{Symbol} <: FT && v isa AbstractVector && return Symbol.(v)
    Float64 <: FT && return _num(v, loc, key)
    Int <: FT && v isa Real && return Int(v)
    String <: FT && v isa AbstractString && return String(v)
    return v
end

function step_from_dict(d, loc)
    t = String(_require(d, "type", loc))
    t == "concurrent" && return Concurrent([step_from_dict(s, "$loc, concurrent step") for s in _require(d, "steps", loc)]...;
        name = get(d, "name", ""))
    haskey(STEP_TYPES, t) || input_error(loc, "type", "unknown step type $t; known types are $(join(sort(vcat(collect(keys(STEP_TYPES)), "concurrent")), ", "))")
    T = STEP_TYPES[t]
    kw = Dict{Symbol, Any}()
    for (k, v) in d
        k == "type" && continue
        f = Symbol(k)
        hasfield(T, f) || input_error("$loc ($t)", k, "unknown field; a $t step has $(join(fieldnames(T), ", "))")
        kw[f] = _from_json(fieldtype(T, f), v, loc, k)
    end
    return _at(() -> T(; kw...), "$loc ($t)")
end

_element_to_dict(e::AbstractStep) = step_to_dict(e)
_element_to_dict(e::Concurrent) = _obj("type" => "concurrent", "name" => isempty(e.name) ? nothing : e.name, "steps" => [_element_to_dict(s) for s in e.steps])
_element_to_dict(e) = input_error("cycle", "", "a $(typeof(e)) in a cycle cannot be written to a file; use step templates")

function cycle_to_dict(c::Cycle)
    if length(c.sequences) == 1 && isempty(first(c.sequences[1]))
        return _obj("steps" => [_element_to_dict(e) for e in last(c.sequences[1])])
    end
    return _obj("sequences" => [_obj("label" => l, "steps" => [_element_to_dict(e) for e in s]) for (l, s) in c.sequences],
        "cycle_time" => c.cycle_time)
end

function cycle_from_dict(d)
    _check_keys(d, ("steps", "sequences", "cycle_time"), "cycle")
    if haskey(d, "steps")
        return Cycle([step_from_dict(s, "step $i of the cycle") for (i, s) in enumerate(d["steps"])])
    end
    seqs = [String(_require(s, "label", "cycle sequence")) => [step_from_dict(e, "step $i of the $(s["label"]) sequence") for (i, e) in enumerate(_require(s, "steps", "cycle sequence"))]
        for s in _require(d, "sequences", "cycle")]
    return Cycle(seqs...; cycle_time = _optnum(d, "cycle_time", "cycle"))
end

# ------------------------------------------------------------------------------
# Processes
# ------------------------------------------------------------------------------

_thermal_to_json(::WithWall) = "with_wall"
_thermal_to_json(::Adiabatic) = "adiabatic"
_thermal_to_json(t::Isothermal) = _obj("type" => "isothermal", "temperature" => t.temperature)
function _thermal_from_json(x, loc)
    x == "with_wall" && return WithWall()
    x == "adiabatic" && return Adiabatic()
    x isa AbstractDict && get(x, "type", "") == "isothermal" && return Isothermal(_num(x["temperature"], loc, "temperature"))
    input_error(loc, "thermal", "use \"with_wall\", \"adiabatic\" or {\"type\": \"isothermal\", \"temperature\": ...}")
end

function _column_to_dict(c::Column)
    return _obj("length" => c.length, "radius" => c.radius, "outer_radius" => c.outer_radius, "porosity" => c.porosity,
        "axial_conductivity" => c.axial_conductivity,
        "wall" => _obj("density" => c.wall.density, "heat_capacity" => c.wall.heat_capacity, "conductivity" => c.wall.conductivity),
        "inner_heat_transfer" => c.inner_heat_transfer, "outer_heat_transfer" => c.outer_heat_transfer,
        "ambient_temperature" => c.ambient_temperature)
end

function _column_from_dict(d, loc)
    keys_ = ("length", "radius", "outer_radius", "porosity", "axial_conductivity", "wall", "inner_heat_transfer", "outer_heat_transfer", "ambient_temperature")
    _check_keys(d, keys_, loc)
    kw = Dict{Symbol, Any}()
    for k in keys_
        haskey(d, k) || continue
        if k == "wall"
            w = d[k]
            _check_keys(w, ("density", "heat_capacity", "conductivity"), "$loc wall")
            kw[:wall] = Wall(; (Symbol(a) => _num(b, "$loc wall", a) for (a, b) in w)...)
        else
            kw[Symbol(k)] = _num(d[k], loc, k)
        end
    end
    return _at(() -> Column(; kw...), loc)
end

function _condition_to_dict(fs, c::PortCondition)
    c.pressure isa Real || input_error("device boundary", "pressure", "only constant boundary pressures can be written")
    return _obj("pressure" => c.pressure, "temperature" => c.temperature, "composition" => _by_name(fs.gas.names, collect(c.composition)))
end

_name_of(x::Symbol) = String(x)

# The transport model of a gas, by name: models other than those Mocca builds
# from the components cannot be written
function _gas_transport_to_dict(gas::Gas)
    t = gas.transport
    t isa ConstantTransport && return nothing
    if t isa SutherlandTransport
        made = sutherland_transport(component.(gas.names), gas.molar_masses; knudsen_pressure = t.P_K, diffusivity_exponent = t.m,
            pressure_dependent = t.pressure_dependent)
        t == made && t.m == 1.75 && !t.pressure_dependent && return "sutherland"
        t == made && return _obj("type" => "sutherland", "diffusivity_exponent" => t.m, "pressure_dependent" => t.pressure_dependent,
            "knudsen_pressure" => t.P_K)
    end
    input_error("gas", "transport", "a transport model of type $(typeof(gas.transport)) not made from the components cannot be written to a file")
end

# The transport of a gas file entry: a name, or a Sutherland model with options
function _transport_from_dict(t)
    t isa AbstractString && t in ("constant", "sutherland") && return Symbol(t)
    t isa AbstractDict || input_error("gas", "transport", "expected \"constant\", \"sutherland\" or an object with a \"type\", got $(repr(t))")
    _check_keys(t, ("type", "knudsen_pressure", "diffusivity_exponent", "pressure_dependent"), "gas transport")
    get(t, "type", "") == "sutherland" || input_error("gas", "transport", "only the \"sutherland\" type takes options")
    P_K = something(_optnum(t, "knudsen_pressure", "gas transport"), 0.0)
    pd = get(t, "pressure_dependent", P_K > 0)
    pd isa Bool || input_error("gas transport", "pressure_dependent", "expected true or false, got $(repr(pd))")
    return (knudsen_pressure = P_K, diffusivity_exponent = something(_optnum(t, "diffusivity_exponent", "gas transport"), 1.75),
        pressure_dependent = pd)
end

# The heat capacities of a gas, as keys of its JSON object
function _gas_thermodynamics_to_dict(gas::Gas)
    th = gas.thermodynamics
    th isa IdealGasConstantCp && return ("heat_capacities" => _by_name(gas.names, gas.heat_capacities),)
    th isa IdealGasShomate && return ("thermodynamics" => "shomate",
        "shomate" => _by_name(gas.names, [collect(c) for c in shomate_coefficients(th)]))
    input_error("gas", "thermodynamics", "thermodynamics of type $(typeof(th)) cannot be written to a file; use constant or Shomate heat capacities")
end

"""
    process_to_dict(p) → Dict

The JSON description of process `p`; see [`write_process`](@ref).
"""
function process_to_dict(p::Process)
    fs = p.flowsheet
    gas = fs.gas
    o = _obj("format" => PROCESS_FORMAT, "version" => FILE_VERSION, "name" => p.name)
    custom = [component(n) for n in gas.names if !(n in BUILTIN_COMPONENTS)]
    isempty(custom) || (o["components"] = [_obj("name" => c.name, "formula" => c.formula, "cas" => c.cas,
        "molar_mass" => c.molar_mass, "heat_capacity" => c.heat_capacity, "viscosity" => c.viscosity,
        "shomate" => isnothing(c.shomate) ? nothing : collect(c.shomate), "sutherland" => c.sutherland,
        "thermal_conductivity" => c.thermal_conductivity, "conductivity_exponent" => c.conductivity_exponent) for c in custom])
    o["gas"] = _obj("components" => gas.names, "molar_masses" => _by_name(gas.names, gas.molar_masses),
        _gas_thermodynamics_to_dict(gas)..., "transport" => _gas_transport_to_dict(gas), "viscosity" => gas.viscosity, "density" => gas.density,
        "diffusivity" => gas.diffusivity isa Number ? gas.diffusivity : _by_name(gas.names, gas.diffusivity))
    o["feed"] = _obj("composition" => _by_name(gas.names, fs.feed.composition), "temperature" => fs.feed.temperature,
        "pressure" => fs.feed.pressure, "light_product" => _by_name(gas.names, fs.feed.light_product))
    # Adsorbents and columns shared by several beds are written once, under
    # their name where they have one
    adsorbents, columns = Pair{String, Any}[], Pair{String, Any}[]
    function key_of(list, x, name, prefix)
        i = findfirst(y -> last(y) === x, list)
        isnothing(i) || return first(list[i])
        k = (isempty(name) || any(y -> first(y) == name, list)) ? "$(prefix)_$(length(list) + 1)" : name
        push!(list, k => x)
        return k
    end
    units = Any[]
    for u in holdup_units(fs)
        spec = get(fs.specs, u, nothing)
        spec isa BedSpec || input_error("unit $u", "", "only beds added with add_bed! can be written to a file")
        ka = key_of(adsorbents, spec.adsorbent, spec.adsorbent.name, "adsorbent")
        kc = key_of(columns, spec.column, "", "column")
        push!(units, _obj("name" => String(u), "type" => "fixed_bed", "adsorbent" => ka, "column" => kc,
            "ncells" => spec.ncells, "thermal" => _thermal_to_json(spec.thermal),
            "dispersion_velocity" => spec.dispersion_velocity isa Symbol ? String(spec.dispersion_velocity) : spec.dispersion_velocity,
            "dispersion" => spec.dispersion, "permeability" => spec.permeability))
    end
    o["adsorbents"] = _obj((k => adsorbent_to_dict(a; gas = gas) for (k, a) in adsorbents)...)
    o["columns"] = _obj((k => _column_to_dict(c) for (k, c) in columns)...)
    o["units"] = units
    devices = Any[]
    for d in unit_names(fs)
        is_device(fs, d) || continue
        role = device_role(fs, d)
        isnothing(role) && input_error("device $d", "", "only devices added with add_device! can be written to a file")
        e = _obj("name" => String(d), "role" => String(role))
        att = attachments(fs, d)
        if role == :link
            from = only(a.unit for a in att if a.device.port == :inlet)
            to = only(a.unit for a in att if a.device.port == :outlet)
            e["from"], e["to"] = [String(from.name), String(from.port)], [String(to.name), String(to.port)]
        else
            a = only(att)
            e["at"] = [String(a.unit.name), String(a.unit.port)]
            port = role in INLET_ROLES ? :inlet : :outlet
            PortRef(d, port) in fs.default_boundaries || (e["boundary"] = _condition_to_dict(fs, boundary_condition(fs, d, port)))
        end
        push!(devices, e)
    end
    o["devices"] = devices
    o["cycle"] = cycle_to_dict(p.cycle)
    ts = p.timestepping
    o["time_stepping"] = _obj("max_dt" => ts.max_dt, "first_dt" => ts.first_dt,
        "keep_step_pattern" => ts.keep_step_pattern isa Symbol ? String(ts.keep_step_pattern) : ts.keep_step_pattern)
    op = p.options
    o["options"] = _obj("key_component" => gas.names[op.key_component], "ambient_pressure" => op.ambient_pressure,
        "pump_efficiency" => op.pump_efficiency, "heat_capacity_ratio" => op.heat_capacity_ratio, "ramp_rate" => op.ramp_rate)
    if !isempty(p.initial_conditions)
        o["initial_state"] = _obj((String(b) => _obj(("$k" => (k == :composition ? _by_name(gas.names, composition(gas, v)) : v) for (k, v) in pairs(ic))...)
            for (b, ic) in p.initial_conditions)...)
    end
    if !isempty(p.parameter_settings)
        o["parameters"] = _obj((String(u) => _obj((String(n) => v for (n, v) in d)...) for (u, d) in p.parameter_settings)...)
    end
    if !isempty(p.variables)
        o["variables"] = [_obj("name" => String(v.name), "places" => [_place_to_dict(x) for x in v.places],
            "bounds" => [_json_float(v.lower), _json_float(v.upper)], "description" => isempty(v.description) ? nothing : v.description)
            for v in p.variables]
    end
    return o
end

_place_to_dict(x::StepPlace) = _obj("step" => x.step, "field" => String(x.field), "sequence" => isempty(x.sequence) ? nothing : x.sequence)
_place_to_dict(x::FeedPlace) = _obj("feed" => String(x.field))
_place_to_dict(x::ParameterPlace) = _obj("unit" => String(x.unit), "parameter" => String(x.name))

function _place_from_dict(d, loc)
    _check_keys(d, ("step", "field", "sequence", "feed", "unit", "parameter"), loc)
    haskey(d, "feed") && return FeedPlace(Symbol(d["feed"]))
    haskey(d, "parameter") && return ParameterPlace(Symbol(_require(d, "unit", loc)), Symbol(d["parameter"]))
    return StepPlace(_require(d, "step", loc), _require(d, "field", loc); sequence = get(d, "sequence", ""))
end

"""
    write_process(path, p)

Write process `p` to a JSON file: its gas, feed, adsorbents, columns, units,
devices, cycle, time stepping, options, initial state, parameter settings and
design variables. Read it back with [`read_process`](@ref).
"""
write_process(path::AbstractString, p::Process) = (JSON.json(path, process_to_dict(p); pretty = true); path)

"""
    read_process(path) → Process

Read a process from a JSON file written by [`write_process`](@ref) or by
hand, or from the library by name (see [`library_processes`](@ref)). Numbers
are SI, or strings with a unit such as `"0.1 bar"`. An adsorbent may be given
by the name of a library adsorbent or the path of an adsorbent file, relative
to the process file.
"""
function read_process(path::AbstractString)
    file = isfile(path) ? path : joinpath(PROCESS_LIBRARY, path * ".json")
    isfile(file) || input_error("process $path", "", "no such file, and not in the library; library processes are $(join(library_processes(), ", "))")
    return process_from_dict(JSON.parsefile(file); base_dir = dirname(abspath(file)))
end

"Names of the processes in Mocca's library."
library_processes() = isdir(PROCESS_LIBRARY) ? sort([splitext(f)[1] for f in readdir(PROCESS_LIBRARY) if endswith(f, ".json")]) : String[]

function process_from_dict(d; base_dir = pwd())
    _check_keys(d, ("format", "version", "name", "components", "gas", "feed", "adsorbents", "columns", "units", "devices",
        "cycle", "time_stepping", "options", "initial_state", "parameters", "variables"), "process")
    get(d, "format", PROCESS_FORMAT) == PROCESS_FORMAT || input_error("process", "format", "expected \"$PROCESS_FORMAT\"")
    get(d, "version", FILE_VERSION) <= FILE_VERSION || input_error("process", "version", "file version $(d["version"]) is newer than this Mocca reads ($FILE_VERSION)")
    for c in get(d, "components", Any[])
        loc = "component $(get(c, "name", "?"))"
        _check_keys(c, ("name", "formula", "cas", "molar_mass", "heat_capacity", "viscosity", "shomate", "sutherland",
            "thermal_conductivity", "conductivity_exponent"), loc)
        register_component!(Component(_require(c, "name", loc); formula = get(c, "formula", c["name"]), cas = get(c, "cas", ""),
            molar_mass = _num(_require(c, "molar_mass", loc), loc, "molar_mass"),
            heat_capacity = _num(_require(c, "heat_capacity", loc), loc, "heat_capacity"),
            viscosity = _num(_require(c, "viscosity", loc), loc, "viscosity"),
            shomate = haskey(c, "shomate") ? [_num(v, loc, "shomate") for v in c["shomate"]] : nothing,
            sutherland = _optnum(c, "sutherland", loc), thermal_conductivity = _optnum(c, "thermal_conductivity", loc),
            conductivity_exponent = _optnum(c, "conductivity_exponent", loc)))
    end
    g = _require(d, "gas", "process")
    _check_keys(g, ("components", "molar_masses", "heat_capacities", "thermodynamics", "shomate", "transport", "viscosity",
        "density", "diffusivity"), "gas")
    pervals(x) = isnothing(x) ? nothing : Dict(String(k) => _num(v, "gas", k) for (k, v) in x)
    diff = get(g, "diffusivity", 1.6e-5)
    thermo = get(g, "thermodynamics", nothing)
    thermo in (nothing, "constant", "shomate") || input_error("gas", "thermodynamics", "expected \"constant\" or \"shomate\", got $(repr(thermo))")
    transport = _transport_from_dict(get(g, "transport", "constant"))
    shomate = haskey(g, "shomate") ? Dict(String(k) => [_num(x, "gas", "shomate") for x in v] for (k, v) in g["shomate"]) : nothing
    gas = Gas(String.(_require(g, "components", "gas")); molar_masses = pervals(get(g, "molar_masses", nothing)),
        heat_capacities = pervals(get(g, "heat_capacities", nothing)),
        thermodynamics = isnothing(thermo) ? nothing : Symbol(thermo), shomate = shomate,
        transport = transport isa Symbol ? transport : sutherland_transport(String.(g["components"]);
            molar_masses = pervals(get(g, "molar_masses", nothing)), transport...),
        viscosity = _optnum(g, "viscosity", "gas"),
        density = _optnum(g, "density", "gas"), diffusivity = diff isa AbstractDict ? pervals(diff) : _num(diff, "gas", "diffusivity"))
    fd = _require(d, "feed", "process")
    _check_keys(fd, ("composition", "temperature", "pressure", "light_product"), "feed")
    comp(x, loc) = x isa AbstractVector ? Float64.(x) : Dict(String(k) => _num(v, loc, k) for (k, v) in x)
    feed = Feed(gas; composition = comp(_require(fd, "composition", "feed"), "feed composition"),
        temperature = something(_optnum(fd, "temperature", "feed"), 298.15), pressure = something(_optnum(fd, "pressure", "feed"), 101325.0),
        light_product = haskey(fd, "light_product") ? comp(fd["light_product"], "feed light product") : nothing)
    adsorbents = Dict{String, Adsorbent}()
    for (k, v) in get(d, "adsorbents", Dict())
        if v isa AbstractString
            rel = joinpath(base_dir, v)
            adsorbents[k] = load_adsorbent(isfile(rel) ? rel : v)
        else
            adsorbents[k] = adsorbent_from_dict(v; name = k)
        end
    end
    columns = Dict(String(k) => _column_from_dict(v, "column $k") for (k, v) in get(d, "columns", Dict()))
    fs = Flowsheet(gas, feed)
    for u in _require(d, "units", "process")
        loc = "unit $(get(u, "name", "?"))"
        _check_keys(u, ("name", "type", "adsorbent", "column", "ncells", "thermal", "dispersion_velocity", "dispersion", "permeability"), loc)
        get(u, "type", "fixed_bed") == "fixed_bed" || input_error(loc, "type", "only \"fixed_bed\" units are known")
        a = String(_require(u, "adsorbent", loc))
        a_ = haskey(adsorbents, a) ? adsorbents[a] : load_adsorbent(a)
        c = String(_require(u, "column", loc))
        haskey(columns, c) || input_error(loc, "column", "no column $c in \"columns\"")
        dv = get(u, "dispersion_velocity", "adsorption")
        add_bed!(fs, Symbol(_require(u, "name", loc)), a_, columns[c]; ncells = Int(get(u, "ncells", 30)),
            thermal = _thermal_from_json(get(u, "thermal", "with_wall"), loc),
            dispersion_velocity = dv isa AbstractString && dv == "adsorption" ? :adsorption : _num(dv, loc, "dispersion_velocity"),
            dispersion = _optnum(u, "dispersion", loc), permeability = _optnum(u, "permeability", loc))
    end
    pair(x, loc, key) = (x isa AbstractVector && length(x) == 2) ? (Symbol(x[1]), Symbol(x[2])) : input_error(loc, key, "expected [unit, port], got $(repr(x))")
    for e in _require(d, "devices", "process")
        loc = "device $(get(e, "name", "?"))"
        _check_keys(e, ("name", "role", "at", "from", "to", "boundary"), loc)
        name = Symbol(_require(e, "name", loc))
        role = Symbol(_require(e, "role", loc))
        if role == :link
            add_device!(fs, name, role; from = pair(_require(e, "from", loc), loc, "from"), to = pair(_require(e, "to", loc), loc, "to"))
        else
            b = get(e, "boundary", nothing)
            bnd = isnothing(b) ? nothing : begin
                _check_keys(b, ("pressure", "temperature", "composition"), "$loc boundary")
                kw = Dict{Symbol, Any}()
                haskey(b, "pressure") && (kw[:pressure] = _num(b["pressure"], loc, "pressure"))
                haskey(b, "temperature") && (kw[:temperature] = _num(b["temperature"], loc, "temperature"))
                haskey(b, "composition") && (kw[:composition] = comp(b["composition"], "$loc boundary"))
                (; kw...)
            end
            add_device!(fs, name, role; at = pair(_require(e, "at", loc), loc, "at"), boundary = bnd)
        end
    end
    cycle = cycle_from_dict(_require(d, "cycle", "process"))
    t = get(d, "time_stepping", Dict())
    _check_keys(t, ("max_dt", "first_dt", "keep_step_pattern"), "time_stepping")
    max_dt = something(_optnum(t, "max_dt", "time_stepping"), 1.0)
    kp = get(t, "keep_step_pattern", "auto")
    timestepping = TimeStepping(; max_dt = max_dt, first_dt = something(_optnum(t, "first_dt", "time_stepping"), min(0.05, max_dt)),
        keep_step_pattern = kp isa Bool ? kp : Symbol(kp))
    op = get(d, "options", Dict())
    _check_keys(op, ("key_component", "ambient_pressure", "pump_efficiency", "heat_capacity_ratio", "ramp_rate"), "options")
    kc = get(op, "key_component", 1)
    ics = Dict{Symbol, Any}()
    for (b, ic) in get(d, "initial_state", Dict())
        loc = "initial state of $b"
        _check_keys(ic, ("pressure", "temperature", "wall_temperature", "composition"), loc)
        kw = Dict{Symbol, Any}()
        for k in ("pressure", "temperature", "wall_temperature")
            haskey(ic, k) && (kw[Symbol(k)] = _num(ic[k], loc, k))
        end
        haskey(ic, "composition") && (kw[:composition] = comp(ic["composition"], loc))
        ics[Symbol(b)] = (; kw...)
    end
    prm = Dict{Symbol, Any}(Symbol(u) => Dict{Symbol, Any}(Symbol(n) => (v isa AbstractVector ? Float64.(v) : _num(v, "parameters of $u", n)) for (n, v) in pd)
        for (u, pd) in get(d, "parameters", Dict()))
    p = Process(fs, cycle; name = String(get(d, "name", "")), timestepping = timestepping,
        key_component = kc isa AbstractString ? kc : Int(kc),
        ambient_pressure = something(_optnum(op, "ambient_pressure", "options"), 101325.0),
        pump_efficiency = something(_optnum(op, "pump_efficiency", "options"), 0.72),
        heat_capacity_ratio = _optnum(op, "heat_capacity_ratio", "options"),
        ramp_rate = something(_optnum(op, "ramp_rate", "options"), 0.5),
        initial_state = ics, parameters = prm)
    for v in get(d, "variables", Any[])
        loc = "variable $(get(v, "name", "?"))"
        _check_keys(v, ("name", "places", "bounds", "description"), loc)
        bnds = get(v, "bounds", [nothing, nothing])
        lo = isnothing(bnds[1]) ? -Inf : _num(bnds[1], loc, "bounds")
        hi = isnothing(bnds[2]) ? Inf : _num(bnds[2], loc, "bounds")
        add_variable!(p, Symbol(_require(v, "name", loc)), [_place_from_dict(x, loc) for x in _require(v, "places", loc)]...;
            bounds = (lo, hi), description = get(v, "description", ""))
    end
    return p
end
