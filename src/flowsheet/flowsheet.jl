"""
    Flowsheet()

A set of named units (beds, devices, ...) and the connections between their
ports. Build it with [`add_unit!`](@ref) and [`connect!`](@ref), then turn it
into a `Jutul.MultiModel` with [`setup_flowsheet_model`](@ref).

Every connection joins a port of a unit with holdup to a port of a
[`FlowDevice`](@ref). Device ports that are not connected are boundaries and
need a [`PortCondition`](@ref), set with [`set_boundary!`](@ref) or per stage.
"""
struct Flowsheet
    units::Dict{Symbol, Any}
    order::Vector{Symbol}
    connections::Vector{Tuple{Symbol, Symbol, Symbol, Symbol}}  # (unit, unit port, device, device port)
    boundaries::Dict{Tuple{Symbol, Symbol}, PortCondition}
end
Flowsheet() = Flowsheet(Dict{Symbol, Any}(), Symbol[], Tuple{Symbol, Symbol, Symbol, Symbol}[], Dict{Tuple{Symbol, Symbol}, PortCondition}())

Base.getindex(fs::Flowsheet, name::Symbol) = fs.units[name]

"""
    add_unit!(fs, name, model)

Add a unit `model` (a `Jutul.SimulationModel` of a [`MoccaSystem`](@ref)) to the
flowsheet under `name`.
"""
function add_unit!(fs::Flowsheet, name::Symbol, model::MoccaModel)
    !haskey(fs.units, name) || error("Flowsheet already has a unit named $name")
    fs.units[name] = model
    push!(fs.order, name)
    return fs
end

_is_device(model) = model isa FlowDeviceModel

"""
    connect!(fs, (unit, port) => (device, device_port))

Connect `port` of a unit with holdup to `device_port` (`:inlet` or `:outlet`)
of a flow device. The pair may be given in either order.
"""
function connect!(fs::Flowsheet, conn::Pair{Tuple{Symbol, Symbol}, Tuple{Symbol, Symbol}})
    (a, pa), (b, pb) = conn
    haskey(fs.units, a) || error("Unknown unit $a")
    haskey(fs.units, b) || error("Unknown unit $b")
    if _is_device(fs.units[a]) && !_is_device(fs.units[b])
        (a, pa), (b, pb) = (b, pb), (a, pa)
    end
    !_is_device(fs.units[a]) || error("Cannot connect two devices ($a and $b) directly; put a unit with holdup between them")
    _is_device(fs.units[b]) || error("Cannot connect two units with holdup ($a and $b) directly; connect them through a device")
    haskey(ports(fs.units[a]), pa) || error("Unit $a has no port $pa, available: $(keys(ports(fs.units[a])))")
    haskey(ports(fs.units[b]), pb) || error("Device $b has no port $pb, available: $(keys(ports(fs.units[b])))")
    for c in fs.connections
        (c[3], c[4]) == (b, pb) && error("Port $pb of device $b is already connected")
    end
    push!(fs.connections, (a, pa, b, pb))
    return fs
end

"""
    set_boundary!(fs, device, port, condition::PortCondition)

Set the default boundary condition for an unconnected device port. Stages can
override it.
"""
function set_boundary!(fs::Flowsheet, device::Symbol, port::Symbol, condition::PortCondition)
    _is_device(fs.units[device]) || error("$device is not a device")
    fs.boundaries[(device, port)] = condition
    return fs
end

"Whether `port` of `device` is connected to a unit."
is_connected(fs::Flowsheet, device::Symbol, port::Symbol) = any(c -> (c[3], c[4]) == (device, port), fs.connections)

"Unconnected device ports, as `(device, port)` pairs."
function boundary_ports(fs::Flowsheet)
    out = Tuple{Symbol, Symbol}[]
    for name in fs.order
        _is_device(fs.units[name]) || continue
        for p in DEVICE_PORTS
            is_connected(fs, name, p) || push!(out, (name, p))
        end
    end
    return out
end

"""
    setup_flowsheet_model(fs) → Jutul.MultiModel

Assemble the flowsheet into a `MultiModel` with a [`PortStateCT`](@ref) and a
[`StreamCT`](@ref) for every connection.
"""
function setup_flowsheet_model(fs::Flowsheet)
    for (d, p) in boundary_ports(fs)
        haskey(fs.boundaries, (d, p)) || error("Port $p of device $d is not connected and has no boundary condition; use set_boundary!")
    end
    models = NamedTuple{Tuple(fs.order)}(Tuple(fs.units[k] for k in fs.order))
    mm = Jutul.MultiModel(models)
    for (u, up, d, dp) in fs.connections
        unit = fs.units[u]
        cell = ports(unit)[up]
        k = port_index(fs.units[d], dp)
        Jutul.add_cross_term!(mm, PortStateCT(cell); target = d, source = u, equation = port_equation(dp))
        Jutul.add_cross_term!(mm, StreamCT(cell, k); target = u, source = d, equation = :mass_conservation)
        if haskey(unit.equations, :energy_column) && unit.equations[:energy_column] isa Jutul.ConservationLaw
            Jutul.add_cross_term!(mm, StreamCT(cell, k); target = u, source = d, equation = :energy_column)
        end
    end
    return mm
end

function _unit_port_condition(fs::Flowsheet, unit_state, u, up)
    cell = ports(fs.units[u])[up]
    return (pressure = unit_state[:Pressure][cell], temperature = unit_state[:Temperature][cell], composition = unit_state[:y][:, cell])
end

"""
    setup_flowsheet_state(fs; kwarg...) → Dict

Initial state for every unit. Pass the initial state of each unit with holdup as
a keyword named after the unit, e.g. `Bed = setup_process_state(fs[:Bed]; ...)`.
Device states are initialised from what their ports are attached to.
"""
function setup_flowsheet_state(fs::Flowsheet; kwarg...)
    states = Dict{Symbol, Any}()
    for name in fs.order
        _is_device(fs.units[name]) && continue
        haskey(kwarg, name) || error("Missing initial state for unit $name")
        states[name] = kwarg[name]
    end
    for name in fs.order
        model = fs.units[name]
        _is_device(model) || continue
        conds = map(DEVICE_PORTS) do p
            idx = findfirst(c -> (c[3], c[4]) == (name, p), fs.connections)
            if isnothing(idx)
                fs.boundaries[(name, p)]
            else
                u, up = fs.connections[idx][1], fs.connections[idx][2]
                _unit_port_condition(fs, states[u], u, up)
            end
        end
        states[name] = setup_device_state(model; inlet = conds[1], outlet = conds[2])
    end
    return states
end

"""
    setup_flowsheet_parameters(fs; kwarg...) → Dict

Parameters for every unit. Beds use [`setup_process_parameters`](@ref); pass a
`Dict` or `NamedTuple` of overrides per unit as a keyword named after it.
"""
function setup_flowsheet_parameters(fs::Flowsheet; kwarg...)
    prm = Dict{Symbol, Any}()
    for name in fs.order
        model = fs.units[name]
        extra = get(kwarg, name, NamedTuple())
        if model isa FixedBedModel
            prm[name] = setup_process_parameters(model; pairs(extra)...)
        else
            prm[name] = Jutul.setup_parameters(model; pairs(extra)...)
        end
    end
    return prm
end
