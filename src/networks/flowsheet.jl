# Flowsheets: named units, the connections between their ports, and the
# boundary conditions of unconnected device ports.
#
# A UnitNetwork is that topology and nothing else; setup_flowsheet_model turns
# it into a Jutul.MultiModel. The Flowsheet of the user interface
# (interface/builder.jl) holds one, and adds the gas, the feed, the roles of
# devices and how each bed was built. Functions here take either, through
# `network`.

"""
    AbstractFlowsheet

Units and the connections between them, which [`setup_flowsheet_model`](@ref)
assembles: a [`UnitNetwork`](@ref), or a [`Flowsheet`](@ref Mocca.Interface.Flowsheet) of the user
interface, which holds one. Subtypes provide [`network`](@ref).
"""
abstract type AbstractFlowsheet end

"""
    network(fs) → UnitNetwork

The units and connections of `fs`.
"""
function network end

"""
    UnitNetwork()

Named units (beds, devices, ...), the connections between their ports, and
the boundary conditions of device ports that are not connected. Build it with
[`add_unit!`](@ref), [`connect!`](@ref) and [`set_boundary!`](@ref), and turn
it into a `Jutul.MultiModel` with [`setup_flowsheet_model`](@ref). For a
process described in terms of materials and device roles, use a
[`Flowsheet`](@ref Mocca.Interface.Flowsheet) instead.

Unconnected device ports are boundaries and need a [`PortCondition`](@ref),
set with [`set_boundary!`](@ref) or per stage.
"""
struct UnitNetwork <: AbstractFlowsheet
    units::Dict{Symbol, MoccaModel}
    order::Vector{Symbol}
    connections::Vector{AbstractConnection}
    boundaries::Dict{PortRef, PortCondition}
end
UnitNetwork() = UnitNetwork(Dict{Symbol, MoccaModel}(), Symbol[], AbstractConnection[], Dict{PortRef, PortCondition}())

network(net::UnitNetwork) = net

Base.getindex(fs::AbstractFlowsheet, name::Symbol) = network(fs).units[name]
Base.haskey(fs::AbstractFlowsheet, name::Symbol) = haskey(network(fs).units, name)

"Names of the units and devices of `fs`, in the order they were added."
unit_names(fs::AbstractFlowsheet) = network(fs).order

"The connections of `fs`; see [`AbstractConnection`](@ref)."
connections(fs::AbstractFlowsheet) = network(fs).connections

_is_device(model) = model isa FlowDeviceModel

"Whether unit `name` of `fs` is a [`FlowDevice`](@ref)."
is_device(fs::AbstractFlowsheet, name::Symbol) = _is_device(fs[name])

"Units with holdup (every unit but the devices), in flowsheet order."
holdup_units(fs::AbstractFlowsheet) = [n for n in unit_names(fs) if !is_device(fs, n)]

"""
    add_unit!(fs, name, model)

Add a unit `model` (a `Jutul.SimulationModel` of a [`MoccaSystem`](@ref)) to the
flowsheet under `name`.
"""
function add_unit!(fs::AbstractFlowsheet, name::Symbol, model::MoccaModel)
    net = network(fs)
    !haskey(net.units, name) || error("Flowsheet already has a unit named $name")
    net.units[name] = model
    push!(net.order, name)
    return fs
end

"""
    connect!(fs, connection::AbstractConnection)
    connect!(fs, (unit, port) => (device, device_port))

Add a connection to the flowsheet, after checking it with
[`check_connection`](@ref). The pair form connects `port` of a unit with
holdup to `device_port` (`:inlet` or `:outlet`) of a flow device, as a
[`StreamConnection`](@ref); the pair may be given in either order.
"""
function connect!(fs::AbstractFlowsheet, c::AbstractConnection)
    net = network(fs)
    check_connection(net, c)
    push!(net.connections, c)
    return fs
end

function connect!(fs::AbstractFlowsheet, conn::Pair{Tuple{Symbol, Symbol}, Tuple{Symbol, Symbol}})
    a, b = PortRef(first(conn)), PortRef(last(conn))
    if haskey(fs, a.name) && haskey(fs, b.name) && is_device(fs, a.name) && !is_device(fs, b.name)
        a, b = b, a
    end
    return connect!(fs, StreamConnection(a, b))
end

"""
    set_boundary!(fs, device, port, condition::PortCondition)

Set the default boundary condition for an unconnected device port. Stages can
override it.
"""
function set_boundary!(fs::UnitNetwork, device::Symbol, port::Symbol, condition::PortCondition)
    is_device(fs, device) || error("$device is not a device")
    fs.boundaries[PortRef(device, port)] = condition
    return fs
end

"The default boundary condition of `port` of `device`; see [`set_boundary!`](@ref)."
boundary_condition(fs::AbstractFlowsheet, device::Symbol, port::Symbol) = network(fs).boundaries[PortRef(device, port)]

"Whether `port` of `device` is joined by a connection."
function is_connected(fs::AbstractFlowsheet, device::Symbol, port::Symbol)
    r = PortRef(device, port)
    return any(c -> r in connection_ports(c), connections(fs))
end

"The [`StreamConnection`](@ref)s of `fs`."
stream_connections(fs::AbstractFlowsheet) = [c for c in connections(fs) if c isa StreamConnection]

"Unconnected device ports, as [`PortRef`](@ref)s."
function boundary_ports(fs::AbstractFlowsheet)
    out = PortRef[]
    for name in unit_names(fs)
        is_device(fs, name) || continue
        for p in DEVICE_PORTS
            is_connected(fs, name, p) || push!(out, PortRef(name, p))
        end
    end
    return out
end

"""
    setup_flowsheet_model(fs) → Jutul.MultiModel

Assemble the flowsheet into a `MultiModel`, with the cross terms of every
connection (see [`add_cross_terms!`](@ref)).
"""
function setup_flowsheet_model(fs::AbstractFlowsheet)
    net = network(fs)
    for r in boundary_ports(net)
        haskey(net.boundaries, r) || error("Port $(r.port) of device $(r.name) is not connected and has no boundary condition; use set_boundary!")
    end
    for name in holdup_units(net)
        # Stage duration gradients stretch the time steps of every unit
        # through its TimeScale; a unit without it would silently keep them
        haskey(net.units[name].parameters, :TimeScale) ||
            error("Unit $name has no TimeScale parameter; select it in Jutul.select_parameters! for its system (see balances/unit_parameters.jl)")
    end
    models = NamedTuple{Tuple(net.order)}(Tuple(net.units[k] for k in net.order))
    mm = Jutul.MultiModel(models)
    for c in net.connections
        add_cross_terms!(mm, net, c)
    end
    return mm
end

function _unit_port_condition(fs, unit_state, r::PortRef)
    cell = ports(fs[r.name])[r.port]
    return (pressure = unit_state[:Pressure][cell], temperature = unit_state[:Temperature][cell], composition = unit_state[:y][:, cell])
end

"""
    setup_flowsheet_state(fs; kwarg...) → Dict

Initial state for every unit. Pass the initial state of each unit with holdup as
a keyword named after the unit, e.g. `Bed = setup_unit_state(fs[:Bed]; ...)`.
Device states are initialised from what their ports are attached to.
"""
function setup_flowsheet_state(fs::AbstractFlowsheet; kwarg...)
    states = Dict{Symbol, Any}()
    for name in holdup_units(fs)
        haskey(kwarg, name) || error("Missing initial state for unit $name")
        states[name] = kwarg[name]
    end
    streams = stream_connections(fs)
    for name in unit_names(fs)
        is_device(fs, name) || continue
        conds = map(DEVICE_PORTS) do p
            i = findfirst(c -> c.device == PortRef(name, p), streams)
            isnothing(i) ? boundary_condition(fs, name, p) : _unit_port_condition(fs, states[streams[i].unit.name], streams[i].unit)
        end
        states[name] = setup_device_state(fs[name]; inlet = conds[1], outlet = conds[2])
    end
    return states
end

"""
    setup_flowsheet_parameters(fs; kwarg...) → Dict

Parameters for every unit, from [`setup_unit_parameters`](@ref); pass a `Dict`
or `NamedTuple` of overrides per unit as a keyword named after it.
"""
function setup_flowsheet_parameters(fs::AbstractFlowsheet; kwarg...)
    prm = Dict{Symbol, Any}()
    for name in unit_names(fs)
        extra = get(kwarg, name, NamedTuple())
        prm[name] = setup_unit_parameters(fs[name]; pairs(extra)...)
    end
    return prm
end
