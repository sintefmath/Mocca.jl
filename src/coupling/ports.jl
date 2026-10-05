# Ports, connections between them, and the cross terms of a gas stream.
#
# A connection is an AbstractConnection: it names the ports it joins, checks
# that it fits the units of a flowsheet, and adds cross terms to the
# MultiModel. The gas stream between a port of a unit with holdup (a bed cell)
# and a port of a FlowDevice, StreamConnection, is a pair of additive cross
# terms:
#
#   PortStateCT  (device ← unit): the device's port copy equals the unit's state
#                                 at the port cell.
#   StreamCT     (unit ← device): the flow through the device enters or leaves
#                                 the unit's mass (and energy) balance at the
#                                 port cell, with upwinded composition and
#                                 temperature.
#
# This is the pattern JutulDarcy uses to couple wells to a reservoir. Other
# kinds of connection, such as sorbent passed between moving beds, are further
# AbstractConnection types.

"""
    PortRef(name, port)

Port `port` of the unit or device `name` in a flowsheet, such as
`PortRef(:Bed, :top)`. A `(name, port)` tuple converts to one.
"""
struct PortRef
    name::Symbol
    port::Symbol
end
PortRef(x::Tuple{Symbol, Symbol}) = PortRef(x...)

"""
    AbstractConnection

A connection between ports of the units of a flowsheet. Add one with
[`connect!`](@ref). An implementation provides:

- [`connection_ports`](@ref)`(c)`: the [`PortRef`](@ref)s it joins;
- [`check_connection`](@ref)`(net, c)`: throw an error if `c` cannot join the
  units of the [`UnitNetwork`](@ref) `net`, given the connections it already
  has;
- [`add_cross_terms!`](@ref)`(mm, net, c)`: add its cross terms to the
  `Jutul.MultiModel` `mm` of `net`.

[`StreamConnection`](@ref), the gas stream through a device port, is the one
built in.
"""
abstract type AbstractConnection end

"The [`PortRef`](@ref)s connection `c` joins."
function connection_ports end

"Throw an error if connection `c` cannot join the units of `net`."
function check_connection end

"Add the cross terms of connection `c` to the `MultiModel` `mm` built from `net`."
function add_cross_terms! end

"""
    StreamConnection(unit::PortRef, device::PortRef)

The gas stream between a port of a unit with holdup and a port of a
[`FlowDevice`](@ref). Each device port takes at most one; a unit port may
take several. Its cross terms are a [`PortStateCT`](@ref) on the device and a
[`StreamCT`](@ref) for each of the unit's [`stream_equations`](@ref).
"""
struct StreamConnection <: AbstractConnection
    unit::PortRef
    device::PortRef
end

connection_ports(c::StreamConnection) = (c.unit, c.device)

function check_connection(net, c::StreamConnection)
    a, b = c.unit, c.device
    for r in (a, b)
        haskey(net.units, r.name) || error("Unknown unit $(r.name)")
    end
    !_is_device(net.units[a.name]) || error("Cannot connect two devices ($(a.name) and $(b.name)) directly; put a unit with holdup between them")
    _is_device(net.units[b.name]) || error("Cannot connect two units with holdup ($(a.name) and $(b.name)) directly; connect them through a device")
    haskey(ports(net.units[a.name]), a.port) || error("Unit $(a.name) has no port $(a.port), available: $(keys(ports(net.units[a.name])))")
    haskey(ports(net.units[b.name]), b.port) || error("Device $(b.name) has no port $(b.port), available: $(keys(ports(net.units[b.name])))")
    any(x -> x isa StreamConnection && x.device == b, net.connections) && error("Port $(b.port) of device $(b.name) is already connected")
    return nothing
end

function add_cross_terms!(mm, net, c::StreamConnection)
    u, d = c.unit.name, c.device.name
    unit = net.units[u]
    cell = ports(unit)[c.unit.port]
    k = port_index(net.units[d], c.device.port)
    Jutul.add_cross_term!(mm, PortStateCT(cell); target = d, source = u, equation = port_equation(c.device.port))
    for eq in stream_equations(unit)
        Jutul.add_cross_term!(mm, StreamCT(cell, k); target = u, source = d, equation = eq)
    end
    return mm
end

"""
    ports(model) → NamedTuple

Named ports of a unit. For a distributed unit each port is the index of the cell
at that end; for a device each port is its index (1 = inlet, 2 = outlet).
"""
function ports end

function ports(model::DistributedUnitModel)
    nc = Jutul.number_of_cells(model.domain)
    return (bottom = 1, top = nc)
end

ports(model::FlowDeviceModel) = (inlet = 1, outlet = 2)

# Jutul adds a cross term's residual through its derivatives with respect to the
# target, and finds which entities (cells, or the unit itself) and, for
# parameter sensitivities, which parameters it depends on by tracing the term
# once, at the current state. Dependencies found then are kept for the whole
# simulation; any others are dropped. So:
#
# - A branch must not change which entities or parameters a term reads, or the
#   trace sees only the branch taken at that state. Use ifelse, which evaluates
#   both sides. Branching between variables of the same entity is fine, such as
#   choosing the upstream port of a device, whose variables all live in its one
#   cell.
# - A term that does not otherwise depend on the target touches it with zero
#   weight through _target_anchor, or its residual is never added.
@inline _target_anchor(x) = zero(x) * x

"""
    PortStateCT(cell)

Cross term setting a device port copy to the state of unit `cell`. Registered
with the device as target and the unit as source, on the device's
`:inlet_state` or `:outlet_state` equation.
"""
struct PortStateCT <: Jutul.AdditiveCrossTerm
    cell::Int
end

Jutul.cross_term_entities(ct::PortStateCT, eq::Jutul.JutulEquation, model) = [1]
Jutul.cross_term_entities_source(ct::PortStateCT, eq::Jutul.JutulEquation, model) = [ct.cell]

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t, state_s, state0_s, model_t, model_s,
        ct::PortStateCT, eq, dt, ldisc = Jutul.local_discretization(ct, i))
    c = ct.cell
    N = number_of_components(model_t.system)
    anchor = _target_anchor(state_t.MolarFlow[1])
    out[1] = -state_s.Pressure[c] + anchor
    out[2] = -state_s.Temperature[c] + anchor
    for k in 1:(N - 1)
        out[2 + k] = -state_s.y[k, c] + anchor
    end
    return out
end

"""
    StreamCT(cell, port)

Cross term adding the flow through device port `port` (1 = inlet, 2 = outlet)
to the balances of unit `cell`. Registered with the unit as target and the
device as source, once for the `:mass_conservation` equation and, if the unit
has an energy balance, once for `:energy_column`.

The stream composition, temperature and pressure are those of the upstream
device port, so component `i` enters or leaves at `F·y_s,i` and every
component is conserved across the connection.
"""
struct StreamCT <: Jutul.AdditiveCrossTerm
    cell::Int
    port::Int
end

Jutul.cross_term_entities(ct::StreamCT, eq::Jutul.JutulEquation, model) = [ct.cell]
Jutul.cross_term_entities_source(ct::StreamCT, eq::Jutul.JutulEquation, model) = [1]

# Molar flow into the unit, and the index of the upstream device port
@inline function _stream_into_unit(ct::StreamCT, device_state)
    F = device_state.MolarFlow[1]
    F_in = ct.port == 2 ? F : -F
    up = F >= 0 ? 1 : 2
    return (F_in, up)
end

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t, state_s, state0_s, model_t, model_s,
        ct::StreamCT, eq::Jutul.ConservationLaw{:ComponentMasses}, dt, ldisc = Jutul.local_discretization(ct, i))
    c = ct.cell
    F_in, up = _stream_into_unit(ct, state_s)
    for k in eachindex(out)
        y_s = port_composition(state_s, up, k)
        # The flow does not depend on the unit's composition, so anchor to it
        # to keep the sparsity pattern
        out[k] = -F_in * y_s + _target_anchor(state_t.y[k, c])
    end
    return out
end

# Molar heat capacity [J/(mol·K)] of the gas at device port `k`, with the
# mixing rule of the energy balance (`C_pg·M̄`, both mole-fraction averages)
@inline function _port_molar_heat_capacity(sys, device_state, k)
    C, M = sys.heat_capacity_gas, sys.molecular_masses
    cp = zero(port_temperature(device_state, k))
    mm = zero(cp)
    for i in eachindex(C)
        y = port_composition(device_state, k, i)
        cp += y * C[i]
        mm += y * M[i]
    end
    return cp * mm
end

# The stream carries the enthalpy of the upstream port, `F·c_p(y)·T`, as the
# interior faces carry that of the upstream cell. Both ends of a device read the
# same upstream port, so what one unit loses the other gains.
function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t, state_s, state0_s, model_t, model_s,
        ct::StreamCT, eq::Jutul.ConservationLaw{:ColumnConservedEnergy}, dt, ldisc = Jutul.local_discretization(ct, i))
    c = ct.cell
    F_in, up = _stream_into_unit(ct, state_s)
    h = _port_molar_heat_capacity(model_t.system, state_s, up) * port_temperature(state_s, up)
    # The enthalpy does not depend on the unit's state, so anchor to it to keep
    # the sparsity pattern
    out[1] = -F_in * h + _target_anchor(state_t.Temperature[c])
    return out
end

"""
    port_conductance(model, cell)

Volumetric conductance [m³/(s·Pa)] of the half-cell connection between `cell`
and the end of a distributed unit, `k·A/(Δx/2)/μ`. A [`LinearValve`](@ref) with this
conductance reproduces the pressure boundary conditions of earlier Mocca
versions.
"""
function port_conductance(model::DistributedUnitModel, cell; parameters = Jutul.setup_parameters(model))
    state = NamedTuple(parameters)
    return calc_bc_trans(model, state, cell) / state.FluidViscosity[1]
end
