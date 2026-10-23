# Ports, connections between them, and the cross terms of a gas stream.
#
# A connection is an AbstractConnection: it names the ports it joins, checks
# that it fits the units of a flowsheet, and adds cross terms to the
# MultiModel. The gas stream between a port of a unit with holdup (the cells at
# one end of a bed) and a port of a FlowDevice, StreamConnection, is a pair of
# additive cross terms:
#
#   PortStateCT  (device ← unit): the device's port copy equals the unit's state
#                                 at the port cells, mixed over them.
#   StreamCT     (unit ← device): the flow through the device enters or leaves
#                                 the unit's mass (and energy) balance at the
#                                 port cells, each taking its share, with
#                                 upwinded composition and temperature.
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
    cells = port_cells(unit, c.unit.port)
    w = port_weights(unit, cells)
    k = port_index(net.units[d], c.device.port)
    Jutul.add_cross_term!(mm, PortStateCT(cells, w); target = d, source = u, equation = port_equation(c.device.port))
    for eq in stream_equations(unit)
        Jutul.add_cross_term!(mm, StreamCT(cells, w, k); target = u, source = d, equation = eq)
    end
    return mm
end

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
    PortStateCT(cells, weights)
    PortStateCT(cell)

Cross term setting a device port copy to the state of a unit at its port
`cells`. Registered with the device as target and the unit as source, on the
device's `:inlet_state` or `:outlet_state` equation.

Over several cells, the state is that of the gas they give when mixed in the
shares `weights` (see [`port_weights`](@ref)): the weighted mean of the
pressure and the mole fractions, and the mixing-cup temperature, at which the
mixed stream carries the enthalpy of its parts. A unit without an energy
balance gives the weighted mean temperature.
"""
struct PortStateCT <: Jutul.AdditiveCrossTerm
    cells::Vector{Int}
    weights::Vector{Float64}
end
PortStateCT(cell::Integer) = PortStateCT([Int(cell)], [1.0])

Jutul.cross_term_entities(ct::PortStateCT, eq::Jutul.JutulEquation, model) = [1]
Jutul.cross_term_entities_source(ct::PortStateCT, eq::Jutul.JutulEquation, model) = ct.cells

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t, state_s, state0_s, model_t, model_s,
        ct::PortStateCT, eq, dt, ldisc = Jutul.local_discretization(ct, i))
    N = number_of_components(model_t.system)
    anchor = _target_anchor(state_t.MolarFlow[1])
    if length(ct.cells) == 1
        c = only(ct.cells)
        out[1] = -state_s.Pressure[c] + anchor
        out[2] = -state_s.Temperature[c] + anchor
        for k in 1:(N - 1)
            out[2 + k] = -state_s.y[k, c] + anchor
        end
    else
        cells, w = ct.cells, ct.weights
        out[1] = -_weighted_mean(state_s.Pressure, cells, w) + anchor
        out[2] = -_port_temperature(model_s, state_s, cells, w) + anchor
        for k in 1:(N - 1)
            out[2 + k] = -_weighted_mean(view(state_s.y, k, :), cells, w) + anchor
        end
    end
    return out
end

@inline function _weighted_mean(v, cells, w)
    s = zero(eltype(v))
    for j in eachindex(cells)
        s += w[j] * v[cells[j]]
    end
    return s
end

# Molar heat capacity [J/(mol·K)] of the gas of unit `cell`, with the mixing
# rule of the energy balance (`C_pg·M̄`, both mole-fraction averages)
@inline function _cell_molar_heat_capacity(sys, y, cell)
    C, M = sys.heat_capacity_gas, sys.molecular_masses
    cp = zero(eltype(y))
    mm = zero(cp)
    for i in eachindex(C)
        cp += y[i, cell] * C[i]
        mm += y[i, cell] * M[i]
    end
    return cp * mm
end

# Temperature of the gas of `cells` mixed in the shares `w`: the mixing-cup
# temperature Σ w_j·c_p(y_j)·T_j / c_p(ȳ) of a unit with an energy balance, so
# that the stream's enthalpy F·c_p(ȳ)·T is what the cells give up
function _port_temperature(model, state, cells, w)
    has_energy_balance(model.system) || return _weighted_mean(state.Temperature, cells, w)
    sys = model.system
    C, M = sys.heat_capacity_gas, sys.molecular_masses
    H = zero(eltype(state.Temperature))
    for j in eachindex(cells)
        H += w[j] * _cell_molar_heat_capacity(sys, state.y, cells[j]) * state.Temperature[cells[j]]
    end
    cp = zero(eltype(state.y))
    mm = zero(cp)
    for i in eachindex(C)
        y_i = _weighted_mean(view(state.y, i, :), cells, w)
        cp += y_i * C[i]
        mm += y_i * M[i]
    end
    return H / (cp * mm)
end

"""
    StreamCT(cells, weights, port)
    StreamCT(cell, port)

Cross term adding the flow through device port `port` (1 = inlet, 2 = outlet)
to the balances of the unit at its port `cells`, each taking its share
`weights` of the flow (see [`port_weights`](@ref)). Registered with the unit
as target and the device as source, once for the `:mass_conservation`
equation and, if the unit has an energy balance, once for `:energy_column`.

A stream into the unit has the composition and temperature of the upstream
device port. A stream out of a port of one cell has those of the device's
copy of the cell, which equal the cell's once converged. Out of a port of
several cells, each cell gives up its own gas, so that a cell does not lose
more of a component than it holds; the device's copy is these mixed (see
[`PortStateCT`](@ref)). Either way every component and the enthalpy are
conserved across the connection.
"""
struct StreamCT <: Jutul.AdditiveCrossTerm
    cells::Vector{Int}
    weights::Vector{Float64}
    port::Int
end
StreamCT(cell::Integer, port::Integer) = StreamCT([Int(cell)], [1.0], Int(port))

Jutul.cross_term_entities(ct::StreamCT, eq::Jutul.JutulEquation, model) = ct.cells
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
    c, w = ct.cells[i], ct.weights[i]
    F_in, up = _stream_into_unit(ct, state_s)
    for k in eachindex(out)
        y_s = port_composition(state_s, up, k)
        if length(ct.cells) > 1
            y_c = state_t.y[k, c]
            y_s = ifelse(F_in >= 0, y_s + zero(y_c), y_c + zero(y_s))
        end
        # The flow does not depend on the unit's composition, so anchor to it
        # to keep the sparsity pattern
        out[k] = -w * F_in * y_s + _target_anchor(state_t.y[k, c])
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
    c, w = ct.cells[i], ct.weights[i]
    F_in, up = _stream_into_unit(ct, state_s)
    h = _port_molar_heat_capacity(model_t.system, state_s, up) * port_temperature(state_s, up)
    if length(ct.cells) > 1
        h_c = _cell_molar_heat_capacity(model_t.system, state_t.y, c) * state_t.Temperature[c]
        h = ifelse(F_in >= 0, h + zero(h_c), h_c + zero(h))
    end
    # The enthalpy does not depend on the unit's state, so anchor to it to keep
    # the sparsity pattern
    out[1] = -w * F_in * h + _target_anchor(state_t.Temperature[c])
    return out
end

"""
    port_conductance(model, cells)
    port_conductance(model, port::Symbol)

Volumetric conductance [m³/(s·Pa)] of the half-cell connections between the
end of a distributed unit and `cells`, the cells of a port, or one of them:
`k·A/(Δx/2)/μ` over the cross-section `A`, each cell counting with its share
of it (see [`port_weights`](@ref)). A [`LinearValve`](@ref) with this
conductance reproduces the pressure boundary conditions of earlier Mocca
versions.
"""
function port_conductance(model::DistributedUnitModel, cells; parameters = Jutul.setup_parameters(model))
    state = NamedTuple(parameters)
    cells = _port_cells(cells)
    w = port_weights(model, cells)
    return sum(w[j] * calc_bc_trans(model, state, cells[j]) for j in eachindex(cells)) / state.FluidViscosity[1]
end

port_conductance(model::DistributedUnitModel, port::Symbol; kwarg...) = port_conductance(model, port_cells(model, port); kwarg...)
