# Ports and the cross terms that connect a unit to a flow device.
#
# A connection joins one port of a unit with holdup (a bed cell) to one port of
# a FlowDevice. It is a pair of additive cross terms:
#
#   PortStateCT  (device ← unit): the device's port copy equals the unit's state
#                                 at the port cell.
#   StreamCT     (unit ← device): the flow through the device enters or leaves
#                                 the unit's mass (and energy) balance at the
#                                 port cell, with upwinded composition and
#                                 temperature.
#
# This is the pattern JutulDarcy uses to couple wells to a reservoir.

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
# target, and finds which target entities (and, for parameter sensitivities,
# which parameters) it depends on by tracing the term once, at the current
# state. So cross terms here avoid branches that change what they read, using
# ifelse on values instead, and a term that does not otherwise depend on the
# target touches it with zero weight through _target_anchor.
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
        out[2 + k] = -DEVICE_FRACTION_SCALE * state_s.y[k, c] + anchor
    end
    return out
end

"""
    StreamCT(cell, port; legacy_inlet = false)

Cross term adding the flow through device port `port` (1 = inlet, 2 = outlet)
to the balances of unit `cell`. Registered with the unit as target and the
device as source, once for the `:mass_conservation` equation and, if the unit
has an energy balance, once for `:energy_column`.

The stream composition, temperature and pressure are those of the upstream
device port, so component `i` enters or leaves at `F·y_s,i` and every
component is conserved across the connection.

`legacy_inlet = true` reproduces the inlet boundary conditions of Mocca 0.1.0,
which add `F·(y_s,i − y_i)` on top of the inflow `F·y_s,i`. That keeps the
total molar flow but not the flow of each component: the unit receives more of
a component than the stream carries when the inlet cell is depleted in it.
"""
struct StreamCT <: Jutul.AdditiveCrossTerm
    cell::Int
    port::Int
    legacy_inlet::Bool
end
StreamCT(cell, port; legacy_inlet = false) = StreamCT(cell, port, legacy_inlet)

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
    F_inflow = ifelse(F_in > 0, F_in, zero(F_in))
    for k in eachindex(out)
        y_s = port_composition(state_s, up, k)
        # Legacy inlet correction, only for flow into the unit
        correction = F_inflow * (y_s - state_t.y[k, c])
        if !ct.legacy_inlet
            correction = _target_anchor(correction)
        end
        out[k] = -F_in * y_s - correction
    end
    return out
end

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t, state_s, state0_s, model_t, model_s,
        ct::StreamCT, eq::Jutul.ConservationLaw{:ColumnConservedEnergy}, dt, ldisc = Jutul.local_discretization(ct, i))
    c = ct.cell
    F_in, up = _stream_into_unit(ct, state_s)
    T_s = port_temperature(state_s, up)
    P_s = port_pressure(state_s, up)
    C_pg = state_t.C_pg[c]
    avm = state_t.AverageMolarMass[c]
    # Enthalpy carried with the pressure-flux term, as in the interior faces
    advective = -F_in * T_s * C_pg * avm
    # Sensible heat of incoming gas relative to the cell, only for inflow
    F_inflow = ifelse(F_in > 0, F_in, zero(F_in))
    q_vol = F_inflow * GAS_CONSTANT * T_s / P_s
    sensible = q_vol * state_t.FluidDensity[1] * C_pg * (T_s - state_t.Temperature[c])
    out[1] = advective - sensible
    return out
end

"""
    port_conductance(model, cell)

Volumetric conductance [m³/(s·Pa)] of the half-cell connection between `cell`
and the end of a bed, `k·A/(Δx/2)/μ`. A [`LinearValve`](@ref) with this
conductance reproduces the pressure boundary conditions of earlier Mocca
versions.
"""
function port_conductance(model::FixedBedModel, cell; parameters = Jutul.setup_parameters(model))
    state = NamedTuple(parameters)
    return calc_bc_trans(model, state, cell) / state.FluidViscosity[1]
end
