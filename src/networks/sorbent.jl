# Connections between units that are not driven by pressure: sorbent moving
# from one moving bed to the next, and heat transfer fluid passing from the
# tubes of one unit to the next. Their flow rates are parameters of the units
# (SorbentMassFlow, HXFluidMassFlow), so no device sits between them.
#
# The receiving unit's own equations take the incoming stream from its inlet
# parameters (SolidInletLoading, SolidInletTemperature and
# SolidInletGasConcentration, or HXFluidInletTemperature). The cross terms swap
# that for the state leaving the source unit: each adds the inlet term with the
# parameter value removed and the source value put in. The inlet parameters
# then have no effect.

"""
    SorbentConnection(upstream, downstream)

Sorbent leaving the bottom of unit `upstream` enters the top of unit
`downstream`, carrying its loading, its temperature and the gas in its pores. Both units must move their
sorbent ([`moves_sorbent`](@ref)), such as a [`MovingBed`](@ref), and have the
same `SorbentMassFlow`; the `SolidInletLoading`, `SolidInletTemperature` and
`SolidInletGasConcentration` of `downstream` are then not used. Each unit passes its sorbent to one other and
takes it from one other. Closed loops, where the sorbent returns to the first
unit, are allowed. [`connect_solids!`](@ref) adds one.
"""
struct SorbentConnection <: AbstractConnection
    upstream::Symbol
    downstream::Symbol
end

connection_ports(c::SorbentConnection) = (PortRef(c.upstream, :sorbent_outlet), PortRef(c.downstream, :sorbent_inlet))

function check_connection(net, c::SorbentConnection)
    a, b = c.upstream, c.downstream
    for u in (a, b)
        haskey(net.units, u) || error("Unknown unit $u")
        moves_sorbent(net.units[u].system) || error("Unit $u does not move its sorbent; only units such as moving beds pass sorbent between them")
    end
    a != b || error("Cannot feed the sorbent of $a back into itself")
    for x in net.connections
        x isa SorbentConnection || continue
        x.upstream == a && error("The sorbent leaving $a already goes to $(x.downstream)")
        x.downstream == b && error("The sorbent entering $b already comes from $(x.upstream)")
    end
    _check_linked_flow(net, a, b, :sorbent_mass_flow, "sorbent")
    return nothing
end

function add_cross_terms!(mm, net, c::SorbentConnection)
    a, b = c.upstream, c.downstream
    ct = SolidStreamCT(ports(net.units[a]).bottom, ports(net.units[b]).top)
    Jutul.add_cross_term!(mm, ct; target = b, source = a, equation = :mass_transfer)
    # The gas in the pores of the sorbent
    Jutul.add_cross_term!(mm, ct; target = b, source = a, equation = :mass_conservation)
    if has_energy_balance(net.units[b].system)
        Jutul.add_cross_term!(mm, ct; target = b, source = a, equation = :energy_column)
    end
    return mm
end

"""
    connect_solids!(fs, upstream => downstream)

Add a [`SorbentConnection`](@ref): the sorbent leaving the bottom of
`upstream` enters the top of `downstream`.
"""
connect_solids!(fs::AbstractFlowsheet, conn::Pair{Symbol, Symbol}) = connect!(fs, SorbentConnection(first(conn), last(conn)))

"""
    sorbent_path(fs[, start]) → Vector{Symbol}

The units the sorbent passes through in turn, following the
[`SorbentConnection`](@ref)s of `fs` from `start`. By default it starts at the
unit that takes its sorbent from no other, or, in a closed loop, at the first
unit of `fs` that moves its sorbent.
"""
function sorbent_path(fs::AbstractFlowsheet, start = nothing)
    links = Dict(c.upstream => c.downstream for c in connections(fs) if c isa SorbentConnection)
    if isnothing(start)
        moving = [u for u in holdup_units(fs) if moves_sorbent(fs[u].system)]
        isempty(moving) && error("The flowsheet has no units whose sorbent moves")
        fed = Set(values(links))
        i = findfirst(u -> !(u in fed), moving)
        start = moving[something(i, 1)]
    end
    path = [Symbol(start)]
    while haskey(links, path[end]) && links[path[end]] != path[1]
        push!(path, links[path[end]])
    end
    return path
end

"""
    HeatTransferFluidConnection(upstream, downstream)

Heat transfer fluid leaving the tubes of unit `upstream` (at its top) enters
the tubes of unit `downstream` (at its bottom). Both units need a
[`HeatExchanger`](@ref) and the same `HXFluidMassFlow`; the
`HXFluidInletTemperature` of `downstream` is then not used. Joining two units
both ways makes a closed loop for heat recovery.
[`connect_heat_transfer_fluid!`](@ref) adds one.
"""
struct HeatTransferFluidConnection <: AbstractConnection
    upstream::Symbol
    downstream::Symbol
end

connection_ports(c::HeatTransferFluidConnection) = (PortRef(c.upstream, :fluid_outlet), PortRef(c.downstream, :fluid_inlet))

function check_connection(net, c::HeatTransferFluidConnection)
    a, b = c.upstream, c.downstream
    for u in (a, b)
        haskey(net.units, u) || error("Unknown unit $u")
        has_heat_exchanger(net.units[u].system) || error("Unit $u has no heat exchanger")
    end
    a != b || error("Cannot feed the heat transfer fluid of $a back into itself")
    for x in net.connections
        x isa HeatTransferFluidConnection || continue
        x.upstream == a && error("The heat transfer fluid leaving $a already goes to $(x.downstream)")
        x.downstream == b && error("The heat transfer fluid entering $b already comes from $(x.upstream)")
    end
    _check_linked_flow(net, a, b, :hx_fluid_mass_flow, "heat transfer fluid")
    return nothing
end

function add_cross_terms!(mm, net, c::HeatTransferFluidConnection)
    a, b = c.upstream, c.downstream
    # The fluid flows up through the tubes, out of the top of a and into the bottom of b
    ct = HXFluidStreamCT(ports(net.units[a]).top, ports(net.units[b]).bottom)
    Jutul.add_cross_term!(mm, ct; target = b, source = a, equation = :energy_hx_fluid)
    return mm
end

"""
    connect_heat_transfer_fluid!(fs, upstream => downstream)

Add a [`HeatTransferFluidConnection`](@ref): the heat transfer fluid leaving
the tubes of `upstream` enters those of `downstream`.
"""
connect_heat_transfer_fluid!(fs::AbstractFlowsheet, conn::Pair{Symbol, Symbol}) =
    connect!(fs, HeatTransferFluidConnection(first(conn), last(conn)))

# Units passing a stream to each other must agree on its flow rate. This checks
# the values the units were built with; parameters overriding them must keep
# them equal.
function _check_linked_flow(net, a, b, key, what)
    Fa, Fb = (first(net.units[u].data_domain[key, Unit()]) for u in (a, b))
    isapprox(Fa, Fb, rtol = 1e-10) || error("$a passes $what to $b, but their mass flows differ ($Fa and $Fb kg/s)")
    return nothing
end

"""
    SolidStreamCT(source_cell, target_cell)

Cross term feeding the sorbent leaving `source_cell` (the bottom of the source
moving bed) into `target_cell` (the top of the target). Registered with the
receiving unit as target, once each for its `:mass_transfer` and
`:mass_conservation` equations (the loading, and the gas in the pores) and,
if it has an energy balance, once for `:energy_column`. See [`SorbentConnection`](@ref).
"""
struct SolidStreamCT <: Jutul.AdditiveCrossTerm
    source_cell::Int
    target_cell::Int
end

Jutul.cross_term_entities(ct::SolidStreamCT, eq::Jutul.JutulEquation, model) = [ct.target_cell]
Jutul.cross_term_entities_source(ct::SolidStreamCT, eq::Jutul.JutulEquation, model) = [ct.source_cell]

# The loading equation has Q_s·(q − q_in); this makes it Q_s·(q − q_source)
function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t, state_s, state0_s, model_t, model_s,
        ct::SolidStreamCT, eq::Jutul.ConservationLaw{:AdsorbedLoading}, dt, ldisc = Jutul.local_discretization(ct, i))
    c, s = ct.target_cell, ct.source_cell
    Q_s = sorbent_volume_flow(state_t)
    for k in eachindex(out)
        q_param = state_t.SolidInletLoading[k, 1]
        out[k] = Q_s * (q_param - state_s.AdsorbedLoading[k, s]) + _target_anchor(state_t.AdsorbedConcentration[k, c])
    end
    return out
end

# The gas balance has −Q_p·(c_in − c) from PoreGasFlow; this makes it −Q_p·(c_source − c)
function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t, state_s, state0_s, model_t, model_s,
        ct::SolidStreamCT, eq::Jutul.ConservationLaw{:ComponentMasses}, dt, ldisc = Jutul.local_discretization(ct, i))
    c, s = ct.target_cell, ct.source_cell
    Q_p = pore_gas_volume_flow(state_t)
    for k in eachindex(out)
        out[k] = Q_p * (state_t.SolidInletGasConcentration[k, 1] - state_s.MolarConcentration[k, s]) + _target_anchor(state_t.y[k, c])
    end
    return out
end

# The energy balance has −C·(T_in − T) from SorbentHeatFlow and −Q_p·e_in from
# PoreGasFlow; this swaps the inlet for the source's sorbent temperature and
# pore gas
function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t, state_s, state0_s, model_t, model_s,
        ct::SolidStreamCT, eq::Jutul.ConservationLaw{:ColumnConservedEnergy}, dt, ldisc = Jutul.local_discretization(ct, i))
    c, s = ct.target_cell, ct.source_cell
    sys = model_t.system
    N = number_of_components(sys)
    C = sorbent_heat_capacity_flow(state_t, c)
    c_in = state_t.SolidInletGasConcentration
    e_in = pore_gas_energy(sys, state_t, c, ntuple(k -> c_in[k, 1], N), state_t.SolidInletTemperature[1])
    e_s = pore_gas_energy(sys, state_t, c, ntuple(k -> state_s.MolarConcentration[k, s], N), state_s.Temperature[s])
    out[1] = C * (state_t.SolidInletTemperature[1] - state_s.Temperature[s]) +
        pore_gas_volume_flow(state_t) * (e_in - e_s) + _target_anchor(state_t.Temperature[c])
    return out
end

"""
    HXFluidStreamCT(source_cell, target_cell)

Cross term feeding the heat transfer fluid leaving `source_cell` (the top of the
source unit's tubes) into `target_cell` (the bottom of the target's). Registered
with the receiving unit as target on its `:energy_hx_fluid` equation. See
[`HeatTransferFluidConnection`](@ref).
"""
struct HXFluidStreamCT <: Jutul.AdditiveCrossTerm
    source_cell::Int
    target_cell::Int
end

Jutul.cross_term_entities(ct::HXFluidStreamCT, eq::Jutul.JutulEquation, model) = [ct.target_cell]
Jutul.cross_term_entities_source(ct::HXFluidStreamCT, eq::Jutul.JutulEquation, model) = [ct.source_cell]

function Jutul.update_cross_term_in_entity!(out, i,
        state_t, state0_t, state_s, state0_s, model_t, model_s,
        ct::HXFluidStreamCT, eq, dt, ldisc = Jutul.local_discretization(ct, i))
    c, s = ct.target_cell, ct.source_cell
    C = hx_fluid_heat_capacity_flow(state_t)
    out[1] = C * (state_t.HXFluidInletTemperature[1] - state_s.HXFluidTemperature[s]) + _target_anchor(state_t.HXFluidTemperature[c])
    return out
end
