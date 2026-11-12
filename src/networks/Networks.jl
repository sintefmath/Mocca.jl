"""
    Mocca.Networks

Ports, connections and flowsheets of units joined into one Jutul MultiModel:
gas streams through devices, and sorbent and heat transfer fluid passed
between units.
"""
module Networks

import Jutul
using ..Framework, ..Physics, ..Balances, ..Devices

export AbstractConnection, AbstractFlowsheet, add_cross_terms!, add_unit!,
    boundary_condition, boundary_ports, check_connection, connect!,
    connect_heat_transfer_fluid!, connect_solids!, connection_ports, connections,
    HeatTransferFluidConnection, holdup_units, HXFluidStreamCT, is_connected, is_device,
    network, port_conductance, PortRef, PortStateCT, set_boundary!, setup_flowsheet_model,
    setup_flowsheet_parameters, setup_flowsheet_state, SolidStreamCT, sorbent_path, SorbentConnection,
    stream_connections, StreamConnection, StreamCT, unit_names, UnitNetwork

include("ports.jl")
include("flowsheet.jl")
include("sorbent.jl")

end
