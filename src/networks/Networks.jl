"""
    Mocca.Networks

Ports, connections and flowsheets of units joined into one Jutul MultiModel.
"""
module Networks

import Jutul
using ..Framework, ..Physics, ..Balances, ..Devices

export AbstractConnection, AbstractFlowsheet, add_cross_terms!, add_unit!,
    boundary_condition, boundary_ports, check_connection, connect!, connection_ports,
    connections, holdup_units, is_connected, is_device, network, port_conductance, PortRef,
    PortStateCT, set_boundary!, setup_flowsheet_model, setup_flowsheet_parameters,
    setup_flowsheet_state, stream_connections, StreamConnection, StreamCT, unit_names,
    UnitNetwork

include("ports.jl")
include("flowsheet.jl")

end
