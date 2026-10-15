"""
    Mocca.Devices

Flow devices: one-cell models that own the flow between two ports.
"""
module Devices

import Jutul
using ..Framework

export AbstractDeviceLaw, Closed, device_flow, DEVICE_PORTS, DeviceMolarFlow,
    evaluate_profile, ExponentialRamp, flow_scale, FlowDevice, FlowDeviceModel,
    law_flow_scale, LinearValve, MIN_FLOW_SCALE, port_composition, port_equation, port_index,
    port_temperature, PortCondition, PortSetpoint, setup_device_model, setup_device_state,
    shift_profile, VolumetricFlow

include("device.jl")

end
