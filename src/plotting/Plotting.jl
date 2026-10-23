"""
    Mocca.Plotting

Plots of states, cells and streams.
"""
module Plotting

import Jutul
using CairoMakie
using ..Framework, ..Physics, ..Devices, ..Interface
using ..Balances: _is_axisymmetric

export plot_bed_section, plot_cell, plot_cell_comparison, plot_flowsheet_streams,
    plot_optimization_history, plot_profiles, plot_radial_profiles, plot_state,
    plot_state_comparison

include("plot.jl")

end
