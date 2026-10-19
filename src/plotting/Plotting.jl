"""
    Mocca.Plotting

Plots of states, cells and streams.
"""
module Plotting

import Jutul
using CairoMakie
using ..Physics, ..Devices, ..Interface

export plot_cell, plot_cell_comparison, plot_flowsheet_streams, plot_optimization_history,
    plot_profiles, plot_state, plot_state_comparison

include("plot.jl")

end
