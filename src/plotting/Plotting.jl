"""
    Mocca.Plotting

Plots of states, cells and streams.
"""
module Plotting

import Jutul
using CairoMakie
using ..Physics, ..Devices, ..Networks, ..Interface
# Makie exports a Unit of its own, so Mocca's is written Framework.Unit
import ..Framework
using ..Framework: component_names

export plot_cell, plot_cell_comparison, plot_flowsheet_streams, plot_moving_bed_profiles,
    plot_optimization_history, plot_profiles, plot_state, plot_state_comparison

include("plot.jl")

end
