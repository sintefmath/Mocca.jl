# What the flowsheet and process layers need from a unit with holdup.
#
# A unit is a Jutul.SimulationModel of a MoccaSystem. To be part of a
# flowsheet, its system provides `component_names`, and its model provides:
#
# - `ports(model)`: named ports, each the cell, or the cells, a stream enters
#   or leaves through (see networks/ports.jl), and, for a port of several
#   cells, [`port_weights`](@ref) if equal shares of the stream do not fit;
# - the variables the port cross terms read at a port cell: `Pressure`,
#   `Temperature` and `y`; and, for streams into an energy balance, the
#   system's `heat_capacity_gas` and `molecular_masses`;
# - the `TimeScale` parameter, used through `_scaled_dt` in every time
#   derivative, for gradients with respect to stage durations;
# - [`stream_equations`](@ref): the equations a stream adds to;
# - [`setup_unit_state`](@ref) and [`setup_unit_parameters`](@ref), if the
#   defaults below do not fit.
#
# A DistributedUnit gets the gas balance for free (balances/gas.jl): its face flux
# reads the `Transmissibilities`, `DispersionTransmissibilities`,
# `AxialDispersion` and `FluidViscosity` parameters. A unit that also
# [`holds_sorbent`](@ref) counts as a bed: steps can act on it and its sorbent
# counts in the productivity of a process.

"""
    ports(model) → NamedTuple

Named ports of a unit. For a distributed unit each port is the vector of cells
at that end, the first and last of its [`axial_layers`](@ref); a unit may also
give a port as the index of one cell. For a device each port is its index
(1 = inlet, 2 = outlet).
"""
function ports end

function ports(model::DistributedUnitModel)
    layers = axial_layers(model)
    return (bottom = first(layers), top = last(layers))
end

"""
    port_cells(model, port) → Vector{Int}

The cells of `port` of a unit (see [`ports`](@ref)).
"""
port_cells(model, port::Symbol) = _port_cells(ports(model)[port])
_port_cells(cell::Integer) = [Int(cell)]
_port_cells(cells) = collect(Int, cells)

"""
    axial_layers(model) → Vector{Vector{Int}}

The cells of a distributed unit in layers across its axis, from the `bottom`
port to the `top`: layer `i` holds the cells at the `i`-th position along the
flow. On a Cartesian mesh the axis is x, so a column of `n` cells has `n`
layers of one cell. On any other mesh each cell is a layer of its own, in
order.
"""
axial_layers(model::DistributedUnitModel) = _axial_layers(Jutul.physical_representation(model.data_domain))

function _axial_layers(g::Jutul.CartesianMesh)
    nx, ny, nz = Jutul.grid_dims_ijk(g)
    return [[i + (j - 1) * nx + (k - 1) * nx * ny for k in 1:nz for j in 1:ny] for i in 1:nx]
end
_axial_layers(g) = [[c] for c in 1:Jutul.number_of_cells(g)]

# The number of axial layers, for use where building them would cost too much,
# such as in an equation; cell `c` is in layer mod1(c, _axial_count(model))
_axial_count(model::DistributedUnitModel) = _axial_count(Jutul.physical_representation(model.data_domain))
_axial_count(g::Jutul.CartesianMesh) = first(Jutul.grid_dims_ijk(g))
_axial_count(g) = Jutul.number_of_cells(g)

"""
    port_weights(model, cells) → Vector{Float64}

The share of each of `cells`, such as the cells of a port, in a stream through
them or in an average over them, adding up to one. In a distributed unit it is
the cell's share of the volume, which for the cells of one layer is its share
of the flow cross-section; in any other unit the shares are equal.
"""
port_weights(model::MoccaModel, cells) = fill(1 / length(cells), length(cells))

function port_weights(model::DistributedUnitModel, cells)
    V = model.data_domain[:volumes][cells]
    return V ./ sum(V)
end

"""
    cross_section_mean(model, values, cells)

The average of a cell variable over `cells` with their [`port_weights`](@ref):
a number for a variable with one value per cell, and a vector for one with a
column per cell, such as the mole fractions. For one cell it is that cell's
value.
"""
function cross_section_mean(model, values::AbstractVector, cells)
    length(cells) == 1 && return values[only(cells)]
    w = port_weights(model, cells)
    return sum(w[j] * values[c] for (j, c) in enumerate(cells))
end

function cross_section_mean(model, values::AbstractMatrix, cells)
    length(cells) == 1 && return values[:, only(cells)]
    w = port_weights(model, cells)
    return sum(w[j] .* values[:, c] for (j, c) in enumerate(cells))
end

"""
    axial_profile(model, values)

The cross-section mean (see [`cross_section_mean`](@ref)) of a cell variable of
a distributed unit in each of its [`axial_layers`](@ref): a vector, or a matrix
with one column per layer for a variable with a column per cell. Plot it
against [`axial_positions`](@ref).
"""
axial_profile(model, values::AbstractVector) = [cross_section_mean(model, values, l) for l in axial_layers(model)]
axial_profile(model, values::AbstractMatrix) = reduce(hcat, [cross_section_mean(model, values, l) for l in axial_layers(model)])

"""
    axial_positions(model) → Vector{Float64}

The position [m] along the axis of each of the [`axial_layers`](@ref) of a
distributed unit, at the centres of its cells.
"""
axial_positions(model) = [model.data_domain[:cell_centroids][1, first(l)] for l in axial_layers(model)]

"""
    port_state(model, state, port) → NamedTuple

The `pressure`, `temperature` and `composition` of a unit at `port`: those of
its port cell, or their [`cross_section_mean`](@ref) over a port of several
cells.
"""
function port_state(model, state, port::Symbol)
    cells = port_cells(model, port)
    return (pressure = cross_section_mean(model, state[:Pressure], cells),
        temperature = cross_section_mean(model, state[:Temperature], cells),
        composition = cross_section_mean(model, state[:y], cells))
end

"""
    stream_equations(model) → Tuple of Symbols

Equations of `model` that a stream through one of its ports adds to: the mass
balance, and the energy balance if the unit solves one.
"""
function stream_equations end

function stream_equations(model::DistributedUnitModel)
    return has_energy_balance(model.system) ? (:mass_conservation, :energy_column) : (:mass_conservation,)
end

function stream_equations(model::MoccaModel)
    error("$(nameof(typeof(model.system))) does not define Mocca.stream_equations, so streams cannot connect to it")
end

# Typed on MoccaModel rather than Jutul.SimulationModel: Julia ranks a plain
# SimulationModel above SimulationModel{<:Any, <:SorbentBed, <:Any, <:Any}, so
# the method for beds would never be called
"""
    setup_unit_parameters(model; overrides...) → Dict

Parameters of a unit, with `overrides` for any of them.
"""
setup_unit_parameters(model::MoccaModel; kwarg...) = Jutul.setup_parameters(model; kwarg...)

"""
    setup_unit_state(model; pressure, temperature, composition,
        wall_temperature = temperature) → Dict

Uniform initial state of a unit with holdup. The default sets the gas state,
and the wall temperature if the unit has a wall; a unit with more primary
variables adds a method.
"""
function setup_unit_state(model::MoccaModel; pressure, temperature, composition, wall_temperature = temperature)
    return Jutul.setup_state(model; _gas_state(model, pressure, temperature, composition, wall_temperature)...)
end

function _gas_state(model, pressure, temperature, composition, wall_temperature)
    kw = Dict{Symbol, Any}(:Pressure => pressure, :Temperature => temperature, :y => collect(composition))
    haskey(model.primary_variables, :WallTemperature) && (kw[:WallTemperature] = wall_temperature)
    return kw
end

"""
    holds_sorbent(model) → Bool

Whether `model` is a bed: a unit holding sorbent, which steps act on and whose
sorbent counts in the productivity of a process.
"""
holds_sorbent(model) = false
holds_sorbent(::SorbentBedModel) = true
