struct WallArea{T} <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::WallArea{T}, symb) where T
    dx = data_domain[:dx]
    if T == :in
        r = first(data_domain[:r_in, Unit()])
    else
        @assert T == :out
        r = first(data_domain[:r_out, Unit()])
    end
    return 2π .* r .* dx  # Lateral surface area of cylindrical shell
end

struct CellDx <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::CellDx, symb)
    return copy(data_domain[:dx])
end

# Cell volumes (solid + pore)
struct BulkVolume <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::BulkVolume, symb)
    if haskey(data_domain, :volumes)
        return copy(data_domain[:volumes])
    else
        error(":volumes must be present in DataDomain to initialize $symb")
    end
end

# Pore volume = volume × porosity
struct FluidVolume <: Jutul.ScalarVariable end
Jutul.minimum_value(::FluidVolume) = eps()

function Jutul.default_parameter_values(data_domain, model, param::FluidVolume, symb)
    if haskey(data_domain, :fluid_volume, Jutul.Cells())
        return copy(data_domain[:fluid_volume])
    elseif haskey(data_domain, :volumes) && haskey(data_domain, :porosity)
        return copy(data_domain[:volumes] .* data_domain[:porosity])
    else
        error(":fluid_volume or (:volumes + :porosity) must be present in DataDomain to initialize $symb")
    end
end

# Cell permeability [m²], used for the half-cell connection at a port
struct Permeability <: Jutul.ScalarVariable end
Jutul.variable_scale(::Permeability) = 1e-10
Jutul.minimum_value(::Permeability) = 0.0

function Jutul.default_parameter_values(data_domain, model, param::Permeability, symb)
    return copy(data_domain[:permeability])
end

# Flow cross-section of the bed [m²]
struct BedCrossSectionArea <: Jutul.ScalarVariable end
Jutul.associated_entity(::BedCrossSectionArea) = Unit()
Jutul.minimum_value(::BedCrossSectionArea) = 0.0

function Jutul.default_parameter_values(data_domain, model, ::BedCrossSectionArea, symb)
    r_in = first(data_domain[:r_in, Unit()])
    return [π * r_in^2]
end

"""
    compute_column_face_area(model, state)

Flow cross-section of a distributed unit [m²]. Read from the `BedCrossSectionArea`
parameter when `state` carries parameters, otherwise from the column radius in
the data domain.
"""
function compute_column_face_area(model::DistributedUnitModel, state)
    if _state_has(state, :BedCrossSectionArea)
        return state.BedCrossSectionArea[1]
    else
        r_in = first(model.data_domain[:r_in, Unit()])
        return π * r_in^2
    end
end

"""
    calc_bc_trans(model, state, cell)

Transmissibility [m³] of the half-cell connection between `cell` and the end of the
unit it touches.
"""
function calc_bc_trans(model::DistributedUnitModel, state, cell)
    k = _state_or_domain(state, model, :Permeability, :permeability)[cell]
    dx = _state_or_domain(state, model, :CellDx, :dx)[cell] / 2.0
    A = compute_column_face_area(model, state)
    return k * A / dx
end

# TPFA transmissibilities (face parameter)
struct Transmissibilities <: Jutul.ScalarVariable end
Jutul.variable_scale(::Transmissibilities) = 1e-10
Jutul.minimum_value(::Transmissibilities) = 0.0
Jutul.associated_entity(::Transmissibilities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::Transmissibilities, symb)
    if haskey(data_domain, :transmissibilities, Jutul.Faces())
        return copy(data_domain[:transmissibilities])
    elseif haskey(data_domain, :permeability, Jutul.Cells())
        g = Jutul.physical_representation(data_domain)
        K = data_domain[:permeability]
        return Jutul.compute_face_trans(g, K)
    else
        error(":permeability or :transmissibilities must be present in DataDomain to initialize $symb")
    end
end

# Inject Transmissibilities when PotentialFlow discretization is present
function Jutul.select_parameters!(S, disc::Jutul.PotentialFlow, model::Jutul.SimulationModel{D, S2}) where {D, S2 <: DistributedUnit}
    S[:Transmissibilities] = Transmissibilities()
end
