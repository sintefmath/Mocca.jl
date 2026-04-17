struct WallArea{T} <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::WallArea{T}, symb) where T
    dx = data_domain[:dx]
    if T == :in
        r = first(data_domain[:r_in, Column()])
    else
        @assert T == :out
        r = first(data_domain[:r_out, Column()])
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
function Jutul.select_parameters!(S, disc::Jutul.PotentialFlow, model::Jutul.SimulationModel{D, S2}) where {D, S2 <: AdsorptionSystem}
    S[:Transmissibilities] = Transmissibilities()
end
