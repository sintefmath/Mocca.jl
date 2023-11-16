#### Secondary variables

struct AverageMolecularMass <: Jutul.ScalarVariable end
Jutul.degrees_of_freedom_per_entity(model::AdsorptionModel, ::AverageMolecularMass) = 1
Jutul.values_per_entity(model::AdsorptionModel, ::AverageMolecularMass) = 1

struct Concentrations <: Jutul.VectorVariables end
Jutul.degrees_of_freedom_per_entity(model::AdsorptionModel, ::Concentrations) = JutulDarcy.number_of_components(model.system)
Jutul.values_per_entity(model::AdsorptionModel, ::Concentrations) = JutulDarcy.number_of_components(model.system)

struct AdsorptionMassTransfer <: Jutul.VectorVariables end
Jutul.degrees_of_freedom_per_entity(model::AdsorptionModel, ::AdsorptionMassTransfer) = JutulDarcy.number_of_components(model.system)
Jutul.values_per_entity(model::AdsorptionModel, ::AdsorptionMassTransfer) = JutulDarcy.number_of_components(model.system)

abstract type Energy <: Jutul.ScalarVariable end
struct ColumnEnergy <: Energy end
struct WallEnergy <: Energy end

struct EnthalpyChange <: Jutul.VectorVariables end
Jutul.degrees_of_freedom_per_entity(model::AdsorptionModel, ::EnthalpyChange) = JutulDarcy.number_of_components(model.system)
Jutul.values_per_entity(model::AdsorptionModel, ::EnthalpyChange) = JutulDarcy.number_of_components(model.system)

abstract type SpecificHeatCapacity <: Jutul.ScalarVariable end
struct SpecificHeatCapacityAdsorbent <: SpecificHeatCapacity end
struct SpecificHeatCapacityFluid <: SpecificHeatCapacity end

struct ThermalConductivities <: Jutul.ScalarVariable end
Jutul.variable_scale(::ThermalConductivities) = 1e-10
Jutul.minimum_value(::ThermalConductivities) = 0.0
Jutul.default_value(model, ::ThermalConductivities) = 1e-3
Jutul.associated_entity(::ThermalConductivities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::ThermalConductivities, symb)
    if haskey(data_domain, :thermal_conductivity, Jutul.Cells())
        U = data_domain[:thermal_conductivity]
        g = Jutul.physical_representation(data_domain)
        T = Jutul.compute_face_trans(g, U)
    else
        error(":thermal_conductivity symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

struct DiffusionTransmissibilities <: Jutul.ScalarVariable end
Jutul.variable_scale(::DiffusionTransmissibilities) = 1e-10
Jutul.minimum_value(::DiffusionTransmissibilities) = 0.0
Jutul.default_value(model, ::DiffusionTransmissibilities) = 1e-3
Jutul.associated_entity(::DiffusionTransmissibilities) = Jutul.Faces()

function Jutul.default_parameter_values(data_domain, model, param::DiffusionTransmissibilities, symb)
    if haskey(data_domain, :diffusion_coefficient, Jutul.Cells())
        U = data_domain[:diffusion_coefficient]
        g = Jutul.physical_representation(data_domain)
        T = Jutul.compute_face_trans(g, U)
    else
        error(":diffusion_coefficient symbol must be present in DataDomain to initialize parameter $symb, had keys: $(keys(data_domain))")
    end
    return T
end

struct WallArea{T} <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, param::WallArea{T}, symb) where T
    T::Symbol
    dx = data_domain[:dx]
    sys = model.system
    if T == :in
        F = x -> area_wall_in(sys, x)
    else
        @assert T == :out
        F = x -> area_wall_out(sys, x)
    end
    return F.(dx)
end
