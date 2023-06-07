struct MySuperStruct end
struct AdsorptionRates <: Jutul.VectorVariables
end

function Jutul.minimum_value(::AdsorptionRates)
    return 1e-10
end

Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::AdsorptionRates) = JutulDarcy.number_of_components(model.system)

Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::AdsorptionRates) = JutulDarcy.number_of_components(model.system)

struct AverageMolecularMass <: Jutul.ScalarVariable
end

Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::AverageMolecularMass) = 1

Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::AverageMolecularMass) = 1

struct Concentrations <: Jutul.VectorVariables end

struct AdsorptionMassTransfer <: Jutul.VectorVariables end

Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::AdsorptionMassTransfer) = JutulDarcy.number_of_components(model.system)

Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::AdsorptionMassTransfer) = JutulDarcy.number_of_components(model.system)

Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::Concentrations) = JutulDarcy.number_of_components(model.system)

Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::Concentrations) = JutulDarcy.number_of_components(model.system)

JutulDarcy.phase_names(::AdsorptionFlowSystem) = ["ouronlyphase"]
JutulDarcy.number_of_phases(::AdsorptionFlowSystem) = 1

struct GasMoleFractions <: JutulDarcy.CompositionalFractions
    dz_max::Float64
    GasMoleFractions(; dz_max=0.2) = new(dz_max)
end

function Jutul.minimum_value(::GasMoleFractions)
    return 1e-12
end

function Jutul.variable_scale(::JutulDarcy.Pressure)
    nothing
end
Jutul.absolute_increment_limit(z::GasMoleFractions) = z.dz_max



# Temperature variables
abstract type Energy <: Jutul.ScalarVariable end
struct ColumnEnergy <: Energy end
struct WallEnergy <: Energy end

struct EnthalpyChange <: Jutul.VectorVariables end
Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::EnthalpyChange) = JutulDarcy.number_of_components(model.system)
Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionFlowSystem}, ::EnthalpyChange) = JutulDarcy.number_of_components(model.system)

abstract type SpecificHeatCapacity <: Jutul.ScalarVariable end
struct SpecificHeatCapasityAdsorbent <: SpecificHeatCapacity end
struct SpecificHeatCapasityFluid <: SpecificHeatCapacity end

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
