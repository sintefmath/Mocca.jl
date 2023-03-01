struct MySuperStruct end
struct AdsorptionRates <: Jutul.VectorVariables
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
    return 1e-20
end

#Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::GasMoleFractions) = JutulDarcy.number_of_components(model.system)

#Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any, AdsorptionFlowSystem}, ::GasMoleFractions) = JutulDarcy.number_of_components(model.system)

Jutul.absolute_increment_limit(z::GasMoleFractions) = z.dz_max

