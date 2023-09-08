# Primary variables

# pressure: p
# component concentration in the gas phase (mol. frac.): yi (i-1 dofs)
# temperature in the column: T
# temperature in the wall: Tw
# adsorbed concentration: qi (i dofs)



# pressure
# JutulDarcy.Pressure

# GasMoleFractions
struct GasMoleFractions <: JutulDarcy.CompositionalFractions
    dz_max::Float64
    GasMoleFractions(; dz_max=0.2) = new(dz_max)
end

const MIN_GAS_MOLEFRACTION = 1e-12

function Jutul.minimum_value(::GasMoleFractions)
    return MIN_GAS_MOLEFRACTION
end

function Jutul.maximum_value(::GasMoleFractions)
    return 1.0 - MIN_GAS_MOLEFRACTION
end

Jutul.absolute_increment_limit(z::GasMoleFractions) = z.dz_max


# Temperature 
# JutulDarcy.Temperature

# WallTemperature
# JutulDarcy.Temperature

# AdsorbedConcentration
struct AdsorbedConcentration <: Jutul.VectorVariables
end

function Jutul.minimum_value(::AdsorbedConcentration)
    return 1e-10
end

function Jutul.degrees_of_freedom_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionSystem}, ::AdsorbedConcentration)
    JutulDarcy.number_of_components(model.system)
end
Jutul.values_per_entity(model::Jutul.SimulationModel{<:Any,AdsorptionSystem}, ::AdsorbedConcentration) = JutulDarcy.number_of_components(model.system)


