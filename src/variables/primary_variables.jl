# Primary variables

# pressure: p
# component concentration in the gas phase (mol. frac.): yi (i-1 dofs)
# temperature in the column: T
# temperature in the wall: Tw
# adsorbed concentration: qi (i dofs)


struct Pressure <: Jutul.ScalarVariable
    minimum::Float64
    max_abs::Float64
    max_rel::Float64
    scale::Float64
end

function Pressure(; minimum = Jutul.si_unit(:bar)*1e-2, max_abs = Inf, max_rel = 0.2, scale = Jutul.si_unit(:bar))
    return Pressure(minimum, max_abs, max_rel, scale)
end

Jutul.minimum_value(p::Pressure) = p.minimum
Jutul.variable_scale(p::Pressure) = p.scale
Jutul.absolute_increment_limit(p::Pressure) = p.max_abs
Jutul.relative_increment_limit(p::Pressure) = p.max_rel

struct Temperature <: Jutul.ScalarVariable
    min::Float64
    max_abs::Float64
    max_rel::Float64
end

function Temperature(; min = 200.0, max_abs = Inf, max_rel = 0.2)
    return Temperature(min, max_abs, max_rel)
end

Jutul.minimum_value(t::Temperature) = t.min
Jutul.absolute_increment_limit(t::Temperature) = t.max_abs
Jutul.relative_increment_limit(t::Temperature) = t.max_rel

struct GasMoleFractions <: Jutul.FractionVariables
    dz_max::Float64
    GasMoleFractions(; dz_max=0.2) = new(dz_max)
end

Jutul.values_per_entity(model::AdsorptionModel, ::GasMoleFractions) = number_of_components(model.system)

const MIN_GAS_MOLEFRACTION = 1e-12

function Jutul.minimum_value(::GasMoleFractions)
    return MIN_GAS_MOLEFRACTION
end

function Jutul.maximum_value(::GasMoleFractions)
    return 1.0 - MIN_GAS_MOLEFRACTION
end

Jutul.absolute_increment_limit(z::GasMoleFractions) = z.dz_max

struct AdsorbedConcentration <: Jutul.VectorVariables
end

function Jutul.minimum_value(::AdsorbedConcentration)
    return 1e-10
end

function Jutul.degrees_of_freedom_per_entity(model::AdsorptionModel, ::AdsorbedConcentration)
    number_of_components(model.system)
end
Jutul.values_per_entity(model::AdsorptionModel, ::AdsorbedConcentration) = number_of_components(model.system)


