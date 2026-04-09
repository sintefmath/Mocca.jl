# Column-entity Jutul parameters
# Each is a single scalar associated with the Column entity (count=1),
# accessed at runtime as state.ParamName[1].

struct AdsorbentDensity <: Jutul.ScalarVariable end
Jutul.associated_entity(::AdsorbentDensity) = Column()
function Jutul.default_parameter_values(data_domain, model, ::AdsorbentDensity, symb)
    return copy(data_domain[:adsorbent_density, Column()])
end

struct AdsorbentHeatCapacity <: Jutul.ScalarVariable end
Jutul.associated_entity(::AdsorbentHeatCapacity) = Column()
function Jutul.default_parameter_values(data_domain, model, ::AdsorbentHeatCapacity, symb)
    return copy(data_domain[:adsorbent_heat_capacity, Column()])
end

struct WallDensity <: Jutul.ScalarVariable end
Jutul.associated_entity(::WallDensity) = Column()
function Jutul.default_parameter_values(data_domain, model, ::WallDensity, symb)
    return copy(data_domain[:wall_density, Column()])
end

struct WallHeatCapacity <: Jutul.ScalarVariable end
Jutul.associated_entity(::WallHeatCapacity) = Column()
function Jutul.default_parameter_values(data_domain, model, ::WallHeatCapacity, symb)
    return copy(data_domain[:wall_heat_capacity, Column()])
end

struct WallConductivity <: Jutul.ScalarVariable end
Jutul.associated_entity(::WallConductivity) = Column()
function Jutul.default_parameter_values(data_domain, model, ::WallConductivity, symb)
    return copy(data_domain[:wall_conductivity, Column()])
end

struct FluidViscosity <: Jutul.ScalarVariable end
Jutul.associated_entity(::FluidViscosity) = Column()
function Jutul.default_parameter_values(data_domain, model, ::FluidViscosity, symb)
    return copy(data_domain[:fluid_viscosity, Column()])
end

struct FluidDensity <: Jutul.ScalarVariable end
Jutul.associated_entity(::FluidDensity) = Column()
function Jutul.default_parameter_values(data_domain, model, ::FluidDensity, symb)
    return copy(data_domain[:fluid_density, Column()])
end

struct InnerHeatTransferCoeff <: Jutul.ScalarVariable end
Jutul.associated_entity(::InnerHeatTransferCoeff) = Column()
function Jutul.default_parameter_values(data_domain, model, ::InnerHeatTransferCoeff, symb)
    return copy(data_domain[:inner_htc, Column()])
end

struct OuterHeatTransferCoeff <: Jutul.ScalarVariable end
Jutul.associated_entity(::OuterHeatTransferCoeff) = Column()
function Jutul.default_parameter_values(data_domain, model, ::OuterHeatTransferCoeff, symb)
    return copy(data_domain[:outer_htc, Column()])
end

struct AmbientTemperature <: Jutul.ScalarVariable end
Jutul.associated_entity(::AmbientTemperature) = Column()
function Jutul.default_parameter_values(data_domain, model, ::AmbientTemperature, symb)
    return copy(data_domain[:ambient_temperature, Column()])
end

struct WallCrossSectionArea <: Jutul.ScalarVariable end
Jutul.associated_entity(::WallCrossSectionArea) = Column()
function Jutul.default_parameter_values(data_domain, model, ::WallCrossSectionArea, symb)
    r_in = first(data_domain[:r_in, Column()])
    r_out = first(data_domain[:r_out, Column()])
    return [π * (r_out^2 - r_in^2)]
end
