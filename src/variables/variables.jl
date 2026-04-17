# Base type for per-component variables (one value per component per entity)
abstract type ComponentVariable <: Jutul.VectorVariables end

Jutul.degrees_of_freedom_per_entity(model::AdsorptionModel, ::ComponentVariable) = number_of_components(model.system)
Jutul.values_per_entity(model::AdsorptionModel, ::ComponentVariable) = number_of_components(model.system)

include("primary_variables.jl")
include("fluid.jl")
include("adsorption.jl")
include("energy.jl")
include("geometry.jl")
include("column_parameters.jl")
