using Parameters

@with_kw struct TwoComponentAdsorptionSystem{T} <: AdsorptionSystem where T<:ConstantsStruct
    component_names::Vector{String} = ["CO2","N2"]
    permeability::Float64
    dispersion::Float64
    p::T
end

JutulDarcy.number_of_components(sys::TwoComponentAdsorptionSystem) = 2
