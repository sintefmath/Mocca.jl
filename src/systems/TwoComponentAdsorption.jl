using Parameters

@with_kw struct TwoComponentAdsorptionSystem <: AdsorptionSystem
   
    number_of_components::Int64 = 2
    component_names::Vector{String} = ["CO2","N2"]
   
    permeability::Float64
    dispersion::Float64
    p::ConstantsStruct 

end
