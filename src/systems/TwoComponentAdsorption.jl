using Parameters

export TwoComponentAdsorptionSystem 
@with_kw struct TwoComponentAdsorptionSystem <: AdsorptionSystem
   
    number_of_components::Int64 = 2
   
    permeability::Float64
    dispersion::Float64
    p::ParameterStruct 

end
