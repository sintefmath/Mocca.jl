using Parameters


@with_kw struct TwoComponentAdsorptionSystem <: AdsorptionSystem
   
    number_of_components::Int64 = 2
   
    permeability::Float64
    dispersion::Float64

    p::HaghpanahParameters = HaghpanahParameters() #TODO Move this out of here
end

