using Parameters


@with_kw struct TwoComponentAdsorptionSystem{vel_model,disp_model} <: AdsorptionSystem
   
    number_of_components::Int64 = 2
   
    velocity_model = vel_model
    dispersion_model = disp_model
   
    permeability = velocity_model.permeability
    dispersion = dispersion_model.dispersion


    p::HaghpanahParameters = HaghpanahParameters() #TODO Move this out of here
end

