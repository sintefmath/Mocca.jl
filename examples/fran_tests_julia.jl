using Parameters

abstract type GenericModel end

@with_kw struct SpecificModel <: GenericModel
   
    number_of_components::Int64 = 2

end

function findnumc(model::GenericModel)
    return model.number_of_components
end



a = SpecificModel()

b = findnumc(a)