using Parameters

@with_kw struct TwoComponentAdsorptionSystem{T, RealT<:Real} <: AdsorptionSystem where T<:ConstantsStruct
    component_names::Vector{String} = ["CO2","N2"]
    permeability::RealT
    dispersion::RealT
    p::T
end

function TwoComponentAdsorptionSystem(constants::ConstantsStruct)
    permeability = compute_permeability(constants)
    axial_dispersion = calc_dispersion(constants)

    return TwoComponentAdsorptionSystem(; permeability = permeability, dispersion = axial_dispersion, p = constants)
end

JutulDarcy.number_of_components(sys::TwoComponentAdsorptionSystem) = 2

JutulDarcy.number_of_phases(sys::TwoComponentAdsorptionSystem) = 1

JutulDarcy.get_reference_phase_index(sys::TwoComponentAdsorptionSystem) = 1

JutulDarcy.eachphase(sys::TwoComponentAdsorptionSystem) = (1,)
