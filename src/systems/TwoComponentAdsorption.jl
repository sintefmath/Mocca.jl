using Parameters

@with_kw struct TwoComponentAdsorptionSystem{T, RealT<:Real, IsoT<:AbstractIsotherm, MtT<:AbstractMassTransfer} <: AdsorptionSystem where T<:ConstantsStruct
    component_names::Vector{String} = ["CO2","N2"]
    molecular_masses::SVector{2, RealT}
    heat_capacity_gas::SVector{2, RealT}
    heat_capacity_adsorbed::SVector{2, RealT}
    permeability::RealT
    dispersion::RealT
    isotherm::IsoT
    mass_transfer::MtT
    p::T
end

function TwoComponentAdsorptionSystem(constants::ConstantsStruct)
    permeability = compute_permeability(constants)
    axial_dispersion = calc_dispersion(constants)
    isotherm = DualSiteLangmuir(constants)
    mass_transfer = LinearDrivingForce(constants.D_m, constants.τ, constants.ϵ_p, constants.d_p)
    molecular_masses = SVector(constants.molecularMassOfCO2, constants.molecularMassOfN2)

    return TwoComponentAdsorptionSystem(;
        molecular_masses = molecular_masses,
        heat_capacity_gas = constants.C_pg,
        heat_capacity_adsorbed = constants.C_pa,
        permeability = permeability,
        dispersion = axial_dispersion,
        isotherm = isotherm,
        mass_transfer = mass_transfer,
        p = constants
    )
end

JutulDarcy.number_of_components(sys::TwoComponentAdsorptionSystem) = 2

JutulDarcy.number_of_phases(sys::TwoComponentAdsorptionSystem) = 1

JutulDarcy.get_reference_phase_index(sys::TwoComponentAdsorptionSystem) = 1

JutulDarcy.eachphase(sys::TwoComponentAdsorptionSystem) = (1,)
