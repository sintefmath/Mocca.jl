struct TwoComponentAdsorptionSystem{RealT<:Real, IsoT<:AbstractIsotherm, MtT<:AbstractMassTransfer} <: AdsorptionSystem
    component_names::Vector{String}
    molecular_masses::SVector{2, RealT}         # Molecular masses per component [kg/mol]
    heat_capacity_gas::SVector{2, RealT}        # Heat capacity of gas per component [J/(kg·K)]
    heat_capacity_adsorbed::SVector{2, RealT}   # Heat capacity of adsorbed phase per component [J/(kg·K)]
    isotherm::IsoT
    mass_transfer::MtT
end

"""
    TwoComponentAdsorptionSystem(; isotherm, mass_transfer, molecular_masses,
        component_names = ["CO2", "N2"], heat_capacity_gas, heat_capacity_adsorbed)

Construct a two-component adsorption system from explicit physics objects.
"""
function TwoComponentAdsorptionSystem(;
    isotherm::AbstractIsotherm,
    mass_transfer::AbstractMassTransfer,
    molecular_masses,
    component_names = ["CO2", "N2"],
    heat_capacity_gas,
    heat_capacity_adsorbed,
)
    return TwoComponentAdsorptionSystem(
        component_names,
        SVector{2}(molecular_masses),
        SVector{2}(heat_capacity_gas),
        SVector{2}(heat_capacity_adsorbed),
        isotherm,
        mass_transfer,
    )
end

function TwoComponentAdsorptionSystem(constants::ConstantsStruct)
    isotherm = DualSiteLangmuir(constants)
    mass_transfer = LinearDrivingForce(constants.D_m, constants.τ, constants.ϵ_p, constants.d_p)

    return TwoComponentAdsorptionSystem(
        isotherm = isotherm,
        mass_transfer = mass_transfer,
        molecular_masses = SVector(constants.molecularMassOfCO2, constants.molecularMassOfN2),
        heat_capacity_gas = constants.C_pg,
        heat_capacity_adsorbed = constants.C_pa,
    )
end

JutulDarcy.number_of_components(sys::TwoComponentAdsorptionSystem) = 2

JutulDarcy.number_of_phases(sys::TwoComponentAdsorptionSystem) = 1

JutulDarcy.get_reference_phase_index(sys::TwoComponentAdsorptionSystem) = 1

JutulDarcy.eachphase(sys::TwoComponentAdsorptionSystem) = (1,)
