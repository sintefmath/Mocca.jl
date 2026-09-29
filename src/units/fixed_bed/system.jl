"""
    FixedBed{N, RealT, IsoT, MtT} <: SorbentBed

A fixed (packed) adsorbent bed: gas flows axially through a stationary sorbent,
optionally exchanging heat with a column wall. `AdsorptionSystem` is kept as an
alias for backwards compatibility.
"""
struct FixedBed{N, RealT<:Real, IsoT<:AbstractIsotherm, MtT<:AbstractMassTransfer} <: SorbentBed
    component_names::Vector{String}
    molecular_masses::SVector{N, RealT}         # Molecular masses per component [kg/mol]
    heat_capacity_gas::SVector{N, RealT}        # Heat capacity of gas per component [J/(kg·K)]
    heat_capacity_adsorbed::SVector{N, RealT}   # Heat capacity of adsorbed phase per component [J/(kg·K)]
    isotherm::IsoT
    mass_transfer::MtT
end

const AdsorptionSystem = FixedBed
const FixedBedModel = Jutul.SimulationModel{<:Any,<:FixedBed,<:Any,<:Any}
const AdsorptionModel = FixedBedModel

"""
    FixedBed(; isotherm, mass_transfer, molecular_masses,
        component_names, heat_capacity_gas, heat_capacity_adsorbed)

Construct an N-component adsorption system from explicit physics objects.
The number of components is inferred from the length of `component_names`.
"""
function FixedBed(;
    isotherm::AbstractIsotherm,
    mass_transfer::AbstractMassTransfer,
    molecular_masses,
    component_names,
    heat_capacity_gas,
    heat_capacity_adsorbed,
)
    N = length(component_names)
    @assert length(molecular_masses) == N "molecular_masses length ($(length(molecular_masses))) must match component_names length ($N)"
    @assert length(heat_capacity_gas) == N "heat_capacity_gas length must match component_names length ($N)"
    @assert length(heat_capacity_adsorbed) == N "heat_capacity_adsorbed length must match component_names length ($N)"
    return FixedBed(
        component_names,
        SVector{N}(molecular_masses),
        SVector{N}(heat_capacity_gas),
        SVector{N}(heat_capacity_adsorbed),
        isotherm,
        mass_transfer,
    )
end

function FixedBed(constants::HaghpanahConstants)
    isotherm = DualSiteLangmuir(constants)
    mass_transfer = LinearDrivingForce(constants.D_m, constants.τ, constants.ϵ_p, constants.d_p)

    return FixedBed(
        isotherm = isotherm,
        mass_transfer = mass_transfer,
        molecular_masses = SVector(constants.molecularMassOfCO2, constants.molecularMassOfN2),
        component_names = ["CO2", "N2"],
        heat_capacity_gas = constants.C_pg,
        heat_capacity_adsorbed = constants.C_pa,
    )
end

function FixedBed(constants::adsorptionConstants)
    isotherm = DualSiteLangmuir(constants)
    mass_transfer = LinearDrivingForce(constants.D_m, constants.τ, constants.ϵ_p, constants.d_p)

    return FixedBed(
        isotherm = isotherm,
        mass_transfer = mass_transfer,
        molecular_masses = constants.molecular_masses,
        component_names = constants.component_names,
        heat_capacity_gas = constants.C_pg,
        heat_capacity_adsorbed = constants.C_pa,
    )
end
