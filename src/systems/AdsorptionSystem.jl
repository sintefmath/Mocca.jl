struct Column <: Jutul.JutulEntity end

struct AdsorptionSystem{N, RealT<:Real, IsoT<:AbstractIsotherm, MtT<:AbstractMassTransfer} <: JutulDarcy.MultiComponentSystem
    component_names::Vector{String}
    molecular_masses::SVector{N, RealT}         # Molecular masses per component [kg/mol]
    heat_capacity_gas::SVector{N, RealT}        # Heat capacity of gas per component [J/(kg·K)]
    heat_capacity_adsorbed::SVector{N, RealT}   # Heat capacity of adsorbed phase per component [J/(kg·K)]
    isotherm::IsoT
    mass_transfer::MtT
end

const AdsorptionModel = Jutul.SimulationModel{<:Any,<:AdsorptionSystem,<:Any,<:Any}

"""
    AdsorptionSystem(; isotherm, mass_transfer, molecular_masses,
        component_names, heat_capacity_gas, heat_capacity_adsorbed)

Construct an N-component adsorption system from explicit physics objects.
The number of components is inferred from the length of `component_names`.
"""
function AdsorptionSystem(;
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
    return AdsorptionSystem(
        component_names,
        SVector{N}(molecular_masses),
        SVector{N}(heat_capacity_gas),
        SVector{N}(heat_capacity_adsorbed),
        isotherm,
        mass_transfer,
    )
end

function AdsorptionSystem(constants::HaghpanahConstants)
    isotherm = DualSiteLangmuir(constants)
    mass_transfer = LinearDrivingForce(constants.D_m, constants.τ, constants.ϵ_p, constants.d_p)

    return AdsorptionSystem(
        isotherm = isotherm,
        mass_transfer = mass_transfer,
        molecular_masses = SVector(constants.molecularMassOfCO2, constants.molecularMassOfN2),
        component_names = ["CO2", "N2"],
        heat_capacity_gas = constants.C_pg,
        heat_capacity_adsorbed = constants.C_pa,
    )
end

function AdsorptionSystem(constants::adsorptionConstants)
    isotherm = DualSiteLangmuir(constants)
    mass_transfer = LinearDrivingForce(constants.D_m, constants.τ, constants.ϵ_p, constants.d_p)

    return AdsorptionSystem(
        isotherm = isotherm,
        mass_transfer = mass_transfer,
        molecular_masses = constants.molecular_masses,
        component_names = constants.component_names,
        heat_capacity_gas = constants.C_pg,
        heat_capacity_adsorbed = constants.C_pa,
    )
end

# Overload JutulDarcy functions

JutulDarcy.component_names(sys::AdsorptionSystem) = sys.component_names
JutulDarcy.number_of_components(sys::AdsorptionSystem) = length(sys.component_names)
JutulDarcy.has_other_phase(::AdsorptionSystem) = false
JutulDarcy.phase_names(::AdsorptionSystem) = ["gas"]
JutulDarcy.number_of_phases(::AdsorptionSystem) = 1
JutulDarcy.get_reference_phase_index(::AdsorptionSystem) = 1
JutulDarcy.eachphase(::AdsorptionSystem) = (1,)

# `AdsorptionSystem` is not a `JutulDarcy.MultiPhaseSystem` in the method table, so without this
# forward the default `Jutul.discretize_domain(::DataDomain, ...)` would wrap the mesh only.
# Delegate to the `MultiPhaseSystem` implementation (TPFA + entity propagation from `DataDomain`).
function Jutul.discretize_domain(d::Jutul.DataDomain, system::AdsorptionSystem, ::Val{:default}; kwarg...)
    return invoke(
        Jutul.discretize_domain,
        Tuple{Jutul.DataDomain, JutulDarcy.MultiPhaseSystem, Val{:default}},
        d, system, Val(:default);
        kwarg...,
    )
end
