struct Column <: Jutul.JutulEntity end

struct AdsorptionSystem{N, RealT<:Real, IsoT<:AbstractIsotherm, MtT<:AbstractMassTransfer} <: Jutul.JutulSystem
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

# System interface functions

component_names(sys::AdsorptionSystem) = sys.component_names
number_of_components(sys::AdsorptionSystem) = length(sys.component_names)

# Set up a discretized domain with TPFA potential flow for the adsorption system.
function Jutul.discretize_domain(d::Jutul.DataDomain, system::AdsorptionSystem, ::Val{:default}; general_ad = true, kwarg...)
    g = Jutul.physical_representation(d)
    N = d[:neighbors]
    nc = Jutul.number_of_cells(g)
    flow_disc = Jutul.PotentialFlow(N, nc)
    disc = (mass_flow = flow_disc,)
    domain = Jutul.DiscretizedDomain(g, disc; kwarg...)
    Jutul.transfer_entities!(domain, d)
    return domain
end
