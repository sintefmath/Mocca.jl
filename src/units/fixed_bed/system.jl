"""
    FixedBed{N, RealT, IsoT, MtT, ThT, AllC} <: SorbentBed

A fixed (packed) adsorbent bed: gas flows axially through a stationary sorbent.
How energy is treated is set by the thermal model (see
[`AbstractThermalModel`](@ref)); the default [`WithWall`](@ref) includes a column
wall exchanging heat with the gas and the surroundings.

`AllC` records whether the heat of sorption and adsorbed-phase heat capacity sum
over all components (`true`) or only the first (`false`, the behaviour of Mocca
0.1.0 and earlier). See [`HeatOfSorption`](@ref).

`AdsorptionSystem` is kept as an alias for backwards compatibility.
"""
struct FixedBed{N, RealT<:Real, IsoT<:AbstractIsotherm, MtT<:AbstractMassTransfer, ThT<:AbstractThermalModel, AllC} <: SorbentBed
    component_names::Vector{String}
    molecular_masses::SVector{N, RealT}         # Molecular masses per component [kg/mol]
    heat_capacity_gas::SVector{N, RealT}        # Heat capacity of gas per component [J/(kg·K)]
    heat_capacity_adsorbed::SVector{N, RealT}   # Heat capacity of adsorbed phase per component [J/(kg·K)]
    isotherm::IsoT
    mass_transfer::MtT
    thermal::ThT
end

const AdsorptionSystem = FixedBed
const FixedBedModel = Jutul.SimulationModel{<:Any,<:FixedBed,<:Any,<:Any}
const AdsorptionModel = FixedBedModel

"""
    FixedBed(; isotherm, mass_transfer, molecular_masses,
        component_names, heat_capacity_gas, heat_capacity_adsorbed,
        thermal = WithWall(), sorption_heat_all_components = false)

Construct an N-component adsorption system from explicit physics objects.
The number of components is inferred from the length of `component_names`.

Set `sorption_heat_all_components = true` to include every component in the
heat of sorption and adsorbed-phase heat capacity. The default `false` only
includes the first component, which reproduces earlier Mocca results.
"""
function FixedBed(;
    isotherm::AbstractIsotherm,
    mass_transfer::AbstractMassTransfer,
    molecular_masses,
    component_names,
    heat_capacity_gas,
    heat_capacity_adsorbed,
    thermal::AbstractThermalModel = WithWall(),
    sorption_heat_all_components::Bool = false,
)
    N = length(component_names)
    @assert length(molecular_masses) == N "molecular_masses length ($(length(molecular_masses))) must match component_names length ($N)"
    @assert length(heat_capacity_gas) == N "heat_capacity_gas length must match component_names length ($N)"
    @assert length(heat_capacity_adsorbed) == N "heat_capacity_adsorbed length must match component_names length ($N)"
    mm = SVector{N}(molecular_masses)
    RealT = promote_type(eltype(mm), eltype(heat_capacity_gas), eltype(heat_capacity_adsorbed))
    IsoT, MtT, ThT = typeof(isotherm), typeof(mass_transfer), typeof(thermal)
    return FixedBed{N, RealT, IsoT, MtT, ThT, sorption_heat_all_components}(
        component_names,
        SVector{N, RealT}(mm),
        SVector{N, RealT}(heat_capacity_gas),
        SVector{N, RealT}(heat_capacity_adsorbed),
        isotherm,
        mass_transfer,
        thermal,
    )
end

function FixedBed(constants::HaghpanahConstants; kwarg...)
    isotherm = DualSiteLangmuir(constants)
    mass_transfer = LinearDrivingForce(constants.D_m, constants.τ, constants.ϵ_p, constants.d_p)

    return FixedBed(;
        isotherm = isotherm,
        mass_transfer = mass_transfer,
        molecular_masses = SVector(constants.molecularMassOfCO2, constants.molecularMassOfN2),
        component_names = ["CO2", "N2"],
        heat_capacity_gas = constants.C_pg,
        heat_capacity_adsorbed = constants.C_pa,
        kwarg...
    )
end

function FixedBed(constants::adsorptionConstants; kwarg...)
    isotherm = DualSiteLangmuir(constants)
    mass_transfer = LinearDrivingForce(constants.D_m, constants.τ, constants.ϵ_p, constants.d_p)

    return FixedBed(;
        isotherm = isotherm,
        mass_transfer = mass_transfer,
        molecular_masses = constants.molecular_masses,
        component_names = constants.component_names,
        heat_capacity_gas = constants.C_pg,
        heat_capacity_adsorbed = constants.C_pa,
        kwarg...
    )
end

# Terms the shared conservation laws pick up for a fixed bed
sorption_heat_all_components(::FixedBed{N, R, I, M, T, AllC}) where {N, R, I, M, T, AllC} = AllC

mass_source_terms(::FixedBed) = (SorptionUptake(),)
function heat_capacity_terms(sys::FixedBed)
    return (SolidHeatCapacity(), AdsorbedPhaseHeatCapacity{sorption_heat_all_components(sys)}())
end
function energy_source_terms(sys::FixedBed)
    hos = HeatOfSorption{sorption_heat_all_components(sys)}()
    return has_wall(sys) ? (hos, WallExchange()) : (hos,)
end
wall_source_terms(::FixedBed) = (WallExchange(), AmbientLoss(), WallEndConduction())
