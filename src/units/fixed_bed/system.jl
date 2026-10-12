"""
    FixedBed{N, RealT, IsoT, MtT, ThT} <: SorbentBed

A fixed (packed) adsorbent bed: gas flows axially through a stationary sorbent.
How energy is treated is set by the thermal model (see
[`AbstractThermalModel`](@ref)); the default [`WithWall`](@ref) includes a column
wall exchanging heat with the gas and the surroundings.

`AdsorptionSystem` is kept as an alias for backwards compatibility.
"""
struct FixedBed{N, RealT<:Real, IsoT<:AbstractIsotherm, MtT<:AbstractMassTransfer, ThT<:AbstractThermalModel} <: SorbentBed
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
        thermal = WithWall())

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
    thermal::AbstractThermalModel = WithWall(),
)
    N = length(component_names)
    @assert length(molecular_masses) == N "molecular_masses length ($(length(molecular_masses))) must match component_names length ($N)"
    @assert length(heat_capacity_gas) == N "heat_capacity_gas length must match component_names length ($N)"
    @assert length(heat_capacity_adsorbed) == N "heat_capacity_adsorbed length must match component_names length ($N)"
    mm = SVector{N}(molecular_masses)
    RealT = promote_type(eltype(mm), eltype(heat_capacity_gas), eltype(heat_capacity_adsorbed))
    IsoT, MtT, ThT = typeof(isotherm), typeof(mass_transfer), typeof(thermal)
    return FixedBed{N, RealT, IsoT, MtT, ThT}(
        component_names,
        SVector{N, RealT}(mm),
        SVector{N, RealT}(heat_capacity_gas),
        SVector{N, RealT}(heat_capacity_adsorbed),
        isotherm,
        mass_transfer,
        thermal,
    )
end

# Terms the shared conservation laws pick up for a fixed bed
mass_source_terms(::FixedBed) = (SorptionUptake(),)
function heat_capacity_terms(sys::FixedBed)
    return (SolidHeatCapacity(), AdsorbedPhaseHeatCapacity())
end
function energy_source_terms(sys::FixedBed)
    return has_wall(sys) ? (HeatOfSorption(), WallExchange()) : (HeatOfSorption(),)
end
wall_source_terms(::FixedBed) = (WallExchange(), AmbientLoss(), WallEndConduction())
