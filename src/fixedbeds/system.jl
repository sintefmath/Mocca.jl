"""
    FixedBed{N, RealT, IsoT, MtT, ThT, GasT, AdsT} <: SorbentBed

A fixed (packed) adsorbent bed: gas flows axially through a stationary sorbent.
How energy is treated is set by the thermal model (see
[`AbstractThermalModel`](@ref)); the default [`WithWall`](@ref) includes a column
wall exchanging heat with the gas and the surroundings. The heat capacities of
the gas and of the adsorbed phase come from two
[`AbstractGasThermodynamics`](@ref) models.

`AdsorptionSystem` is kept as an alias for backwards compatibility.
"""
struct FixedBed{N, RealT<:Real, IsoT<:AbstractIsotherm, MtT<:AbstractMassTransfer, ThT<:AbstractThermalModel,
        GasT<:AbstractGasThermodynamics, AdsT<:AbstractGasThermodynamics} <: SorbentBed
    component_names::Vector{String}
    molecular_masses::SVector{N, RealT}         # Molecular masses per component [kg/mol]
    thermodynamics::GasT                        # Heat capacities and equation of state of the gas
    adsorbed_thermodynamics::AdsT               # Heat capacities of the adsorbed phase
    isotherm::IsoT
    mass_transfer::MtT
    thermal::ThT
end

const AdsorptionSystem = FixedBed
const FixedBedModel = Jutul.SimulationModel{<:Any,<:FixedBed,<:Any,<:Any}
const AdsorptionModel = FixedBedModel

"""
    FixedBed(; isotherm, mass_transfer, molecular_masses, component_names,
        thermodynamics, adsorbed_thermodynamics = thermodynamics,
        thermal = WithWall())

Construct an N-component adsorption system from explicit physics objects.
The number of components is inferred from the length of `component_names`.

`thermodynamics` and `adsorbed_thermodynamics` are
[`AbstractGasThermodynamics`](@ref) models of the gas and the adsorbed phase.
Instead of them, `heat_capacity_gas` and `heat_capacity_adsorbed` give
constant heat capacities per mass [J/(kg·K)] for each component, as
[`IdealGasConstantCp`](@ref) models.
"""
function FixedBed(;
    isotherm::AbstractIsotherm,
    mass_transfer::AbstractMassTransfer,
    molecular_masses,
    component_names,
    thermodynamics = nothing,
    adsorbed_thermodynamics = nothing,
    heat_capacity_gas = nothing,
    heat_capacity_adsorbed = nothing,
    thermal::AbstractThermalModel = WithWall(),
)
    N = length(component_names)
    @assert length(molecular_masses) == N "molecular_masses length ($(length(molecular_masses))) must match component_names length ($N)"
    gas = _thermodynamics(thermodynamics, heat_capacity_gas, molecular_masses, "gas")
    isnothing(gas) && error("FixedBed needs the thermodynamics of the gas, or its heat_capacity_gas")
    adsorbed = something(_thermodynamics(adsorbed_thermodynamics, heat_capacity_adsorbed, molecular_masses, "adsorbed"), gas)
    @assert number_of_components(gas) == N "the gas thermodynamics must have $N components"
    @assert number_of_components(adsorbed) == N "the adsorbed thermodynamics must have $N components"
    mm = SVector{N}(molecular_masses)
    RealT = eltype(mm)
    IsoT, MtT, ThT = typeof(isotherm), typeof(mass_transfer), typeof(thermal)
    return FixedBed{N, RealT, IsoT, MtT, ThT, typeof(gas), typeof(adsorbed)}(
        component_names,
        mm,
        gas,
        adsorbed,
        isotherm,
        mass_transfer,
        thermal,
    )
end

# A thermodynamics model, or constant heat capacities per mass that make one
function _thermodynamics(model, heat_capacities, molar_masses, phase)
    isnothing(model) || isnothing(heat_capacities) || error("give the $phase thermodynamics or its heat capacities, not both")
    isnothing(heat_capacities) || return IdealGasConstantCp(; molar_masses = molar_masses, heat_capacities = heat_capacities)
    return model
end

# Terms the shared conservation laws pick up for a fixed bed
Framework.mass_source_terms(::FixedBed) = (SorptionUptake(),)
function Framework.heat_capacity_terms(sys::FixedBed)
    return (SolidHeatCapacity(), AdsorbedPhaseHeatCapacity())
end
function Framework.energy_source_terms(sys::FixedBed)
    return has_wall(sys) ? (HeatOfSorption(), WallExchange()) : (HeatOfSorption(),)
end
Framework.wall_source_terms(::FixedBed) = (WallExchange(), AmbientLoss(), WallEndConduction())
