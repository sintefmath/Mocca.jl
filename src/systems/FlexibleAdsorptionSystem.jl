using Parameters

"""
    FlexibleAdsorptionSystem{T,I}

A flexible two-component adsorption system that supports pluggable isotherm models.

# Fields  
- `component_names`: Names of the components
- `permeability`: Permeability of the packed bed [m²]
- `dispersion`: Axial dispersion coefficient [m²/s] 
- `p`: Physical parameters struct (e.g., HaghpanahConstants)
- `isotherm_model`: Isotherm model implementing AbstractIsothermModel interface
"""
struct FlexibleAdsorptionSystem{T,I} <: AdsorptionSystem where {T<:ConstantsStruct,I<:AbstractIsothermModel}
    component_names::Vector{String}
    permeability::Float64
    dispersion::Float64
    p::T
    isotherm_model::I
end

JutulDarcy.number_of_components(sys::FlexibleAdsorptionSystem) = length(sys.component_names)
JutulDarcy.number_of_phases(sys::FlexibleAdsorptionSystem) = 1
JutulDarcy.get_reference_phase_index(sys::FlexibleAdsorptionSystem) = 1
JutulDarcy.eachphase(sys::FlexibleAdsorptionSystem) = (1,)

# Convenience constructors with different names to avoid ambiguity

# Constructor that creates a system with Haghpanah dual-site isotherm (original behavior)
function FlexibleAdsorptionSystemWithHaghpanah(; permeability, dispersion, p::T,
    component_names::Vector{String}=["CO2", "N2"]) where T<:ConstantsStruct
    isotherm_model = HaghpanahDualSiteIsotherm(p)
    return FlexibleAdsorptionSystem(
        component_names,
        permeability,
        dispersion,
        p,
        isotherm_model
    )
end

# Constructor with Langmuir.jl backend - keyword arguments version
function FlexibleAdsorptionSystemWithLangmuir(; permeability, dispersion, p::T,
    use_iast::Bool=true,
    component_names::Vector{String}=["CO2", "N2"]) where T<:ConstantsStruct
    isotherm_model = LangmuirIsothermModel(p; use_iast=use_iast)
    return FlexibleAdsorptionSystem(
        component_names,
        permeability,
        dispersion,
        p,
        isotherm_model
    )
end

# Simplified constructor with Langmuir.jl backend - takes only constants
function FlexibleAdsorptionSystemWithLangmuir(constants::ConstantsStruct; use_iast::Bool=true)
    permeability = compute_permeability(constants)
    axial_dispersion = calc_dispersion(constants)
    return FlexibleAdsorptionSystemWithLangmuir(
        permeability = permeability,
        dispersion = axial_dispersion,
        p = constants,
        use_iast = use_iast
    )
end

# Generic constructor with explicit isotherm model
function FlexibleAdsorptionSystem(permeability, dispersion, p::T, isotherm_model::I;
    component_names::Vector{String}=["CO2", "N2"]) where {T<:ConstantsStruct,I<:AbstractIsothermModel}
    return FlexibleAdsorptionSystem(
        component_names,
        permeability,
        dispersion,
        p,
        isotherm_model
    )
end
