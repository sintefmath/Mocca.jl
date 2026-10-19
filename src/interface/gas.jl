# Components, gas mixtures and feeds: the first layer of the user interface.
#
# Everything here is plain SI: molar masses in kg/mol, heat capacities in
# J/(kg·K), viscosities in Pa·s, densities in kg/m³, diffusivities in m²/s,
# pressures in Pa and temperatures in K. See Mocca.Units for conversions.

"""
    Component(name; formula = name, cas = "", molar_mass, heat_capacity, viscosity)

A chemical component: its molar mass [kg/mol], gas heat capacity [J/(kg·K)]
and gas viscosity [Pa·s], both near 25 °C. Mocca knows common gases by name
(see [`component`](@ref)); add others with [`register_component!`](@ref).
"""
struct Component
    name::String
    formula::String
    cas::String
    molar_mass::Float64
    heat_capacity::Float64
    viscosity::Float64
end
Component(name; formula = name, cas = "", molar_mass, heat_capacity, viscosity) =
    Component(String(name), String(formula), String(cas), molar_mass, heat_capacity, viscosity)

const COMPONENTS = Dict{String, Component}()

"""
    register_component!(c::Component)

Make `c` known by name, so that gases, compositions and isotherms can use it.
Replaces a component of the same name.
"""
register_component!(c::Component) = (COMPONENTS[c.name] = c; c)

"""
    component(name) → Component

The component called `name`, for example `"CO2"`. Known: $(join(["CO2", "N2", "O2", "Ar", "H2O", "CH4", "H2", "CO", "He"], ", ")),
and any added with [`register_component!`](@ref).
"""
function component(name)
    n = String(name)
    haskey(COMPONENTS, n) || input_error("component $n", "", "unknown component; known components are $(join(sort(collect(keys(COMPONENTS))), ", ")). Add it with register_component!")
    return COMPONENTS[n]
end

# Molar masses from IUPAC atomic weights; heat capacities and viscosities of the
# gas at 25 °C and 1 atm
for c in [
        Component("CO2"; formula = "CO2", cas = "124-38-9", molar_mass = 44.01e-3, heat_capacity = 846.0, viscosity = 1.49e-5),
        Component("N2"; formula = "N2", cas = "7727-37-9", molar_mass = 28.0134e-3, heat_capacity = 1040.0, viscosity = 1.78e-5),
        Component("O2"; formula = "O2", cas = "7782-44-7", molar_mass = 31.998e-3, heat_capacity = 918.0, viscosity = 2.06e-5),
        Component("Ar"; formula = "Ar", cas = "7440-37-1", molar_mass = 39.948e-3, heat_capacity = 520.0, viscosity = 2.26e-5),
        Component("H2O"; formula = "H2O", cas = "7732-18-5", molar_mass = 18.015e-3, heat_capacity = 1864.0, viscosity = 0.99e-5),
        Component("CH4"; formula = "CH4", cas = "74-82-8", molar_mass = 16.043e-3, heat_capacity = 2226.0, viscosity = 1.11e-5),
        Component("H2"; formula = "H2", cas = "1333-74-0", molar_mass = 2.016e-3, heat_capacity = 14300.0, viscosity = 0.89e-5),
        Component("CO"; formula = "CO", cas = "630-08-0", molar_mass = 28.010e-3, heat_capacity = 1040.0, viscosity = 1.77e-5),
        Component("He"; formula = "He", cas = "7440-59-7", molar_mass = 4.0026e-3, heat_capacity = 5193.0, viscosity = 1.99e-5),
    ]
    register_component!(c)
end

"""
    Gas(names; molar_masses, heat_capacities, viscosity, density, diffusivity = 1.6e-5)

A gas mixture of the named components, in order. The order sets the order of
compositions and loadings everywhere; by convention the first component is the
one to capture (it is the key component of purity and recovery).

Properties default to those of the components (see [`component`](@ref)), and
can be overridden, all at once as a vector in gas order or for some
components by name, such as `heat_capacities = (CO2 = 697.6,)`:

- `molar_masses` [kg/mol] and `heat_capacities` [J/(kg·K)] per component;
- `viscosity` [Pa·s] and `density` [kg/m³] of the mixture. By default they
  are those of the feed: the Herning–Zipperer mixture of the component
  viscosities, and the ideal gas density at the feed pressure and temperature;
- `diffusivity` [m²/s]: molecular diffusivity, one value or one per
  component. It sets the macropore diffusion of the default linear driving
  force and, through the first component, the axial dispersion. The default is
  that of CO2 in N2 at 25 °C.
"""
struct Gas
    names::Vector{String}
    molar_masses::Vector{Float64}
    heat_capacities::Vector{Float64}
    viscosity::Union{Nothing, Float64}
    density::Union{Nothing, Float64}
    diffusivity::Union{Float64, Vector{Float64}}
end

function Gas(names::AbstractVector; molar_masses = nothing, heat_capacities = nothing, viscosity = nothing, density = nothing, diffusivity = 1.6e-5)
    names = String.(names)
    isempty(names) && input_error("gas", "names", "a gas needs at least one component")
    allunique(names) || input_error("gas", "names", "components appear more than once: $(join(names, ", "))")
    comps = component.(names)
    M = _by_component(names, molar_masses, [c.molar_mass for c in comps], "molar_masses")
    C = _by_component(names, heat_capacities, [c.heat_capacity for c in comps], "heat_capacities")
    D = diffusivity isa Number ? Float64(diffusivity) : _by_component(names, diffusivity, fill(NaN, length(names)), "diffusivity")
    any(isnan, D) && input_error("gas", "diffusivity", "give one diffusivity, or one for every component")
    return Gas(names, M, C, _maybe_float(viscosity), _maybe_float(density), D)
end
Gas(names::Union{AbstractString, Symbol}...; kwarg...) = Gas(collect(names); kwarg...)

_maybe_float(x) = isnothing(x) ? nothing : Float64(x)

component_names(gas::Gas) = gas.names
number_of_components(gas::Gas) = length(gas.names)

function component_index(gas::Gas, name)
    i = findfirst(==(String(name)), gas.names)
    isnothing(i) && input_error("gas", "", "$name is not a component of the gas ($(join(gas.names, ", ")))")
    return i
end

# A per-component setting: `nothing` (the defaults), a vector in gas order, or
# values for some components by name
function _by_component(names, x, default, what)
    isnothing(x) && return Vector{Float64}(default)
    v = Vector{Float64}(default)
    if x isa AbstractVector && !(eltype(x) <: Pair)
        length(x) == length(names) || input_error("gas", what, "expected $(length(names)) values, one per component ($(join(names, ", "))), got $(length(x))")
        v .= x
    else
        for (k, val) in pairs_of(x)
            i = findfirst(==(String(k)), names)
            isnothing(i) && input_error("gas", what, "$k is not a component of the gas ($(join(names, ", ")))")
            v[i] = val
        end
    end
    return v
end

pairs_of(x::NamedTuple) = pairs(x)
pairs_of(x::AbstractDict) = pairs(x)
pairs_of(x::AbstractVector{<:Pair}) = x
pairs_of(x::Pair) = (x,)
pairs_of(x) = pairs(x)

"""
    composition(gas, spec; trace = 1e-10) → Vector{Float64}

Mole fractions in gas order. `spec` is a vector in gas order, or fractions by
name, such as `(CO2 = 0.15, N2 = 0.85)` or `Dict("CO2" => 0.15, "N2" => 0.85)`.
Components not named get `trace`, to keep every fraction positive, and the
named ones are scaled down to make room. Fractions must sum to one.
"""
function composition(gas::Gas, spec; trace = 1e-10, location = "composition")
    N = number_of_components(gas)
    if spec isa AbstractVector && !(eltype(spec) <: Pair)
        length(spec) == N || input_error(location, "", "expected $N mole fractions ($(join(gas.names, ", "))), got $(length(spec))")
        y = Vector{Float64}(spec)
    else
        y = fill(NaN, N)
        for (k, v) in pairs_of(spec)
            i = findfirst(==(String(k)), gas.names)
            isnothing(i) && input_error(location, "", "$k is not a component of the gas ($(join(gas.names, ", ")))")
            y[i] = v
        end
        named = .!isnan.(y)
        total = sum(y[named])
        abs(total - 1) < 1e-6 || input_error(location, "", "mole fractions sum to $total, not 1")
        missing_ = count(!, named)
        y[named] .*= (1 - missing_ * trace)
        y[.!named] .= trace
    end
    any(<(0), y) && input_error(location, "", "mole fractions must not be negative: $y")
    abs(sum(y) - 1) < 1e-6 || input_error(location, "", "mole fractions sum to $(sum(y)), not 1")
    return y
end

"""
    Feed(gas; composition, temperature = 298.15, pressure = 101325.0, light_product)

The feed gas of a flowsheet: its composition (see [`composition`](@ref)),
temperature [K] and pressure [Pa]. `light_product` is the composition of the
gas leaving the top of a bed during adsorption, used where light product flows
back into a bed (for example in light product pressurisation); by default the
feed with its first component reduced to a trace.
"""
mutable struct Feed
    composition::Vector{Float64}
    temperature::Float64
    pressure::Float64
    light_product::Vector{Float64}
end

function Feed(gas::Gas; composition, temperature = 298.15, pressure = 101325.0, light_product = nothing)
    y = Mocca.composition(gas, composition; location = "feed")
    yl = isnothing(light_product) ? light_product_composition(y) : Mocca.composition(gas, light_product; location = "feed light product")
    temperature > 0 || input_error("feed", "temperature", "must be positive (in K), got $temperature")
    pressure > 0 || input_error("feed", "pressure", "must be positive (in Pa), got $pressure")
    return Feed(y, Float64(temperature), Float64(pressure), yl)
end

"Mixture viscosity [Pa·s] of `gas` at composition `y`: the gas value if set, else Herning–Zipperer."
function gas_viscosity(gas::Gas, y)
    isnothing(gas.viscosity) || return gas.viscosity
    μ = [component(n).viscosity for n in gas.names]
    w = y .* sqrt.(gas.molar_masses)
    return sum(w .* μ) / sum(w)
end

"Gas density [kg/m³]: the gas value if set, else the ideal gas at pressure `P` and temperature `T`."
function gas_density(gas::Gas, y, P, T)
    isnothing(gas.density) || return gas.density
    return P * sum(y .* gas.molar_masses) / (GAS_CONSTANT * T)
end

"Diffusivity [m²/s] for the axial dispersion: the gas value, or that of its first component."
axial_diffusivity(gas::Gas) = gas.diffusivity isa Number ? gas.diffusivity : gas.diffusivity[1]

function Base.show(io::IO, ::MIME"text/plain", gas::Gas)
    println(io, "Gas of ", join(gas.names, ", "))
    for (i, n) in enumerate(gas.names)
        println(io, "  ", rpad(n, 6), "M = ", round(1e3 * gas.molar_masses[i], sigdigits = 5), " g/mol, cp = ", round(gas.heat_capacities[i], sigdigits = 5), " J/(kg·K)")
    end
    print(io, "  viscosity ", isnothing(gas.viscosity) ? "from the feed" : "$(gas.viscosity) Pa·s",
        ", density ", isnothing(gas.density) ? "from the feed" : "$(gas.density) kg/m³",
        ", diffusivity ", gas.diffusivity, " m²/s")
end
