# Components, gas mixtures and feeds: the first layer of the user interface.
#
# Everything here is plain SI: molar masses in kg/mol, heat capacities in
# J/(kg·K) (J/(mol·K) in the thermodynamics models), viscosities in Pa·s,
# thermal conductivities in W/(m·K), densities in kg/m³, diffusivities in
# m²/s, pressures in Pa and temperatures in K. See Mocca.Units for conversions.

"""
    Component(name; formula = name, cas = "", molar_mass, heat_capacity, viscosity,
        shomate = nothing, sutherland = nothing, thermal_conductivity = nothing,
        conductivity_exponent = nothing)

A chemical component: its molar mass [kg/mol], gas heat capacity [J/(kg·K)]
and gas viscosity [Pa·s], both near 25 °C, and optionally:

- `shomate`: the coefficients `(A, B, C, D, E)` of the Shomate equation for
  its molar heat capacity (see [`IdealGasShomate`](@ref));
- `sutherland` [K], `thermal_conductivity` [W/(m·K)] at 25 °C and
  `conductivity_exponent`, for viscosity and conductivity that change with
  temperature (see [`SutherlandTransport`](@ref)).

Mocca knows common gases by name (see [`component`](@ref)); add others with
[`register_component!`](@ref).
"""
struct Component
    name::String
    formula::String
    cas::String
    molar_mass::Float64
    heat_capacity::Float64
    viscosity::Float64
    shomate::Union{Nothing, NTuple{5, Float64}}
    sutherland::Union{Nothing, Float64}
    thermal_conductivity::Union{Nothing, Float64}
    conductivity_exponent::Union{Nothing, Float64}
end
Component(name; formula = name, cas = "", molar_mass, heat_capacity, viscosity, shomate = nothing, sutherland = nothing,
        thermal_conductivity = nothing, conductivity_exponent = nothing) =
    Component(String(name), String(formula), String(cas), molar_mass, heat_capacity, viscosity, _shomate_tuple(shomate),
        _maybe_float(sutherland), _maybe_float(thermal_conductivity), _maybe_float(conductivity_exponent))

_maybe_float(x) = isnothing(x) ? nothing : Float64(x)

_shomate_tuple(::Nothing) = nothing
function _shomate_tuple(x)
    length(x) == 5 || input_error("component", "shomate", "expected the five coefficients (A, B, C, D, E), got $(length(x))")
    return NTuple{5, Float64}(Tuple(Float64.(x)))
end

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
# gas at 25 °C and 1 atm. Shomate coefficients from the NIST Chemistry WebBook,
# for the range from 298 K (100 K for N2 and O2) upwards; that of H2O is fitted
# from 500 K and is within 0.1% of the tabulated values down to 298 K.
# Sutherland constants and conductivity exponents fitted to the viscosities and
# thermal conductivities of the NIST Chemistry WebBook at 300 and 500 K and
# 0.01 bar: within 0.6% at 400 K (water 1.8%), and within 1.6% for the
# conductivities; conductivities at 25 °C from the same fit.
for c in [
        Component("CO2"; formula = "CO2", cas = "124-38-9", molar_mass = 44.01e-3, heat_capacity = 846.0, viscosity = 1.49e-5,
            shomate = (24.99735, 55.18696, -33.69137, 7.948387, -0.136638), sutherland = 273.0,
            thermal_conductivity = 0.01658, conductivity_exponent = 1.32),
        Component("N2"; formula = "N2", cas = "7727-37-9", molar_mass = 28.0134e-3, heat_capacity = 1040.0, viscosity = 1.78e-5,
            shomate = (28.98641, 1.853978, -9.647459, 16.63537, 0.000117), sutherland = 120.0,
            thermal_conductivity = 0.02581, conductivity_exponent = 0.8),
        Component("O2"; formula = "O2", cas = "7782-44-7", molar_mass = 31.998e-3, heat_capacity = 918.0, viscosity = 2.06e-5,
            shomate = (31.32234, -20.23531, 57.86644, -36.50624, -0.007374), sutherland = 138.0,
            thermal_conductivity = 0.0263, conductivity_exponent = 0.86),
        Component("Ar"; formula = "Ar", cas = "7440-37-1", molar_mass = 39.948e-3, heat_capacity = 520.0, viscosity = 2.26e-5,
            shomate = (20.786, 2.825911e-7, -1.464191e-7, 1.092131e-8, -3.661371e-8), sutherland = 160.0,
            thermal_conductivity = 0.01772, conductivity_exponent = 0.79),
        Component("H2O"; formula = "H2O", cas = "7732-18-5", molar_mass = 18.015e-3, heat_capacity = 1864.0, viscosity = 0.99e-5,
            shomate = (30.092, 6.832514, 6.793435, -2.53448, 0.082139), sutherland = 640.0,
            thermal_conductivity = 0.01842, conductivity_exponent = 1.29),
        Component("CH4"; formula = "CH4", cas = "74-82-8", molar_mass = 16.043e-3, heat_capacity = 2226.0, viscosity = 1.11e-5,
            shomate = (-0.703029, 108.4773, -42.52157, 5.862788, 0.678565), sutherland = 180.0,
            thermal_conductivity = 0.03386, conductivity_exponent = 1.36),
        Component("H2"; formula = "H2", cas = "1333-74-0", molar_mass = 2.016e-3, heat_capacity = 14300.0, viscosity = 0.89e-5,
            shomate = (33.066178, -11.363417, 11.432816, -2.772874, -0.158558), sutherland = 92.0,
            thermal_conductivity = 0.1857, conductivity_exponent = 0.73),
        Component("CO"; formula = "CO", cas = "630-08-0", molar_mass = 28.010e-3, heat_capacity = 1040.0, viscosity = 1.77e-5,
            shomate = (25.56759, 6.09613, 4.054656, -2.671301, 0.131021), sutherland = 117.0,
            thermal_conductivity = 0.02487, conductivity_exponent = 0.81),
        Component("He"; formula = "He", cas = "7440-59-7", molar_mass = 4.0026e-3, heat_capacity = 5193.0, viscosity = 1.99e-5,
            shomate = (20.78603, 4.850638e-10, -1.582916e-10, 1.525102e-11, 3.196347e-11), sutherland = 91.0,
            thermal_conductivity = 0.1552, conductivity_exponent = 0.69),
    ]
    register_component!(c)
end

"""
    Gas(names; molar_masses, heat_capacities, thermodynamics, shomate, viscosity,
        density, diffusivity = 1.6e-5, transport = :constant)

A gas mixture of the named components, in order. The order sets the order of
compositions and loadings everywhere; by convention the first component is the
one to capture (it is the key component of purity and recovery).

Properties default to those of the components (see [`component`](@ref)), and
can be overridden, all at once as a vector in gas order or for some
components by name, such as `heat_capacities = (CO2 = 697.6,)`:

- `molar_masses` [kg/mol] per component;
- the heat capacities, through the `thermodynamics` of the gas (see
  [`AbstractGasThermodynamics`](@ref)):
  - by default constant, `heat_capacities` [J/(kg·K)] per component
    ([`IdealGasConstantCp`](@ref));
  - `thermodynamics = :shomate`: changing with temperature by the Shomate
    equation ([`IdealGasShomate`](@ref)), with the coefficients of the
    components, or `shomate`, coefficients `(A, B, C, D, E)` per component.
    Giving `shomate` chooses this;
  - or an [`AbstractGasThermodynamics`](@ref) model for the components in
    gas order.

  The heat capacity of the mixture is the mole average of the molar heat
  capacities, ``\\sum_i y_i M_i c_{p,i}``;
- `viscosity` [Pa·s] of the mixture at the feed. By default that of the feed,
  the Herning–Zipperer mixture of the component viscosities (at the feed
  temperature with a `transport` model);
- `density` [kg/m³] of the mixture, by default the ideal gas at the feed
  pressure and temperature. The built-in units do not use it;
- `diffusivity` [m²/s]: molecular diffusivity at the feed, one value or one per
  component. It sets the macropore diffusion of the default linear driving
  force and, through the first component, the axial dispersion. The default is
  that of CO2 in N2 at 25 °C and 1 atm;
- `transport`: how the viscosity, the diffusivity and the thermal
  conductivity of a bed (its column's `axial_conductivity`) change from their
  values at the feed with the temperature, pressure and composition of each
  cell (see [`AbstractGasTransport`](@ref)). `:constant` (the default) keeps
  them; `:sutherland` uses the viscosities, Sutherland constants and thermal
  conductivities of the components ([`SutherlandTransport`](@ref)); or give an
  [`AbstractGasTransport`](@ref) model for the components in gas order.

`gas.heat_capacities` holds the heat capacities [J/(kg·K)] per component at
25 °C, `gas.thermodynamics` the thermodynamics model and `gas.transport` the
transport model.
"""
struct Gas
    names::Vector{String}
    molar_masses::Vector{Float64}
    heat_capacities::Vector{Float64}
    thermodynamics::AbstractGasThermodynamics
    transport::AbstractGasTransport
    viscosity::Union{Nothing, Float64}
    density::Union{Nothing, Float64}
    diffusivity::Union{Float64, Vector{Float64}}
end

function Gas(names::AbstractVector; molar_masses = nothing, heat_capacities = nothing, thermodynamics = nothing, shomate = nothing,
        viscosity = nothing, density = nothing, diffusivity = 1.6e-5, transport = :constant)
    names = String.(names)
    isempty(names) && input_error("gas", "names", "a gas needs at least one component")
    allunique(names) || input_error("gas", "names", "components appear more than once: $(join(names, ", "))")
    comps = component.(names)
    M = _by_component(names, molar_masses, [c.molar_mass for c in comps], "molar_masses")
    thermo = _gas_thermodynamics(names, comps, M, thermodynamics, heat_capacities, shomate)
    # Per mass at 25 °C: as given for constant heat capacities, else from the model
    C = thermo isa IdealGasConstantCp && !(thermodynamics isa AbstractGasThermodynamics) ?
        _by_component(names, heat_capacities, [c.heat_capacity for c in comps], "heat_capacities") :
        collect(component_heat_capacities(thermo, THERMODYNAMIC_REFERENCE_TEMPERATURE)) ./ M
    D = diffusivity isa Number ? Float64(diffusivity) : _by_component(names, diffusivity, fill(NaN, length(names)), "diffusivity")
    any(isnan, D) && input_error("gas", "diffusivity", "give one diffusivity, or one for every component")
    tr = _gas_transport(names, comps, M, transport)
    return Gas(names, M, C, thermo, tr, _maybe_float(viscosity), _maybe_float(density), D)
end
Gas(names::Union{AbstractString, Symbol}...; kwarg...) = Gas(collect(names); kwarg...)

# The transport model of a gas from the keyword of Gas
function _gas_transport(names, comps, M, transport)
    transport isa AbstractGasTransport && return transport
    transport == :constant && return ConstantTransport()
    transport == :sutherland || input_error("gas", "transport", "use :constant, :sutherland or an AbstractGasTransport, got $(repr(transport))")
    for (field, what) in ((:sutherland, "Sutherland constant"), (:thermal_conductivity, "thermal conductivity"),
            (:conductivity_exponent, "conductivity exponent"))
        missing_ = [c.name for c in comps if isnothing(getfield(c, field))]
        isempty(missing_) || input_error("gas", "transport", "no $what for $(join(missing_, ", ")); register the component with its `$field`")
    end
    return SutherlandTransport(; viscosity = [c.viscosity for c in comps], sutherland = [c.sutherland for c in comps],
        thermal_conductivity = [c.thermal_conductivity for c in comps], conductivity_exponent = [c.conductivity_exponent for c in comps],
        molar_masses = M)
end

# The thermodynamics model of a gas from the keywords of Gas
function _gas_thermodynamics(names, comps, M, model, heat_capacities, shomate)
    N = length(names)
    if model isa AbstractGasThermodynamics
        isnothing(heat_capacities) && isnothing(shomate) ||
            input_error("gas", "thermodynamics", "a thermodynamics model sets the heat capacities; do not also give heat_capacities or shomate")
        number_of_components(model) == N || input_error("gas", "thermodynamics", "the model has $(number_of_components(model)) components, the gas $N")
        return model
    end
    kind = isnothing(model) ? (isnothing(shomate) ? :constant : :shomate) : model
    if kind == :constant
        isnothing(shomate) || input_error("gas", "shomate", "Shomate coefficients need thermodynamics = :shomate")
        C = _by_component(names, heat_capacities, [c.heat_capacity for c in comps], "heat_capacities")
        return IdealGasConstantCp(; molar_masses = M, heat_capacities = C)
    elseif kind == :shomate
        isnothing(heat_capacities) ||
            input_error("gas", "heat_capacities", "sets constant heat capacities; with Shomate heat capacities give the coefficients through shomate")
        coefficients = Vector{Union{Nothing, NTuple{5, Float64}}}([c.shomate for c in comps])
        if !isnothing(shomate)
            if shomate isa AbstractVector && !(eltype(shomate) <: Pair)
                length(shomate) == N || input_error("gas", "shomate", "expected $N sets of coefficients, one per component ($(join(names, ", "))), got $(length(shomate))")
                coefficients .= _shomate_tuple.(shomate)
            else
                for (k, v) in pairs_of(shomate)
                    i = findfirst(==(String(k)), names)
                    isnothing(i) && input_error("gas", "shomate", "$k is not a component of the gas ($(join(names, ", ")))")
                    coefficients[i] = _shomate_tuple(v)
                end
            end
        end
        missing_ = names[isnothing.(coefficients)]
        isempty(missing_) || input_error("gas", "shomate", "no Shomate coefficients for $(join(missing_, ", ")); give them through shomate")
        return IdealGasShomate(collect(NTuple{5, Float64}, coefficients))
    end
    input_error("gas", "thermodynamics", "use :constant, :shomate or an AbstractGasThermodynamics, got $model")
end

Framework.component_names(gas::Gas) = gas.names
Framework.number_of_components(gas::Gas) = length(gas.names)

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
    light_product_composition(y_feed; component = 1, trace = 1e-10)

The feed composition with `component` (the one captured) reduced to `trace`
and the rest renormalised: a stand-in for the light product of a separation.
"""
function light_product_composition(y_feed; component = 1, trace = 1e-10)
    rest = sum(y_feed) - y_feed[component]
    y = collect(Float64, y_feed) .* ((1 - trace) / rest)
    y[component] = trace
    return y
end

"""
    Feed(gas; composition, temperature = 298.15, pressure = 101325.0, light_product)

The feed gas of a flowsheet: its composition (see [`composition`](@ref)),
temperature [K] and pressure [Pa]. `light_product` is the composition of the
gas leaving the top of a bed during adsorption, used where light product flows
back into a bed (for example in light product pressurisation). By default it
is the feed with the component captured reduced to a trace: the first
component, or the `key_component` of the [`Process`](@ref) the feed is used in.
"""
mutable struct Feed
    composition::Vector{Float64}
    temperature::Float64
    pressure::Float64
    light_product::Vector{Float64}
    # Whether `light_product` was given, rather than derived from the key
    # component, which a Process may change
    light_product_given::Bool
end

function Feed(gas::Gas; composition, temperature = 298.15, pressure = 101325.0, light_product = nothing)
    y = Interface.composition(gas, composition; location = "feed")
    yl = isnothing(light_product) ? light_product_composition(y) : Interface.composition(gas, light_product; location = "feed light product")
    temperature > 0 || input_error("feed", "temperature", "must be positive (in K), got $temperature")
    pressure > 0 || input_error("feed", "pressure", "must be positive (in Pa), got $pressure")
    return Feed(y, Float64(temperature), Float64(pressure), yl, !isnothing(light_product))
end

"""
    gas_viscosity(gas, y, T = 298.15)

Mixture viscosity [Pa·s] of `gas` at composition `y` and temperature `T` [K]:
the gas value if set; else that of its transport model, or with constant
transport the Herning–Zipperer mixture of the component viscosities at 25 °C.
"""
function gas_viscosity(gas::Gas, y, T = THERMODYNAMIC_REFERENCE_TEMPERATURE)
    isnothing(gas.viscosity) || return gas.viscosity
    gas.transport isa ConstantTransport || return mixture_viscosity(gas.transport, y, T)
    μ = [component(n).viscosity for n in gas.names]
    w = y .* sqrt.(gas.molar_masses)
    return sum(w .* μ) / sum(w)
end

"Gas density [kg/m³]: the gas value if set, else from the equation of state at pressure `P` and temperature `T`."
function gas_density(gas::Gas, y, P, T)
    isnothing(gas.density) || return gas.density
    return molar_concentration(gas.thermodynamics, P, T) * sum(y .* gas.molar_masses)
end

"""
    gas_heat_capacity_ratio(gas, y, T)

Ratio of heat capacities ``\\gamma = c_p/c_v`` of `gas` at composition `y` and
temperature `T` [K], from its thermodynamics (see [`heat_capacity_ratio`](@ref)).
"""
gas_heat_capacity_ratio(gas::Gas, y, T) = heat_capacity_ratio(gas.thermodynamics, y, T)

"Diffusivity [m²/s] for the axial dispersion: the gas value, or that of its first component."
axial_diffusivity(gas::Gas) = gas.diffusivity isa Number ? gas.diffusivity : gas.diffusivity[1]

function Base.show(io::IO, ::MIME"text/plain", gas::Gas)
    constant = gas.thermodynamics isa IdealGasConstantCp
    println(io, "Gas of ", join(gas.names, ", "), constant ? "" : ", heat capacities from $(nameof(typeof(gas.thermodynamics)))")
    for (i, n) in enumerate(gas.names)
        println(io, "  ", rpad(n, 6), "M = ", round(1e3 * gas.molar_masses[i], sigdigits = 5), " g/mol, cp", constant ? "" : " at 25 °C", " = ",
            round(gas.heat_capacities[i], sigdigits = 5), " J/(kg·K)")
    end
    gas.transport isa ConstantTransport || println(io, "  viscosity, conductivity and diffusivity from the feed by $(nameof(typeof(gas.transport)))")
    print(io, "  viscosity ", isnothing(gas.viscosity) ? "from the feed" : "$(gas.viscosity) Pa·s",
        ", density ", isnothing(gas.density) ? "from the feed" : "$(gas.density) kg/m³",
        ", diffusivity ", gas.diffusivity, " m²/s")
end
