# Isotherms by component, adsorbents and columns: the material and geometry
# layer of the user interface. Plain SI throughout (see Mocca.Units).

# ------------------------------------------------------------------------------
# Isotherm data by component
# ------------------------------------------------------------------------------

"""
    AbstractIsothermData

Isotherm parameters given by component name, independent of the order of the
components in a gas. [`build_isotherm`](@ref) turns them into an
[`AbstractIsotherm`](@ref) for a particular [`Gas`](@ref). Components of the
gas without parameters do not adsorb; parameters for components the gas does
not have are ignored, so one adsorbent can serve gases with different
components.

Implement `build_isotherm(data, gas; density)` for a new kind, and register
it with [`register_isotherm!`](@ref) to read and write it in files.
"""
abstract type AbstractIsothermData end

"""
    build_isotherm(data, gas; density) → AbstractIsotherm

The isotherm for the components of `gas`, in its order. `density` [kg/m³] is
the adsorbent density, which converts capacities given per kilogram. An
[`AbstractIsotherm`](@ref) is used as it is.
"""
function build_isotherm end

function build_isotherm(iso::AbstractIsotherm, gas::Gas; density = nothing)
    if iso isa DualSiteLangmuir
        length(iso.qsb) == number_of_components(gas) || input_error("isotherm", "", "it has $(length(iso.qsb)) components, the gas has $(number_of_components(gas)) ($(join(gas.names, ", ")))")
    end
    return iso
end

"""
    DualSiteLangmuirData

Parameters of the dual-site Langmuir isotherm (see [`DualSiteLangmuir`](@ref))
by component. Build it with [`dual_site_langmuir`](@ref) or
[`langmuir`](@ref).
"""
struct DualSiteLangmuirData <: AbstractIsothermData
    components::Vector{String}
    qsb::Vector{Float64}
    b0::Vector{Float64}
    ΔUb::Vector{Float64}
    qsd::Vector{Float64}
    d0::Vector{Float64}
    ΔUd::Vector{Float64}
    reference_temperature::Float64
    basis::Symbol
end

const _DSL_KEYS = (:qsb, :b0, :ΔUb, :qsd, :d0, :ΔUd)
const _DSL_ALIASES = Dict(:dUb => :ΔUb, :dUd => :ΔUd)

"""
    dual_site_langmuir(name => parameters, ...; reference_temperature = 298.15, basis = :volume)

Dual-site Langmuir isotherm data, one entry per component:

```julia
dual_site_langmuir(
    "CO2" => (qsb = 3489.44, b0 = 8.65e-7, ΔUb = -36.64kJ_per_mol, qsd = 2872.35, d0 = 2.63e-8, ΔUd = -35.69kJ_per_mol),
    "N2" => (qsb = 6613.55, b0 = 2.50e-6, ΔUb = -15.8kJ_per_mol))
```

For each site, b (strong) and d (weak): saturation capacity `qsb`, `qsd`;
affinity at infinite temperature `b0`, `d0` [m³/mol], on a concentration
basis; and internal energy of adsorption `ΔUb`, `ΔUd` [J/mol], negative.
Omitted d-site values are zero. `dUb` and `dUd` may be written for `ΔUb` and
`ΔUd`.

Capacities are in mol per m³ of adsorbent with `basis = :volume`, or mol/kg
with `basis = :mass`. `reference_temperature` [K] is where the isosteric heat
of adsorption is evaluated.
"""
function dual_site_langmuir(entries::Pair...; reference_temperature = 298.15, basis = :volume)
    basis in (:volume, :mass) || input_error("isotherm", "basis", "use :volume (mol per m³ of adsorbent) or :mass (mol/kg), got $basis")
    names = String[]
    cols = Dict(k => Float64[] for k in _DSL_KEYS)
    for (name, p) in entries
        n = String(name)
        component(n)  # catch misspelt names
        n in names && input_error("isotherm", n, "given more than once")
        push!(names, n)
        given = Dict{Symbol, Float64}()
        for (k, v) in pairs_of(p)
            key = get(_DSL_ALIASES, Symbol(k), Symbol(k))
            key in _DSL_KEYS || input_error("isotherm for $n", String(k), "unknown parameter; dual-site Langmuir takes $(join(_DSL_KEYS, ", "))")
            given[key] = v
        end
        haskey(given, :qsb) && haskey(given, :b0) && haskey(given, :ΔUb) || input_error("isotherm for $n", "", "qsb, b0 and ΔUb are required")
        for k in _DSL_KEYS
            push!(cols[k], get(given, k, 0.0))
        end
        cols[:qsb][end] >= 0 && cols[:qsd][end] >= 0 || input_error("isotherm for $n", "", "capacities must not be negative")
        cols[:b0][end] >= 0 && cols[:d0][end] >= 0 || input_error("isotherm for $n", "", "affinities must not be negative")
    end
    return DualSiteLangmuirData(names, cols[:qsb], cols[:b0], cols[:ΔUb], cols[:qsd], cols[:d0], cols[:ΔUd],
        Float64(reference_temperature), basis)
end

"""
    langmuir(name => (qs = ..., b0 = ..., ΔU = ...), ...; reference_temperature = 298.15, basis = :volume)

Single-site Langmuir isotherm data by component: a [`dual_site_langmuir`](@ref)
with an empty second site.
"""
function langmuir(entries::Pair...; kwarg...)
    converted = map(entries) do (name, p)
        d = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in pairs_of(p))
        keys_ok = issubset(keys(d), (:qs, :b0, :ΔU, :dU))
        keys_ok || input_error("isotherm for $name", "", "Langmuir takes qs, b0 and ΔU")
        name => (qsb = d[:qs], b0 = d[:b0], ΔUb = haskey(d, :ΔU) ? d[:ΔU] : d[:dU])
    end
    return dual_site_langmuir(converted...; kwarg...)
end

function build_isotherm(data::DualSiteLangmuirData, gas::Gas; density = nothing)
    N = number_of_components(gas)
    pick(v) = [(i = findfirst(==(n), data.components); isnothing(i) ? 0.0 : v[i]) for n in gas.names]
    scale = 1.0
    if data.basis == :mass
        isnothing(density) && input_error("isotherm", "basis", "capacities per kg need the adsorbent density")
        scale = density
    end
    any(n -> n in data.components, gas.names) || input_error("isotherm", "", "it has none of the gas components ($(join(gas.names, ", "))); it has $(join(data.components, ", "))")
    return DualSiteLangmuir(qsb = scale .* pick(data.qsb), b0 = pick(data.b0), ΔUb = pick(data.ΔUb),
        qsd = scale .* pick(data.qsd), d0 = pick(data.d0), ΔUd = pick(data.ΔUd), T0 = data.reference_temperature)
end

"""
    VirialData

Parameters of the multicomponent Virial isotherm (see [`Virial`](@ref)) by
component. Build it with [`virial`](@ref).
"""
struct VirialData <: AbstractIsothermData
    components::Vector{String}
    K∞::Vector{Float64}
    ΔH::Vector{Float64}
    A0::Vector{Float64}
    A1::Vector{Float64}
    B0::Vector{Float64}
    B1::Vector{Float64}
end

const _VIRIAL_KEYS = (:K∞, :ΔH, :A0, :A1, :B0, :B1)
const _VIRIAL_ALIASES = Dict(:Kinf => :K∞, :dH => :ΔH)

"""
    virial(name => parameters, ...)

Virial isotherm data, one entry per component, on a mass basis:

```julia
virial(
    "CO2" => (K∞ = 2.6969e-10, ΔH = 30.006kJ_per_mol, A0 = -4.3235, A1 = 1474.0, B0 = 1.4239, B1 = -465.20),
    "N2" => (K∞ = 5.5486e-10, ΔH = 21.934kJ_per_mol, A0 = -22.982, A1 = 7121.1, B0 = 51.644, B1 = -15756.0))
```

The partial pressure is `Pᵢ = qᵢ/K_H,ᵢ·exp(Σⱼ Aᵢⱼqⱼ + Σⱼₖ Bᵢⱼₖqⱼqₖ)`, with
loadings in mol/kg and:

- `K∞` [mol/(kg·Pa)], the Henry constant at infinite temperature, and `ΔH`
  [J/mol], the heat of adsorption at zero coverage, positive:
  `K_H = K∞·exp(ΔH/(R·T))`;
- `A0` [kg/mol] and `A1` [K·kg/mol]: `A = A0 + A1/T`;
- `B0` [kg²/mol²] and `B1` [K·kg²/mol²]: `B = B0 + B1/T`.

Omitted A and B values are zero. `Kinf` and `dH` may be written for `K∞` and
`ΔH`. Henry constants are often given per kPa; divide them by 1000.
"""
function virial(entries::Pair...)
    names = String[]
    cols = Dict(k => Float64[] for k in _VIRIAL_KEYS)
    for (name, p) in entries
        n = String(name)
        component(n)  # catch misspelt names
        n in names && input_error("isotherm", n, "given more than once")
        push!(names, n)
        given = Dict{Symbol, Float64}()
        for (k, v) in pairs_of(p)
            key = get(_VIRIAL_ALIASES, Symbol(k), Symbol(k))
            key in _VIRIAL_KEYS || input_error("isotherm for $n", String(k), "unknown parameter; the Virial isotherm takes $(join(_VIRIAL_KEYS, ", "))")
            given[key] = v
        end
        haskey(given, :K∞) && haskey(given, :ΔH) || input_error("isotherm for $n", "", "K∞ and ΔH are required")
        given[:K∞] > 0 || input_error("isotherm for $n", "K∞", "must be positive")
        for k in _VIRIAL_KEYS
            push!(cols[k], get(given, k, 0.0))
        end
    end
    return VirialData(names, (cols[k] for k in _VIRIAL_KEYS)...)
end

function build_isotherm(data::VirialData, gas::Gas; density = nothing)
    isnothing(density) && input_error("isotherm", "", "a Virial isotherm needs the adsorbent density, since its loadings are per kg")
    any(n -> n in data.components, gas.names) || input_error("isotherm", "", "it has none of the gas components ($(join(gas.names, ", "))); it has $(join(data.components, ", "))")
    missing_ = [n for n in gas.names if !(n in data.components)]
    isempty(missing_) || input_error("isotherm", "", "a Virial isotherm needs parameters for every component of the gas; it has none for $(join(missing_, ", "))")
    pick(v) = [v[findfirst(==(n), data.components)] for n in gas.names]
    return Virial(K∞ = pick(data.K∞), ΔH = pick(data.ΔH), A0 = pick(data.A0), A1 = pick(data.A1),
        B0 = pick(data.B0), B1 = pick(data.B1), density = density)
end

"""
    TothGABData

Parameters of the CO₂–water co-adsorption isotherm [`TothGAB`](@ref) by
component name. Build it with [`toth_gab`](@ref).
"""
struct TothGABData <: AbstractIsothermData
    toth::String
    water::String
    ns0::Float64
    b0::Float64
    t0::Float64
    α::Float64
    χ::Float64
    ΔH0::Float64
    reference_temperature::Float64
    γ::Float64
    β::Float64
    cm::Float64
    cG::Float64
    K::Float64
    heat_components::Vector{String}
    heats::Vector{Float64}
    basis::Symbol
end

"""
    toth_gab(; component = "CO2", water = "H2O", ns0, b0, t0, α = 0, χ = 0, ΔH0,
        reference_temperature, γ = 0, β = 0, cm, cG, K, heats, basis = :mass)

CO₂–water co-adsorption on an amine sorbent (see [`TothGAB`](@ref)): water on a
GAB isotherm with monolayer capacity `cm` and constants `cG` and `K`, and
`component` on a Toth isotherm with capacity `ns0`, affinity `b0` [1/Pa] and
exponent `t0` at the `reference_temperature` [K], their temperature
dependence `χ`, `ΔH0` [J/mol, negative] and `α`, and their enhancement by the
water loading, `γ` and `β`. `heats` gives the heat of adsorption [J/mol,
negative] by component, such as `(CO2 = -57e3, H2O = -49e3)`.

Capacities are in mol/kg and `γ`, `β` in kg/mol with `basis = :mass`, or per
m³ of adsorbent and in m³/mol with `basis = :volume`. Stampi-Bombelli et al.
(2020) give b0 in 1/kPa: divide it by 1000.
"""
function toth_gab(; component = "CO2", water = "H2O", ns0, b0, t0, α = 0.0, χ = 0.0, ΔH0, reference_temperature,
        γ = 0.0, β = 0.0, cm, cG, K, heats, basis = :mass)
    loc = "toth_gab isotherm"
    basis in (:volume, :mass) || input_error(loc, "basis", "use :mass (mol/kg) or :volume (mol per m³ of adsorbent), got $basis")
    c, w = String(component), String(water)
    Interface.component(c); Interface.component(w)  # catch misspelt names
    c == w && input_error(loc, "water", "must differ from the Toth component $c")
    for (k, v) in (:ns0 => ns0, :b0 => b0, :t0 => t0, :reference_temperature => reference_temperature, :cm => cm, :cG => cG, :K => K)
        v > 0 || input_error(loc, k, "must be positive, got $v")
    end
    γ >= 0 && β >= 0 || input_error(loc, "", "the water enhancements γ and β must not be negative")
    names = [String(k) for (k, _) in pairs_of(heats)]
    foreach(Interface.component, names)
    return TothGABData(c, w, ns0, b0, t0, α, χ, ΔH0, reference_temperature, γ, β, cm, cG, K,
        names, [Float64(v) for (_, v) in pairs_of(heats)], basis)
end

function build_isotherm(data::TothGABData, gas::Gas; density = nothing)
    i = findfirst(==(data.toth), gas.names)
    isnothing(i) && input_error("isotherm", "", "the gas has no $(data.toth), the component of the Toth isotherm; it has $(join(gas.names, ", "))")
    j = something(findfirst(==(data.water), gas.names), 0)
    scale = 1.0
    if data.basis == :mass
        isnothing(density) && input_error("isotherm", "basis", "capacities per kg need the adsorbent density")
        scale = density
    end
    ΔH = [(k = findfirst(==(n), data.heat_components); isnothing(k) ? 0.0 : data.heats[k]) for n in gas.names]
    return TothGAB(toth = i, water = j, ns0 = scale * data.ns0, b0 = data.b0, t0 = data.t0, α = data.α, χ = data.χ,
        ΔH0 = data.ΔH0, T0 = data.reference_temperature, γ = data.γ / scale, β = data.β / scale,
        cm = scale * data.cm, cG = data.cG, K = data.K, ΔH = ΔH)
end

# ------------------------------------------------------------------------------
# Mass transfer by component
# ------------------------------------------------------------------------------

"""
    LumpedLDFData

Rate coefficients of a [`LumpedLinearDrivingForce`](@ref) by component name.
Build it with [`lumped_ldf`](@ref).
"""
struct LumpedLDFData
    components::Vector{String}
    k::Vector{Float64}
    default::Float64
end

"""
    lumped_ldf(name => k, ...; default = 1.0)

Linear driving force mass transfer with a fixed rate coefficient `k` [1/s]
per component, such as `lumped_ldf("CO2" => 2e-4, "H2O" => 2e-3)`, for an
adsorbent whose rates were fitted to breakthrough curves (see
[`LumpedLinearDrivingForce`](@ref)). Components not named take `default`,
fast enough that a component the isotherm does not take up stays unloaded.
A rate of zero would leave its loading as it starts, and so the cyclic
steady state undetermined.
"""
function lumped_ldf(entries::Pair...; default = 1.0)
    names = [String(first(e)) for e in entries]
    foreach(component, names)
    allunique(names) || input_error("mass transfer", "", "components appear more than once: $(join(names, ", "))")
    k = [Float64(last(e)) for e in entries]
    all(>(0), k) && default > 0 || input_error("mass transfer", "", "rate coefficients must be positive")
    return LumpedLDFData(names, k, Float64(default))
end

"The mass transfer model of `data` for the components of `gas`, in its order."
build_mass_transfer(mt::AbstractMassTransfer, gas::Gas) = mt
function build_mass_transfer(data::LumpedLDFData, gas::Gas)
    k = [(i = findfirst(==(n), data.components); isnothing(i) ? data.default : data.k[i]) for n in gas.names]
    return LumpedLinearDrivingForce(SVector{length(k)}(k))
end

# ------------------------------------------------------------------------------
# Adsorbents
# ------------------------------------------------------------------------------

"""
    Adsorbent(name = ""; isotherm, density, heat_capacity, particle_diameter,
        particle_porosity, tortuosity = 3.0, macropore_radius = nothing,
        adsorbed_heat_capacity = nothing, mass_transfer = nothing, reference = "")

An adsorbent material:

- `isotherm`: isotherm data by component, such as [`dual_site_langmuir`](@ref),
  or an [`AbstractIsotherm`](@ref) in the order of the gas it is used with;
- `density` [kg/m³] and `heat_capacity` [J/(kg·K)] of the adsorbent;
- `particle_diameter` [m], `particle_porosity` and `tortuosity`, which set the
  macropore diffusion of the default linear driving force;
- `macropore_radius` [m]: if given, Knudsen diffusion in the macropores adds
  to the molecular diffusion of the default linear driving force;
- `adsorbed_heat_capacity` [J/(kg·K)] of the adsorbed phase, constant, per
  component (a vector in gas order or by name; components not named take the
  gas heat capacity at 25 °C); by default the heat capacities of the gas,
  following its thermodynamics;
- `mass_transfer`: an [`AbstractMassTransfer`](@ref) in gas order, or rates
  by component from [`lumped_ldf`](@ref); by default a
  [`LinearDrivingForce`](@ref) from the particle properties, the gas
  diffusivity and, with a macropore radius, the gas molar masses;
- `reference`: where the data come from.

Load one from the library with [`load_adsorbent`](@ref).
"""
struct Adsorbent
    name::String
    isotherm::Any
    density::Float64
    heat_capacity::Float64
    particle_diameter::Float64
    particle_porosity::Float64
    tortuosity::Float64
    macropore_radius::Union{Nothing, Float64}
    adsorbed_heat_capacity::Any
    mass_transfer::Any
    reference::String
end

function Adsorbent(name = ""; isotherm, density, heat_capacity, particle_diameter, particle_porosity,
        tortuosity = 3.0, macropore_radius = nothing, adsorbed_heat_capacity = nothing, mass_transfer = nothing, reference = "")
    loc = isempty(name) ? "adsorbent" : "adsorbent $name"
    isotherm isa Union{AbstractIsothermData, AbstractIsotherm} || input_error(loc, "isotherm", "expected isotherm data such as dual_site_langmuir(...), got $(typeof(isotherm))")
    for (k, v) in (:density => density, :heat_capacity => heat_capacity, :particle_diameter => particle_diameter, :tortuosity => tortuosity)
        v > 0 || input_error(loc, k, "must be positive, got $v")
    end
    0 < particle_porosity < 1 || input_error(loc, "particle_porosity", "must be between 0 and 1, got $particle_porosity")
    isnothing(macropore_radius) || macropore_radius > 0 || input_error(loc, "macropore_radius", "must be positive, got $macropore_radius")
    isnothing(macropore_radius) || isnothing(mass_transfer) ||
        input_error(loc, "macropore_radius", "sets the default mass transfer, so it cannot be given with mass_transfer")
    isnothing(mass_transfer) || mass_transfer isa Union{AbstractMassTransfer, LumpedLDFData} || input_error(loc, "mass_transfer", "expected an AbstractMassTransfer or lumped_ldf(...), got $(typeof(mass_transfer))")
    return Adsorbent(String(name), isotherm, density, heat_capacity, particle_diameter, particle_porosity, tortuosity,
        isnothing(macropore_radius) ? nothing : Float64(macropore_radius), adsorbed_heat_capacity, mass_transfer, String(reference))
end

"Isotherm of `adsorbent` for the components of `gas`."
isotherm(adsorbent::Adsorbent, gas::Gas) =
    _at(() -> build_isotherm(adsorbent.isotherm, gas; density = adsorbent.density), isempty(adsorbent.name) ? "adsorbent" : "adsorbent $(adsorbent.name)")

function mass_transfer(adsorbent::Adsorbent, gas::Gas, feed)
    isnothing(adsorbent.mass_transfer) || return build_mass_transfer(adsorbent.mass_transfer, gas)
    D = gas.diffusivity isa Number ? gas.diffusivity : gas_diffusivities(gas, feed.composition, feed.temperature, feed.pressure)
    D_m = D isa Number ? D : SVector{length(D)}(D)
    isnothing(adsorbent.macropore_radius) && return LinearDrivingForce(D_m, adsorbent.tortuosity, adsorbent.particle_porosity, adsorbent.particle_diameter)
    return LinearDrivingForce(D_m, adsorbent.tortuosity, adsorbent.particle_porosity, adsorbent.particle_diameter;
        macropore_radius = adsorbent.macropore_radius, molar_masses = gas.molar_masses)
end
# Without a feed, for a gas whose diffusivities do not depend on one
mass_transfer(adsorbent::Adsorbent, gas::Gas) = mass_transfer(adsorbent, gas, _no_feed(gas, "mass_transfer"))

"Thermodynamics of the phase adsorbed on `adsorbent` from `gas`, for its heat capacities."
function adsorbed_thermodynamics(adsorbent::Adsorbent, gas::Gas)
    isnothing(adsorbent.adsorbed_heat_capacity) && return gas.thermodynamics
    C = _by_component(gas.names, adsorbent.adsorbed_heat_capacity, gas.heat_capacities, "adsorbed_heat_capacity")
    return IdealGasConstantCp(; molar_masses = gas.molar_masses, heat_capacities = C)
end

"""
    equilibrium_loading(adsorbent, gas; pressure, temperature, composition) → Vector

Equilibrium loading [mol/kg] of each component of `gas` on `adsorbent`, at a
total `pressure` [Pa], `temperature` [K] and gas `composition` (see
[`composition`](@ref)). Use it to check isotherm data and units.
"""
function equilibrium_loading(adsorbent::Adsorbent, gas::Gas; pressure, temperature, composition)
    y = Interface.composition(gas, composition)
    c = y .* molar_concentration(gas.thermodynamics, pressure, temperature)
    q = compute_equilibrium(isotherm(adsorbent, gas), c, temperature)
    return collect(q) ./ adsorbent.density
end

"""
    isotherm_table(adsorbent, gas, name; pressures, temperatures = [298.15]) → Matrix

Pure-component loading [mol/kg] of component `name` on `adsorbent`, one row
per partial pressure [Pa] in `pressures` and one column per temperature [K]
in `temperatures`. The other components are at trace levels.
"""
function isotherm_table(adsorbent::Adsorbent, gas::Gas, name; pressures, temperatures = [298.15])
    i = component_index(gas, name)
    iso = isotherm(adsorbent, gas)
    N = number_of_components(gas)
    return [begin
        c = zeros(N)
        c[i] = molar_concentration(gas.thermodynamics, P, T)
        compute_equilibrium(iso, c, T)[i] / adsorbent.density
    end for P in pressures, T in temperatures]
end

function Base.show(io::IO, ::MIME"text/plain", a::Adsorbent)
    println(io, "Adsorbent ", isempty(a.name) ? "(unnamed)" : a.name)
    a.isotherm isa DualSiteLangmuirData && println(io, "  dual-site Langmuir isotherm for ", join(a.isotherm.components, ", "))
    a.isotherm isa VirialData && println(io, "  Virial isotherm for ", join(a.isotherm.components, ", "))
    a.isotherm isa TothGABData && println(io, "  Toth isotherm for ", a.isotherm.toth, " enhanced by ", a.isotherm.water, " on a GAB isotherm")
    println(io, "  density ", a.density, " kg/m³, heat capacity ", a.heat_capacity, " J/(kg·K)")
    print(io, "  particles ", 1e3 * a.particle_diameter, " mm, porosity ", a.particle_porosity, ", tortuosity ", a.tortuosity)
    isnothing(a.macropore_radius) || print(io, ", macropores ", round(1e6 * a.macropore_radius, sigdigits = 4), " µm in radius")
    isempty(a.reference) || print(io, "\n  ", a.reference)
end

# ------------------------------------------------------------------------------
# Columns
# ------------------------------------------------------------------------------

"""
    Wall(; density = 7800.0, heat_capacity = 502.0, conductivity = 16.0)

Column wall material, by default stainless steel: density [kg/m³], heat
capacity [J/(kg·K)] and conductivity [W/(m·K)].
"""
struct Wall
    density::Float64
    heat_capacity::Float64
    conductivity::Float64
end
Wall(; density = 7800.0, heat_capacity = 502.0, conductivity = 16.0) = Wall(density, heat_capacity, conductivity)

"""
    Column(; length, radius, porosity, outer_radius = radius, axial_conductivity = 0.09,
        wall = Wall(), inner_heat_transfer = 8.6, outer_heat_transfer = 2.5,
        ambient_temperature = 298.15)

The geometry and heat transfer of a packed column:

- `length` [m], inner `radius` [m] and bed `porosity` (void fraction between
  particles);
- `outer_radius` [m] of the wall, needed when the bed has a wall (the default
  thermal model);
- `axial_conductivity` [W/(m·K)] of the packed bed;
- `wall`, a [`Wall`](@ref);
- `inner_heat_transfer` and `outer_heat_transfer` [W/(m²·K)], bed to wall and
  wall to surroundings, and `ambient_temperature` [K].

The heat transfer defaults are those of Haghpanah et al. (2013).
"""
struct Column
    length::Float64
    radius::Float64
    outer_radius::Float64
    porosity::Float64
    axial_conductivity::Float64
    wall::Wall
    inner_heat_transfer::Float64
    outer_heat_transfer::Float64
    ambient_temperature::Float64
end

function Column(; length = nothing, radius = nothing, porosity = nothing, outer_radius = radius, axial_conductivity = 0.09, wall = Wall(),
        inner_heat_transfer = 8.6, outer_heat_transfer = 2.5, ambient_temperature = 298.15)
    for (k, v) in (:length => length, :radius => radius, :porosity => porosity)
        isnothing(v) && input_error("column", k, "is required")
    end
    length > 0 || input_error("column", "length", "must be positive, got $length")
    radius > 0 || input_error("column", "radius", "must be positive, got $radius")
    outer_radius >= radius || input_error("column", "outer_radius", "must be at least the inner radius $radius, got $outer_radius")
    0 < porosity < 1 || input_error("column", "porosity", "must be between 0 and 1, got $porosity")
    return Column(length, radius, outer_radius, porosity, axial_conductivity, wall, inner_heat_transfer, outer_heat_transfer, ambient_temperature)
end

function Base.show(io::IO, ::MIME"text/plain", c::Column)
    print(io, "Column ", c.length, " m long, radius ", c.radius, " m (wall to ", c.outer_radius, " m), bed porosity ", c.porosity)
end
