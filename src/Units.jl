"""
    Mocca.Units

Mocca works in plain SI units throughout: Pa, K, m, s, mol, kg and J. These
constants convert other units to SI by multiplication, so `0.1bar` is `1e4`
(Pa) and `-36.6kJ_per_mol` is `-36600.0` (J/mol):

```julia
using Mocca.Units
Evacuate(40.0; to = 0.1bar)
```

[`Units.celsius`](@ref) converts a temperature in °C to K. Files read with
[`read_process`](@ref Mocca.Interface.read_process) may give a number with one of these units as a string,
such as `"0.1 bar"`.
"""
module Units

export Pa, kPa, MPa, bar, mbar, atm
export kJ, MJ, kWh, J_per_mol, kJ_per_mol
export mmol, kmol, tonne
export mm, cm, μm, minute, hour, day
export celsius

const Pa = 1.0
const kPa = 1e3
const MPa = 1e6
const bar = 1e5
const mbar = 1e2
const atm = 101325.0

const J = 1.0
const kJ = 1e3
const MJ = 1e6
const kWh = 3.6e6
const J_per_mol = 1.0
const kJ_per_mol = 1e3

const mol = 1.0
const mmol = 1e-3
const kmol = 1e3

const g = 1e-3
const kg = 1.0
const tonne = 1e3

const m = 1.0
const cm = 1e-2
const mm = 1e-3
const μm = 1e-6

const s = 1.0
const minute = 60.0
const hour = 3600.0
const day = 86400.0

"Temperature in K of `T` °C."
celsius(T) = T + 273.15

# Units that may follow a number in an input file, with their factor to SI.
# Temperatures in °C are handled separately, since they have an offset.
const FILE_UNITS = Dict{String, Float64}(
    "Pa" => Pa, "kPa" => kPa, "MPa" => MPa, "bar" => bar, "mbar" => mbar, "atm" => atm,
    "J" => J, "kJ" => kJ, "MJ" => MJ, "kWh" => kWh,
    "J/mol" => J_per_mol, "kJ/mol" => kJ_per_mol,
    "mol" => mol, "mmol" => mmol, "kmol" => kmol,
    "g" => g, "kg" => kg, "t" => tonne,
    "m" => m, "cm" => cm, "mm" => mm, "um" => μm, "μm" => μm,
    "s" => s, "min" => minute, "h" => hour, "d" => day,
    "m/s" => 1.0, "cm/s" => cm, "m3/s" => 1.0, "m^3/s" => 1.0,
    "K" => 1.0, "kg/m3" => 1.0, "kg/m^3" => 1.0, "J/kg/K" => 1.0, "W/m/K" => 1.0, "W/m2/K" => 1.0,
    "mol/m3" => 1.0, "mol/m^3" => 1.0, "m3/mol" => 1.0, "m^3/mol" => 1.0, "m2/s" => 1.0, "m^2/s" => 1.0,
    "Pa.s" => 1.0, "Pa s" => 1.0, "1/s" => 1.0,
)

"""
    parse_quantity(s::AbstractString) → Float64

The SI value of a number followed by a unit, such as `"0.1 bar"` or `"25 °C"`.
A bare number is returned as it is.
"""
function parse_quantity(str::AbstractString)
    m = match(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(.*?)\s*$", str)
    isnothing(m) && throw(ArgumentError("Expected a number with an optional unit, got \"$str\""))
    x = parse(Float64, m[1])
    u = m[2]
    isempty(u) && return x
    u in ("°C", "degC", "C") && return celsius(x)
    haskey(FILE_UNITS, u) || throw(ArgumentError("Unknown unit \"$u\" in \"$str\"; known units: $(join(sort(collect(keys(FILE_UNITS))), ", ")), °C"))
    return x * FILE_UNITS[u]
end

end
