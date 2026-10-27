```@meta
CurrentModule = Mocca
```

# Gas properties

The gas of a bed is an ideal gas: its density and concentration follow
``P/(RT)`` in each cell. Its heat capacity, viscosity and diffusivities come
from two models held by the [`Gas`](@ref), and both are constant by default:

| Property | Default | Changing with the state of each cell |
|---|---|---|
| Heat capacity | Constant per component, as given (`heat_capacities`) | `thermodynamics = :polynomial`: cubic in temperature; `:shomate`: the Shomate equation |
| Viscosity | The value at the feed (`viscosity`) | `transport = :chapman_enskog`: kinetic theory with Wilke's mixing rule |
| Molecular diffusivities | The values at the feed (`diffusivity`) | `transport = :chapman_enskog`: kinetic theory with Wilke's mixture rule |
| Axial dispersion | Fixed, ``0.7 D_m + 0.5 v_0 d_p`` | With a transport model: Wakao and Funazkri, per component and face |
| Axial thermal conductivity | The `Column`'s `axial_conductivity` | Always constant |

The equations are in [Model equations](@ref): see [Energy balance of the
column](@ref), [Axial dispersion](@ref) and [Transport properties](@ref).

## Choosing the models

```julia
using Mocca
gas = Gas(["CO2", "N2", "H2O"]; thermodynamics = :polynomial, transport = :chapman_enskog)
```

The two choices are independent. A transport model does not replace the
viscosity and diffusivities you give: they are the values at the feed, and
each cell scales them by the ratio of the model's value at its state to that at
the feed. Give none and they are the model's own values at the feed. So a
preset keeps the properties of its paper at the feed, and a design variable or
history match on the viscosity or the diffusivities keeps its meaning:

```julia
gas = Mocca.krishnamurthy_2014_gas(transport = :chapman_enskog)
process = Mocca.lpp_vsa_process(gas = gas)
```

In a process file, set `"thermodynamics": "polynomial"` and
`"transport": "chapman_enskog"` in the `gas` (see [Process files](@ref)).

The adsorbed phase takes the heat capacities of the gas, weighted by the
loadings. The heat capacity ratio of pumps and compressors is that of the feed,
unless the process sets `heat_capacity_ratio` (the presets set 1.4, as their
papers do).

## Component data

Mocca knows these components by name (see [`component`](@ref)). Heat
capacities and viscosities at 25 °C are the constant defaults. The cubic heat
capacities, ``\tilde c_p = A + BT + CT^2 + DT^3`` in J/(mol·K) with ``T`` in K,
are from Appendix A of Reid, Prausnitz and Poling (1987); the Lennard-Jones
parameters are from Table B.1 of Poling, Prausnitz and O'Connell (2001). Both
are the sets of the MATLAB thermodynamics library Mocca's rules follow. The
Shomate coefficients are the lowest-temperature sets of the NIST Chemistry
WebBook.

```@eval
using Mocca, Markdown
num(x) = (y = round(x, sigdigits = 4); abs(y) >= 100 && isinteger(y) ? string(Int(y)) : string(y))
opt(x, f) = isnothing(x) ? "–" : f(x)
rows = ["| | ``M`` [g/mol] | ``c_p`` [J/(kg·K)] | ``\\mu`` [µPa·s] | ``A`` | ``B`` | ``C`` | ``D`` | ``\\sigma`` [Å] | ``\\varepsilon/k`` [K] |",
    "|---|---|---|---|---|---|---|---|---|---|"]
for name in ["CO2", "N2", "O2", "Ar", "H2O", "CH4", "H2", "CO", "He"]
    c = Mocca.component(name)
    p = isnothing(c.cp_polynomial) ? fill("–", 4) : [x == 0 ? "0" : num(x) for x in c.cp_polynomial]
    lj = isnothing(c.lennard_jones) ? ["–", "–"] : [num(x) for x in c.lennard_jones]
    push!(rows, "| $name | $(num(1e3c.molar_mass)) | $(num(c.heat_capacity)) | $(num(1e6c.viscosity)) | $(join(p, " | ")) | $(join(lj, " | ")) |")
end
Markdown.parse(join(rows, "\n"))
```

CO and He have no cubic heat capacities or Lennard-Jones parameters, so they
can take only constant or Shomate heat capacities and constant transport. Add
the data for a component, or a new component, with [`register_component!`](@ref):

```julia
register_component!(Component("CO"; formula = "CO", molar_mass = 28.010e-3, heat_capacity = 1040.0,
    viscosity = 1.77e-5, cp_polynomial = (A, B, C, D), lennard_jones = (σ, ε_k)))
```

The figure shows the three heat capacity models, and the Chapman–Enskog
viscosities and binary diffusivities with N₂ at 1 atm (for N₂ its
self-diffusivity), from 280 to 500 K.

```@example gas_properties
using Mocca, CairoMakie # hide
names = ["CO2", "N2", "H2O"] # hide
T = range(280.0, 500.0, length = 45) # hide
models = [(Gas(names).thermodynamics, "constant", :dot), (Gas(names; thermodynamics = :polynomial).thermodynamics, "cubic", :solid), # hide
    (Gas(names; thermodynamics = :shomate).thermodynamics, "Shomate", :dash)] # hide
ce = Gas(names; transport = :chapman_enskog).transport # hide
colors = Makie.wong_colors() # hide
fig = Figure(size = (900, 300)) # hide
ax1 = Axis(fig[1, 1], xlabel = "T [K]", ylabel = "heat capacity [J/(mol·K)]") # hide
ax2 = Axis(fig[1, 2], xlabel = "T [K]", ylabel = "viscosity [µPa·s]") # hide
ax3 = Axis(fig[1, 3], xlabel = "T [K]", ylabel = "diffusivity with N₂ [cm²/s]") # hide
for (i, name) in enumerate(names) # hide
    for (m, _, style) in models # hide
        lines!(ax1, T, [Mocca.component_heat_capacities(m, t)[i] for t in T], color = colors[i], linestyle = style) # hide
    end # hide
    lines!(ax2, T, [1e6 * Mocca.component_viscosities(ce, t)[i] for t in T], color = colors[i], label = name) # hide
    lines!(ax3, T, [1e4 * Mocca.binary_diffusivity(ce, i, 2, t, 101325.0) for t in T], color = colors[i]) # hide
end # hide
Legend(fig[1, 4], [[LineElement(color = colors[i]) for i in eachindex(names)], [LineElement(color = :black, linestyle = s) for (_, _, s) in models]], # hide
    [names, [l for (_, l, _) in models]], ["component", "heat capacity"], framevisible = false) # hide
fig # hide
```

## Accuracy

Against the NIST Chemistry WebBook, from 300 to 500 K: the largest deviation of
the cubic heat capacities, and the deviation of the Chapman–Enskog viscosity
from that of the gas at 0.01 bar.

| | Cubic ``c_p`` | ``\mu`` at 300 K | 400 K | 500 K |
|---|---|---|---|---|
| CO2 | 0.7% | +1.3% | −0.3% | −1.9% |
| N2 | 0.7% | −1.0% | −2.6% | −3.6% |
| O2 | 1.3% | −0.3% | −1.1% | −1.6% |
| H2O | 0.6% | +10.1% | +6.0% | +2.1% |
| CH4 | 2.0% | +0.6% | −0.2% | −0.7% |
| H2 | 0.9% | −0.6% | −1.0% | −1.8% |
| Ar | 0.0% | +0.4% | −1.8% | −3.3% |

Water is polar, and a Lennard-Jones potential describes it poorly: its
viscosity is 10% high at 300 K. Since each cell takes the ratio to the feed,
that offset cancels, but the temperature trend of water's viscosity is too
steep. The Shomate coefficients match the NIST sets exactly. That of water is
fitted from 500 K, but is within 0.05% of the NIST heat capacity of the vapour
from 300 to 400 K. Those of N₂ and O₂ are valid only to 500 and 700 K.

At 25 °C and 1 atm the binary diffusivities are 1.55×10⁻⁵ m²/s for CO₂–N₂
(Mocca's constant default is 1.6×10⁻⁵), 2.14×10⁻⁵ for N₂–H₂O and 1.48×10⁻⁵
for CO₂–H₂O.

## Limits

- **The diffusivity of the linear driving force** is taken at the feed
  pressure; only its temperature and composition change. Molecular diffusion
  goes as ``1/P``, but at the 0.03 bar of an evacuation that makes mass
  transfer about 30 times faster than at 1 bar, and the cyclic steady state of
  the wet flue gas processes does not converge. In the macropores of a pellet,
  Knudsen diffusion would limit that rise; Mocca does not model it. The axial
  dispersion takes the diffusivities at the pressure of each cell.
- **Wilke's mixture diffusivity** has a floor of 10⁻⁶ on the mole fractions of
  the other components (see [Transport properties](@ref)). Without it, near a
  pure component in a gas of three or more, Newton's method fails. The floor
  changes no weight by more than 1% for mole fractions above 10⁻⁴, and a binary
  gas gets its binary diffusivity exactly.
- **Wakao and Funazkri's molecular term**, ``20 D_m/\varepsilon``, is large:
  at 1 bar it is comparable to the velocity term, and at evacuation pressures
  it dominates the dispersion. Scale it with the bed's `DispersionScale`
  parameter, or give `dispersion` to [`add_bed!`](@ref) to fix the dispersion.
- **Thermal conductivity** stays as given: the bed's axial conductivity is one
  effective value for gas and solid, and the wall's conductivity and the heat
  transfer coefficients are constant.

## Effect on the presets

The recovery and the energy of the presets at cyclic steady state (30 cells),
with the properties of their papers at the feed. "Chapman–Enskog" is
`transport = :chapman_enskog`, with Wakao–Funazkri dispersion; "+ cubic" adds
the cubic heat capacities in place of the papers' constant ones.

| Preset | Recovery: constant | Chapman–Enskog | + cubic | Energy [kWh/t]: constant | Chapman–Enskog | + cubic |
|---|---|---|---|---|---|---|
| Four-step VSA, Haghpanah et al. (2013) | 0.367 | 0.350 | 0.349 | 141.0 | 144.7 | 144.7 |
| Two-bed VSA, Haghpanah et al. (2013) | 0.612 | 0.599 | 0.599 | 92.3 | 92.7 | 92.6 |
| LPP VSA, Krishnamurthy et al. (2014) | 0.788 | 0.708 | 0.714 | 169.2 | 182.3 | 181.3 |
| Dual-adsorbent VSA, Krishnamurthy et al. (2014) | 0.655 | 0.577 | 0.581 | 196.0 | 215.7 | 214.6 |

The heat capacities change little. Nearly all of the drop of the wet flue gas
processes is the mass transfer of water: its diffusivity in the CO₂-rich gas
of an evacuation is 0.75 of that in the feed, so water is taken up more
slowly, its front spreads further into the zeolite and takes CO₂ capacity.
With only the viscosity following the state, the LPP process recovers 0.793.

## References

- Neufeld, P. D., Janzen, A. R. and Aziz, R. A. (1972). Empirical equations to
  calculate 16 of the transport collision integrals Ω(l, s)* for the
  Lennard-Jones (12–6) potential. *J. Chem. Phys.* 57, 1100–1102.
- Poling, B. E., Prausnitz, J. M. and O'Connell, J. P. (2001). *The Properties
  of Gases and Liquids*, 5th ed. McGraw-Hill.
- Reid, R. C., Prausnitz, J. M. and Poling, B. E. (1987). *The Properties of
  Gases and Liquids*, 4th ed. McGraw-Hill.
- Wakao, N. and Funazkri, T. (1978). Effect of fluid dispersion coefficients
  on particle-to-fluid mass transfer coefficients in packed beds. *Chem. Eng.
  Sci.* 33, 1375–1384.
- Wilke, C. R. (1950). A viscosity equation for gas mixtures. *J. Chem. Phys.*
  18, 517–519.
- Wilke, C. R. (1950). Diffusional properties of multicomponent gases. *Chem.
  Eng. Prog.* 46, 95–104.
- Linstrom, P. J. and Mallard, W. G. (eds.). *NIST Chemistry WebBook*, NIST
  Standard Reference Database 69. https://doi.org/10.18434/T4D303
