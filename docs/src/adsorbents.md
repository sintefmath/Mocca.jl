```@meta
CurrentModule = Mocca
```

# Adsorbents

Mocca's library holds the adsorbents of published studies. Load one by name
with [`load_adsorbent`](@ref), and list them with `library_adsorbents()`:

```julia
using Mocca
library_adsorbents()
adsorbent = load_adsorbent("zeolite_13x_haghpanah_2013")
```

The zeolites and the silica gel have dual-site Langmuir isotherms (see
[Model equations](@ref)) and the particle properties that set the mass
transfer rate. The figures show the
pure-component isotherms in mol/kg at three temperatures, computed with
[`isotherm_table`](@ref); the tables give the parameters as they are in the
file, with capacities in mol per m³ of adsorbent particles.

```@setup adsorbents
using Mocca, CairoMakie
function isotherm_figure(ads, gas)
    names = ads.isotherm.components
    T = [298.15, 323.15, 348.15]
    fig = Figure(size = (260 * length(names) + 90, 280))
    for (j, name) in enumerate(names)
        pmax = name == "H2O" ? 4000.0 : 1e5
        p = range(0, pmax, length = 101)
        ax = Axis(fig[1, j], title = name, xlabel = name == "H2O" ? "partial pressure [kPa]" : "partial pressure [bar]",
            ylabel = j == 1 ? "loading [mol/kg]" : "")
        q = Mocca.isotherm_table(ads, gas, name; pressures = p, temperatures = T)
        for k in eachindex(T)
            lines!(ax, p ./ (name == "H2O" ? 1e3 : 1e5), q[:, k], label = "$(round(Int, T[k] - 273.15)) °C")
        end
        j == 1 && Legend(fig[1, length(names) + 1], ax, framevisible = false)
    end
    return fig
end
```

## Zeolite 13X (Haghpanah et al. 2013)

The adsorbent of the four-step and two-bed VSA processes, for dry flue gas. Library name `"zeolite_13x_haghpanah_2013"`.

```@example adsorbents
ads = Mocca.load_adsorbent("zeolite_13x_haghpanah_2013") # hide
isotherm_figure(ads, Mocca.haghpanah_2013_gas()) # hide
```

```@eval
using Mocca, Markdown
ads = Mocca.load_adsorbent("zeolite_13x_haghpanah_2013")
iso = ads.isotherm
num(x) = x == 0 ? "–" : (y = round(x, sigdigits = 4); abs(y) >= 100 && isinteger(y) ? string(Int(y)) : string(y))
rows = ["| | ``q_{sb}`` [mol/m³] | ``b_0`` [m³/mol] | ``\\Delta U_b`` [kJ/mol] | ``q_{sd}`` [mol/m³] | ``d_0`` [m³/mol] | ``\\Delta U_d`` [kJ/mol] |",
    "|---|---|---|---|---|---|---|"]
for (i, c) in enumerate(iso.components)
    push!(rows, "| $c | $(num(iso.qsb[i])) | $(num(iso.b0[i])) | $(num(iso.ΔUb[i] / 1e3)) | $(num(iso.qsd[i])) | $(num(iso.d0[i])) | $(num(iso.ΔUd[i] / 1e3)) |")
end
push!(rows, "", "Density $(num(ads.density)) kg/m³, particle diameter $(num(1e3ads.particle_diameter)) mm, particle porosity $(num(ads.particle_porosity)), tortuosity $(num(ads.tortuosity)), heat capacity $(num(ads.heat_capacity)) J/(kg·K).")
Markdown.parse(join(rows, "\n"))
```

## Zeolite 13X (Krishnamurthy et al. 2014)

The adsorbent of the wet flue gas processes. Its CO₂ and N₂ parameters are close to those above, and it adds water. Water is held very strongly: at 25 °C the zeolite is nearly saturated with it below 0.5 kPa. Library name `"zeolite_13x_krishnamurthy_2014"`.

```@example adsorbents
ads = Mocca.load_adsorbent("zeolite_13x_krishnamurthy_2014") # hide
isotherm_figure(ads, Mocca.krishnamurthy_2014_gas()) # hide
```

```@eval
using Mocca, Markdown
ads = Mocca.load_adsorbent("zeolite_13x_krishnamurthy_2014")
iso = ads.isotherm
num(x) = x == 0 ? "–" : (y = round(x, sigdigits = 4); abs(y) >= 100 && isinteger(y) ? string(Int(y)) : string(y))
rows = ["| | ``q_{sb}`` [mol/m³] | ``b_0`` [m³/mol] | ``\\Delta U_b`` [kJ/mol] | ``q_{sd}`` [mol/m³] | ``d_0`` [m³/mol] | ``\\Delta U_d`` [kJ/mol] |",
    "|---|---|---|---|---|---|---|"]
for (i, c) in enumerate(iso.components)
    push!(rows, "| $c | $(num(iso.qsb[i])) | $(num(iso.b0[i])) | $(num(iso.ΔUb[i] / 1e3)) | $(num(iso.qsd[i])) | $(num(iso.d0[i])) | $(num(iso.ΔUd[i] / 1e3)) |")
end
push!(rows, "", "Density $(num(ads.density)) kg/m³, particle diameter $(num(1e3ads.particle_diameter)) mm, particle porosity $(num(ads.particle_porosity)), tortuosity $(num(ads.tortuosity)), heat capacity $(num(ads.heat_capacity)) J/(kg·K).")
Markdown.parse(join(rows, "\n"))
```

## Silica gel (Krishnamurthy et al. 2014)

The water-removal bed of the dual-adsorbent process. It holds about a fifth as much CO₂ as zeolite 13X at 1 bar, and takes up water more gradually. Library name `"silica_gel_krishnamurthy_2014"`.

```@example adsorbents
ads = Mocca.load_adsorbent("silica_gel_krishnamurthy_2014") # hide
isotherm_figure(ads, Mocca.krishnamurthy_2014_gas()) # hide
```

```@eval
using Mocca, Markdown
ads = Mocca.load_adsorbent("silica_gel_krishnamurthy_2014")
iso = ads.isotherm
num(x) = x == 0 ? "–" : (y = round(x, sigdigits = 4); abs(y) >= 100 && isinteger(y) ? string(Int(y)) : string(y))
rows = ["| | ``q_{sb}`` [mol/m³] | ``b_0`` [m³/mol] | ``\\Delta U_b`` [kJ/mol] | ``q_{sd}`` [mol/m³] | ``d_0`` [m³/mol] | ``\\Delta U_d`` [kJ/mol] |",
    "|---|---|---|---|---|---|---|"]
for (i, c) in enumerate(iso.components)
    push!(rows, "| $c | $(num(iso.qsb[i])) | $(num(iso.b0[i])) | $(num(iso.ΔUb[i] / 1e3)) | $(num(iso.qsd[i])) | $(num(iso.d0[i])) | $(num(iso.ΔUd[i] / 1e3)) |")
end
push!(rows, "", "Density $(num(ads.density)) kg/m³, particle diameter $(num(1e3ads.particle_diameter)) mm, particle porosity $(num(ads.particle_porosity)), tortuosity $(num(ads.tortuosity)), heat capacity $(num(ads.heat_capacity)) J/(kg·K).")
Markdown.parse(join(rows, "\n"))
```

## APDES-NFC (Stampi-Bombelli et al. 2020)

An amine-functionalised nanofibrillated cellulose for direct air capture,
the adsorbent of the S-TVSA process. CO₂ is held on a Toth isotherm whose
affinity and capacity grow with the water loading, and water on a GAB
isotherm in the relative humidity (see [Model equations](@ref)); the mass
transfer coefficients were fitted to breakthrough curves. The figure shows the
pure-component isotherms; at the 0.4 mbar of CO₂ in air the sorbent is close to
its capacity at 20 °C, and with 50% humidity it holds about 4% more. Library
name `"apdes_nfc_stampi_bombelli_2020"`.

```@example adsorbents
ads = Mocca.load_adsorbent("apdes_nfc_stampi_bombelli_2020") # hide
gas = Mocca.stampi_bombelli_2020_gas() # hide
T = [293.15, 323.15, 368.15] # hide
fig = Figure(size = (700, 280)) # hide
for (j, (name, pmax)) in enumerate((("CO2", 2e3), ("H2O", 2e4))) # hide
    p = range(0, pmax, length = 201) # hide
    ax = Axis(fig[1, j], title = name, xlabel = "partial pressure [kPa]", ylabel = j == 1 ? "loading [mol/kg]" : "") # hide
    q = Mocca.isotherm_table(ads, gas, name; pressures = p, temperatures = T) # hide
    for k in eachindex(T) # hide
        ok = name == "CO2" ? trues(length(p)) : p .< 0.9 .* Mocca.water_vapour_pressure(T[k]) # hide
        lines!(ax, p[ok] ./ 1e3, q[ok, k], label = "$(round(Int, T[k] - 273.15)) °C") # hide
    end # hide
    j == 1 && Legend(fig[1, 3], ax, framevisible = false) # hide
end # hide
fig # hide
```

| Toth, CO₂ | | GAB, H₂O | |
|---|---|---|---|
| ``n_{s0}`` | 2.38 mol/kg | ``c_m`` | 36.48 mol/kg |
| ``b_0`` | 707.4 1/Pa (70.74×10⁴ 1/kPa) | ``c_G`` | 0.1489 |
| ``t_0``, ``\alpha`` | 0.4148, −1.606 | ``K`` | 0.5751 |
| ``\Delta H_0``, ``\chi`` | −57.047 kJ/mol, 0 | ``\Delta H`` | −49 kJ/mol |
| ``\gamma``, ``\beta`` | 0.0061, 28.907 kg/mol | ``k`` | 2×10⁻³ 1/s |
| ``T_0`` | 296 K | | |
| ``\Delta H`` | −57 kJ/mol | | |
| ``k`` | 2×10⁻⁴ 1/s | | |

Density 1590 kg/m³ (the material; the bed holds 55.4 kg/m³ at a porosity of
0.965), particle diameter 7.5 mm, heat capacity 2070 J/(kg·K).

## Your own adsorbent

Write the isotherm with [`dual_site_langmuir`](@ref), [`langmuir`](@ref) or
[`toth_gab`](@ref),
make an [`Adsorbent`](@ref) and save it with [`save_adsorbent`](@ref) to use it
again; the [Building a process from its parts](@ref) example does this, and
[Process files](@ref) shows the file format. Check it against the source with
[`equilibrium_loading`](@ref) and [`isotherm_table`](@ref), which report
mol/kg. Other isotherm models can be added with [`register_isotherm!`](@ref).
