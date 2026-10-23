```@meta
CurrentModule = Mocca
```

# Processes and validation

Mocca's library holds processes from published studies, set up at the papers'
operating conditions. Read one by name with [`read_process`](@ref), change a
setting and run it:

```julia
using Mocca
library_processes()
process = read_process("krishnamurthy_2014_lpp_vsa")
```

This page compares Mocca's results for each of them with the paper's, and
says where they differ and why, as far as it is known.

## How Mocca compares

Mocca's results against the published ones. A point on the diagonal is an
exact match; the dotted lines are ±10%.

```@example validation
using CairoMakie # hide
# Published values and Mocca's, as computed by the examples (update these when # hide
# the model changes): purity and recovery [%], pump energy [kWh/t CO₂] and # hide
# productivity [t CO₂ per m³ of 13X per day] # hide
cases = [ # hide
    (name = "Four-step VSA (Haghpanah 2013)", paper = (88.7, 35.8, NaN, NaN), mocca = (79.5, 36.7, NaN, NaN), marker = :circle), # hide
    (name = "Wet flue gas, LPP (Krishnamurthy 2014)", paper = (78.8, 88.6, 188.8, 1.73), mocca = (82.6, 79.2, 185.0, 1.56), marker = :utriangle), # hide
    (name = "Wet flue gas, dual adsorbent (Krishnamurthy 2014)", paper = (95.0, 90.0, 177.0, 1.82), mocca = (95.6, 85.6, 192.0, 1.74), marker = :diamond), # hide
] # hide
titles = ["purity [%]", "recovery [%]", "pump energy [kWh/t]", "productivity [t/(m³·day)]"] # hide
limits = [(70, 100), (30, 100), (150, 210), (1.4, 2.0)] # hide
colours = ["#5b2f73", "#c4362e", "#3d8c47"] # hide
fig = Figure(size = (960, 290)) # hide
for (j, title) in enumerate(titles) # hide
    lo, hi = limits[j] # hide
    ax = Axis(fig[1, j], title = title, xlabel = "paper", ylabel = j == 1 ? "Mocca" : "", limits = (lo, hi, lo, hi), aspect = 1, xticks = WilkinsonTicks(4), yticks = WilkinsonTicks(4)) # hide
    lines!(ax, [lo, hi], [lo, hi], color = :gray) # hide
    lines!(ax, [lo, hi], [1.1lo, 1.1hi], color = :gray, linestyle = :dot) # hide
    lines!(ax, [lo, hi], [0.9lo, 0.9hi], color = :gray, linestyle = :dot) # hide
    for (k, c) in enumerate(cases) # hide
        isnan(c.paper[j]) && continue # hide
        scatter!(ax, [c.paper[j]], [c.mocca[j]], marker = c.marker, color = colours[k], markersize = 14, label = c.name) # hide
    end # hide
end # hide
Legend(fig[2, 1:4], [MarkerElement(marker = c.marker, color = colours[k], markersize = 14) for (k, c) in enumerate(cases)], # hide
    [c.name for c in cases], orientation = :horizontal, framevisible = false) # hide
fig # hide
```

| Process | | Purity | Recovery | Pump energy | Productivity |
|---|---|---|---|---|---|
| Four-step VSA | paper | 88.7% | 35.8% | | |
| | Mocca | 79.5% | 36.7% | | |
| Wet flue gas, LPP | paper | 78.8% (95.2% dry) | 88.6% | 188.8 kWh/t | 1.73 t/(m³·day) |
| | Mocca | 82.6% (97.5% dry) | 79.2% | 185 kWh/t | 1.56 t/(m³·day) |
| Wet flue gas, dual adsorbent | paper | 95% | 90% | 177 kWh/t | 1.82 t/(m³·day) |
| | Mocca | 95.6% | 85.6% | 192 kWh/t | 1.74 t/(m³·day) |

The energy counts the vacuum pumps only, as the papers do;
[`performance`](@ref) also counts the feed compressor. Productivity is per m³
of 13X. Purities with water count the water in the CO₂ product.

What is known about the differences:

- **Four-step VSA.** Mocca's purity is about 9 points below the paper's at the
  same recovery, and the difference does not change with the grid, the time
  step or the heat model. Its cause has not been found yet.
- **Wet flue gas, LPP.** Water builds up in the bed over hundreds of cycles,
  and after the 300 cycles of the example the CO₂ balance is still 0.7% open.
  Mocca pressurises with nitrogen where the paper uses stored light product.
  Recovery is 9 points below the paper's.
- **Wet flue gas, dual adsorbent.** The CO₂ balance closes after 131 cycles,
  but the silica gel is still filling with water, and in Mocca a few hundred
  ppm of water reach the 13X bed, against 20 ppm in the paper. Recovery is 4
  points lower and the energy 8% higher.

## Four-step VSA

[Haghpanah et al. (2013)](https://doi.org/10.1021/ie302658y): dry flue gas,
one bed of zeolite 13X. Library name `"haghpanah_2013_four_stage_vsa"`.

```@raw html
<figure class="mocca-schematic mocca-process">
<div class="mocca-process-figs"><div><svg viewBox="0 0 400 290" role="img" aria-label="Flowsheet: one bed with a feed valve and a vacuum valve at the bottom and a product valve at the top"><defs><marker id="fs-ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="fs-head" d="M0,0 L10,5 L0,10 z"/></marker></defs><rect class="fs-bed" x="178" y="75" width="44" height="140" rx="6"/><text class="fs-name" x="170" y="149.0" text-anchor="end">13X</text><path class="fs-pipe" d="M200 75 V40"/><path class="fs-pipe" d="M200 215 V250"/><path class="fs-pipe" d="M200 40 H360" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M283 33 L297 47 L297 33 L283 47 Z"/><text class="fs-dev" x="290" y="28" text-anchor="middle">V_product</text><text class="fs-text" x="360" y="62" text-anchor="end">light product, N₂</text><path class="fs-pipe" d="M40 250 H196" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M103 243 L117 257 L117 243 L103 257 Z"/><text class="fs-dev" x="110" y="270" text-anchor="middle">V_feed</text><text class="fs-text" x="40" y="238" text-anchor="start">feed, 15% CO₂</text><path class="fs-pipe" d="M200 250 H360" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M283 243 L297 257 L297 243 L283 257 Z"/><text class="fs-dev" x="290" y="270" text-anchor="middle">V_vacuum</text><text class="fs-text" x="360" y="238" text-anchor="end">CO₂ product</text></svg></div><div><svg viewBox="0 0 440 90" role="img" aria-label="Cycle: Bed: pressurisation 15 s, adsorption 15 s, blowdown 30 s, evacuation 40 s"><text class="fs-dev" x="44" y="22" text-anchor="end">Bed</text><rect class="fs-seg " x="52.5" y="10" width="55.7" height="16" style="fill:var(--mc-step-pressurise)"><title>pressurisation, 15 s</title></rect><rect class="fs-seg " x="109.2" y="10" width="55.7" height="16" style="fill:var(--mc-step-adsorb)"><title>adsorption, 15 s</title></rect><rect class="fs-seg " x="165.9" y="10" width="112.4" height="16" style="fill:var(--mc-step-blowdown)"><title>blowdown, 30 s</title></rect><rect class="fs-seg " x="279.3" y="10" width="150.2" height="16" style="fill:var(--mc-step-evacuate)"><title>evacuation, 40 s</title></rect><text class="fs-dev" x="52.0" y="46" text-anchor="middle">0</text><text class="fs-dev" x="146.5" y="46" text-anchor="middle">25</text><text class="fs-dev" x="241.0" y="46" text-anchor="middle">50</text><text class="fs-dev" x="335.5" y="46" text-anchor="middle">75</text><text class="fs-dev" x="430.0" y="46" text-anchor="middle">100</text><text class="fs-dev" x="430" y="60" text-anchor="end">time in the cycle [s]</text><rect x="8" y="70" width="12" height="8" rx="1" style="fill:var(--mc-step-pressurise)"/><text class="fs-dev" x="24" y="78">pressurisation</text><rect x="125.8" y="70" width="12" height="8" rx="1" style="fill:var(--mc-step-adsorb)"/><text class="fs-dev" x="141.8" y="78">adsorption</text><rect x="216.8" y="70" width="12" height="8" rx="1" style="fill:var(--mc-step-blowdown)"/><text class="fs-dev" x="232.8" y="78">blowdown</text><rect x="294.4" y="70" width="12" height="8" rx="1" style="fill:var(--mc-step-evacuate)"/><text class="fs-dev" x="310.4" y="78">evacuation</text></svg></div></div>
<figcaption>The flowsheet and one cycle, to scale. Valve names are those of the process.</figcaption>
</figure>
```

| Setting | Value |
|---|---|
| Feed | 15% CO₂ in N₂, 25 °C, 1 bar, interstitial velocity 1 m/s |
| Column | 1 m long, 0.289 m inner diameter, with the steel wall |
| Pressures | 1 bar adsorption, 0.2 bar blowdown, 0.1 bar evacuation |
| Steps | pressurisation 15 s, adsorption 15 s, blowdown 30 s, evacuation 40 s |

The [Getting started](@ref) example runs it, and
[How a VSA cycle works](@ref) follows it through one cycle.

## Two-bed VSA

Built in Mocca from the beds of Haghpanah et al. (2013), to show a multi-bed
cycle; it has no published counterpart. Make it with
`two_bed_vsa_process()`.

```@raw html
<figure class="mocca-schematic mocca-process">
<div class="mocca-process-figs"><div><svg viewBox="0 0 440 300" role="img" aria-label="Flowsheet: two beds A and B, each with a feed valve, a vacuum valve and a product valve, and a link valve between their tops"><defs><marker id="fs-ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="fs-head" d="M0,0 L10,5 L0,10 z"/></marker></defs><rect class="fs-bed" x="128" y="75" width="44" height="140" rx="6"/><text class="fs-name" x="120" y="149.0" text-anchor="end">A</text><rect class="fs-bed" x="268" y="75" width="44" height="140" rx="6"/><text class="fs-name" x="320" y="149.0">B</text><path class="fs-pipe" d="M150 75 V40"/><path class="fs-pipe" d="M150 215 V245"/><path class="fs-pipe" d="M290 75 V40"/><path class="fs-pipe" d="M290 215 V245"/><path class="fs-pipe" d="M150 40 H20" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M73 33 L87 47 L87 33 L73 47 Z"/><text class="fs-dev" x="80" y="28" text-anchor="middle">V_product_A</text><path class="fs-pipe" d="M290 40 H420" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M353 33 L367 47 L367 33 L353 47 Z"/><text class="fs-dev" x="360" y="28" text-anchor="middle">V_product_B</text><path class="fs-pipe" d="M150 40 H290"/><path class="fs-valve" d="M213.0 33 L227.0 47 L227.0 33 L213.0 47 Z"/><text class="fs-dev" x="220.0" y="28" text-anchor="middle">V_eq</text><path class="fs-pipe" d="M20 245 H146" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M78 238 L92 252 L92 238 L78 252 Z"/><text class="fs-dev" x="85" y="233" text-anchor="middle">V_feed_A</text><text class="fs-text" x="20" y="262" text-anchor="start">feed</text><path class="fs-pipe" d="M420 245 H294" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M348 238 L362 252 L362 238 L348 252 Z"/><text class="fs-dev" x="355" y="233" text-anchor="middle">V_feed_B</text><text class="fs-text" x="420" y="262" text-anchor="end">feed</text><path class="fs-pipe" d="M150 245 V292" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M143 263 L157 263 L143 277 L157 277 Z"/><text class="fs-dev" x="138" y="274" text-anchor="end">V_vacuum_A</text><path class="fs-pipe" d="M290 245 V292" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M283 263 L297 263 L283 277 L297 277 Z"/><text class="fs-dev" x="302" y="274">V_vacuum_B</text><text class="fs-text" x="220.0" y="286" text-anchor="middle">CO₂ product</text></svg></div><div><svg viewBox="0 0 440 112" role="img" aria-label="Cycle: A: pressurisation 15 s, adsorption 15 s, equalisation 10 s, evacuation 30 s, equalisation 10 s; B: evacuation 30 s, equalisation 10 s, pressurisation 15 s, adsorption 15 s, equalisation 10 s"><text class="fs-dev" x="44" y="22" text-anchor="end">A</text><rect class="fs-seg " x="52.5" y="10" width="69.9" height="16" style="fill:var(--mc-step-pressurise)"><title>pressurisation, 15 s</title></rect><rect class="fs-seg " x="123.4" y="10" width="69.9" height="16" style="fill:var(--mc-step-adsorb)"><title>adsorption, 15 s</title></rect><rect class="fs-seg " x="194.2" y="10" width="46.2" height="16" style="fill:var(--mc-step-equalise)"><title>equalisation, 10 s</title></rect><rect class="fs-seg " x="241.5" y="10" width="140.8" height="16" style="fill:var(--mc-step-evacuate)"><title>evacuation, 30 s</title></rect><rect class="fs-seg " x="383.2" y="10" width="46.2" height="16" style="fill:var(--mc-step-equalise)"><title>equalisation, 10 s</title></rect><text class="fs-dev" x="44" y="44" text-anchor="end">B</text><rect class="fs-seg " x="52.5" y="32" width="140.8" height="16" style="fill:var(--mc-step-evacuate)"><title>evacuation, 30 s</title></rect><rect class="fs-seg " x="194.2" y="32" width="46.2" height="16" style="fill:var(--mc-step-equalise)"><title>equalisation, 10 s</title></rect><rect class="fs-seg " x="241.5" y="32" width="69.9" height="16" style="fill:var(--mc-step-pressurise)"><title>pressurisation, 15 s</title></rect><rect class="fs-seg " x="312.4" y="32" width="69.9" height="16" style="fill:var(--mc-step-adsorb)"><title>adsorption, 15 s</title></rect><rect class="fs-seg " x="383.2" y="32" width="46.2" height="16" style="fill:var(--mc-step-equalise)"><title>equalisation, 10 s</title></rect><text class="fs-dev" x="52.0" y="68" text-anchor="middle">0</text><text class="fs-dev" x="170.1" y="68" text-anchor="middle">25</text><text class="fs-dev" x="288.2" y="68" text-anchor="middle">50</text><text class="fs-dev" x="406.4" y="68" text-anchor="middle">75</text><text class="fs-dev" x="430" y="82" text-anchor="end">time in the cycle [s]</text><rect x="8" y="92" width="12" height="8" rx="1" style="fill:var(--mc-step-pressurise)"/><text class="fs-dev" x="24" y="100">pressurisation</text><rect x="125.8" y="92" width="12" height="8" rx="1" style="fill:var(--mc-step-adsorb)"/><text class="fs-dev" x="141.8" y="100">adsorption</text><rect x="216.8" y="92" width="12" height="8" rx="1" style="fill:var(--mc-step-equalise)"/><text class="fs-dev" x="232.8" y="100">equalisation</text><rect x="321.20000000000005" y="92" width="12" height="8" rx="1" style="fill:var(--mc-step-evacuate)"/><text class="fs-dev" x="337.20000000000005" y="100">evacuation</text></svg></div></div>
<figcaption>The beds run half a cycle apart and pass gas through the link between their tops.</figcaption>
</figure>
```

Evacuating to 0.1 bar, it gives 57.4% purity, 60.6% recovery and 93.4 kWh/t,
pumps and compressor together. The
[Two-bed VSA with pressure equalisation](@ref) example builds it step by step
and shows what the equalisation is worth.

## Wet flue gas, light product pressurisation

Section 3 of [Krishnamurthy et al. (2014)](https://doi.org/10.1021/ie5024723):
wet flue gas, one bed of zeolite 13X. Library name
`"krishnamurthy_2014_lpp_vsa"`.

```@raw html
<figure class="mocca-schematic mocca-process">
<div class="mocca-process-figs"><div><svg viewBox="0 0 400 290" role="img" aria-label="Flowsheet: one bed with a feed valve and a vacuum valve at the bottom and a product valve at the top and a light product inlet at the top"><defs><marker id="fs-ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="fs-head" d="M0,0 L10,5 L0,10 z"/></marker></defs><rect class="fs-bed" x="178" y="75" width="44" height="140" rx="6"/><text class="fs-name" x="170" y="149.0" text-anchor="end">13X</text><path class="fs-pipe" d="M200 75 V40"/><path class="fs-pipe" d="M200 215 V250"/><path class="fs-pipe" d="M200 40 H360" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M283 33 L297 47 L297 33 L283 47 Z"/><text class="fs-dev" x="290" y="28" text-anchor="middle">V_product</text><text class="fs-text" x="360" y="62" text-anchor="end">light product, N₂</text><path class="fs-pipe" d="M40 40 H196" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M103 33 L117 47 L117 33 L103 47 Z"/><text class="fs-dev" x="110" y="28" text-anchor="middle">V_lpp</text><text class="fs-text" x="40" y="62" text-anchor="start">N₂ in (LPP)</text><path class="fs-pipe" d="M40 250 H196" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M103 243 L117 257 L117 243 L103 257 Z"/><text class="fs-dev" x="110" y="270" text-anchor="middle">V_feed</text><text class="fs-text" x="40" y="238" text-anchor="start">feed, 15% CO₂, 3% H₂O</text><path class="fs-pipe" d="M200 250 H360" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M283 243 L297 257 L297 243 L283 257 Z"/><text class="fs-dev" x="290" y="270" text-anchor="middle">V_vacuum</text><text class="fs-text" x="360" y="238" text-anchor="end">CO₂ product</text></svg></div><div><svg viewBox="0 0 440 90" role="img" aria-label="Cycle: Bed: LPP 20 s, adsorption 78.4 s, blowdown 37.9 s, evacuation 111.6 s"><text class="fs-dev" x="44" y="22" text-anchor="end">Bed</text><rect class="fs-seg " x="52.5" y="10" width="29.5" height="16" style="fill:var(--mc-step-lpp)"><title>LPP, 20 s</title></rect><rect class="fs-seg " x="83.0" y="10" width="118.5" height="16" style="fill:var(--mc-step-adsorb)"><title>adsorption, 78.4 s</title></rect><rect class="fs-seg " x="202.5" y="10" width="56.8" height="16" style="fill:var(--mc-step-blowdown)"><title>blowdown, 37.9 s</title></rect><rect class="fs-seg " x="260.3" y="10" width="169.2" height="16" style="fill:var(--mc-step-evacuate)"><title>evacuation, 111.6 s</title></rect><text class="fs-dev" x="52.0" y="46" text-anchor="middle">0</text><text class="fs-dev" x="128.2" y="46" text-anchor="middle">50</text><text class="fs-dev" x="204.5" y="46" text-anchor="middle">100</text><text class="fs-dev" x="280.7" y="46" text-anchor="middle">150</text><text class="fs-dev" x="357.0" y="46" text-anchor="middle">200</text><text class="fs-dev" x="430" y="60" text-anchor="end">time in the cycle [s]</text><rect x="8" y="70" width="12" height="8" rx="1" style="fill:var(--mc-step-lpp)"/><text class="fs-dev" x="24" y="78">LPP</text><rect x="52.1" y="70" width="12" height="8" rx="1" style="fill:var(--mc-step-adsorb)"/><text class="fs-dev" x="68.1" y="78">adsorption</text><rect x="143.1" y="70" width="12" height="8" rx="1" style="fill:var(--mc-step-blowdown)"/><text class="fs-dev" x="159.1" y="78">blowdown</text><rect x="220.7" y="70" width="12" height="8" rx="1" style="fill:var(--mc-step-evacuate)"/><text class="fs-dev" x="236.7" y="78">evacuation</text></svg></div></div>
<figcaption>Light product pressurisation fills the bed from the top with nitrogen before the feed starts.</figcaption>
</figure>
```

| Setting | Value |
|---|---|
| Feed | 15% CO₂, 82% N₂, 3% H₂O, 25 °C, 1.013 bar, interstitial velocity 0.45 m/s |
| Column | 1 m long, 0.3 m inner diameter |
| Pressures | 1.013 bar adsorption, 0.08 bar blowdown, 0.03 bar evacuation |
| Steps | LPP 20 s, adsorption 78.4 s, blowdown 37.9 s, evacuation 111.6 s |

The [Wet flue gas: single-bed VSA with light product pressurisation](@ref)
example runs it and compares its bed profiles with the paper's.

## Wet flue gas, dual adsorbent

Section 4 of [Krishnamurthy et al. (2014)](https://doi.org/10.1021/ie5024723):
a silica gel bed takes the water out of the wet flue gas before a 13X bed
captures the CO₂. Library name `"krishnamurthy_2014_dual_adsorbent_vsa"`.

```@raw html
<figure class="mocca-schematic mocca-process">
<div class="mocca-process-figs"><div><svg viewBox="0 0 460 300" role="img" aria-label="Flowsheet: a silica gel bed with feed and waste valves at the bottom, linked from its top to the bottom of a 13X bed, which has a vacuum valve at the bottom and product and light product valves at the top"><defs><marker id="fs-ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="fs-head" d="M0,0 L10,5 L0,10 z"/></marker></defs><rect class="fs-bed" x="108" y="150" width="44" height="62" rx="6"/><text class="fs-name" x="100" y="185.0" text-anchor="end">silica gel</text><rect class="fs-bed" x="278" y="75" width="44" height="140" rx="6"/><text class="fs-name" x="330" y="149.0">13X</text><path class="fs-pipe" d="M130 212 V245"/><path class="fs-pipe" d="M300 215 V245"/><path class="fs-pipe" d="M300 75 V40"/><path class="fs-pipe" d="M130 150 V115 H215 V245 H296" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M208 173 L222 173 L208 187 L222 187 Z"/><text class="fs-dev" x="227" y="184">V_link</text><path class="fs-pipe" d="M20 245 H126" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M63 238 L77 252 L77 238 L63 252 Z"/><text class="fs-dev" x="70" y="265" text-anchor="middle">V_feed</text><text class="fs-text" x="20" y="233" text-anchor="start">wet feed</text><path class="fs-pipe" d="M130 245 V292" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M123 263 L137 263 L123 277 L137 277 Z"/><text class="fs-dev" x="142" y="274">V_waste</text><text class="fs-text" x="118" y="290" text-anchor="end">waste</text><path class="fs-pipe" d="M300 245 H440" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M373 238 L387 252 L387 238 L373 252 Z"/><text class="fs-dev" x="380" y="265" text-anchor="middle">V_vacuum</text><text class="fs-text" x="440" y="233" text-anchor="end">CO₂ product</text><path class="fs-pipe" d="M300 40 H440" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M373 33 L387 47 L387 33 L373 47 Z"/><text class="fs-dev" x="380" y="28" text-anchor="middle">V_product</text><text class="fs-text" x="440" y="62" text-anchor="end">light product</text><path class="fs-pipe" d="M215 40 H296" marker-end="url(#fs-ah)"/><path class="fs-valve" d="M243 33 L257 47 L257 33 L243 47 Z"/><text class="fs-dev" x="250" y="28" text-anchor="middle">V_lpp</text><text class="fs-text" x="215" y="62" text-anchor="start">N₂ in</text></svg></div><div><svg viewBox="0 0 440 112" role="img" aria-label="Cycle: SG: pressurisation 20 s, adsorption 46.2 s, blowdown 42.21 s, evacuation 46.02 s, idle 69.27 s; 13X: LPP 20 s, adsorption 46.2 s, blowdown 56.3 s, evacuation 101.2 s"><text class="fs-dev" x="44" y="22" text-anchor="end">SG</text><rect class="fs-seg " x="52.5" y="10" width="32.8" height="16" style="fill:var(--mc-step-pressurise)"><title>pressurisation, 20 s</title></rect><rect class="fs-seg " x="86.3" y="10" width="77.1" height="16" style="fill:var(--mc-step-adsorb)"><title>adsorption, 46.2 s</title></rect><rect class="fs-seg " x="164.4" y="10" width="70.3" height="16" style="fill:var(--mc-step-blowdown)"><title>blowdown, 42.21 s</title></rect><rect class="fs-seg " x="235.7" y="10" width="76.8" height="16" style="fill:var(--mc-step-evacuate)"><title>evacuation, 46.02 s</title></rect><rect class="fs-seg fs-idle" x="313.5" y="10" width="116.0" height="16"><title>idle, 69.27 s</title></rect><text class="fs-dev" x="44" y="44" text-anchor="end">13X</text><rect class="fs-seg " x="52.5" y="32" width="32.8" height="16" style="fill:var(--mc-step-lpp)"><title>LPP, 20 s</title></rect><rect class="fs-seg " x="86.3" y="32" width="77.1" height="16" style="fill:var(--mc-step-adsorb)"><title>adsorption, 46.2 s</title></rect><rect class="fs-seg " x="164.4" y="32" width="94.1" height="16" style="fill:var(--mc-step-blowdown)"><title>blowdown, 56.3 s</title></rect><rect class="fs-seg " x="259.5" y="32" width="170.0" height="16" style="fill:var(--mc-step-evacuate)"><title>evacuation, 101.2 s</title></rect><text class="fs-dev" x="52.0" y="68" text-anchor="middle">0</text><text class="fs-dev" x="136.5" y="68" text-anchor="middle">50</text><text class="fs-dev" x="221.0" y="68" text-anchor="middle">100</text><text class="fs-dev" x="305.5" y="68" text-anchor="middle">150</text><text class="fs-dev" x="390.0" y="68" text-anchor="middle">200</text><text class="fs-dev" x="430" y="82" text-anchor="end">time in the cycle [s]</text><rect x="8" y="92" width="12" height="8" rx="1" style="fill:var(--mc-step-pressurise)"/><text class="fs-dev" x="24" y="100">pressurisation</text><rect x="125.8" y="92" width="12" height="8" rx="1" style="fill:var(--mc-step-adsorb)"/><text class="fs-dev" x="141.8" y="100">adsorption</text><rect x="216.8" y="92" width="12" height="8" rx="1" style="fill:var(--mc-step-blowdown)"/><text class="fs-dev" x="232.8" y="100">blowdown</text><rect x="294.4" y="92" width="12" height="8" rx="1" style="fill:var(--mc-step-evacuate)"/><text class="fs-dev" x="310.4" y="100">evacuation</text><rect x="385.4" y="92" width="12" height="8" rx="1" class="fs-idle"/><text class="fs-dev" x="401.4" y="100">idle</text><rect x="436.2" y="92" width="12" height="8" rx="1" style="fill:var(--mc-step-lpp)"/><text class="fs-dev" x="452.2" y="100">LPP</text></svg></div></div>
<figcaption>The silica gel bed takes the wet feed and passes dry gas to the 13X bed during adsorption; it then waits, closed, until the 13X bed finishes its cycle.</figcaption>
</figure>
```

| Setting | Value |
|---|---|
| Feed | 15% CO₂, 82% N₂, 3% H₂O, 25 °C, 1.013 bar, interstitial velocity 0.7 m/s |
| Columns | silica gel 0.41 m long, 13X 1 m long, both 0.3 m inner diameter |
| Pressures | silica gel down to 0.48 bar and then 0.3 bar through its bottom; 13X down to 0.07 bar through its top and then 0.03 bar through its bottom |
| Steps | silica gel 154 s and then closed until the 13X bed finishes its 224 s |

The [Wet flue gas: dual-adsorbent, two-bed VSA](@ref) example runs it and
follows the water into the 13X bed.

## Conduction across the radius

No process in the library is meshed across its radius, so a bed that is
(see [Beds meshed across their radius](@ref)) is checked against an exact
solution instead. A bed of argon, which the sorbent does not take up, starts
at 320 K inside a wall held at 300 K, with no flow, and cools by conduction
across its radius alone. Its temperature then follows the series solution of
transient conduction in a cylinder,

```math
\frac{T - T_w}{T_0 - T_w} = \sum_n \frac{2}{\beta_n J_1(\beta_n)}\,
J_0\!\left(\beta_n \frac{r}{r_{in}}\right) e^{-\beta_n^2\, \mathrm{Fo}},
\qquad J_0(\beta_n) = 0,
\qquad \mathrm{Fo} = \frac{K_r\, t}{\rho_s c_{ps} (1-\varepsilon)\, r_{in}^2}.
```

The largest difference from it, relative to the initial 20 K, at
``\mathrm{Fo} = 0.094`` (``K_r`` = 5 W/(m·K), 300 s, time steps of 2 s):

| Rings | 4 | 8 | 16 |
|---|---|---|---|
| Mocca | 3.0% | 0.78% | 0.29% |
| Without the half ring at the wall | 18% | 9.0% | 4.5% |

Each halving of the rings divides the error by about 4, as for a second-order
method, until the error of the time steps (about 0.2%) is reached. Without
the conduction across the half of the outer ring between its centre and the
wall ([`RingWallResistance`](@ref)), the wall would take the heat as if the
outer ring were at the wall, and the error would be first order: this is why
Mocca includes it. The test suite runs this comparison.

In a process the fronts are sharper. For a breakthrough of the four-step VSA's
feed (15% CO₂ at 1 bar and 0.5 m/s) into a 1 m column of 10 cells along its
length, with ``K_r`` = 0.5 W/(m·K) and ``h_w`` = 30 W/(m²·K), the peak mean
temperature 0.15 m up the bed over 300 s changes as the rings are refined:

| Rings | 4 | 8 | 16 | 32 |
|---|---|---|---|---|
| Peak temperature [K] | 351.40 | 351.00 | 350.79 | 350.71 |
| Without the half ring [K] | 350.04 | 350.20 | 350.38 | 350.51 |

Both converge to about 350.7 K, from either side. The changes from one grid to
the next shrink by 1.9 and then 2.9 times, approaching second order as the
rings become thin against the layer the wall cools; 8 rings are about 0.3 K
from it, in a front 52 K above the feed.
