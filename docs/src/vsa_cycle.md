```@meta
CurrentModule = Mocca
```

# How a VSA cycle works

A vacuum swing adsorption (VSA) process separates CO₂ from flue gas with a bed
of adsorbent that holds CO₂ much more strongly than N₂. The bed is fed at about
atmospheric pressure, so CO₂ is taken up while N₂ passes through. The bed is
then evacuated, and the CO₂ that comes off is the product. The bed repeats the
same steps over and over.

This page follows the four-step cycle of
[Haghpanah et al. (2013)](https://doi.org/10.1021/ie302658y), which is in
Mocca's library, through one steady cycle.

```@raw html
<figure class="mocca-schematic">
<svg viewBox="0 0 660 330" role="img" aria-label="The four steps of the VSA cycle: pressurisation with feed, adsorption with feed in and nitrogen-rich product out, blowdown through the top, and evacuation through the bottom">
<style>
.cy-bed{fill:var(--vp-c-bg-soft);stroke:var(--vp-c-text-1);stroke-width:1.4}
.cy-load{fill:var(--mc-co2);opacity:.35}
.cy-line{stroke:var(--vp-c-text-1);stroke-width:1.6;fill:none}
.cy-head{fill:var(--vp-c-text-1)}
.cy-text{fill:var(--vp-c-text-2);font:12.5px var(--vp-font-family-base)}
.cy-name{fill:var(--vp-c-text-1);font:600 14px var(--vp-font-family-base)}
.cy-sym{font-style:italic}
.cy-sub{font-size:10px}
.cy-num{fill:var(--vp-c-text-2);font:600 12px var(--vp-font-family-mono)}
</style>
<defs><marker id="cy-ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="cy-head" d="M0,0 L10,5 L0,10 z"/></marker></defs>
<rect class="cy-bed" x="40" y="70" width="40" height="160" rx="6"/>
<rect class="cy-load" x="41" y="195" width="38" height="34"/>
<line class="cy-line" x1="60" y1="270" x2="60" y2="234" marker-end="url(#cy-ah)"/>
<text class="cy-text" x="70" y="256">feed</text>
<text class="cy-num" x="8" y="20">1</text>
<text class="cy-name" x="22" y="20">Pressurisation</text>
<rect x="8" y="28" width="141" height="5" rx="1" style="fill:var(--mc-step-pressurise)"/>
<text class="cy-text" x="60" y="296" text-anchor="middle"><tspan class="cy-sym">P</tspan><tspan class="cy-sub" dy="3">L</tspan><tspan dy="-3">&#8203; </tspan>→ <tspan class="cy-sym">P</tspan><tspan class="cy-sub" dy="3">H</tspan><tspan dy="-3">&#8203; </tspan></text>
<rect class="cy-bed" x="205" y="70" width="40" height="160" rx="6"/>
<rect class="cy-load" x="206" y="131" width="38" height="98"/>
<line class="cy-line" x1="225" y1="66" x2="225" y2="30" marker-end="url(#cy-ah)"/>
<text class="cy-text" x="235" y="48">N₂-rich product</text>
<line class="cy-line" x1="225" y1="270" x2="225" y2="234" marker-end="url(#cy-ah)"/>
<text class="cy-text" x="235" y="256">feed</text>
<text class="cy-num" x="173" y="20">2</text>
<text class="cy-name" x="187" y="20">Adsorption</text>
<rect x="173" y="28" width="141" height="5" rx="1" style="fill:var(--mc-step-adsorb)"/>
<text class="cy-text" x="225" y="296" text-anchor="middle">at <tspan class="cy-sym">P</tspan><tspan class="cy-sub" dy="3">H</tspan><tspan dy="-3">&#8203; </tspan></text>
<rect class="cy-bed" x="370" y="70" width="40" height="160" rx="6"/>
<rect class="cy-load" x="371" y="131" width="38" height="98"/>
<line class="cy-line" x1="390" y1="66" x2="390" y2="30" marker-end="url(#cy-ah)"/>
<text class="cy-text" x="400" y="48">N₂ and some CO₂</text>
<text class="cy-num" x="338" y="20">3</text>
<text class="cy-name" x="352" y="20">Blowdown</text>
<rect x="338" y="28" width="141" height="5" rx="1" style="fill:var(--mc-step-blowdown)"/>
<text class="cy-text" x="390" y="296" text-anchor="middle"><tspan class="cy-sym">P</tspan><tspan class="cy-sub" dy="3">H</tspan><tspan dy="-3">&#8203; </tspan>→ <tspan class="cy-sym">P</tspan><tspan class="cy-sub" dy="3">I</tspan><tspan dy="-3">&#8203; </tspan></text>
<rect class="cy-bed" x="535" y="70" width="40" height="160" rx="6"/>
<rect class="cy-load" x="536" y="131" width="38" height="98" style="opacity:.15"/>
<line class="cy-line" x1="555" y1="234" x2="555" y2="270" marker-end="url(#cy-ah)"/>
<text class="cy-text" x="565" y="256">CO₂ product</text>
<text class="cy-num" x="503" y="20">4</text>
<text class="cy-name" x="517" y="20">Evacuation</text>
<rect x="503" y="28" width="141" height="5" rx="1" style="fill:var(--mc-step-evacuate)"/>
<text class="cy-text" x="555" y="296" text-anchor="middle"><tspan class="cy-sym">P</tspan><tspan class="cy-sub" dy="3">I</tspan><tspan dy="-3">&#8203; </tspan>→ <tspan class="cy-sym">P</tspan><tspan class="cy-sub" dy="3">L</tspan><tspan dy="-3">&#8203; </tspan></text>
</svg>
<figcaption>One bed through the four steps. The shading is the CO₂ held on the adsorbent; it is drawn as a sharp front for clarity, which the real profile is not.</figcaption>
</figure>

```

| Step | Duration | Open | Pressure | What happens |
|---|---|---|---|---|
| Pressurisation | 15 s | feed, at the bottom | 0.1 → 1 bar | Feed gas fills the evacuated bed. |
| Adsorption | 15 s | feed at the bottom, product at the top | 1 bar | CO₂ is taken up near the feed end; N₂-rich gas leaves at the top. |
| Blowdown | 30 s | top | 1 → 0.2 bar | Most of the N₂ left in the bed leaves through the top, with some CO₂. |
| Evacuation | 40 s | bottom | 0.2 → 0.1 bar | The vacuum pump draws the CO₂ off the adsorbent: the product. |

The bottom of the bed is the feed end, `z = 0`, and the top the light product
end, `z = L`. The steps are described on the [Cycle steps](@ref) page.

## Running the cycle

Read the process, run it to cyclic steady state and measure it:

```@example cycle
using Mocca, Mocca.Units, CairoMakie
process = Mocca.read_process("haghpanah_2013_four_stage_vsa")
result = Mocca.cyclic_steady_state(process; progress = false)
perf = Mocca.performance(result)
println("purity $(round(100perf.purity, digits = 1))%, recovery $(round(100perf.recovery, digits = 1))%, ",
    "$(round(perf.specific_energy, digits = 1)) kWh per tonne of CO₂")
```

At this evacuation pressure the cycle recovers only about a third of the CO₂
fed. The figures below show where the rest goes.

## The pressure

The pressure at the two ends of the bed through the cycle:

```@example cycle
stages = Mocca.stages(process) # hide
step_colour = Dict("pressurisation" => "#8a63b8", "adsorption" => "#3d8c47", "blowdown" => "#c48a1f", "evacuation" => "#3f63c9") # hide
t = cumsum(result.timesteps) # hide
fig = Figure(size = (760, 320)) # hide
ax = Axis(fig[1, 1], xlabel = "time in the cycle [s]", ylabel = "pressure [bar]") # hide
t0 = 0.0 # hide
for st in stages # hide
    vspan!(ax, t0, t0 + st.duration, color = (step_colour[st.name], 0.15)) # hide
    text!(ax, t0 + st.duration / 2, 1.08, text = st.name, align = (:center, :bottom), fontsize = 13) # hide
    global t0 += st.duration # hide
end # hide
lines!(ax, t, [s[:Bed][:Pressure][1] / bar for s in result.states], label = "bottom (feed end)", linewidth = 2) # hide
lines!(ax, t, [s[:Bed][:Pressure][end] / bar for s in result.states], label = "top", linewidth = 2, linestyle = :dash) # hide
ylims!(ax, 0, 1.25) # hide
axislegend(ax, position = :rc) # hide
fig # hide
```

Pressurisation takes the bed from 0.1 to 1 bar in about 10 s. During
adsorption the feed end is about 0.02 bar above the top: the pressure drop of
the gas flowing through the bed. Otherwise the two ends move together, as
blowdown takes the bed to 0.2 bar and evacuation to 0.1 bar.

## Where the CO₂ is held

Mole fractions and loadings along the bed at the end of each step:

```@example cycle
bed = process.flowsheet[:Bed]
step_ends = cumsum([st.duration for st in stages])
Mocca.plot_profiles(result.states, result.timesteps, bed, step_ends;
    labels = [st.name for st in stages], unit = :Bed)
```

Evacuated only to 0.1 bar, the adsorbent stays more than half loaded with CO₂
all along the bed. Each cycle uses just the difference between the loading at
the end of adsorption and at the end of evacuation, under 2% of its saturation
capacity. The bed is so full that at the end of adsorption the gas at the top
already holds about 9% CO₂. Blowdown removes most of the N₂, and at the end of
evacuation the gas in the bed is 77 to 97% CO₂.

## What leaves the bed

The flow through each valve, and the CO₂ mole fraction in it, over the cycle:

```@example cycle
Mocca.plot_flowsheet_streams(result.states, process.model, result.timesteps; stages = stages)
```

The feed valve fills the bed during pressurisation and then delivers a steady
1 mol/s during adsorption. The product line carries about 9% CO₂ during
adsorption and up to 40% in blowdown. That CO₂ is lost with the light
product, which is why the recovery is low. The vacuum line delivers 60% CO₂
at first, rising to 97% as the bed empties.

How the purity, recovery and energy are worked out from these streams is on
the [Model equations](@ref) page, under Performance.

## The cyclic steady state

A bed that starts from some initial state ends each cycle a little different
from how it began. After enough cycles the state at the end of a cycle is the
state it started from: the cyclic steady state, where the performance is
measured. [`cyclic_steady_state`](@ref) gets there faster than running cycle
after cycle, by Anderson acceleration and then Newton's method. The largest
relative change of any variable over a cycle, at each iteration:

```@example cycle
fig = Figure(size = (620, 300)) # hide
ax = Axis(fig[1, 1], yscale = log10, xlabel = "iteration", ylabel = "change over a cycle") # hide
scatterlines!(ax, 1:length(result.history), result.history) # hide
fig # hide
```

From the initial state, a bed at 1 bar, Anderson acceleration brings the
change below 10⁻³ in about 270 cycles. Newton's method then takes it below
10⁻¹⁰ in three iterations.

## Changing a setting

Each run from the last steady state of the process starts close to the new
one, so a sweep of a setting is quick. The evacuation pressure from 0.03 to
0.15 bar:

```@example cycle
P_L = [0.03, 0.05, 0.07, 0.1, 0.15] .* bar
runs = map(P_L) do p
    process[:P_L] = p
    Mocca.performance(Mocca.cyclic_steady_state(process; progress = false))
end
fig = Figure(size = (760, 300)) # hide
ax1 = Axis(fig[1, 1], xlabel = "evacuation pressure [bar]", ylabel = "%") # hide
scatterlines!(ax1, P_L ./ bar, [100r.purity for r in runs], label = "purity") # hide
scatterlines!(ax1, P_L ./ bar, [100r.recovery for r in runs], label = "recovery") # hide
axislegend(ax1, position = :lb) # hide
ax2 = Axis(fig[1, 2], xlabel = "evacuation pressure [bar]", ylabel = "specific energy [kWh/t CO₂]") # hide
scatterlines!(ax2, P_L ./ bar, [r.specific_energy for r in runs], color = :black) # hide
fig # hide
```

A lower evacuation pressure takes more CO₂ off the adsorbent each cycle, so
recovery rises steeply: from 12% at 0.15 bar to 82% at 0.03 bar. Purity
changes much less and is highest near 0.07 bar. The energy per tonne is
lowest between 0.07 and 0.1 bar. Below that range the pump does more work for
each mole; above it, too little CO₂ is recovered to share the work.

To go further, the [Gradients at cyclic steady state](@ref) example gets the
sensitivity to every setting at once, and [Optimization of a cyclic
process](@ref) searches for the best settings.
