```@meta
CurrentModule = Mocca
```

# Cycle steps

A [`Cycle`](@ref) is a list of steps. Each step opens the devices of one bed
that it needs, found by their roles (see [`add_device!`](@ref)), and keeps the
others closed. Pressures carry over from step to step: a step that changes the
pressure starts where the step before left the bed, unless given `from`.

```@raw html
<figure class="mocca-schematic">
<svg viewBox="0 0 640 500" role="img" aria-label="The step types and the ends of the bed each one opens">
<style>
.st-bed{fill:var(--vp-c-bg-soft);stroke:var(--vp-c-text-1);stroke-width:1.4}
.st-line{stroke:var(--vp-c-text-1);stroke-width:1.6;fill:none}
.st-head{fill:var(--vp-c-text-1)}
.st-text{fill:var(--vp-c-text-2);font:12.5px var(--vp-font-family-base)}
.st-name{fill:var(--vp-c-text-1);font:600 14px var(--vp-font-family-base)}
</style>
<defs><marker id="st-ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path class="st-head" d="M0,0 L10,5 L0,10 z"/></marker></defs>
<text class="st-name" x="10" y="22">Pressurise</text>
<rect x="10" y="30" width="130" height="5" rx="1" style="fill:var(--mc-step-pressurise)"/>
<rect class="st-bed" x="40" y="85" width="32" height="110" rx="5"/>
<line class="st-line" x1="56.0" y1="235" x2="56.0" y2="199" marker-end="url(#st-ah)"/>
<text class="st-text" x="68.0" y="221">feed</text>
<text class="st-name" x="170" y="22">Pressurise, LPP</text>
<rect x="170" y="30" width="130" height="5" rx="1" style="fill:var(--mc-step-lpp)"/>
<rect class="st-bed" x="200" y="85" width="32" height="110" rx="5"/>
<line class="st-line" x1="216.0" y1="45" x2="216.0" y2="81" marker-end="url(#st-ah)"/>
<text class="st-text" x="228.0" y="67">light product</text>
<text class="st-name" x="330" y="22">Adsorb</text>
<rect x="330" y="30" width="130" height="5" rx="1" style="fill:var(--mc-step-adsorb)"/>
<rect class="st-bed" x="360" y="85" width="32" height="110" rx="5"/>
<line class="st-line" x1="376.0" y1="81" x2="376.0" y2="45" marker-end="url(#st-ah)"/>
<text class="st-text" x="388.0" y="67">product</text>
<line class="st-line" x1="376.0" y1="235" x2="376.0" y2="199" marker-end="url(#st-ah)"/>
<text class="st-text" x="388.0" y="221">feed</text>
<text class="st-name" x="490" y="22">Blowdown</text>
<rect x="490" y="30" width="130" height="5" rx="1" style="fill:var(--mc-step-blowdown)"/>
<rect class="st-bed" x="520" y="85" width="32" height="110" rx="5"/>
<line class="st-line" x1="536.0" y1="81" x2="536.0" y2="45" marker-end="url(#st-ah)"/>
<text class="st-text" x="548.0" y="67">product line</text>
<text class="st-name" x="10" y="272">Evacuate</text>
<rect x="10" y="280" width="130" height="5" rx="1" style="fill:var(--mc-step-evacuate)"/>
<rect class="st-bed" x="40" y="335" width="32" height="110" rx="5"/>
<line class="st-line" x1="56.0" y1="449" x2="56.0" y2="485" marker-end="url(#st-ah)"/>
<text class="st-text" x="68.0" y="471">extract</text>
<text class="st-name" x="170" y="272">Purge</text>
<rect x="170" y="280" width="130" height="5" rx="1" style="fill:var(--mc-step-lpp)"/>
<rect class="st-bed" x="200" y="335" width="32" height="110" rx="5"/>
<line class="st-line" x1="216.0" y1="295" x2="216.0" y2="331" marker-end="url(#st-ah)"/>
<text class="st-text" x="228.0" y="317">light product</text>
<line class="st-line" x1="216.0" y1="449" x2="216.0" y2="485" marker-end="url(#st-ah)"/>
<text class="st-text" x="228.0" y="471">extract</text>
<text class="st-name" x="330" y="272">Equalise</text>
<rect x="330" y="280" width="130" height="5" rx="1" style="fill:var(--mc-step-equalise)"/>
<rect class="st-bed" x="348" y="335" width="32" height="110" rx="5"/>
<rect class="st-bed" x="412" y="335" width="32" height="110" rx="5"/>
<path class="st-line" d="M364.0 331 V309 H428.0 V327" marker-end="url(#st-ah)"/>
<text class="st-text" x="396.0" y="301" text-anchor="middle">link, high to low</text>
<text class="st-name" x="490" y="272">Idle</text>
<rect x="490" y="280" width="130" height="5" rx="1" style="fill:var(--mc-step-equalise)"/>
<rect class="st-bed" x="520" y="335" width="32" height="110" rx="5"/>
<text class="st-text" x="548.0" y="317">all closed</text>
</svg>
<figcaption>The step types, with the ends of the bed each one opens by default. The bottom of each bed is its feed end.</figcaption>
</figure>

```

| Step | Opens | Pressure | Default name |
|---|---|---|---|
| [`Pressurise`](@ref) | the feed device, or the light inlet with `source = :light_product` | rises to `to`, by default the feed pressure | `pressurisation`, or `LPP` |
| [`Adsorb`](@ref) | the feed device and the product device at the other end | held at `pressure`, by default the feed pressure | `adsorption` |
| [`Blowdown`](@ref) | an outlet at its `through` end, by default the product device at the top | falls to `to` | `blowdown` |
| [`Evacuate`](@ref) | an outlet at its `through` end, by default the extract device at the bottom | falls to `to` | `evacuation` |
| [`Purge`](@ref) | the light inlet and the extract device at the other end | held at `pressure`, by default where the step before left it | `purge` |
| [`Equalise`](@ref) | the link between two beds | both meet between their start pressures | `equalisation A to B` |
| [`Idle`](@ref) | nothing | unchanged | `idle` |

Pressures in steps are in Pa; `Mocca.Units` has `bar`, so `0.2bar` is
20000 Pa. While the pressure changes, the pressure on the far side of the
valve follows an exponential ramp (see [Model equations](@ref)).

## Pressurise

Raises the bed from where the step before left it to the feed pressure, or to
`to`, through the feed device. With `source = :light_product` it takes light
product through the light inlet at the top instead: light product
pressurisation (LPP), as in the wet flue gas process of Krishnamurthy et al.
(2014).

```julia
Pressurise(15.0)                                  # feed, to the feed pressure
Pressurise(20.0; source = :light_product)         # LPP
```

```json
{"type": "pressurise", "duration": 15}
{"type": "pressurise", "duration": 20, "source": "light_product"}
```

## Adsorb

Feeds the bed at constant flow, with the light product leaving the other end
at the feed pressure, or at `pressure`. Give the flow as one of an interstitial
`velocity` [m/s], a `superficial_velocity` [m/s] or a `flow_rate` [m³/s at the
bed inlet].

```julia
Adsorb(15.0; velocity = 1.0)
```

```json
{"type": "adsorb", "duration": 15, "velocity": 1.0}
```

A bed without a feed device of its own is fed through a link from the bed
before it, as in the dual-adsorbent process, and takes no flow setting.

## Blowdown

Lowers the bed to `to` through an outlet: by default the product device at the
top (co-current). With `through = :bottom` it goes through a waste or extract
device at the bottom instead (counter-current).

```julia
Blowdown(30.0; to = 0.2bar)
```

```json
{"type": "blowdown", "duration": 30, "to": "0.2 bar"}
```

## Evacuate

Lowers the bed to `to` through the extract (vacuum) device at the bottom,
which takes the CO₂ product. Like blowdown, `through` picks the other end.

```julia
Evacuate(40.0; to = 0.1bar)
```

```json
{"type": "evacuate", "duration": 40, "to": "0.1 bar"}
```

## Purge

Light product enters through the light inlet at a constant flow and leaves
through the extract device at the other end, sweeping the CO₂ out. The flow is
given as in `Adsorb`.

```julia
Purge(20.0; velocity = 0.2)
```

```json
{"type": "purge", "duration": 20, "velocity": 0.2}
```

## Equalise

Opens the link between two beds, such as their tops, so that gas flows from
the bed at the higher pressure into the other. The size of the link valve, its
`conductance`, sets how far the pressures come together in the time given.
Afterwards both beds are taken to be at the mean of their start pressures,
where the next ramps start.

```julia
Equalise(10.0; beds = (:A, :B), conductance = 1e-6)
```

```json
{"type": "equalise", "duration": 10, "beds": ["A", "B"], "conductance": 1e-6}
```

The [Two-bed VSA with pressure equalisation](@ref) example builds a process
with it.

## Idle

Keeps every device of the bed closed, for example while another bed finishes
its steps.

```julia
Idle(10.0; bed = :B)
```

```json
{"type": "idle", "duration": 10, "bed": "B"}
```

## Several beds

With several beds, each step names its `bed`. Steps run one after the other;
[`Concurrent`](@ref) runs steps on different beds at the same time, as one
stage:

```julia
Concurrent(Pressurise(15.0; bed = :A), Evacuate(15.0; bed = :B, to = 0.1bar))
```

Beds can also follow sequences of their own, side by side, each a list of
steps with a label. A shorter sequence waits with its devices closed until the
longest one ends:

```julia
Cycle("silica gel" => [...], "13X" => [...])
```

```json
{"sequences": [{"label": "silica gel", "steps": [...]}, {"label": "13X", "steps": [...]}]}
```

The [Wet flue gas: dual-adsorbent, two-bed VSA](@ref) example runs a cycle like this.
New kinds of step can be added with [`register_step!`](@ref); see
[Extending Mocca](@ref).
