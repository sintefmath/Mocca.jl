```@meta
CurrentModule = Mocca
```

# Optimisation

Mocca optimises a process at the state it settles to: the cyclic steady state
of a cyclic process, or the steady state of a moving bed. Each design the
optimiser tries is settled, measured and differentiated, and the gradients
come from adjoints, so a design of a few variables needs tens of evaluations
rather than the thousands of a genetic algorithm. The same steps serve every
kind of [`Process`](@ref):

1. name the inputs to change, with their bounds, as design variables;
2. pick the measure to maximise or minimise and the targets the others must
   meet, from the built-in measures or your own;
3. call [`optimise!`](@ref), or first [`explore`](@ref) the bounds to see
   where the targets can be met.

```julia
using Mocca, Mocca.Units

p = four_stage_vsa_process()
add_variable!(p, :P_L, ("evacuation", :to); bounds = (0.05bar, 0.5bar))
add_variable!(p, :P_I, ("blowdown", :to); bounds = (0.05bar, 0.5bar))
add_variable!(p, :t_ads, ("adsorption", :duration); bounds = (10.0, 100.0))

best = optimise!(p; maximise = :recovery, at_least = (purity = 0.9,),
    linear_constraints = [(Dict(:P_L => 1.0, :P_I => -1.0), -0.01bar)])
best.values, best.metrics, best.feasible
```

The [Optimization](examples/optimization.md) example runs this with the
feed velocity as well, and then minimises the energy from the optimum it
finds.

## Design variables

A design variable is a named number that sets one or more places in the
description of a process: a field of a step, the feed, a parameter of a
unit, a device setting or the geometry of a bed. [`add_variable!`](@ref)
adds one, and `p[name]` reads and sets it:

```julia
add_variable!(p, :P_L, ("evacuation", :to); bounds = (0.05bar, 0.5bar),
    description = "evacuation pressure")
p[:P_L] = 0.1bar
```

Setting a variable recompiles the cycle, or rebuilds the unit it changes, so
the process is consistent at once; nothing is simulated until the process
is run. A process read from the library with [`read_process`](@ref) comes
with the design variables of its study and their bounds, such as `P_L`,
`P_I`, `v_feed` and the step durations of the four-stage VSA;
`Mocca.variable_names(p)` lists them, and the optimisers use all of them
unless given `variables`. The places, with the tuple each can be written as:

| Place | Written as | Sets | Process |
|---|---|---|---|
| [`StepPlace`](@ref) | `("evacuation", :to)`, `(3, :duration)`, `("A", "evacuation", :to)` | a field of a step, by its name or position, in the sequence labelled `"A"` if the cycle has several | cyclic |
| [`FeedPlace`](@ref) | `(:feed, :pressure)` | the feed's `:pressure` or `:temperature` | cyclic |
| [`ParameterPlace`](@ref) | `(:parameter, :Bed, :AdsorbentDensity)` | a parameter of a unit, one value throughout it | both |
| [`SettingPlace`](@ref) | `(:setting, :V_purge, :rate)`, `(:setting, :V_product, :outlet, :pressure)` | a field of a device's flow law, or the pressure or temperature at a boundary port | steady |
| [`GeometryPlace`](@ref) | `(:geometry, :Bed, :length)` | a field of a bed's `Column`, or a moving bed's `:packing_fraction`; the unit is rebuilt on the same number of cells | both |
| [`ScaledPlace`](@ref) | `ScaledPlace(place, factor)` | `place` at `factor` times the variable | both |

A few things to know when choosing them:

- **One variable can set several places**, which then move together, such
  as the sorbent flow of every section of a moving bed:
  `add_variable!(p, :sorbent_flow, ((:parameter, s, :SorbentMassFlow) for s in sections)...)`.
  A `ScaledPlace` sets a place in proportion instead, so that a variable can
  be in other units than its place, such as a purge as a percentage of the
  feed.
- **Pressures carry over from step to step**, so one place is usually enough
  for a pressure: the evacuation pressure is also where the next
  pressurisation starts (see [Cycle steps](cycle_steps.md)). A design
  that a step cannot follow, such as a blowdown that would raise the
  pressure, cannot be built; keep the optimiser away from it with a linear
  constraint, as in the example above.
- **A stage duration keeps its time steps.** While a duration is a design
  variable, each stage keeps the number of time steps it has and stretches
  them, so that the performance changes smoothly with the duration and the
  gradients describe the model the optimiser sees (`keep_step_pattern` of
  [`TimeStepping`](@ref)). [`optimise!`](@ref) takes the durations it starts
  from as the pattern; [`reset_step_pattern!`](@ref) does so by hand.
- **A variable changes values, not structure.** The number of cells or of
  stages, the isotherm, or which devices a step opens are not places. To
  compare such choices, build a process for each and optimise each.

## Measures

The measures of the performance are what [`performance`](@ref) reports, and
any of them can be the objective or a target:

| Measure | Unit | |
|---|---|---|
| `purity` | – | mole fraction of the key component in the extract |
| `recovery` | – | fraction of the key component fed that leaves in the extract |
| `specific_energy` | kWh/t | work of the vacuum pumps, compressors and blowers per tonne of key component in the extract |
| `productivity` | mol/(m³·s) | key component in the extract per m³ of adsorbent per second |
| `specific_heat_duty` | J/kg | heat given by heat transfer fluids and heating jackets per kg of key component; with heating only |
| `specific_equivalent_work` | kWh/t | the energy plus the work equivalent to the heat, each bit weighted by its Carnot factor; with heating only |

The key component is the one the process captures, `key_component` of the
process, by default the first component of the gas. Each measure is a
function of a few integrals over the cycle, such as `key_extract`, the moles
of the key component through the extract devices, and of quantities that do
not depend on the state, such as the cycle time `duration`. A measure of
your own is a function of those, and gets its gradient the same way:

```julia
add_integral!(p, :n2_extract, StreamObjective(:V_vacuum; component = 2))
add_measure!(p, :nitrogen_fraction, s -> s.n2_extract / s.total_extract)
```

[`performance_integrals`](@ref), [`performance_quantities`](@ref) and
[`performance_measures`](@ref) list what `s` holds. The
[Optimising direct air capture](examples/dac_optimisation.md) example adds
the production rate and the heat per kg of CO₂ of its paper this way.

## Gradients

[`performance_gradient`](@ref) settles the process, measures it and
differentiates every measure with respect to the design variables:

```julia
g = performance_gradient(p; variables = [:P_L, :t_ads])
g.metrics.recovery          # the performance at the steady state
g.gradients.recovery.P_L    # ∂ recovery / ∂ P_L, per Pa
```

The gradient goes through two steps. First, how the inputs of the
simulation change with each variable: the device settings and duration of
every stage, the parameters of every unit and quantities such as the
adsorbent volume. Setting the variable a little up and down and recompiling
gives these without simulating, so a new kind of place needs no gradient
code of its own. Second, how the measures change with those inputs, at the
steady state. The cyclic steady state solves ``x = \Phi(x, \theta)``, where
``\Phi`` maps the state at the start of a cycle to the state at its end, so
it moves with the inputs ``\theta``; the implicit function theorem gives the
gradient of a measure from adjoint solves over one cycle (see
[`cyclic_steady_state_gradient`](@ref)). A steady state is the fixed point
of one long backward-Euler step, and its gradient follows the same way.

The cost is one steady state and an adjoint solve per integral the measures
depend on (or per measure, when there are fewer measures), and hardly grows
with the number of variables. Two things to keep in mind:

- **The gradients are only as accurate as the steady state.** Newton's
  method converges it far enough by default (see
  [Cyclic steady state](steady_state.md)). The cycle map must also be smooth:
  each stage should start with short time steps, as the default
  `TimeStepping(first_dt = 0.05)` does.
- **Some inputs stay fixed.** The gas properties of a bed stay as they were
  when the process was built, and so does its axial dispersion unless a
  variable sets its `AxialDispersion` parameter.

[`check_gradients`](@ref) compares the gradients with central differences at
two step sizes, and marks each pair of variable and measure as passing,
inconclusive (the performance has a kink or a jump there, such as a pump
switching on) or failing. Run it when adding a kind of place, unit or
measure (see [Extending Mocca](extending.md)).

## Optimising

[`optimise!`](@ref) maximises or minimises one measure over the box of the
variables' bounds, with targets on others:

```julia
optimise!(p; minimise = :specific_energy,
    at_least = (purity = 0.9, recovery = 0.6),
    scale = (specific_energy = 100.0,),
    linear_constraints = [(Dict(:P_L => 1.0, :P_I => -1.0), -0.01bar)],
    max_it = 15, log = (values, r) -> println(values, " ", r.metrics))
```

- `maximise` or `minimise`: the measure to optimise, one of the two.
- `at_least` and `at_most`: targets on other measures.
- `linear_constraints`: pairs `(Dict(name => coefficient, ...), bound)`,
  each meaning ``\sum \text{coefficient} \cdot \text{variable} \le
  \text{bound}``.
- `scale`: the unit each measure is compared in, 1 by default. Choose it so
  that the objective and the targets are of order one, such as the energy in
  units of 100 kWh/t.
- `variables`: the variables to change, by default all. Every one needs
  finite bounds; one whose bounds are equal stays at that value.

### How it works

The variables are scaled to the unit box of their bounds, and BFGS within
the box (`Jutul.unit_box_bfgs`) follows the gradients for up to `max_it`
iterations per run. The targets are met with an augmented Lagrangian: each
run maximises the objective less a quadratic penalty, of `weight`, on the
shortfall from each target, shifted by a multiplier per target. After each
run the multipliers grow by the shortfall left, and the penalty grows
tenfold where the shortfall did not shrink enough. Runs follow one another,
each from the best design of the last, until every target is met to
`constraint_tol` (in units of its scale) and a run converges within
`max_it` iterations, or after `outer_iterations` runs. With
`constraints = :penalty` only the first run is made, the plain penalty
method, which may stop a little short of the targets.

Each evaluation finds the steady state starting from the last one found, so
it takes a few Newton iterations once the optimiser takes small steps. A
design that cannot be built, fails to simulate or does not settle is left
out with a warning, and the optimiser steps back from it. A design the
optimiser asks for twice is evaluated once.

The process is left at the best design: the best objective among the
designs that meet the targets, or, if none does, the one closest to meeting
them. The result holds:

| Field | |
|---|---|
| `values` | the best design, a NamedTuple by variable |
| `metrics` | its performance |
| `feasible` | whether it meets every target to `constraint_tol` |
| `multipliers` | the multiplier of each target |
| `evaluations` | the number of designs evaluated |
| `history` | every design evaluated, with its `values`, `metrics`, whether it `settled`, and its `violation` of the targets |

### Without gradients

With `method = :pattern`, a compass search in the unit box needs no
gradients: it steps along each variable in turn, takes the first step that
improves the objective less a penalty of `weight` on any shortfall, and
halves the steps when none does, for up to `max_it` evaluations. It suits
an objective with kinks or jumps, such as near a change of regime, where
the gradients point the wrong way.

### When it stops early

BFGS stops when the objective barely changes from one iteration to the
next. Along a flat ridge it can stop partway: in the
[Optimization](examples/optimization.md) example the recovery keeps rising
slowly as the feed slows and the adsorption lengthens, and a longer run
goes further than the default. Raise `max_it`, or run [`optimise!`](@ref)
again from where it stopped, when the history is still improving at the
end.

## Trade-offs

[`pareto_front!`](@ref) traces the trade-off between two measures by
holding one at each of a series of targets and optimising the other:

```julia
front = pareto_front!(p; minimise = :specific_energy,
    sweep = :recovery => [0.6, 0.7, 0.8], at_least = (purity = 0.9,))
[(pt.target, pt.metrics.specific_energy, pt.feasible) for pt in front]
```

The swept measure is held `at_least` each target, or at most with
`sense = :at_most`, and the other keywords go to [`optimise!`](@ref). Each
optimisation starts from the optimum of the one before, so order the
targets so that they tighten. The process is left at the last optimum.

## Looking around first

A process can behave quite differently in different parts of the box of its
bounds: a purge too small to strip the sorbent, or a target that only part of
the box can meet. [`explore`](@ref) evaluates designs spread over the box
(points of a Halton sequence), each settled from the initial state so that
one regime does not carry over to the next, and leaves the process as it
was:

```julia
points = explore(p; samples = 10)
[(pt.values, pt.metrics.recovery) for pt in points if pt.settled]
```

Use it to see where the targets can be met, and start [`optimise!`](@ref)
from a promising design by setting the variables to it.

## Examples

- [Gradients at cyclic steady state](examples/flowsheet_vsa_gradients.md):
  the gradients of a VSA cycle, checked against finite differences.
- [Optimization](examples/optimization.md): the four-stage VSA, for the most
  recovery and then the least energy, with targets and a linear constraint.
- [Optimising direct air capture](examples/dac_optimisation.md): measures of
  your own, and the desorption time and steam flow of a TVSA.
- [Optimising the moving bed TSA](examples/mbtsa_optimisation.md): a steady
  process, with scaled, tied and geometry variables, [`explore`](@ref) and a
  Pareto front.
- [History matching](examples/history_matching.md): fitting the parameters
  of a bed to a measured breakthrough curve, with Jutul's optimisation and
  adjoints directly rather than [`optimise!`](@ref).
