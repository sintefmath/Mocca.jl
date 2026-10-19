```@meta
CurrentModule = Mocca
```

# Mocca

Mocca simulates cyclic adsorption processes for CO₂ capture, such as pressure
and vacuum swing adsorption (PSA/VSA). It finds the cyclic steady state of a
process, measures its purity, recovery, energy use and productivity, and
computes how these change with every operating setting, for optimisation. It is
written in [Julia](https://julialang.org/) and built on the
[Jutul](https://github.com/sintefmath/Jutul.jl) simulator.

## Quick start

Read the four-step VSA process of
[Haghpanah et al. (2013)](https://doi.org/10.1021/ie302658y) from Mocca's
library, lower its evacuation pressure, run it to cyclic steady state and
measure it:

```julia
using Mocca
using Mocca.Units

process = read_process("haghpanah_2013_four_stage_vsa")
process[:P_L] = 0.05bar                  # evacuation pressure
result = cyclic_steady_state(process)
perf = performance(result)
```

```julia-repl
julia> perf.purity, perf.recovery, perf.specific_energy
(0.7900223577190472, 0.7088768107066176, 149.85070321200155)
```

That is 79.0% CO₂ purity, 70.9% recovery and 150 kWh per tonne of CO₂. The
first run in a session takes a few minutes while Julia compiles; after that,
changing a setting and finding the new steady state takes seconds.

## What Mocca models

- **Beds**: one-dimensional, axially dispersed plug flow through a packed
  column, with Darcy's law for the pressure drop, a linear driving force for
  mass transfer, and an energy balance including the column wall (or an
  adiabatic or isothermal bed). See [Model equations](@ref).
- **Isotherms**: multicomponent dual-site and single-site Langmuir, with a
  library of published adsorbents (zeolite 13X, silica gel).
- **Flowsheets**: any number of beds connected to the feed, the products and
  each other through valves and flow controllers.
- **Cycles**: pressurisation (with feed or light product), adsorption,
  blowdown, evacuation, purge, pressure equalisation and idle steps, run in
  sequence or side by side in several beds.
- **Design**: gradients of the performance at cyclic steady state with
  respect to step durations, pressures, flow rates and model parameters, and
  gradient-based optimisation.
- **Ready-made processes**: the four-step VSA of Haghpanah et al. (2013), a
  two-bed VSA with pressure equalisation, and the wet flue gas processes of
  [Krishnamurthy et al. (2014)](https://doi.org/10.1021/ie5024723).

Not yet included: temperature swing (heated or cooled steps), the inertial
(Ergun) term of the pressure drop, and other isotherm models (new ones can be
added with [`register_isotherm!`](@ref)).

## Where next

| To | Read |
|---|---|
| Install Julia and Mocca | [Installation](@ref) |
| Run, change and save a process step by step | [Getting started](@ref "Quick start: processes from files") |
| See the equations solved and how performance is defined | [Model equations](@ref) |
| Get numbers and plots out of a simulation | [Results and plots](@ref) |
| Set up your own gas, adsorbent, column and cycle | [Building a process from its parts](@ref) |
| Optimise a process | [Optimization of a cyclic process](@ref) |
| Look up a function, or the process file format | [Reference](@ref), [Process files](@ref) |
