```@meta
CurrentModule = Mocca
```


# Mocca

[Mocca.jl](https://github.com/sintefmath/Mocca.jl) provides a [Julia](https://julialang.org/) based framework for the simulating pressure / temperature swing adsorption processes for CO2 capture.

A process is described in layers — a gas and its feed, adsorbents and columns, a flowsheet of beds and devices, and a cycle of steps — each of which can be edited, saved to a file and reused on its own; see [Setting up processes](@ref) and [Process files](@ref). Mocca finds the cyclic steady state of a process, measures its performance, and computes the gradients of its performance with respect to its operating settings and parameters by adjoints, for optimisation.

Mocca's library has the 4-stage vacuum swing adsorption process for CO2 capture of Haghpanah et al. (2013), from a two-component flue gas on zeolite 13X with a dual-site Langmuir isotherm, and the wet flue gas processes of Krishnamurthy et al. (2014). See [Haghpanah DCB](https://github.com/sintefmath/Mocca.jl/blob/main/examples/dcb_haghpanah_2013_co2_n2.jl) and [Haghpanah Cyclic VSA](https://github.com/sintefmath/Mocca.jl/blob/main/examples/cyclic_vsa_haghpanah_2013_co2_n2.jl). For an example of setting up a simulation through defining your own parameters and potentially custom models in a script, see [Custom setup cyclic VSA](https://github.com/sintefmath/Mocca.jl/blob/main/examples/custom_setup_cyclic_vsa.jl). Additionally, we have made examples demonstrating capabilities for doing [Optimization](https://github.com/sintefmath/Mocca.jl/blob/main/examples/optimization.jl) and [History matching](https://github.com/sintefmath/Mocca.jl/blob/main/examples/history_matching.jl) in Mocca.jl.

In the future we hope to implement examples of other systems and isotherms e.g. temperature swing adsorption for Direct Air Capture (DAC).


```@index
```
