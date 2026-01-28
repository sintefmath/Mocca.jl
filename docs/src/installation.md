```@meta
CurrentModule = Mocca
```

# Installation

First install Julia from [here](https://julialang.org/downloads/). 

Mocca can be downloaded by cloning the repository ....

```bash
git clone mocca ...
```

We recommend running in a specific environment (similar to a virtual environment in python). More information on environments in Julia can be found [here](https://pkgdocs.julialang.org/v1/environments/).

To create an environment in Mocca.jl navigate to the Mocca.jl folder, start the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and type the following at the Julia prompt:

```julia
Pkg.activate(".")
Pkg.instantiate()
```

This will activate the environment in the current directory and install all necessary dependencies. Mocca is now installed and ready to use.

A good starting example to try is [Direct Column Breakthrough simulation](@ref). Bear in mind that the first time you run the code in the Julia REPL it may take several minutes to run as Julia needs to compile all the necessary code. As long as you do not close the REPL, the second time you run the code will be much quicker!