```@meta
CurrentModule = Mocca
```

The latest stable version of Mocca can be installed directly from Julia. First install Julia from [here](https://julialang.org/downloads/).

## Working in an environment (optional)
We recommend running in a specific environment (similar to a virtual environment in python). More information on environments in Julia can be found [here](https://pkgdocs.julialang.org/v1/environments/).

To create an environment in Julia, navigate to the folder where you want the environment to be, start the [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and type the following at the Julia prompt:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Installing Mocca in your environment

To install Mocca just add the package to your current environment in the Julia REPL:

```julia
Pkg.add("Mocca")
```

This will add Mocca to the current environment and install all necessary dependencies. Mocca is now installed and ready to use.

To get started try the [Quick start example](@ref) or [Direct Column Breakthrough simulation](@ref) examples. Bear in mind that the first time you run the code in the Julia REPL it may take several minutes to run as Julia needs to compile all the necessary code. As long as you do not close the REPL, the second time you run the code will be much quicker!
