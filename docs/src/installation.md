```@meta
CurrentModule = Mocca
```

# Installation

## 1. Install Julia

Install Julia 1.10 or later with [juliaup](https://julialang.org/downloads/),
the official installer. For writing and running scripts, [VS Code with the
Julia extension](https://www.julia-vscode.org/) works well: it runs a script
line by line and shows plots in a side pane.

## 2. Make a project folder

Julia keeps the packages of a project in an *environment*, much like a Python
virtual environment. Make a folder for your work, open a terminal in it and
start Julia with that folder as the environment:

```
julia --project=.
```

## 3. Add Mocca and CairoMakie

At the `julia>` prompt:

```julia
using Pkg
Pkg.add(url = "https://github.com/sintefmath/Mocca.jl", rev = "main")
Pkg.add("CairoMakie")
```

!!! warning "Install from GitHub for now"
    These docs describe the development version of Mocca. The version in
    Julia's package registry (0.1.0) has the earlier interface, without
    `read_process`, `Process` or `cyclic_steady_state`, so install from GitHub
    as above until the next release.

CairoMakie is only needed to save plots to files and to run the examples that
draw their own figures.

The next time you start Julia in the folder with `julia --project=.`, the
packages are there.

## 4. Check that it works

```julia
using Mocca
library_processes()
```

The first `using Mocca` in a new session compiles the code, and the first
simulation compiles the solver, which can take a few minutes. Later runs in
the same session take seconds, so keep the session open while you work.

## 5. Run an example

Each example in these docs is a Julia script in the [`examples` folder](https://github.com/sintefmath/Mocca.jl/tree/main/examples)
on GitHub. Download one into your project folder and run it from the Julia
prompt:

```julia
include("json_input.jl")
```

or open it in VS Code and run it a block at a time with Shift+Enter.

Next: [Getting started](@ref "Quick start: processes from files") reads a
process from Mocca's library, runs it to cyclic steady state and measures it.
