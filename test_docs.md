```@meta
EditURL = "test_docs.jl"
```

# **8.** Example

[![](https://mybinder.org/badge_logo.svg)](<unknown>/generated/example.ipynb)
[![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](<unknown>/generated/example.ipynb)

This is an example generated with Literate based on this
source file: [`example.jl`](<unknown>/examples/example.jl).
You are seeing the
HTML-output which Documenter has generated based on a markdown
file generated with Literate. The corresponding notebook
can be viewed in [nbviewer](http://nbviewer.jupyter.org/) here:
[`example.ipynb`](<unknown>/generated/example.ipynb),
and opened in [binder](https://mybinder.org/) here:
[`example.ipynb`](<unknown>/generated/example.ipynb),
and the plain script output can be found here: [`example.jl`](./example.jl).

It is recommended to have the [source file](<unknown>/examples/example.jl)
available when reading this, to better understand how the syntax in the source file
corresponds to the output you are seeing.

### Basic syntax
The basic syntax for Literate is simple, lines starting with `# ` is interpreted
as markdown, and all the other lines are interpreted as code. Here is some code:

````@example test_docs
x = 1//3
y = 2//5
````

In markdown sections we can use markdown syntax. For example, we can
write *text in italic font*, **text in bold font** and use
[links](https://www.youtube.com/watch?v=dQw4w9WgXcQ).

It is possible to filter out lines depending on the output using the
`#md`, `#nb`, `#jl` and `#src` tags (see [Filtering lines](@ref)):
- This line starts with `#md` and is thus only visible in the markdown output.

The source file is parsed in chunks of markdown and code. Starting a line
with `#-` manually inserts a chunk break. For example, if we want to
display the output of the following operations we may insert `#-` in
between. These two code blocks will now end up in different
`@example`-blocks in the markdown output, and two different notebook cells
in the notebook output.

````@example test_docs
x + y
````

````@example test_docs
x * y
````

### Output capturing
Code chunks are by default placed in Documenter `@example` blocks in the generated
markdown. This means that the output will be captured in a block when Documenter is
building the docs. In notebooks the output is captured in output cells, if the
`execute` keyword argument is set to true. Output to `stdout`/`stderr` is also
captured.

!!! note
    Note that Documenter currently only displays output to `stdout`/`stderr`
    if there is no other result to show. Since the vector `[1, 2, 3, 4]` is
    returned from `foo`, the printing of `"This string is printed to stdout."`
    is hidden.

````@example test_docs
function foo()
    println("This string is printed to stdout.")
    return [1, 2, 3, 4]
end

foo()
````

Just like in the REPL, outputs ending with a semicolon hides the output:

````@example test_docs
1 + 1;
nothing #hide
````

Both Documenter's `@example` block and notebooks can display images. Here is an example
where we generate a simple plot using the
[Plots.jl](https://github.com/JuliaPlots/Plots.jl) package

````@example test_docs
using Plots
x = range(0, stop=6π, length=1000)
y1 = sin.(x)
y2 = cos.(x)
plot(x, [y1, y2])
````

### Custom processing

It is possible to give Literate custom pre- and post-processing functions.
For example, here we insert a placeholder value `x = 123` in the source, and use a
preprocessing function that replaces it with `y = 321` in the rendered output.

````@example test_docs
x = 123
````

In this case the preprocessing function is defined by

````@example test_docs
function pre(s::String)
    s = replace(s, "x = 123" => "y = 321")
    return s
end
````

### [Documenter.jl interaction](@id documenter-interaction)

In the source file it is possible to use Documenter.jl style references,
such as `@ref` and `@id`. These will be filtered out in the notebook output.
For example, [here is a link](@ref documenter-interaction), but it is only
visible as a link if you are reading the markdown output. We can also
use equations:

```math
\int_\Omega \nabla v \cdot \nabla u\ \mathrm{d}\Omega = \int_\Omega v f\ \mathrm{d}\Omega
```

using Documenters math syntax. Documenters syntax is automatically changed to
`\begin{equation} ... \end{equation}` in the notebook output to display correctly.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

