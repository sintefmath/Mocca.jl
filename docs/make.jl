using Documenter
using Mocca
using Literate


function build_mocca_docs(; build_examples = true)

    function update_footer(content, pth)
        return content*"\n\n # ## Example on GitHub\n "*
        "# If you would like to run this example yourself, it can be downloaded from "*
        "[the Mocca.jl GitHub repository](https://github.com/sintefmath/Mocca.jl/blob/main/examples/$pth.jl)."
    end

    mocca_dir = joinpath(dirname(pathof(Mocca)),"..")

    ## Build examples
    # <example name> => <example path>
    examples = [
        "Simulate DCB" => "simulate_DCB"
        "Simulate cyclic" => "simulate_cyclic"
        "History matching" => "history_matching"
        "Optimization" => "optimization"
    ]

    examples_markdown = []


    for (ex, pth) in examples
        in_pth = joinpath(mocca_dir, "examples", "$pth.jl")
        out_dir = joinpath(mocca_dir, "docs", "src", "examples")
        out_dir_notebooks = joinpath(mocca_dir, "docs", "src", "notebooks")
        push!(examples_markdown, ex => joinpath("examples", "$pth.md"))
        if build_examples
            upd(content) = update_footer(content, pth)
            Literate.markdown(in_pth, out_dir, preprocess = upd, flavor = Literate.DocumenterFlavor())
            Literate.notebook(in_pth, out_dir_notebooks, preprocess = upd, flavor = Literate.DocumenterFlavor())
        end
    end

    ## Make docs

    makedocs(;
        modules = [Mocca],
        warnonly = [:missing_docs],
        sitename="Mocca.jl",
        pages=[
            "Home" => "index.md",
            "Installation" => "installation.md",
            "Quick start" => "quick_start.md",
            "Examples" => examples_markdown
        ],
        format=Documenter.HTML(;
            edit_link="main",
            # size_threshold_ignore = [
            #     "examples/quick_start.md",
            #     "examples/simulate_DCB.md",
            #     "examples/history_matching.md",
            #     "examples/optimization.md",                
            # ]        
            size_threshold=300_000
        )

    )

    # ## Deploy docs

    deploydocs(;
        repo="github.com/sintefmath/Mocca.jl",
        devbranch="main",
        push_preview = true
    )
end

build_mocca_docs(build_examples=true)
