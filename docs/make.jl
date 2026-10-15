using Documenter
using DocumenterVitepress
using Mocca
using Literate


function build_mocca_docs(; build_examples = true, build_notebooks = true)

    # The examples show the progress of their main cyclic steady state when
    # run as scripts; the pages built from them leave it out
    hide_progress(content) = replace(content, r"[;,] progress = true(?=\))" => "")

    function update_footer(content, pth)
        content = hide_progress(content)
        return content*"\n\n # ## Example on GitHub\n "*
        "# If you would like to run this example yourself, it can be downloaded from "*
        "[the Mocca.jl GitHub repository](https://github.com/sintefmath/Mocca.jl/blob/main/examples/$pth.jl)."
    end

    mocca_dir = joinpath(dirname(pathof(Mocca)),"..")

    ## Build examples
    # <example name> => <example path>
    examples = [
        "Getting started" => "getting_started"
        "Building a process from its parts" => "process_from_parts"
        "Breakthrough curve" => "breakthrough"
        "Two-bed VSA with pressure equalisation" => "two_bed_vsa"
        "Wet flue gas: single-bed VSA" => "wet_flue_gas_lpp_vsa"
        "Wet flue gas: dual-adsorbent VSA" => "wet_flue_gas_dual_adsorbent_vsa"
        "Gradients at cyclic steady state" => "flowsheet_vsa_gradients"
        "Optimization" => "optimization"
        "History matching" => "history_matching"
        "A flowsheet by hand" => "unit_network"
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
            if build_notebooks
                Literate.notebook(in_pth, out_dir_notebooks; preprocess = upd, flavor = Literate.DocumenterFlavor(), execute = false)
            end
        end
    end

    ## Make docs

    makedocs(;
        modules = [Mocca],
        warnonly = [:missing_docs, :cross_references],
        sitename = "Mocca.jl",
        authors = "SINTEF Digital",
        repo = "https://github.com/sintefmath/Mocca.jl",
        pages = [
            "Home" => "index.md",
            "Getting started" => [
                "Installation" => "installation.md",
                "How a VSA cycle works" => "vsa_cycle.md",
                "Getting started" => last(first(examples_markdown)),
            ],
            "Manual" => [
                "Model equations" => "model.md",
                "Cycle steps" => "cycle_steps.md",
                "Adsorbents" => "adsorbents.md",
                "Results and plots" => "results.md",
            ],
            "Examples" => examples_markdown[2:end],
            "Reference" => [
                "Functions" => "user_interface.md",
                "Process files" => "json_input.md",
                "Extending Mocca" => "extending.md",
                "Migrating from 0.1.0" => "migration.md",
            ],
        ],
        format = DocumenterVitepress.MarkdownVitepress(
            repo = "https://github.com/sintefmath/Mocca.jl",
            devbranch = "main",
            devurl = "dev",
            deploy_url = "sintefmath.github.io/Mocca.jl",
        ),
    )

    # ## Deploy docs

    DocumenterVitepress.deploydocs(;
        repo = "github.com/sintefmath/Mocca.jl",
        target = joinpath(@__DIR__, "build"),
        branch = "gh-pages",
        devbranch = "main",
        push_preview = true,
    )
end

build_mocca_docs(build_examples = true, build_notebooks = true)
