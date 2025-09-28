using MoleculeScreen
using Documenter

DocMeta.setdocmeta!(MoleculeScreen, :DocTestSetup, :(using MoleculeScreen, MoleculeFlow); recursive=true)

makedocs(;
    modules=[MoleculeScreen],
    authors="Renee Gil",
    sitename="MoleculeScreen.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MoleculeHub.github.io/MoleculeScreen",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/MoleculeHub/MoleculeScreen.git",
    devbranch="main",
)