using PlantRayTracer
using Documenter

DocMeta.setdocmeta!(PlantRayTracer, :DocTestSetup, :(using PlantRayTracer); recursive=true)

makedocs(;
    modules=[PlantRayTracer],
    authors="Alejandro Morales Sierra <alejandro.moralessierra@wur.nl> and contributors",
    repo="https://github.com/AleMorales/PlantRayTracer.jl/blob/{commit}{path}#{line}",
    sitename="PlantRayTracer.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
