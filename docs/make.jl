using VPLRayTracer
using Documenter

DocMeta.setdocmeta!(VPLRayTracer, :DocTestSetup, :(using VPLRayTracer); recursive=true)

makedocs(;
    modules=[VPLRayTracer],
    authors="Alejandro Morales Sierra <alejandro.moralessierra@wur.nl> and contributors",
    repo="https://github.com/AleMorales/VPLRayTracer.jl/blob/{commit}{path}#{line}",
    sitename="VPLRayTracer.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
