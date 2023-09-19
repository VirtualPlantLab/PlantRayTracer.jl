using PlantRayTracer
using Documenter

makedocs(;
    doctest = false,
    modules = [PlantRayTracer],
    authors = "Alejandro Morales Sierra <alejandro.moralessierra@wur.nl> and contributors",
    repo = "https://github.com/VirtualPlantLab/PlantRayTracer.jl/blob/{commit}{path}#{line}",
    sitename = "PlantRayTracer.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "master",
        assets = String[]),
    pages = [
        "Home" => "index.md",
    ])

deploydocs(;
    repo = "github.com/VirtualPlantLab/PlantRayTracer.jl.git",
    devbranch = "master")
