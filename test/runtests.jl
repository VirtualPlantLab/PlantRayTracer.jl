using PlantRayTracer
using Test
using Documenter
import Aqua

# Test examples on documentation (jldoctest blocks)
DocMeta.setdocmeta!(PlantRayTracer,
    :DocTestSetup,
    :(using PlantRayTracer);
    recursive = true)
doctest(PlantRayTracer)

# Aqua
@testset "Aqua" begin
    Aqua.test_all(PlantRayTracer, ambiguities = false)
    Aqua.test_ambiguities([PlantRayTracer])
end

# Test individual elements of the ray tracer
include("elements_raytracer.jl")

# Test ray tracing of scenes
include("raytracer.jl")
