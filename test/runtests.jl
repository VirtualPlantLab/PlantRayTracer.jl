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
@testset "Elements of the ray tracer" begin
    include("elements_raytracer.jl")
end

# Test ray tracing of tiles
@testset "Ray tracing scenes" begin
    include("raytracer_tiles.jl")
end

# Test ray tracing of scenes
@testset "Ray tracing scenes" begin
    include("raytracer.jl")
end
