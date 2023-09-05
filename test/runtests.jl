using PlantRayTracer
using Test
import Aqua

@testset "PlantRayTracer.jl" begin

    # Aqua
    @testset "Aqua" begin
        Aqua.test_all(PlantRayTracer, ambiguities = false)
        Aqua.test_ambiguities([PlantRayTracer])
    end

    include("elements_raytracer.jl")
    #include("scenes_raytracer.jl")
end
