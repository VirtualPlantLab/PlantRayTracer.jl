using VPLRayTracer
using Test
import Aqua

@testset "VPLRayTracer.jl" begin

    # Aqua
    @testset "Aqua" begin
        Aqua.test_all(VPLRayTracer, ambiguities = false)
        Aqua.test_ambiguities([VPLRayTracer])
    end


    include("elements_raytracer.jl")
    #include("scenes_raytracer.jl")
end
