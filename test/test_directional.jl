using Test
using LinearAlgebra
import PlantRayTracer as RT
import PlantGeomPrimitives as PGP

@testset "rotate_coordinates" begin

    FT = Float64
    tol = 1e-10

    @testset "Sun directly overhead (θ=0)" begin
        # At zenith = 0, azimuth doesn't matter - sun is straight up
        # The z axis of the local system should point down (-Z world)
        axes = RT.rotate_coordinates(FT(0), FT(0))
        @test all(abs.(axes.z .- PGP.Vec(0, 0, -1)) .< tol)
    end

    @testset "Sun on horizon from North (θ=π/2, Φ=0)" begin
        # Sun on horizon, azimuth=0 means North
        # The z axis of the new coordinated system (pointing from sun) should be along X
        axes = RT.rotate_coordinates(FT(π/2), FT(0))
        @test all(abs.(axes.z .- PGP.Vec(1, 0.0, 0)) .< tol)
    end

    @testset "Sun on horizon from East (θ=π/2, Φ=π/2)" begin
        # Sun on horizon, azimuth=π/2 means East (Y axis direction)
        axes = RT.rotate_coordinates(FT(π/2), FT(π/2))
        @test all(abs.(axes.z .- PGP.Vec(0, -1, 0)) .< tol)
    end

    @testset "Axes are orthonormal" begin
        # For arbitrary angles, the three axes should be orthogonal unit vectors
        for (θ, Φ) in [(π/4, π/4), (π/3, π/6), (π/6, 2π/3)]
            axes = RT.rotate_coordinates(FT(θ), FT(Φ))
            @test abs(norm(axes.x) - 1) < tol
            @test abs(norm(axes.y) - 1) < tol
            @test abs(norm(axes.z) - 1) < tol
            @test abs(axes.x ⋅ axes.y) < tol
            @test abs(axes.x ⋅ axes.z) < tol
            @test abs(axes.y ⋅ axes.z) < tol
        end
    end

    @testset "Orientation angle α=pi matches original function" begin
        for (θ, Φ) in [(π/4, π/4), (π/3, π/6), (0.0, 0.0)]
            axes_orig = RT.rotate_coordinates(FT(θ), FT(Φ))
            axes_new  = RT.rotate_coordinates(FT(θ), FT(Φ), FT(pi))
            @test all(abs.(axes_orig.x .- axes_new.x) .< tol)
            @test all(abs.(axes_orig.y .- axes_new.y) .< tol)
            @test all(abs.(axes_orig.z .- axes_new.z) .< tol)
        end
    end

    @testset "X points East (α=π/2): sun on horizon from North becomes East" begin
        # If X points East (α=π/2), then geographic North is now -Y.
        # Sun on horizon at Φ=0 (geographic North) should come from -Y direction.
        axes = RT.rotate_coordinates(FT(π/2), FT(0), FT(π/2))
        @test all(abs.(axes.z .- PGP.Vec(0, -1, 0)) .< tol)
    end

end
