using Test
using LinearAlgebra
import PlantRayTracer as RT
import PlantGeomPrimitives as PGP

@testset "rotate_coordinates" begin

    FT = Float64
    tol = 1e-10

    # Auxilliary function to compute the cosine of the angle between PGP.

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
        # If X points East (α=π/2), then geographic North is +Y.
        # Sun on horizon at Φ=0 (geographic North) should come from +Y, rays toward -Y.
        axes = RT.rotate_coordinates(FT(π/2), FT(0), FT(π/2))
        @test all(abs.(axes.z .- PGP.Vec(0, -1, 0)) .< tol)
    end

    # --- Tilted surface (alpha_soil, beta_soil) ---

    @testset "Flat surface (alpha_soil=0) is independent of beta_soil" begin
        # When alpha_soil=0 the slope term vanishes; beta_soil should have no effect
        axes_ref = RT.rotate_coordinates(FT(π/4), FT(π/3))
        for beta in [FT(0), FT(π/2), FT(π), FT(3π/2)]
            axes = RT.rotate_coordinates(FT(π/4), FT(π/3), FT(π), FT(0), beta)
            @test all(abs.(axes_ref.z .- axes.z) .< tol)
        end
    end

    @testset "South-facing 45° slope, sun overhead" begin
        # Slope normal tilted 45° toward South (+X). Overhead sun (θ=0) hits at 45°
        # from the slope normal → equal South (+X) and downward (-Z) components in SCS.
        axes = RT.rotate_coordinates(FT(0), FT(0), FT(π), FT(π/4), FT(π))
        @test all(abs.(axes.z .- PGP.Vec(1/sqrt(2), 0, -1/sqrt(2))) .< tol)
    end

    @testset "Sun perpendicular to south-facing 45° slope" begin
        # Sun at zenith θ=π/4 from South (Φ=π) faces the 45° south-facing slope directly.
        # Ray should be straight into the slope normal: z = (0, 0, -1) in SCS.
        axes = RT.rotate_coordinates(FT(π/4), FT(π), FT(π), FT(π/4), FT(π))
        @test all(abs.(axes.z .- PGP.Vec(0, 0, -1)) .< tol)
    end

    @testset "East-horizon sun on NS-tilted slope equals flat result" begin
        # East direction is perpendicular to the North-South tilt axis, so a slope
        # tilted toward South does not change how an East-horizon ray looks in SCS.
        axes_flat   = RT.rotate_coordinates(FT(π/2), FT(π/2))
        axes_tilted = RT.rotate_coordinates(FT(π/2), FT(π/2), FT(π), FT(π/4), FT(π))
        @test all(abs.(axes_flat.z .- axes_tilted.z) .< tol)
    end

    @testset "North-facing 45° slope, sun overhead" begin
        # Slope normal tilted 45° toward North (-X). Overhead sun hits from the South
        # side → equal North (-X) and downward (-Z) components in SCS.
        axes = RT.rotate_coordinates(FT(0), FT(0), FT(π), FT(π/4), FT(0))
        @test all(abs.(axes.z .- PGP.Vec(-1/sqrt(2), 0, -1/sqrt(2))) .< tol)
    end

    @testset "Orthonormality for tilted surfaces" begin
        for (θ, Φ, α, α_s, β_s) in [
            (π/4, π/4,  π,   π/6, π),
            (π/3, π/6,  π/2, π/4, π/2),
            (π/6, 2π/3, π,   π/3, 3π/4),
        ]
            axes = RT.rotate_coordinates(FT(θ), FT(Φ), FT(α), FT(α_s), FT(β_s))
            @test abs(norm(axes.x) - 1) < tol
            @test abs(norm(axes.y) - 1) < tol
            @test abs(norm(axes.z) - 1) < tol
            @test abs(axes.x ⋅ axes.y) < tol
            @test abs(axes.x ⋅ axes.z) < tol
            @test abs(axes.y ⋅ axes.z) < tol
        end
    end

    @testset "α rotation combined with tilted surface" begin
        # X points East (α=π/2) on a south-facing 45° slope. Overhead sun should
        # have ray components along -Y (South in this frame) and -Z in SCS.
        axes = RT.rotate_coordinates(FT(0), FT(0), FT(π/2), FT(π/4), FT(π))
        @test all(abs.(axes.z .- PGP.Vec(0, -1/sqrt(2), -1/sqrt(2))) .< tol)
    end

end
