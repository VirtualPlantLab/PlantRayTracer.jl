using Test
using LinearAlgebra
import PlantRayTracer as RT
import PlantGeomPrimitives as PGP

@testset "rotate_coordinates" begin

    FT = Float64
    tol = 1e-10

    # Auxilliary function to compute the cosine of the angle of incidence
    function check_cos(dir, θ::FT, Φ::FT, alpha_soil = zero(FT), beta_soil = FT(180)) where {FT}
        actual_cos = dot(dir, PGP.Z())
        expected_cos = cosd(alpha_soil)*cosd(θ) + sind(alpha_soil)*sind(θ)*cosd(Φ - beta_soil)
        abs(actual_cos + expected_cos)
    end

    @testset "Sun directly overhead (θ=0)" begin
        # At zenith = 0, azimuth doesn't matter - sun is straight up
        # The z axis of the local system should point down (-Z world)
        axes = RT.rotate_coordinates(FT(0), FT(0))
        @test all(abs.(axes.z .- PGP.Vec(0, 0, -1)) .< tol)
        @test check_cos(axes.z, FT(0), FT(0)) < tol
    end

    @testset "Sun on horizon from North (θ=90, Φ=0)" begin
        # Sun on horizon, azimuth=0 means North
        # The z axis of the new coordinated system (pointing from sun) should be along X
        axes = RT.rotate_coordinates(FT(90), FT(0))
        @test all(abs.(axes.z .- PGP.Vec(1, 0.0, 0)) .< tol)
        @test check_cos(axes.z, FT(90), FT(0)) < tol
    end

    @testset "Sun on horizon from East (θ=90, Φ=90)" begin
        # Sun on horizon, azimuth=90 means East (Y axis direction)
        axes = RT.rotate_coordinates(FT(90), FT(90))
        @test all(abs.(axes.z .- PGP.Vec(0, -1, 0)) .< tol)
        @test check_cos(axes.z, FT(90), FT(90)) < tol
    end

    @testset "Axes are orthonormal" begin
        # For arbitrary angles, the three axes should be orthogonal unit vectors
        for (θ, Φ) in [(45, 45), (60, 30), (30, 120)]
            axes = RT.rotate_coordinates(FT(θ), FT(Φ))
            @test abs(norm(axes.x) - 1) < tol
            @test abs(norm(axes.y) - 1) < tol
            @test abs(norm(axes.z) - 1) < tol
            @test abs(axes.x ⋅ axes.y) < tol
            @test abs(axes.x ⋅ axes.z) < tol
            @test abs(axes.y ⋅ axes.z) < tol
        end
    end

    @testset "Orientation angle α=180 matches original function" begin
        for (θ, Φ) in [(45, 45), (60, 30), (0.0, 0.0)]
            axes_orig = RT.rotate_coordinates(FT(θ), FT(Φ))
            axes_new  = RT.rotate_coordinates(FT(θ), FT(Φ), FT(180))
            @test all(abs.(axes_orig.x .- axes_new.x) .< tol)
            @test all(abs.(axes_orig.y .- axes_new.y) .< tol)
            @test all(abs.(axes_orig.z .- axes_new.z) .< tol)
            @test check_cos(axes_new.z, FT(θ), FT(Φ)) < tol
        end
    end

    @testset "X points East (α=90): sun on horizon from North becomes East" begin
        # If X points East (α=90), then geographic North is +Y.
        # Sun on horizon at Φ=0 (geographic North) should come from +Y, rays toward -Y.
        axes = RT.rotate_coordinates(FT(90), FT(0), FT(90))
        @test all(abs.(axes.z .- PGP.Vec(0, -1, 0)) .< tol)
        @test check_cos(axes.z, FT(90), FT(0)) < tol
    end

    # --- Tilted surface (alpha_soil, beta_soil) ---

    @testset "Flat surface (alpha_soil=0) is independent of beta_soil" begin
        # When alpha_soil=0 the slope term vanishes; beta_soil should have no effect
        axes_ref = RT.rotate_coordinates(FT(45), FT(60))
        for beta in [FT(0), FT(90), FT(180), FT(270)]
            axes = RT.rotate_coordinates(FT(45), FT(60), FT(180), FT(0), beta)
            @test all(abs.(axes_ref.z .- axes.z) .< tol)
            @test check_cos(axes.z, FT(45), FT(60), FT(0), beta) < tol
        end
    end

    @testset "South-facing 45° slope, sun overhead" begin
        # Slope normal tilted 45° toward South (+X). Overhead sun (θ=0) hits at 45°
        # from the slope normal → equal South (+X) and downward (-Z) components in SCS.
        axes = RT.rotate_coordinates(FT(0), FT(0), FT(180), FT(45), FT(180))
        @test all(abs.(axes.z .- PGP.Vec(1/sqrt(2), 0, -1/sqrt(2))) .< tol)
        @test check_cos(axes.z, FT(0), FT(0), FT(45), FT(180)) < tol
    end

    @testset "Sun perpendicular to south-facing 45° slope" begin
        # Sun at zenith θ=45 from South (Φ=180) faces the 45° south-facing slope directly.
        # Ray should be straight into the slope normal: z = (0, 0, -1) in SCS.
        axes = RT.rotate_coordinates(FT(45), FT(180), FT(180), FT(45), FT(180))
        @test all(abs.(axes.z .- PGP.Vec(0, 0, -1)) .< tol)
        @test check_cos(axes.z, FT(45), FT(180), FT(45), FT(180)) < tol
    end

    @testset "East-horizon sun on NS-tilted slope equals flat result" begin
        # East direction is perpendicular to the North-South tilt axis, so a slope
        # tilted toward South does not change how an East-horizon ray looks in SCS.
        axes_flat   = RT.rotate_coordinates(FT(90), FT(90))
        axes_tilted = RT.rotate_coordinates(FT(90), FT(90), FT(180), FT(45), FT(180))
        @test all(abs.(axes_flat.z .- axes_tilted.z) .< tol)
        @test check_cos(axes_tilted.z, FT(90), FT(90), FT(45), FT(180)) < tol
    end

    @testset "North-facing 45° slope, sun overhead" begin
        # Slope normal tilted 45° toward North (-X). Overhead sun hits from the South
        # side → equal North (-X) and downward (-Z) components in SCS.
        axes = RT.rotate_coordinates(FT(0), FT(0), FT(180), FT(45), FT(0))
        @test all(abs.(axes.z .- PGP.Vec(-1/sqrt(2), 0, -1/sqrt(2))) .< tol)
        @test check_cos(axes.z, FT(0), FT(0), FT(45), FT(0)) < tol
    end

    @testset "Orthonormality for tilted surfaces" begin
        for (θ, Φ, α, α_s, β_s) in [
            (45, 45,  180,  30, 180),
            (60, 30,   90,  45,  90),
            (30, 120, 180,  60, 135),
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
        # X points East (α=90) on a south-facing 45° slope. Overhead sun should
        # have ray components along -Y (South in this frame) and -Z in SCS.
        axes = RT.rotate_coordinates(FT(0), FT(0), FT(90), FT(45), FT(180))
        @test all(abs.(axes.z .- PGP.Vec(0, -1/sqrt(2), -1/sqrt(2))) .< tol)
        @test check_cos(axes.z, FT(0), FT(0), FT(45), FT(180)) < tol
    end

    @testset "Float32 and mixed-precision inputs" begin
        tol32 = 1f-5

        # All five arguments Float32 → output is Vec{Float32}
        axes = RT.rotate_coordinates(Float32(0), Float32(0), Float32(180), Float32(0), Float32(180))
        @test eltype(axes.z) == Float32
        @test all(abs.(axes.z .- PGP.Vec(0f0, 0f0, -1f0)) .< tol32)

        axes = RT.rotate_coordinates(Float32(90), Float32(0), Float32(180), Float32(0), Float32(180))
        @test eltype(axes.z) == Float32
        @test all(abs.(axes.z .- PGP.Vec(1f0, 0f0, 0f0)) .< tol32)

        # Float32 axes are orthonormal (with relaxed tolerance)
        axes = RT.rotate_coordinates(Float32(45), Float32(60), Float32(180), Float32(30), Float32(180))
        @test eltype(axes.z) == Float32
        @test abs(norm(axes.x) - 1) < tol32
        @test abs(norm(axes.y) - 1) < tol32
        @test abs(norm(axes.z) - 1) < tol32
        @test abs(axes.x ⋅ axes.y) < tol32
        @test abs(axes.x ⋅ axes.z) < tol32
        @test abs(axes.y ⋅ axes.z) < tol32

        # Float32 θ/Φ with Float64 α or beta_soil → Float64
        axes = RT.rotate_coordinates(Float32(0), Float32(0), 180.0, Float32(0), 180.0)
        @test eltype(axes.z) == Float64

        # Float32 θ/Φ with the default alpha_soil=0.0 (Float64) → Float64
        axes = RT.rotate_coordinates(Float32(0), Float32(0))
        @test eltype(axes.z) == Float64

        # Mixed Float32/Float64 → Float64
        axes = RT.rotate_coordinates(Float32(90), 0.0, Float32(180), Float32(0), Float32(180))
        @test eltype(axes.z) == Float64
        @test all(abs.(axes.z .- PGP.Vec(1.0, 0.0, 0.0)) .< 1e-7)
    end

end
