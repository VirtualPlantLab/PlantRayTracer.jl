import PlantRayTracer as RT
import PlantRayTracer
using PlantRayTracer
import PlantGeomPrimitives as PGP
using Test
import LinearAlgebra: ⋅, norm
import Random
using StaticArrays
import ColorTypes: RGB, RGBA
import StatsBase

let

    ##### Test Ray, Triangle and AABB #####

    # Simple ray construction (check FT)
    o64 = PGP.Z(Float64)
    dir64 = .-PGP.Z(Float64)
    r64 = RT.Ray(o64, dir64)
    @test r64.o == o64
    @test r64.dir == dir64

    o32 = PGP.Z(Float32)
    dir32 = .-PGP.Z(Float32)
    r32 = RT.Ray(o32, dir32)
    @test r32.o == o32
    @test r32.dir == dir32

    # Simple Triangle construction (check FT)
    FT = Float64
    t64 = RT.Triangle(PGP.O(FT), PGP.X(FT), PGP.Y(FT))
    @test t64.p  == PGP.O(FT)
    @test t64.e1 == PGP.X(FT)
    @test t64.e2 == PGP.Y(FT)
    @test t64.n  == -PGP.Z(FT)

    FT = Float32
    t32 = RT.Triangle(PGP.O(FT), PGP.X(FT), PGP.Y(FT))
    @test t32.p  == PGP.O(FT)
    @test t32.e1 == PGP.X(FT)
    @test t32.e2 == PGP.Y(FT)
    @test t32.n  == -PGP.Z(FT)

    # Simple AABB construction (check FT)
    FT = Float64
    t64 = RT.Triangle(PGP.Z(FT), PGP.X(FT), PGP.Y(FT))
    aabb64 = RT.AABB(t64)
    @test all(aabb64.min .≈ zeros(FT, 3))
    @test all(aabb64.max .≈ ones(FT, 3))

    FT = Float32
    t32 = RT.Triangle(PGP.Z(FT), PGP.X(FT), PGP.Y(FT))
    aabb32 = RT.AABB(t32)
    @test all(aabb32.min .≈ zeros(FT, 3))
    @test all(aabb32.max .≈ ones(FT, 3))

    # Ray triangle intersection - 1
    r1 = RT.Ray(PGP.Vec(0.1, 0.1, 1), .-PGP.Z())
    t1 = RT.Triangle(PGP.Vec(0, 0, 0.5), PGP.Vec(1, 0, 0.5), PGP.Vec(0, 1, 0.5))
    i1 = intersect(r1, t1)
    @test i1[1]
    @test i1[2] ≈ 0.5
    @test !i1[3]

    # Ray triangle intersection - 2
    r2 = RT.Ray(PGP.Vec(0.1, 0.1, 1), .-PGP.Z())
    t2 = RT.Triangle(PGP.Vec(0, 0, 0.5), PGP.Vec(0, 1, 0.5), PGP.Vec(1, 0, 0.5))
    i2 = intersect(r2, t2)
    @test i2[1]
    @test i2[2] ≈ 0.5
    @test i2[3]

    # Ray triangle intersection - 3
    r3 = RT.Ray(PGP.Vec(0.0, 0, 1), .-PGP.Z())
    t3 = RT.Triangle(PGP.Vec(0, 0, 0.5), PGP.Vec(0, 1, 0.5), PGP.Vec(1, 0, 0.5))
    i3 = intersect(r3, t3)
    @test !i3[1]

    # Ray triangle intersection - 4
    r4 = RT.Ray(PGP.Vec(-0.1, -0.1, 1), .-PGP.Z())
    t4 = RT.Triangle(PGP.Vec(0, 0, 0.5), PGP.Vec(0, 1, 0.5), PGP.Vec(1, 0, 0.5))
    i4 = intersect(r4, t4)
    @test !i4[1]

    # Ray AABB intersection - 1
    r1 = RT.Ray(PGP.Vec(0.1, 0.1, 1), .-PGP.Z())
    a1 = RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(1, 1, 1.0))
    i1 = intersect(r1, a1)
    @test i1[1]
    @test i1[2] ≈ 0.0

    # Ray AABB intersection - 2
    r2 = RT.Ray(PGP.Vec(0.1, 0.1, 2), .-PGP.Z())
    a2 = RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(1, 1, 1.0))
    i2 = intersect(r2, a2)
    @test i2[1]
    @test i2[2] ≈ 1.0

    # Ray AABB intersection - 3
    r3 = RT.Ray(PGP.Vec(0.1, 0.1, 0.5), .-PGP.Z())
    a3 = RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(1, 1, 1.0))
    i3 = intersect(r3, a3)
    @test i3[1]
    @test i3[2] ≈ -0.5

    # Ray AABB intersection - 4
    r4 = RT.Ray(PGP.Vec(2.0, 0.1, 0.5), .-PGP.Z())
    a4 = RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(1, 1, 1.0))
    i4 = intersect(r4, a4)
    @test !i4[1]

    # Area of a triangle
    @test RT.area(t4) ≈ 0.5

    # Axes of a triangle
    axs = RT.axes(t4)
    @test axs[1] ⋅ axs[2] ≈ 0.0
    @test axs[1] ⋅ axs[3] ≈ 0.0
    @test axs[2] ⋅ axs[3] ≈ 0.0
    @test all(norm.(axs) .≈ [1, 1, 1.0])

    # Generate a random point from within a triangle
    rng = Random.MersenneTwister(123456789)
    p = RT.generate_point(t4, rng)
    @test p isa PGP.Vec

    # Area of an AABB
    a1 = RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(1, 1, 1.0))
    a2 = RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(2, 2, 2.0))
    @test RT.area(a1) ≈ 6.0
    @test RT.area(a2) ≈ 24.0

    # Center of an AABB
    @test RT.center(a1) ≈ PGP.Vec(0.5, 0.5, 0.5)
    @test RT.center(a2) ≈ PGP.Vec(1, 1, 1.0)

    # Longest axis of an AABB
    @test RT.longest(RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(3, 2, 1.0))) == 1
    @test RT.longest(RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(2, 3, 1.0))) == 2
    @test RT.longest(RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(1, 2, 3.0))) == 3

    # Height of an AABB
    @test RT.height(RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(3, 2, 1.0))) ≈ 1.0
    @test RT.height(RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(1, 2, 3.0))) ≈ 3.0

    # Center point on the top face of an AABB
    @test RT.topcenter(a1) ≈ PGP.Vec(0.5, 0.5, 1)

    # Extract 8 vertices from an AABB
    verts = PGP.vertices(a1)
    @test length(verts) == 8
    @test verts[1] ≈ PGP.Vec(0, 0, 0)

    # Test for an empty AABB (basically a point)
    a0 = RT.AABB(PGP.Vec(0.0, 0, 0), PGP.Vec(0, 0, 0.0))
    isempty(a0)

    # AABB around vector of triangles
    ta = RT.Triangle(PGP.Vec(0, 0, 0.5), PGP.Vec(0, 1, 0.5), PGP.Vec(1, 0, 0.5))
    tb = RT.Triangle(PGP.Vec(1, 0, 0.5), PGP.Vec(0, 1, 0.5), PGP.Vec(-1, 0, 0.5))
    aabb = RT.AABB([ta, tb])
    @test aabb.min ≈ PGP.Vec(-1.0, 0.0, 0.5)
    @test aabb.max ≈ PGP.Vec(1.0, 1.0, 0.5)

    # Merge a list of AABBs
    aabb1 = RT.AABB([ta, tb])
    aabb2 = RT.AABB([RT.AABB(ta), RT.AABB(tb)])
    @test aabb1 == aabb2

    ##### Test light source geometry #####

    # Point source
    source64 = RT.PointSource(PGP.Vec(0.0, 0.0, 1.0))
    @test source64 isa RT.SourceGeometry
    @test source64.loc == PGP.Vec(0.0, 0.0, 1.0)
    @test eltype(source64.loc) == Float64

    source32 = RT.PointSource(PGP.Vec(0.0f0, 0.0f0, 1.0f0))
    @test source32 isa RT.SourceGeometry
    @test source32.loc == PGP.Vec(0.0f0, 0.0f0, 1.0f0)
    @test eltype(source32.loc) == Float32

    rng = Random.MersenneTwister(123456789)
    @test RT.generate_point(source32, rng) == source32.loc

    # Line source
    source64 = RT.LineSource(PGP.Vec(0.0, 0.0, 1.0), PGP.X(Float64))
    @test source64 isa RT.SourceGeometry
    @test source64.p == PGP.Vec(0.0, 0.0, 1.0)
    @test source64.line == PGP.X(Float64)
    @test eltype(source64.p) == Float64

    source32 = RT.LineSource(PGP.Vec(0.0f0, 0.0f0, 1.0f0), PGP.X(Float32))
    @test source32 isa RT.SourceGeometry
    @test source32.p == PGP.Vec(0.0f0, 0.0f0, 1.0f0)
    @test source32.line == PGP.X(Float32)
    @test eltype(source32.p) == Float32

    rng = Random.MersenneTwister(123456789)
    for i in 1:10
        local p = RT.generate_point(source32, rng)
        dist = p .- source32.p
        norm(dist) <= norm(source32.line)
        @test dist ⋅ source32.line ≈ norm(dist) * norm(source32.line)
    end

    # Area source
    r = PGP.Rectangle(length = 2.0, width = 2.0)
    source = RT.AreaSource(r)
    source.tvec isa Vector
    @test first(source.tvec) isa RT.Triangle
    @test source.areas isa StatsBase.Weights
    @test sum(source.areas) ≈ 2 * 2

    for i in 1:10
        local p = RT.generate_point(source, rng)
        @test p[1] ≈ 0.0
        @test -1.0 <= p[2] <= 1.0
        @test 0.0 <= p[3] <= 2.0
    end

    r = PGP.Rectangle(length = 3.0, width = 2.5)
    source = RT.AreaSource(r)
    @test sum(source.areas) ≈ 3 * 2.5
    for i in 1:10
        local p = RT.generate_point(source, rng)
        @test p[1] ≈ 0.0
        @test -1.25 <= p[2] <= 1.25
        @test 0.0 <= p[3] <= 3.0
    end

    ##### Test light source angle #####

    # Fixed angle source
    source = RT.FixedSource(0.0, 0.0)
    @test source.dir ≈ PGP.Vec(0, 0, -1.0)
    @test RT.generate_direction(source, rng) == source.dir

    # Lambertian source
    source1 = RT.LambertianSource(PGP.X(), PGP.Y(), PGP.Z())
    source2 = RT.LambertianSource((PGP.X(), PGP.Y(), PGP.Z()))
    @test source1 == source2

    for i in 1:10
        local p = RT.generate_direction(source1, rng)
        @test p[3] > 0
    end

    ##### Test light source - Point & Fixed #####
    geom = RT.PointSource(PGP.Vec(0.0, 0.0, 1.0))
    angle = RT.FixedSource(0.0, 0.0)

    # Single wavelength
    source = RT.Source(geom, angle, 1.0, 1_000)
    @test source isa RT.Source
    @test RT.get_nw(source) == 1
    @test RT.generate_point(source, rng) == geom.loc
    @test RT.generate_direction(source, rng) == angle.dir
    source_power = [0.0]
    ray = RT.shoot_ray!(source, source_power, rng)
    @test source_power == [1.0]
    @test ray.o == geom.loc
    @test ray.dir == angle.dir

    # Multiple wavelengths
    source = RT.Source(geom, angle, (1.0, 2.0), 1_000)
    @test source isa RT.Source
    @test RT.get_nw(source) == 2
    @test RT.generate_point(source, rng) == geom.loc
    @test RT.generate_direction(source, rng) == angle.dir
    source_power = [0.0, 1.0]
    ray = RT.shoot_ray!(source, source_power, rng)
    @test source_power == [1.0, 2.0]
    @test ray.o == geom.loc
    @test ray.dir == angle.dir

    ##### Test directional light source #####
    gbox = RT.AABB(PGP.O(), PGP.Vec(1.0, 1.0, 1.0))
    dsource = RT.DirectionalSource(gbox, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000)
    @test dsource isa RT.Source
    @test dsource.geom isa RT.Directional

    ##### Test materials #####

    # Lambertian material
    mat = RT.Lambertian(τ = 0.0, ρ = 0.3)
    @test mat isa Material
    @test length(power(mat)) == 1
    @test all(power(mat) .== zeros(1))
    source_power = [1.0]
    mode, coef = RT.choose_outcome(mat, source_power, rng)
    @test coef == [0.3]
    @test mode == :ρ
    mat.power[1] = 8.0
    RT.reset!(mat)
    @test all(power(mat) .== [0.0])

    mat = RT.Lambertian(τ = 0.3, ρ = 0.0)
    @test mat isa Material
    @test length(power(mat)) == 1
    source_power = [1.0]
    mode, coef = RT.choose_outcome(mat, source_power, rng)
    @test coef == [0.3]
    @test mode == :τ

    mat = RT.Lambertian(τ = 0.3, ρ = 0.7)
    @test mat isa Material
    @test length(power(mat)) == 1
    source_power = [1.0]
    mode, coef = RT.choose_outcome(mat, source_power, rng)
    @test coef == [1.0]
    @test mode == :ρ || mode == :τ

    mat = RT.Lambertian(τ = (0.0, 0.0), ρ = (0.3, 0.3))
    @test mat isa Material
    @test length(power(mat)) == 2
    @test all(power(mat) .== zeros(2))
    source_power = [1.0, 1.0]
    mode, coef = RT.choose_outcome(mat, source_power, rng)
    @test coef == [0.3, 0.3]
    @test mode == :ρ
    mat.power[1] = 8.0
    mat.power[2] = 2.0
    RT.reset!(mat)
    @test all(power(mat) .== zeros(2))

    mat = RT.Lambertian(τ = (0.3, 0.3), ρ = (0.0, 0.0))
    @test mat isa Material
    @test length(power(mat)) == 2
    source_power = [1.0, 1.0]
    mode, coef = RT.choose_outcome(mat, source_power, rng)
    @test coef == [0.3, 0.3]
    @test mode == :τ

    mat = RT.Lambertian(τ = (0.3, 0.7), ρ = (0.7, 0.3))
    @test mat isa Material
    @test length(power(mat)) == 2
    source_power = [1.0, 1.0]
    mode, coef = RT.choose_outcome(mat, source_power, rng)
    @test coef == 2 * [0.3, 0.7] || coef == 2 * [0.7, 0.3]
    @test mode == :ρ || mode == :τ

    # Phong material
    mat = RT.Phong(τ = 0.2, ρd = 0.3, ρsmax = 0.7, n = 2)
    @test mat isa Material
    mat = RT.Phong(τ = (0.2, 0.2), ρd = (0.3, 0.3), ρsmax = (0.7, 0.7), n = 2)
    @test length(mat.τ) == 2

    # Black material
    mat = RT.Black(1)
    @test all(power(mat) .== 0.0)
    mat = RT.Black(3)
    @test all(power(mat) .== zeros(3))

    # Sensor material
    mat = RT.Sensor(1)
    @test all(power(mat) .== 0.0)
    mat = RT.Sensor(3)
    @test all(power(mat) .== zeros(3))

    # Two-sided Sensor material
    mat = RT.TwoSidedSensor(1)
    @test all(power(mat) .== 0.0)
    @test all(power(mat, front = false) .== 0.0)
    mat = RT.TwoSidedSensor(3)
    @test all(power(mat) .== zeros(3))
    @test all(power(mat, front = false) .== zeros(3))

    ##### Triangle conversion #####
    # Make sure normal vectors are not inverted!
    r = PGP.Rectangle()
    ts = RT.Triangle(r)
    @test all(PGP.normals(r)[1] .== ts[1].n)
    @test all(PGP.normals(r)[2] .== ts[2].n)
end
