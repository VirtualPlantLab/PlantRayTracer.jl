module PlantRayTracer

using PrecompileTools

# Dependencies from Julia's standard library
using LinearAlgebra
import Statistics
import Base
import Random
import Base.Threads as T
import Atomix

# External dependencies
import StaticArrays as SA
import Unrolled: @unroll
import CoordinateTransformations as CT
import Rotations
import StatsBase

# VPL dependencies
import PlantGeomPrimitives as PG
import PlantGeomPrimitives as PGP

# Raytracing API
export RayTracer, RTSettings, trace!, Naive, accelerate, Directional, materials, Material,
    Source, LambertianSource, DirectionalSource, PointSource, LineSource, AreaSource, TwoSidedSensor,
    tau, rho, Lambertian, Phong, Sensor, Black,
    get_nw, FixedSource, reset!, power, BVH, SAH, AvgSplit

# Helpers and auxilliary functions
include("utils.jl")

# Basic geometry
include("Geometry/Ray.jl")
include("Geometry/Triangle.jl")
include("Geometry/AABB.jl")

# Materials
include("Materials/Material.jl")

# Acceleration structures
include("Acceleration/Acceleration.jl")

# Sources
include("Sources/Source.jl")

# RayTracer
include("RayTracer/RayTracer.jl")

# @compile_workload begin

#     rng = Random.MersenneTwister(123456789)

#     # Simple ray construction (check FT)
#     for FT in (Float32, Float64)
#         o = O(FT)
#         dir = .-PGP.Z(FT)
#         r = Ray(o, dir)
#     end
#     # Simple Triangle construction (check FT)
#     for FT in (Float32, Float64)
#         t = Triangle(O(FT), X(FT), Y(FT))
#         a = area(t)
#         axs = axes(t)
#         p = generate_point(t, rng)
#     end

#     # Simple AABB construction (check FT)
#     for FT in (Float32, Float64)
#         t = Triangle(PGP.Z(FT), X(FT), Y(FT))
#         a = AABB(t)
#         longest(a)
#         base_area(a)
#         height(a)
#         topcenter(a)
#         vertices(a)
#         isempty(a)
#     end

#     # Ray triangle intersection
#     for FT in (Float32, Float64)
#         r = Ray(PGP.Vec(FT(0.1), FT(0.1), FT(1)), .-PGP.Z(FT))
#         t = Triangle(PGP.Vec(FT(0), FT(0), FT(0.5)), PGP.Vec(FT(1), FT(0), FT(0.5)), PGP.Vec(FT(0), FT(1), FT(0.5)))
#         i = intersect(r, t)
#     end

#     # Ray AABB intersection
#     for FT in (Float32, Float64)
#         r = Ray(PGP.Vec(FT(0.1), FT(0.1), FT(1)), .-PGP.Z(FT))
#         a = AABB(PGP.Vec(FT(0), FT(0), FT(0)), PGP.Vec(FT(1), FT(1), FT(1)))
#         i = intersect(r, a)
#     end

#     # AABB around vector of triangles
#     for FT in (Float32, Float64)
#         ta = Triangle(PGP.Vec(FT(0), FT(0), FT(0.5)), PGP.Vec(FT(1), FT(0), FT(0.5)), PGP.Vec(FT(0), FT(1), FT(0.5)))
#         tb = Triangle(PGP.Vec(FT(1), FT(0), FT(0.5)), PGP.Vec(FT(0), FT(1), FT(0.5)), PGP.Vec(FT(-1), FT(0), FT(0.5)))
#         a = AABB([ta, tb])
#         aabb = AABB(AABB(ta), AABB(tb))
#     end

#     # Point source
#     for FT in (Float32, Float64)
#         source = PointSource(PGP.Vec(FT(0), FT(0), FT(1)))
#         generate_point(source, rng)
#     end

#     # Line source
#     for FT in (Float32, Float64)
#         source = LineSource(PGP.Vec(FT(0), FT(0), FT(1)), X(FT))
#         p = generate_point(source, rng)
#     end

#     # Area source
#     for FT in (Float32, Float64)
#         r = PG.Rectangle(length = FT(2), width = FT(2))
#         source = AreaSource(r)
#         p = generate_point(source, rng)
#     end

#     # Fixed angle source
#     for FT in (Float32, Float64)
#         source = FixedSource(FT(0), FT(0))
#         generate_direction(source, rng)
#     end

#     # Lambertian source
#     for FT in (Float32, Float64)
#         source = LambertianSource(X(FT), Y(FT), PGP.Z(FT))
#         generate_direction(source, rng)
#     end

#     # Point + fixed light source up to 5 wavelengths
#     for FT in (Float32, Float64)
#         geom = PointSource(PGP.Vec(FT(0), FT(0), FT(1)))
#         angle = FixedSource(FT(0), FT(0))
#         for i in 1:5
#             source = Source(geom, angle, Tuple(ones(i)), 1_000)
#             generate_point(source, rng)
#             generate_direction(source, rng)
#             source_power = ones(i)
#             ray = shoot_ray!(source, source_power, rng)
#         end
#     end
#     # Point + fixed light source up to 5 wavelengths
#     for FT in (Float32, Float64)
#         for i in 1:5
#             gbox = AABB(O(FT), PGP.Vec(FT(1), FT(1), FT(1)))
#             rad = SA.SVector{i, Float64}(ones(i)...)
#             source = DirectionalSource(gbox, θ = FT(0), Φ = FT(0), radiosity = rad, nrays = 1_000)
#             generate_point(source, rng)
#             generate_direction(source, rng)
#             source_power = ones(i)
#             ray = shoot_ray!(source, source_power, rng)
#         end
#     end

#     ##### Test materials #####

#     # Lambertian material up to 5 wavelengths
#     for i in 1:5
#         mat = Lambertian(τ = tau(zeros(i)...), ρ = rho((0.3 .+ zeros(i))...))
#         pow = power(mat)
#         source_power = ones(i)
#         mode, coef = choose_outcome(mat, source_power, rng)
#         reset!(mat)
#     end

#     # Phong material up to 5 wavelengths
#     for i in 1:5
#         mat = Phong(τ = tau((0.2 .+ zeros(i))...), ρd =  rho((0.3 .+ zeros(i))...),
#                     ρsmax =  rho((0.7 .+ zeros(i))...), n =  2)
#         pow = power(mat)
#         source_power = ones(i)
#         mode, coef = choose_outcome(mat, 0.2, source_power, rng)
#         reset!(mat)
#     end

#     # Black material up to 5 wavelengths
#     for i in 1:5
#         mat = Black(i)
#         pow = power(mat)
#     end

#     # Sensor up to 5 wavelengths
#     for i in 1:5
#         mat = Sensor(i)
#         pow = power(mat)
#     end
# end

end
