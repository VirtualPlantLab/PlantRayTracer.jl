module PlantRayTracer

# Dependencies from Julia's standard library
using LinearAlgebra
import Statistics: quantile, mean
import Base: intersect, length
import Random
import Base.Threads: @threads, nthreads
import Atomix: @atomic

# External dependencies
import StaticArrays: SVector,
    SArray, SizedVector, SMatrix, @SVector, MVector, @MVector, SDiagonal
import Unrolled: @unroll
import CoordinateTransformations: compose, Translation, LinearMap, AbstractAffineMap
import Rotations: RotX, RotY, RotZ
import StatsBase: sample, Weights
import ColorTypes: RGBA

# VPL dependencies
import PlantGeomPrimitives: Vec, O, X, Y, Z, Mesh, areas,
    ntriangles, Ellipse, rotate!, translate!, BBox,
    Scene, vertices, mesh, faces, material_ids, materials,
    Material, add!

# Raytracing API
export RayTracer, RTSettings, trace!, Naive, accelerate, Directional,
    Source, LambertianSource, DirectionalSource, PointSource, LineSource, AreaSource,
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

end
