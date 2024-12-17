import PlantGraphs as PG
import PlantGeomPrimitives as PGP
import PlantGeomTurtle as PGT
import PlantRayTracer as PRT
#import GLMakie
#import PlantViz as PV
using Test
import ColorTypes: RGB, RGBA, Colorant
import StaticArrays: @SVector

# Simple graph that creates tiles
module Tiles
    using PlantGraphs
    using PlantRayTracer
    @kwdef struct Shade <: Node
        height::Float64 = 1.0
        mat::Lambertian{1} = Lambertian(τ = 0.1, ρ = 0.0)
    end
    @kwdef struct Tile <: Node
        height::Float64 = 1.0
        mat::Black{1} = Black(1)
    end
    @kwdef struct LightSensor <: Node
        height::Float64 = 1.0
        mat::Sensor{1} = Sensor(1)
    end
end

let
import .Tiles


function PGT.feed!(turtle::PGT.Turtle, tile::Tiles.Tile, data)
    PGT.t!(turtle, to = PGP.Vec(0.0, 0.0, tile.height))
    PGP.Rectangle!(turtle, length = 1.0, width = 1.0,
                   materials = tile.mat, colors = rand(RGB))
    return nothing
end
function PGT.feed!(turtle::PGT.Turtle, shade::Tiles.Shade, data)
    PGT.t!(turtle, to = PGP.Vec(0.0, 0.0, shade.height))
    PGP.Rectangle!(turtle, length = 1.0, width = 1.0,
                   materials = shade.mat, colors = rand(RGBA))
    return nothing
end
function PGT.feed!(turtle::PGT.Turtle, sensor::Tiles.LightSensor, data)
    PGT.t!(turtle, to = PGP.Vec(0.0, 0.0, sensor.height))
    PGP.Rectangle!(turtle, length = 1.0, width = 1.0,
                   materials = sensor.mat, colors = rand(RGBA))
    return nothing
end

# Create the shade cloth and soil tile
axiom = PGT.RA(90.0)
for i in 1:2:20
    axiom += Tiles.Shade(height = Float64(i)) + Tiles.LightSensor(height = Float64(i + 1))
end
g     = PG.Graph(axiom = axiom);
scene = PGP.Mesh(g);
#PV.render(scene, normals = true)

# Run ray tracer
settings = PRT.RTSettings(pkill = 0.9, maxiter = 10, nx = 10, ny = 10, parallel = true)
source   = PRT.DirectionalSource(scene, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000)
rtobj    = PRT.RayTracer(scene, source, settings = settings, acceleration = PRT.BVH,
                         rule = PRT.SAH{1}(2, 10));
nrays, trays = PRT.trace!(rtobj)

# Extract the sensors
sensors = PG.apply(g, PG.Query(Tiles.LightSensor))
@test all(Int.(round.(log10.([PRT.power(sensor.mat)[1] for sensor in sensors]))) .== -9:0)
shades  = PG.apply(g, PG.Query(Tiles.Shade))
@test all(Int.(round.(log10.([PRT.power(shade.mat)[1] for shade in shades]))) .== -9:0)

end
