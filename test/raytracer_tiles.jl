import PlantGraphs as PG
import PlantGeomPrimitives as PGP
import PlantGeomTurtle as PGT
import PlantRayTracer as PRT
#import GLMakie
#import PlantViz as PV
using Test
import ColorTypes: RGB
import StaticArrays: @SVector

# Simple graph that creates tiles
module Tiles
using PlantGraphs
using PlantRayTracer
struct Tile{N, M} <: Node
    length::Float64
    mat::Vector{Lambertian{M}}
end
end
import .Tiles

let

    # Power and absorbed radiance are extracted with a query
    function get_power(graph, N, M)
        tiles = PG.apply(graph, PG.Query(Tiles.Tile{N, M}))
        power = sum(sum(sum(mat.power) for mat in tile.mat)
                    for tile in tiles)
        area = length(tiles) * tiles[1].length^2
        return power, power ./ area
    end

    # Auxilliary function for relative differences
    rel(x, y) = abs(x - y) / x

    ##################### Common elements ####################

    # Length of unit tile
    L = 0.5

    # Construct the tile with N triangles
    function PGT.feed!(turtle::PGT.Turtle, tile::Tiles.Tile{N, M}, data) where {N, M}
        for _ in 1:N
            col = rand(RGB)
            PGP.Rectangle!(turtle, length = tile.length / N, width = tile.length, move = true,
                materials = [tile.mat[1] for _ in 1:2], colors = [col for _ in 1:2])
        end
        return nothing
    end

    # Common settings for the ray tracer
    settings = PRT.RTSettings(pkill = 0.3, maxiter = 3, nx = 3, ny = 3)
    psettings = PRT.RTSettings(pkill = 0.3, maxiter = 3, nx = 3, ny = 3, parallel = true)

    # Create the ray tracing scene and run it
    function ray_trace!(scene, settings, acceleration = PRT.Naive; radiosity = 1.0,
                        nrays = 1000, θ = 0.0)
        source = PRT.DirectionalSource(scene, θ = θ, Φ = 0.0, radiosity = radiosity, nrays = nrays)
        if acceleration == PRT.Naive
            rtobj = PRT.RayTracer(scene, source, settings = settings, acceleration = acceleration)
        else
            rtobj = PRT.RayTracer(scene, source, settings = settings, acceleration = PRT.BVH,
                                  rule = PRT.SAH{1}(2, 5))
        end
        nrays = PRT.trace!(rtobj)
        return nothing
    end

    # Test results
    function test_results(irradiance, target = radiosity * 0.7)
        @test rel(irradiance, target) < 1e-2
        return nothing
    end

    ##################### One tile with two triangles and single wavelength ####################

    # Construct the scene
    N = 1
    nw = 1
    radiosity = 1.0
    nrays = 1_000
    axiom = PGT.RA(90.0) + Tiles.Tile{N, nw}(L, [PRT.Lambertian(τ = 0.1, ρ = 0.2)])
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Mesh(graph)
    #render(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Source is at an angle
    θ = π/3
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity*cos(θ), nrays = nrays, θ = θ)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7 * cos(θ))

    # Tile is also at an angle (because of cloner only a fraction is exposed!)
    N = 1
    nw = 1
    radiosity = 1.0
    nrays = 1_000
    axiom = PGT.RA(90.0) + PGT.RH(60.0) + Tiles.Tile{N, nw}(L, [PRT.Lambertian(τ = 1e-6, ρ = 1e-6)])
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Mesh(graph)
    #render(scene)
    θ = π/3
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity*cos(θ), nrays = nrays, θ = θ)
    pow, irradiance = get_power(graph, N, nw)
    d = 0.5
    sin_alpha = sin(pi/3)
    h = d*sin_alpha
    dp = d - h*h/d
    test_results(irradiance, sum(radiosity) * dp/d)

    ##################### One tile with multiple triangles and single wavelength ####################
    N = 10
    nw = 1
    radiosity = 1.0
    nrays = 1_000
    axiom = PGT.RA(90.0) + Tiles.Tile{N, nw}(L, [PRT.Lambertian(τ = 0.1, ρ = 0.2)])
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Mesh(graph)
    #render_scene(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Source is at an angle
    θ = π/3
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity*cos(θ), nrays = nrays, θ = θ)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7 * cos(θ))

    ##################### One tile with multiple triangles and multiple wavelength ####################
    N = 10
    nw = 2
    radiosity = @SVector [0.5, 0.5]
    nrays = 1_000
    axiom = PGT.RA(90.0) +
            Tiles.Tile{N, nw}(L, [PRT.Lambertian(τ = (0.1, 0.1), ρ = (0.2, 0.2))])
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Mesh(graph)
    #render_scene(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Source is at an angle
    θ = π/3
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity*cos(θ), nrays = nrays, θ = θ)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7 * cos(θ))

    ##################### Many tiles with multiple triangles and single wavelength ####################
    N = 10
    nw = 1
    radiosity = 1.0
    nrays = 1_000
    make_tile() = Tiles.Tile{N, nw}(L, [PRT.Lambertian(τ = 0.1, ρ = 0.2)])
    axiom = PGT.RA(90.0) + make_tile() + make_tile()
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Mesh(graph)
    #render_scene(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Source is at an angle
    θ = π/3
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity*cos(θ), nrays = nrays, θ = θ)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7 * cos(θ))


    ##################### Many tiles with multiple triangles and multiple wavelengths ####################
    N = 10
    nw = 2
    radiosity = @SVector [0.5, 0.5]
    nrays = 1_000
    make_tile() = Tiles.Tile{N, nw}(L, [PRT.Lambertian(τ = (0.1, 0.1), ρ = (0.2, 0.2))])
    axiom = PGT.RA(90.0) + make_tile() + make_tile()
    graph = PG.Graph(axiom = axiom)
    scene = PGP.Mesh(graph)
    #render_scene(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Source is at an angle
    θ = π/3
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity*cos(θ), nrays = nrays, θ = θ)
    pow, irradiance = get_power(graph, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7 * cos(θ))


    ##################### Two graphs - many tiles with multiple triangles and multiple wavelengths ####################
    N = 10
    nw = 2
    radiosity = @SVector [0.5, 0.5]
    nrays = 1_000_000
    make_tile() = Tiles.Tile{N, nw}(L, [PRT.Lambertian(τ = (0.1, 0.1), ρ = (0.2, 0.2))])
    axiom1 = PGT.RA(90.0) + make_tile() + make_tile()
    graph1 = PG.Graph(axiom = axiom1)
    axiom2 = PGT.RA(90.0) + PGT.F(2.0) + make_tile() + make_tile()
    graph2 = PG.Graph(axiom = axiom2)
    scene = PGP.Mesh([graph1, graph2])
    # #render_scene(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Source is at an angle
    θ = π/3
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity*cos(θ), nrays = nrays, θ = θ)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7 * cos(θ))
    pow, irradiance = get_power(graph2, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7 * cos(θ))


    ##################### Two graphs - many tiles with multiple triangles and multiple wavelengths - use add! ####################
    N = 10
    nw = 2
    radiosity = @SVector [0.5, 0.5]
    nrays = 500_000
    make_tile() = Tiles.Tile{N, nw}(L, [PRT.Lambertian(τ = (0.1, 0.1), ρ = (0.2, 0.2))])
    axiom1 = PGT.RA(90.0) + make_tile() + make_tile()
    graph1 = PG.Graph(axiom = axiom1)
    axiom2 = PGT.RA(90.0) + PGT.F(2.0) + make_tile() + make_tile()
    graph2 = PG.Graph(axiom = axiom2)
    scene = PGP.Mesh(graph1)
    scene_extra = PGP.Mesh(graph2)
    mat = PRT.Lambertian(τ = (0.1, 0.1), ρ = (0.2, 0.2))
    PGP.add!(scene, scene_extra, materials = mat, colors = rand(RGB))
    #render_scene(scene)

    # Naive + Serial
    ray_trace!(scene, settings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Naive + Parallel
    ray_trace!(scene, psettings, PRT.Naive, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Serial
    ray_trace!(scene, settings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra)
    test_results(irradiance, sum(radiosity) * 0.7)

    # BVH + Parallel
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity, nrays = nrays)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7)
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra)
    test_results(irradiance, sum(radiosity) * 0.7)

    # Source is at an angle
    θ = π/3
    ray_trace!(scene, psettings, PRT.BVH, radiosity = radiosity*cos(θ), nrays = nrays, θ = θ)
    pow, irradiance = get_power(graph1, N, nw)
    test_results(irradiance, sum(radiosity) * 0.7 * cos(θ))
    pow = sum(mat.power)
    irradiance = pow / PGP.area(scene_extra)
    test_results(irradiance, sum(radiosity) * 0.7 * cos(θ))
end
