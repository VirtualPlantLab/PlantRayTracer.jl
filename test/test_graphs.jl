using PlantGraphs
using PlantGeomPrimitives
import PlantGeomPrimitives as Geom
using PlantGeomTurtle
using PlantRayTracer
import PlantRayTracer as RT
using Test
import ColorTypes: RGB

##### Test turtle-based construction of a scene #####

# Modules needed to test ray tracing of graphs
module sn
   using PlantGraphs
   using PlantRayTracer
   struct E64 <: Node
       length::Float64
       mat::Black{1}
   end
   struct E32 <: Node
       length::Float32
       mat::Black{1}
   end
   struct E2 <: Node
       length::Float64
       mat::Vector{Black{1}}
   end
end
import .sn

module btree
    import PlantGraphs as G
    using PlantRayTracer
    # Meristem
    struct Meristem <: G.Node end
    # Node
    struct Node <: G.Node end
    # Internode
    mutable struct Internode <: G.Node
        length::Float64
        mat::Lambertian{1}
    end
    # Graph-level variables
    struct treeparams
        growth::Float64
    end
end
import .btree

let

    # Koch curve @ 64 bits
    L = 1.0
    axiom = sn.E64(L, Black(1)) + RU(120.0) + sn.E64(L, Black(1)) +
            RU(120.0) + sn.E64(L, Black(1))
    function Kochsnowflake(x)
        L = data(x).length
        sn.E64(L / 3, Black(1)) + RU(-60.0) + sn.E64(L / 3, Black(1)) +
        RU(120.0) + sn.E64(L / 3, Black(1)) + RU(-60.0) +
        sn.E64(L / 3, Black(1))
    end
    function PlantGeomTurtle.feed!(turtle::Turtle, e::sn.E64, data)
        HollowCylinder!(turtle, length = e.length, width = e.length / 10,
            height = e.length / 10, move = true, materials = e.mat,
            colors = rand(RGB))
        return nothing
    end
    rule = Rule(sn.E64, rhs = Kochsnowflake)
    Koch = Graph(axiom = axiom, rules = Tuple(rule))
    mesh = Mesh(Koch)
    @test ntriangles(mesh) == 120
    @test length(materials(mesh)) == ntriangles(mesh)
    @test eltype(materials(mesh)) == Black{1}

    ##### Test intersection of specific rays with Naive acc #####

    # Intersection of a rectangle from a directional light source downwards
    nrays = 100_000
    radiosity = 1.0
    mesh = Rectangle(length = 1.0, width = 1.0)
    rotatey!(mesh, -π / 2) # To put it in the XY plane
    material = Black()
    add_property!(mesh, :materials, material)
    #mesh = Mesh(mesh = rect, materials = material)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(gbox, θ = 0.0, Φ = 0.0, radiosity = radiosity, nrays = nrays)
    settings = RTSettings(pkill = 1.0, maxiter = 1)
    rtobj = RayTracer(mesh, source, settings = settings, acceleration = Naive)
    nrays, _ = trace!(rtobj)

    @test nrays == nrays
    pow_abs = material.power[1]
    pow_gen = source.power[1] * source.nrays
    @test pow_abs ≈ pow_gen

    ##### Test intersection code of specific rays with Naive acc + grid cloner #####

    ##### Using black materials #####

    # Intersection of a rectangle from a directional light source downwards
    nrays = 100_000
    radiosity = 1.0
    mesh = Rectangle(length = 2.0, width = 1.0)
    rotatey!(mesh, -π / 2) # To put it in the XY plane
    material = Black()
    add_property!(mesh, :materials, material)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(gbox, θ = 0.0, Φ = 0.0, radiosity = radiosity, nrays = nrays)
    settings = RTSettings(pkill = 1.0, maxiter = 1, nx = 3, ny = 3)
    rtobj = RayTracer(mesh, source, settings = settings, acceleration = Naive)
    nrays_traced, _ = trace!(rtobj)
    @test nrays_traced == nrays
    pow_abs = material.power[1]
    pow_gen = source.power[1] * source.nrays
    @test pow_abs ≈ pow_gen

    # Intersection of a rectangle from a directional light source that is horizontal
    nrays = 100_000
    radiosity = 1.0
    mesh = Rectangle(length = 2.0, width = 1.0)
    rotatey!(mesh, -π / 2) # To put it in the XY plane
    material = Black()
    add_property!(mesh, :materials, material)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(gbox,
        θ = π / 2,
        Φ = 0.0,
        radiosity = radiosity,
        nrays = nrays)
    settings = RTSettings(pkill = 1.0, maxiter = 1, nx = 3, ny = 3, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, source, settings = settings, acceleration = Naive)
    nrays_traced, _ = trace!(rtobj)
    @test nrays_traced == nrays
    pow_abs = material.power[1]
    pow_gen = source.power[1] * source.nrays
    @test pow_abs == 0.0

    # Intersection of a rectangle from a directional light source that is at an angle
    nrays = 100_000
    radiosity = 1.0
    mesh = Rectangle(length = 2.0, width = 1.0)
    rotatey!(mesh, -π / 2) # To put it in the XY plane
    material = Black()
    add_property!(mesh, :materials, material)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(gbox,
        θ = π / 4,
        Φ = 0.0,
        radiosity = radiosity*cos(π / 4),
        nrays = nrays)
    settings = RTSettings(pkill = 1.0, maxiter = 1, nx = 3, ny = 3, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, [source], settings = settings, acceleration = Naive)
    nrays_traced, _ = trace!(rtobj)
    @test nrays_traced == nrays
    pow_abs = material.power[1]
    pow_gen = source.power[1] * source.nrays
    @test pow_abs ≈ pow_gen
    @test pow_abs / area(mesh) ≈ radiosity*cos(π/4)

    ##### Using sensors #####

    # Intersection of a rectangle from a directional light source downwards
    nrays = 100_000
    radiosity = 1.0
    rect1, rect2, rect3 = collect(begin
        r = Rectangle(length = 1.5, width = 1.0)
        rotatey!(r, -π / 2) # To put it in the XY plane
        r
      end
    for i in 1:3)
    translate!(rect2, Z())
    translate!(rect3, 2.0 * Z())
    rectangles = Geom.Mesh([rect1, rect2, rect3])
    mats = [Sensor() for _ in 1:3]
    ext_mats = mats[[1,1,2,2,3,3]]
    mesh = add_property!(rectangles, :materials,ext_mats)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(mesh,
        θ = 0.0,
        Φ = 0.0,
        radiosity = radiosity,
        nrays = nrays)

    # Need to make sure maxiter > 1 or it will stop after the first sensor
    settings = RTSettings(pkill = 1.0, maxiter = 2, nx = 1, ny = 1, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, source, settings = settings, acceleration = Naive)
    nrays_traced, _ = trace!(rtobj)
    @test nrays_traced == nrays
    pow_abs = [material.power[1] for material in mats]
    pow_gen = source.power[1] * source.nrays
    @test all(pow_abs .≈ pow_gen)
    @test all(pow_abs ./ area(rect1) .≈ radiosity)

    # Intersection of a rectangle from a directional light source at an angle
    source = DirectionalSource(gbox,
        θ = π / 4,
        Φ = 0.0,
        radiosity = radiosity*cos(π / 4),
        nrays = nrays)
    settings = RTSettings(pkill = 1.0, maxiter = 3, nx = 2, ny = 2, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, source, settings = settings, acceleration = Naive)
    nrays_traced, _ = trace!(rtobj)

    @test nrays_traced == nrays
    pow_abs = [material.power[1] for material in mats]
    pow_gen = source.power[1] * source.nrays
    @test all(pow_abs .≈ pow_gen)
    @test all(pow_abs ./ area(rect1) .≈ radiosity .* cos.(π/4))

    ##### Using Lambertian #####

    # Intersection of a rectangle from a directional light source downwards
    nrays = 1_000_000
    radiosity = 1.0
    rect1, rect2, rect3 = collect(begin
        r = Rectangle(length = 1.0, width = 1.5)
        rotatey!(r, -π / 2) # To put it in the XY plane
        r
    end for i in 1:3)
    translate!(rect2, Z())
    translate!(rect3, 2.0 * Z())
    rectangles = Geom.Mesh([rect1, rect2, rect3])
    mats = [Lambertian(τ = 0.3, ρ = 0.3) for i in 1:3]
    ext_mats = mats[[1, 1, 2, 2, 3, 3]]
    mesh = add_property!(rectangles, :materials,ext_mats)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(gbox, θ = 0.0, Φ = 0.0, radiosity = radiosity, nrays = nrays)
    # Need to make sure maxiter > 1 or it will stop after the first sensor
    settings = RTSettings(pkill = 0.9, maxiter = 4, nx = 1, ny = 1, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, source, settings = settings, acceleration = Naive)
    nrays_traced, _ = trace!(rtobj)

    @test nrays_traced > nrays
    pow_abs = [material.power[1] for material in mats]
    pow_gen = source.power[1] * source.nrays
    @test sum(pow_abs) < pow_gen
    @test pow_abs[1] < pow_abs[2] < pow_abs[3]

    # Intersection of a rectangle from a directional light source at an angle
    source = DirectionalSource(gbox,
        θ = π / 4,
        Φ = 0.0,
        radiosity = radiosity*cos(π / 4),
        nrays = nrays)
    settings = RTSettings(pkill = 0.9, maxiter = 4, nx = 2, ny = 2, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, source, settings = settings, acceleration = Naive)
    nrays_traced, _ = trace!(rtobj)

    @test nrays_traced > nrays
    RTirrs = [mats[i].power[1] / area(rect1) for i in 1:3]
    @test all(RTirrs .< [cos(π / 4) for i in 1:3])
    @test RTirrs[1] < RTirrs[2] < RTirrs[3]
    RTirrs_naive = deepcopy(RTirrs) # for comparison with BVH later

    ##### Test intersection of specific rays with BVH acc #####

    # Intersection of a rectangle from a directional light source downwards
    nrays = 100_000
    radiosity = 1.0
    rect = Rectangle(length = 1.0, width = 1.5)
    rotatey!(rect, -π / 2) # To put it in the XY plane
    mat = Black()
    ex_mat = [mat,mat]
    mesh = add_property!(rect, :materials, ex_mat)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(mesh,
        θ = 0.0,
        Φ = 0.0,
        radiosity = radiosity,
        nrays = nrays)
    settings = RTSettings(pkill = 1.0, maxiter = 1)
    rtobj = RayTracer(mesh, source, settings = settings, acceleration = BVH,
        rule = SAH{1}(2, 5))
    nrays, _ = trace!(rtobj)

    @test nrays == nrays
    RTirr = mat.power[1] ./ area(rect)
    @test RTirr ≈ radiosity

    ##### Test intersection code of specific rays with BVH acc + grid cloner #####

    # Intersection of a rectangle from a directional light source downwards (black)
    nrays = 100_000
    radiosity = 1.0
    rect = Rectangle(length = 3.0, width = 1.0)
    rotatey!(rect, -π / 2) # To put it in the XY plane
    mat = Black()
    ex_mat = [mat,mat]
    mesh = add_property!(rect, :materials, ex_mat)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(mesh,
        θ = 0.0,
        Φ = 0.0,
        radiosity = radiosity,
        nrays = nrays)
    settings = RTSettings(pkill = 1.0, maxiter = 1, nx = 3, ny = 3, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, [source], settings = settings, acceleration = BVH,
        rule = SAH{1}(2, 5))
    nrays_traced, _ = trace!(rtobj)
    @test nrays_traced == nrays
    RTirr = mat.power[1] / area(rect)
    @test RTirr ≈ radiosity

    # Intersection of a rectangle from a directional light source downwards (sensor)
    nrays = 100_000
    radiosity = 1.0
    rect1, rect2, rect3 = collect(begin
        r = Rectangle(length = 1.5, width = 1.0)
        rotatey!(r, -π / 2) # To put it in the XY plane
        r
    end
                                  for i in 1:3)
    translate!(rect2, Z())
    translate!(rect3, 2.0 * Z())
    mesh = Geom.Mesh([rect1, rect2, rect3])
    mats = [Sensor() for i in 1:3]
    ext_mats = mats[[1, 1, 2, 2, 3, 3]]
    add_property!(mesh, :materials, ext_mats)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(gbox, θ = 0.0, Φ = 0.0, radiosity = radiosity, nrays = nrays)
    # Need to make sure maxiter > 1 or it will stop after the first sensor
    settings = RTSettings(pkill = 1.0, maxiter = 3, nx = 1, ny = 1, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, [source], settings = settings, acceleration = BVH,
        rule = SAH{1}(2, 5))
    nrays_traced, _ = trace!(rtobj)
    @test nrays_traced == nrays
    RTirrs = [mats[i].power[1] / area(rect1) for i in 1:3]
    @test RTirrs ≈ [1.0 for i in 1:3]

    # Intersection of a rectangle from a directional light source at an angle (sensor)
    source = DirectionalSource(gbox, θ = π / 4, Φ = 0.0, radiosity = cos(π / 4), nrays = nrays)
    settings = RTSettings(pkill = 1.0, maxiter = 3, nx = 2, ny = 2, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, [source], settings = settings, acceleration = BVH,
        rule = SAH{1}(2, 5))
    nrays_traced, _ = trace!(rtobj)

    @test nrays_traced == nrays
    RTirrs = [mats[i].power[1] / area(rect1) for i in 1:3]
    @test RTirrs ≈ [cos(π / 4) for i in 1:3]

    # Intersection of a rectangle from a directional light source at an angle (Lambertian)
    nrays = 1_000_000
    radiosity = 1.0
    rect1, rect2, rect3 = collect(begin
        r = Rectangle(length = 2.0, width = 1.0)
        rotatey!(r, -π / 2) # To put it in the XY plane
        r
    end
                                  for i in 1:3)
    translate!(rect2, Z())
    translate!(rect3, 2.0 * Z())
    mesh = Geom.Mesh([rect1, rect2, rect3])
    mats = [Lambertian(τ = 0.3, ρ = 0.3) for i in 1:3]
    ext_mats = mats[[1, 1, 2, 2, 3, 3]]
    add_property!(mesh, :materials, ext_mats)
    gbox = RT.AABB(mesh)
    source = DirectionalSource(gbox,
        θ = π / 4,
        Φ = 0.0,
        radiosity = radiosity*cos(π / 4),
        nrays = nrays)
    settings = RTSettings(pkill = 0.9, maxiter = 4, nx = 4, ny = 4, dx = 1.0, dy = 1.0)
    rtobj = RayTracer(mesh, [source], settings = settings, acceleration = BVH,
        rule = SAH{1}(2, 5))
    nrays_traced, _ = trace!(rtobj)

    @test nrays_traced > nrays
    RTirrs = [mats[i].power[1] / area(rect1) for i in 1:3]
    @test all(RTirrs .< [cos(π / 4) for i in 1:3])
    @test RTirrs[1] < RTirrs[2] < RTirrs[3]

    # Should yield the same results as the naive acceleration structure!
    @test all(abs.(RTirrs .- RTirrs_naive) .< 0.02)

    ##### Ray trace binary tree #####

    function PlantGeomTurtle.feed!(turtle::Turtle, i::btree.Internode, vars)
        HollowCube!(turtle, length = i.length, height = i.length / 10,
            width = i.length / 10, move = true, colors = RGB(0, 1, 0),
            materials = i.mat)
        return nothing
    end
    rule = Rule(btree.Meristem,
        rhs = mer -> btree.Node() +
                     (RU(-60.0) + btree.Internode(0.1, Lambertian(τ = 0.0, ρ = 0.3)) +
                      RH(90.0) + btree.Meristem(),
            RU(60.0) + btree.Internode(0.1, Lambertian(τ = 0.0, ρ = 0.3)) +
            RH(90.0) + btree.Meristem()))
    axiom = btree.Internode(0.1, Lambertian(τ = 0.0, ρ = 0.3)) + btree.Meristem()
    tree = Graph(axiom, Tuple(rule), btree.treeparams(0.5))
    getInternode = Query(btree.Internode)
    function elongate!(tree, query)
        for x in apply(tree, query)
            x.length = x.length * (1.0 + data(tree).growth)
        end
    end
    function growth!(tree, query)
        elongate!(tree, query)
        rewrite!(tree)
    end
    function simulate(tree, query, nsteps)
        new_tree = deepcopy(tree)
        for i in 1:nsteps
            growth!(new_tree, query)
        end
        return new_tree
    end
    function getpower(tree, query)
        powers = Float64[]
        for x in apply(tree, query)
            push!(powers, x.mat.power[1])
        end
        return powers
    end

    newtree = simulate(tree, getInternode, 4)

    # Ray trace the tree with a single directional light source
    nrays = 1_000_000
    mesh = Mesh(newtree)
    source = DirectionalSource(mesh, θ = π / 4, Φ = 0.0, radiosity = cos(π / 4), nrays = nrays)
    # Tracing with BVH acceleration structure
    settings = RTSettings(pkill = 0.9, maxiter = 4, nx = 5, ny = 5, dx = 1.0,
        dy = 1.0, parallel = true)
    rtobj = RayTracer(mesh, source, settings = settings, acceleration = BVH,
        rule = SAH{6}(5, 10))
    nrays_traced, _ = trace!(rtobj)
    powers_bvh = getpower(newtree, getInternode)
    # Tracing with Naive acceleration structure
    settings = RTSettings(pkill = 0.9, maxiter = 4, nx = 5, ny = 5, dx = 1.0,
        dy = 1.0, parallel = true)
    rtobj = RayTracer(mesh, [source], settings = settings, acceleration = Naive)
    nrays_traced, _ = trace!(rtobj)
    powers_naive = getpower(newtree, getInternode)

    # For low number of rays the results are the same for the Naive and BVH
    # acceleration structures, but for large number of rays the results diverge
    # This divergence does not seem to decrease with the number of rays. It may
    # depend on the mesh itself though?
    @test maximum(abs.((powers_bvh .- powers_naive) ./ (powers_naive .+ eps(Float64)))) <
          0.008
    @test abs(sum(powers_bvh) - sum(powers_naive)) / sum(powers_bvh) < 6e-5

    # Simple test to make sure that rays are always generated from above the mesh
    r = Rectangle(length = 2.0, width = 1.0)
    rotatey!(r, -π / 2) # To put it in the XY plane
    translate!(r, Vec(0.0, 0.5, 0.0))
    r2 = deepcopy(r)
    translate!(r2, Vec(0.0, 0.0, -1.0))
    mats = [Black(), Black()]
    ex_mats = mats[[1, 1, 2, 2]]
    mesh = Geom.Mesh([r, r2])
    Geom.add_property!(mesh, :materials, ex_mats)
    sources = DirectionalSource(mesh,
        θ = π / 2 * 0.99,
        Φ = 0.0,
        radiosity = cos(π / 2 * 0.99),
        nrays = nrays)
    power_out = sources.power * sources.nrays
    settings = RTSettings(nx = 15, ny = 15, dx = 2.0, dy = 1.0, parallel = true)
    rtobj = RayTracer(mesh, sources, settings = settings)
    nrays_traced, _ = trace!(rtobj)
    @test mats[1].power[1] / power_out[1] ≈ 1.0
    @test mats[2].power[1] / power_out[1] ≈ 0.0

    # Simple test of having multiple materials per geometry
    r = Rectangle(length = 2.0, width = 1.0)
    rotatey!(r, -π / 2) # To put it in the XY plane
    translate!(r, Vec(1.0, 0.5, 0.0))
    L = 1.0
    m() = [Black(1) for i in 1:8] # Number of triangles in the hollow cube
    axiom = sn.E2(L, m()) + RU(120.0) + sn.E2(L, m()) +
            RU(120.0) + sn.E2(L, m())
    function Kochsnowflake(x)
        L = data(x).length
        sn.E2(L / 3, m()) + RU(-60.0) + sn.E2(L / 3, m()) +
        RU(120.0) + sn.E2(L / 3, m()) + RU(-60.0) +
        sn.E2(L / 3, m())
    end
    function PlantGeomTurtle.feed!(turtle::Turtle, e::sn.E2, vars)
        if turtle.message == "raytracer"
            HollowCube!(turtle, length = e.length, width = e.length / 10,
                height = e.length / 10, move = true, materials = e.mat)
        elseif turtle.message == "render"
            powers = getindex.(power.(e.mat), 1)
            color = any(powers .> 0.0) ? RGB(1.0, 0.0, 0.0) : RGB(0.0, 0.0, 0.0)
            HollowCube!(turtle, length = e.length, width = e.length / 10,
                height = e.length / 10, move = true,
                colors = color)
        end
        return nothing
    end
    rule = Rule(sn.E2, rhs = Kochsnowflake)
    Koch = Graph(axiom = axiom, rules = Tuple(rule))
    mesh = Mesh(Koch, message = "raytracer")
    Geom.add!(mesh, r, materials = Black(1), colors = RGB(0.5, 0.5, 0.0))
    sources = DirectionalSource(mesh, θ = π / 4, Φ = 0.0, radiosity = 1.0, nrays = nrays)
    settings = RTSettings(parallel = true)
    rtobj = RayTracer(mesh, sources, settings = settings)
    nrays_traced, _ = trace!(rtobj)
    @test length(filter(x -> x.power[1] > 0.0, materials(mesh))) == 6 # only 4 faces are seen (+ soil)
    mesh = Mesh(Koch, message = "render")
    Geom.add!(mesh, r, materials = Black(1), colors = RGB(0.5, 0.5, 0.0))

    mesh = Mesh(Koch, message = "raytracer")
    Geom.add!(mesh, r, materials = Black(1), colors = RGB(0.5, 0.5, 0.0))
    sources = DirectionalSource(mesh, θ = π / 4, Φ = π / 2, radiosity = cos(π / 4), nrays = nrays)
    settings = RTSettings(parallel = true)
    rtobj = RayTracer(mesh, sources, settings = settings)
    nrays_traced, _ = trace!(rtobj)
    @test length(filter(x -> x.power[1] > 0.0, materials(mesh))) == 30 # 8 faces seen (+ soil)
    mesh = Mesh(Koch, message = "render")
    Geom.add!(mesh, r, materials = Black(1), colors = RGB(0.5, 0.5, 0.0))

    # Check that we can processes an array of meshs properly
    mesh = Mesh([Koch, deepcopy(Koch)], message = "raytracer")

end
