### This file contains public API ###
# RTSettings
# RayTracer
# get_nw
# trace!
# accelerate

"""
    RTSettings(;verbose = true, parallel = false, pkill = 0.2, maxiter = 2,
                sampler = Random.Xoshiro(123456789), nx = 3, ny = 3, nz = 0, dx = 0.0,
                dy = 0.0, dz = 0.0)

Settings for the ray tracer:

- `verbose` indicates if the ray tracer will print warnings and other information.

- `parallel` indicates if the raytracer will run on a single core or make use of multiple
cores in the machine based on Julia's multithreading support. `pkill` is the probably that a
ray is terminated by the Russian roulette after it has been scattered a `maxiter` number of
times (this excludes interactions with materials of type `Sensor`).

- `sampler` is the pseudo-random number generator to be used by the ray tracer.

- `nx` and `ny` are the number of times the mesh will be cloned by the grid cloner in each
direction along the x and y axis (e.g., setting `nx = 1` and `ny = 1` will generate a grid
of 3 x 3 clones of the original mesh), whereas `dx` and `dy` will be distance at which each
new clone will be generated (along the axis).

- `nz` is the number of times the mesh will be cloned in the vertical direction. Unlike
horizontal cloning, the vertical cloning is always done in the positive direction of `z`
axis and the number of clones will be exactly `nz`.

See VPL documentation for more details on the ray  tracer.

## Examples
```jldoctest
julia> RTSettings(parallel = true, maxiter = 3);
```
"""
Base.@kwdef struct RTSettings{RNG, dt <: Union{Real, Nothing}}
    verbose::Bool = true
    parallel::Bool = false
    pkill::Float64 = 0.2
    maxiter::Int = 2
    sampler::RNG = Random.Xoshiro(123456789)
    nx::Int = 3
    ny::Int = 3
    nz::Int = 0
    dx::dt = nothing
    dy::dt = nothing
    dz::dt = dy isa Nothing ? nothing : zero(dy)
end

# Wrap the acceleration structure and the grid cloner in a single object
struct AccMesh{A, G}
    acc::A
    grid::G
end

"""
    accelerate(mesh::Mesh; settings = RTSettings(), acceleration = Naive, rule = nothing)

Create an `AccMesh` object from a mesh, settings and acceleration function
(choose from `Naive` or  `BVH`). The argument `rule` is only required for the
accelerator `BVH` and it must be an object of type `SAH` or `AvgSplit`
(it is ignored for the `Naive` accelerator). The `AccMesh` object contains
the acceleration structure and the grid cloner structure built on top of the
original 3D meshes in `mesh`. See VPL documentation for details.

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> mesh = Ellipse();

julia> acc = accelerate(mesh);
```
"""
function accelerate(mesh::PGP.Mesh; settings = RTSettings(), acceleration = Naive,
                    rule = nothing)
    triangles = Triangle(mesh)
    acc = acceleration(triangles, collect(1:length(triangles)), rule)
    grid = GridCloner(acc; nx = settings.nx, ny = settings.ny, dx = settings.dx,
        dy = settings.dy)
    AccMesh(acc, grid)
end

"""
    RayTracer(mesh, materials, sources, settings)

Create a ray tracer object from an acceleration structure built around a 3D mesh, a grid cloner structure around
the acceleration structure (`mesh`), a vector of materials associated to the mesh (`materials`), a vector of sources of irradiance
(`sources`) and settings. (as generated by `RTSettings()`). See VPL documentation for more details on the ray tracer.
"""
struct RayTracer{A, M, S, RT}
    mesh::A
    materials::M
    sources::S
    settings::RT
end

"""
    RayTracer(mesh, sources; settings = RTSettings(), acceleration = Naive, rule = nothing)

Create a `RayTracer` object from a mesh, a tuple of sources (objects that inherit from `Source`),
or a single source, settings and acceleration function (choose from `Naive` or  `BVH`). The argument `rule` is only required for the
accelerator `BVH` and it must be an object of type `SAH` or `AvgSplit` (it is ignored for the `Naive` accelerator).

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> mesh = Ellipse();

julia> mat = Lambertian(τ = 0.1, ρ = 0.2);

julia> add_property!(mesh, :materials, [mat for _ in 1:ntriangles(mesh)]);

julia> source = DirectionalSource(mesh, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);

julia> rt = RayTracer(mesh, source);

julia> sources = (DirectionalSource(mesh, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000),
                  DirectionalSource(mesh, θ = 45.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000));

julia> rt = RayTracer(mesh, sources);
```
"""
function RayTracer(mesh::PGP.Mesh, sources; settings = RTSettings(),
    acceleration = Naive, rule = nothing)
    # Construct the acceleration structure and grid cloner around the mesh
    acc_mesh = accelerate(mesh, settings = settings, acceleration = acceleration,
        rule = rule)
    # Check for use directional light sources in the absence of the grid cloner
    acc_mesh.grid.nleaves == 1 && any(source.geom isa Directional for source in sources) &&
    settings.verbose &&
        println("Using Directional sources in the absence of a grid cloner may lead to incorrect results. See VPL documentation for details.")
    RayTracer(acc_mesh, materials(mesh), sources, settings)
end

function RayTracer(mesh::PGP.Mesh, source::Source; kwargs...)
    RayTracer(mesh, (source,); kwargs...)
end

###################################
########## Shooting rays ##########
###################################

# Retrieve how many rays are going to be generated from each source
# in order to facilitate choosing the right source later on
function distribute_rays(sources)
    @inbounds begin
        nsources = length(sources)
        if nsources > 1
            nrays = cumsum(source.nrays for source in sources)
        else
            nrays = [sources[1].nrays]
        end
        return nrays, nrays[end]
    end
end

# Choose the correct light source and generate a ray from it knowing
# how many rays need to be shot from each source
function generate_ray!(sources, nrays, iray, power, rng)
    @inbounds begin
        nsources = length(sources)
        if nsources > 1
            block = 1
            for i in 1:nsources
                nrays[i] > iray && (block = i; break)
            end
            ray = shoot_ray!(sources[block], power, rng)
        else
            ray = shoot_ray!(sources[1], power, rng)
        end
        return ray
    end
end

#########################################################
###### Recursive ray tracing with Russian roulette ######
#########################################################

"""
    get_nw(rt::RayTracer)

Retrieve the number of wavelengths being simulated by the ray tracer.
See VPL documentation for more details on the ray tracer.

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> mesh = Ellipse();

julia> mat = Lambertian(τ = 0.1, ρ = 0.2);

julia> add_property!(mesh, :materials, [mat for _ in 1:ntriangles(mesh)]);

julia> source = DirectionalSource(mesh, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);

julia> rt = RayTracer(mesh, source);

julia> get_nw(rt)
1
```
"""
function get_nw(rt::RayTracer)
    @inbounds nw = get_nw(rt.sources[1])
    return nw
end

# Random termination of rays after maximum iteration with asymptotic bias correction
function roulette!(power, settings, iteration, rng)
    all(power .== 0.0) && return true
    if iteration >= settings.maxiter
        terminate = rand(rng) < settings.pkill
        if terminate
            return true
        else
            power .= power ./ (1.0 .- settings.pkill)
            return false
        end
    else
        return false
    end
end

"""
    trace!(rt)

Run the ray tracing simulations. This function will overwrite the `power` component of any material object
that is included in the mesh. It returns the total number of rays being traced (primary and secondary) without
including interactions with `Sensor` objects (first value returned) or including interactions with `Sensor`
objects (second value returned).
See VPL documentation for more details on the ray tracer.

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> mesh = Ellipse();

julia> mat = Lambertian(τ = 0.1, ρ = 0.2);

julia> add_property!(mesh, :materials, [mat for _ in 1:ntriangles(mesh)]);

julia> source = DirectionalSource(mesh, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);

julia> rt = RayTracer(mesh, source);

julia> trace!(rt);
```
"""
function trace!(rt::RayTracer)
    # Reset all materials
    reset!(rt.materials)
    # Distribute primary rays across sources (nrays is per source, total_rays is aggregated)
    nrays, total_rays = distribute_rays(rt.sources)
    # Run ray tracer across multiple threads
    if rt.settings.parallel
        # Each thread gets a deep copy of materials to avoid data races
        #materials = [deepcopy(rt.materials) for _ in 1:T.nthreads()]
        # Each thread gets a different pseudo-random number generator which seed depends on the original one
        samplers = [deepcopy(rt.settings.sampler) for _ in 1:T.nthreads()]
        for i in 1:T.nthreads()
            Random.seed!(samplers[i], rand(rt.settings.sampler, UInt))
        end
        # Allocated primary rays evenly across threads
        Δrays = div(total_rays, T.nthreads())
        irays = collect(0:(T.nthreads() - 1)) .* Δrays .+ 1
        erays = collect(1:T.nthreads()) .* Δrays
        erays[end] = total_rays
        # Loop over threads while updating materials and collecting total number of rays traced
        net_rays_thread = zeros(Int, T.nthreads())
        total_rays_thread = zeros(Int, T.nthreads())
        T.@threads for i in 1:T.nthreads()
            @inbounds net_rays_thread[i], total_rays_thread[i]  = trace_thread!(rt,
                rt.materials,
                nrays,
                irays[i],
                erays[i],
                samplers[i])
        end
        net_rays = sum(net_rays_thread)
        total_rays = sum(total_rays_thread)
        # Run ray tracer in a single multiple thread
    else
        net_rays, total_rays = trace_thread!(rt,
            rt.materials,
            nrays,
            1,
            total_rays,
            rt.settings.sampler)
    end
    return net_rays, total_rays
end

# Trace a subset of the rays on a given PC thread
function trace_thread!(rt, materials, nrays, i_ray, e_ray, rng)
    # Get the current floating point precision out of the ray tracer
    FT = eltype(rt.mesh.acc.gbox.min)
    # Number of wavelengths
    nw = get_nw(rt)
    # Keep track of all the rays being traced (including secondary)
    net_rays_traced = 0
    total_rays_traced = 0
    # Setup nodestack and power payload and trace all primary rays
    nodestack = Int[]
    tnodestack = Int[]
    dstack = FT[]
    tdstack = FT[]
    power = SA.MVector{nw, Float64}(Tuple(0.0 for _ in 1:nw))
    for i in i_ray:e_ray
        r = generate_ray!(rt.sources, nrays, i, power, rng)
        temp = trace!(r,
            rt,
            materials,
            tnodestack,
            tdstack,
            nodestack,
            dstack,
            power,
            rng)
        net_rays_traced += temp[1]
        total_rays_traced += temp[2]
    end
    return net_rays_traced, total_rays_traced
end

# Trace a ray through a mesh with Russian roulette
function trace!(r::Ray{FT}, rt::RayTracer, materials, tnodestack::Vector{Int},
    tdstack::Vector{FT},
    nodestack::Vector{Int}, dstack::Vector{FT}, power, rng) where {FT}
    # Recursive Monte Carlo ray tracing with Russian roulette criterion
    iteration = 1
    nrays_traced = 1
    while true
        # Check intersection against mesh
        hit, intersection, disp = intersect(r,
            rt.mesh.grid,
            rt.mesh.acc,
            tnodestack,
            tdstack,
            nodestack,
            dstack)
        !hit && (return (iteration, nrays_traced))
        material = materials[intersection.id]
        # Interaction with surface material
        interaction = calculate_interaction(material, power, r, intersection, rng)
        absorb_power!(material, power, interaction)
        # Russian roulette
        roulette!(power, rt.settings, iteration, rng) && (return (iteration, nrays_traced))
        # Generate new ray
        r = generate_ray(material, r, disp, intersection, interaction, rng)
        # Increase iteration counter (unless it is a sensor-like material)
        if interaction.mode != :sensor
            iteration += 1
        end
        # Do keep track of all the rays tracer (because of sensors)
        nrays_traced += 1
    end
    return iteration, nrays_traced
end
