### This file contains public API ###
# RTSettings
# RayTracer
# get_nw
# trace!
# accelerate

"""
    RTSettings(;parallel = false, pkill = 0.2, maxiter = 2, sampler = Random.Xoshiro(123456789),
                nx = 3, ny = 3, nz = 0, dx = 0.0, dy = 0.0, dz = 0.0)

Settings for the ray tracer: `parallel` indicates if the raytracer will run on a single core or make use of
multiple cores in the machine based on Julia's multithreading support. `pkill` is the probably that a ray is
terminated by the Russian roulette after it has been scattered a `maxiter` number of times. `sampler` is the
pseudo-random number generator to be used by the ray tracer. `nx` and `ny` are the number of times the scene
will be clone by the grid cloner in each direction along the x and y axis (e.g., setting `nx = 1` and `ny = 1`
will generate a grid of 3 x 3 clones of the original scene), whereas `dx` and `dy` will be distance at which
each new clone will be generated (along the axis). `nz` is the number of times the scene will be cloned in the
vertical direction. Unlike horizontal cloning, the vertical cloning is always done in the positive direction of
`z` axis and the number of clones will be exactly `nz`. See VPL documentation for more details on the ray  tracer.

## Examples
```jldoctest
julia> RTSettings(parallel = true, maxiter = 3);
```
"""
Base.@kwdef struct RTSettings{RNG, dt <: Union{Real, Nothing}}
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
struct AccScene{A, G}
    acc::A
    grid::G
end

"""
    accelerate(scene::Scene; settings = RTSettings(), acceleration = Naive, rule = nothing)

Create an `AccScene` object from a scene, settings and acceleration function
(choose from `Naive` or  `BVH`). The argument `rule` is only required for the
accelerator `BVH` and it must be an object of type `SAH` or `AvgSplit`
(it is ignored for the `Naive` accelerator). The `AccScene` object contains
the acceleration structure and the grid cloner structure built on top of the
original 3D meshes in `scene`. See VPL documentation for details.

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> sc = Scene(mesh = Ellipse());

julia> acc = accelerate(sc);
```
"""
function accelerate(scene::Scene; settings = RTSettings(), acceleration = Naive,
    rule = nothing)
    triangles = Triangle(scene)
    acc = acceleration(triangles, material_ids(scene), rule)
    grid = GridCloner(acc; nx = settings.nx, ny = settings.ny, dx = settings.dx,
        dy = settings.dy)
    AccScene(acc, grid)
end

"""
    RayTracer(scene, materials, sources, settings)

Create a ray tracer object from an acceleration structure built around a 3D mesh, a grid cloner structure around
the acceleration structure (`scene`), a vector of materials associated to the mesh (`materials`), a vector of sources of irradiance
(`sources`) and settings. (as generated by `RTSettings()`). See VPL documentation for more details on the ray tracer.
"""
struct RayTracer{A, M, S, RT}
    scene::A
    materials::M
    sources::S
    settings::RT
end

"""
    RayTracer(scene, sources; settings = RTSettings(), acceleration = Naive, rule = nothing)

Create a `RayTracer` object from a scene (as generated by `Scene`), a tuple of sources (objects that inherit from `Source`),
or a single source, settings and acceleration function (choose from `Naive` or  `BVH`). The argument `rule` is only required for the
accelerator `BVH` and it must be an object of type `SAH` or `AvgSplit` (it is ignored for the `Naive` accelerator).

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> sc = Scene(mesh = Ellipse());

julia> source = DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);

julia> rt = RayTracer(sc, source);

julia> sources = (DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000),
                  DirectionalSource(sc, θ = 45.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000));

julia> rt = RayTracer(sc, sources);
```
"""
function RayTracer(scene::Scene, sources; settings = RTSettings(),
    acceleration = Naive, rule = nothing)
    # Construct the acceleration structure and grid cloner around the scene
    acc_scene = accelerate(scene, settings = settings, acceleration = acceleration,
        rule = rule)
    # Check for use directional light sources in the absence of the grid cloner
    acc_scene.grid.nleaves == 1 && any(source.geom isa Directional for source in sources) &&
        println("Using Directional sources in the absence of a grid cloner may lead to incorrect results. See VPL documentation for details.")
    RayTracer(acc_scene, materials(scene), sources, settings)
end

function RayTracer(scene::Scene, source::Source; kwargs...)
    RayTracer(scene, (source,); kwargs...)
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

julia> sc = Scene(mesh = Ellipse());

julia> source = DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);

julia> rt = RayTracer(sc, source);

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
that is included in the scene. It returns the total number of rays being traced (primary and secondary).
See VPL documentation for more details on the ray tracer.

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> sc = Scene(mesh = Ellipse());

julia> source = DirectionalSource(sc, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);

julia> rt = RayTracer(sc, source);

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
        #materials = [deepcopy(rt.materials) for _ in 1:nthreads()]
        # Each thread gets a different pseudo-random number generator which seed depends on the original one
        samplers = [deepcopy(rt.settings.sampler) for _ in 1:Threads.nthreads()]
        for i in 1:Threads.nthreads()
            Random.seed!(samplers[i], rand(rt.settings.sampler, UInt))
        end
        # Allocated primary rays evenly across threads
        Δrays = div(total_rays, nthreads())
        irays = collect(0:(nthreads() - 1)) .* Δrays .+ 1
        erays = collect(1:nthreads()) .* Δrays
        erays[end] = total_rays
        # Loop over threads while updating materials and collecting total number of rays traced
        nrays_thread = zeros(Int, nthreads())
        @threads for i in 1:nthreads()
            @inbounds nrays_thread[i] = trace_thread!(rt,
                rt.materials,
                #materials[i],
                nrays,
                irays[i],
                erays[i],
                samplers[i])
        end
        nrays_traced = sum(nrays_thread)
        # Copy the power stored in each material back to the original
        # for it in 1:nthreads()
        #     for im in 1:length(rt.materials)
        #         @inbounds rt.materials[im].power .+= materials[it][im].power
        #     end
        # end
        # Run ray tracer in a single multiple thread
    else
        nrays_traced = trace_thread!(rt,
            rt.materials,
            nrays,
            1,
            total_rays,
            rt.settings.sampler)
    end
    return nrays_traced
end

# Trace a subset of the rays on a given PC thread
function trace_thread!(rt, materials, nrays, i_ray, e_ray, rng)
    # Get the current floating point precision out of the ray tracer
    FT = eltype(rt.scene.acc.gbox.min)
    # Number of wavelengths
    nw = get_nw(rt)
    # Keep track of all the rays being traced (including secondary)
    nrays_traced = 0
    # Setup nodestack and power payload and trace all primary rays
    nodestack = Int[]
    tnodestack = Int[]
    dstack = FT[]
    tdstack = FT[]
    power = MVector{nw, Float64}(Tuple(0.0 for _ in 1:nw))
    for i in i_ray:e_ray
        r = generate_ray!(rt.sources, nrays, i, power, rng)
        nrays_traced += trace!(r,
            rt,
            materials,
            tnodestack,
            tdstack,
            nodestack,
            dstack,
            power,
            rng)
    end
    return nrays_traced
end

# Trace a ray through a scene with Russian roulette
function trace!(r::Ray{FT}, rt::RayTracer, materials, tnodestack::Vector{Int},
    tdstack::Vector{FT},
    nodestack::Vector{Int}, dstack::Vector{FT}, power, rng) where {FT}
    # Recursive Monte Carlo ray tracing with Russian roulette criterion
    iteration = 1
    while true
        # Check intersection against scene
        hit, intersection, disp = intersect(r,
            rt.scene.grid,
            rt.scene.acc,
            tnodestack,
            tdstack,
            nodestack,
            dstack)
        !hit && (return iteration)#(@show loops; return iteration)
        material = materials[intersection.id]
        # Interaction with surface material
        interaction = calculate_interaction(material, power, r, intersection, rng)
        absorb_power!(material, power, interaction)
        # Russian roulette
        roulette!(power, rt.settings, iteration, rng) && (return iteration)#(@show loops; return iteration)
        # Generate new ray
        r = generate_ray(material, r, disp, intersection, interaction, rng)
        # Increase iteration counter (unless it is a sensor-like material)
        # TODO: if sensors do not increase iteration counter the ray tracer seems to "get stuck" in some cases
        #if interaction.mode != :sensor
        iteration += 1
        #end
    end
    return iteration
end
