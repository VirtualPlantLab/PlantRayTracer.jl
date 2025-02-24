### This file contains public API ###
# DirectionalSource

"""
    DirectionalSource(box::AABB; θ, Φ, radiosity, nrays)
    DirectionalSource(mesh::Mesh; θ, Φ, radiosity, nrays)

Create a Directional source (including geometry and angle components) by providing an axis-aligned
bounding box (`box`) or an `Mesh` object  as well as the zenith (`θ`) and azimuth (`Φ`)
angles, the radiosity of the source projected on the horizontal plane and the number of
rays to be generated. Directional sources may generate incorrect results in the absence of
a grid cloner that extends the mesh. This is because the rays are generated from the upper
face of the mesh's bounding box. See VPL documentation for details on light sources.

## Examples
```jldoctest
julia> using PlantGeomPrimitives;

julia> mesh = Ellipse();

julia> source = DirectionalSource(mesh, θ = 0.0, Φ = 0.0, radiosity = 1.0, nrays = 1_000);
```
"""
function DirectionalSource(box::AABB; θ, Φ, radiosity, nrays)
    dir_geom = create_directional(box)
    # Radiosity is projected onto horizontal plane as we sample from the top of the bounding box
    # The code below ensures that we get the right irradiance onto the mesh
    power = radiosity .* area(dir_geom)
    out = Source(dir_geom, FixedSource(θ, Φ), power./nrays, nrays)
end
function DirectionalSource(mesh::PGP.Mesh; θ, Φ, radiosity, nrays)
    box = AABB(mesh)
    DirectionalSource(box, θ = θ, Φ = Φ, radiosity = radiosity, nrays = nrays)
end

# Upper face of the mesh's AABB
struct Directional{FT} <: SourceGeometry
    xmin::FT
    xmax::FT
    ymin::FT
    ymax::FT
    zmax::FT
end

# Calculate the area of the directional light source
area(d::Directional) = (d.xmax - d.xmin) * (d.ymax - d.ymin)

# Create geometry for a directional light source
function create_directional(box::AABB{FT}) where {FT}
    @inbounds Directional(box.min[1], box.max[1], box.min[2], box.max[2], box.max[3])
end

# Select a random point from the upper face of the mesh's AABB
function generate_point(d::Directional{FT}, rng) where {FT}
    x = d.xmin + rand(rng, FT) * (d.xmax - d.xmin)
    y = d.ymin + rand(rng, FT) * (d.ymax - d.ymin)
    PGP.Vec(x, y, d.zmax + 100eps(FT))
end
