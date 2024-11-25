### This file contains public API ###
# PointSource
# LineSource
# AreaSource

###############################################################################
############################## Implementations ################################
###############################################################################

"""
    PointSource(vec)

Create a point irradiance source geometry at given 3D location `vec`, defined as vector
of Cartesian coordinates (`Vec(x, y, z)`).

## Examples
```jldoctest
julia> import PlantGeomPrimitives as PG

julia> source_geom = PointSource(PG.Vec(1.0, 1.0, 1.0));
```
"""
struct PointSource{FT} <: SourceGeometry
    loc::PGP.Vec{FT}
end

generate_point(g::PointSource, rng) = g.loc

"""
    LineSource(p, line)

Create a line irradiance source geometry given an origin (`p`) and a segment (`line`) both
specified as vector of Cartesian coordinates (`Vec(x, y, z)`). This will create a
line source between the points `p` and `p .+ line`.

## Examples
```jldoctest
julia> import PlantGeomPrimitives as PG

julia> source_geom = LineSource(PG.Vec(1.0, 1.0, 1.0), PG.Y());
```
"""
struct LineSource{FT} <: SourceGeometry
    p::PGP.Vec{FT}
    line::PGP.Vec{FT} # not an unit vector!
end

function generate_point(g::LineSource{FT}, rng) where {FT}
    @. g.p + rand(rng, FT) * g.line
end

struct AreaSource{FT} <: SourceGeometry
    tvec::Vector{Triangle{FT}}
    areas::StatsBase.Weights{FT, FT, Vector{FT}}
end

"""
    AreaSource(mesh)

Create an area irradiance source geometry given a triangular mesh.

## Examples
```jldoctest
julia> using PlantGeomPrimitives

julia> e = Ellipse();

julia> source_geom = AreaSource(e);
```
"""
AreaSource(mesh::PGP.Mesh) = AreaSource(Triangle(mesh), StatsBase.Weights(PGP.areas(mesh)))

# Select randomly the triangle within the mesh and select randomly a point within the triangle
function generate_point(g::AreaSource, rng)
    t = StatsBase.sample(rng, g.tvec, g.areas)
    generate_point(t, rng)
end
