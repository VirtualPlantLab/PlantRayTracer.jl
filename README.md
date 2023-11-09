# PlantRayTracer

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://virtualplantlab.com/stable/api/raytracer/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://virtualplantlab.com/dev/api/raytracer/)
[![CI](https://github.com/VirtualPlantLab/PlantRayTracer.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/VirtualPlantLab/PlantRayTracer.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/VirtualPlantLab/PlantRayTracer.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/VirtualPlantLab/PlantRayTracer.jl)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This package is a component of the VPL ecosystem. It implements and ray tracer algorithm for
simulating light interception by plants. This package is a component
of the [Virtual Plant Lab](http://virtualplantlab.com/). Users should install instead the
interface package [VirtualPlantLab.jl](https://github.com/VirtualPlantLab/VirtualPlantLab.jl).

# 1. Installation

You can install the latest stable version of PlantRayTracer.jl with the Julia package manager:

```julia
] add PlantRayTracer
```

Or the development version directly from here:

```julia
import Pkg
Pkg.add(url="https://github.com/VirtualPlantLab/PlantRayTracer.jl", rev = "master")
```
