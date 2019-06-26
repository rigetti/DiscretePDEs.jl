using UniqueVectors: UniqueVector

abstract type Geo end

struct GeoPoint <: Geo
    x::Float64
    y::Float64
    z::Float64
    mesh_size::Float64
end

struct GeoLine <: Geo
    p1::GeoPoint
    p2::GeoPoint
end

struct GeoCurve <: Geo
    lines::Vector{GeoLine}
end

struct GeoPlane <: Geo
    curve::GeoCurve
end

struct GeoSurfaceDifference <: Geo
    p1::GeoPlane
    p2::GeoPlane
    keep_tool::Bool
end

struct GeoSurfaceExtrusion <: Geo
    x::Float64
    y::Float64
    z::Float64
    plane::GeoPlane
end

struct GeoVolumeDifference <: Geo
    v1::GeoSurfaceExtrusion
    v2::GeoSurfaceExtrusion
    keep_tool::Bool
end

struct GeoSurfaceGroup <: Geo
    name::String
    planes::Vector{GeoPlane}
end

struct GeoVolumeGroup <: Geo
    name::String
    extrusions::Vector{GeoSurfaceExtrusion}
end

export GeoCode
"""
    GeoCode(points::UniqueVector{GeoPoint}, lines::UniqueVector{GeoLine}
        curves::UniqueVector{GeoCurve}, planes::UniqueVector{GeoPlane},
        surface_differences::UniqueVector{GeoSurfaceDifference},
        extrusions::UniqueVector{GeoSurfaceExtrusion},
        volume_differences::UniqueVector{GeoVolumeDifference},
        surface_groups::UniqueVector{GeoSurfaceGroup},
        volume_groups::UniqueVector{GeoVolumeGroup})
    GeoCode()

A representation of the .geo format for constructive solid geometry accepted by Gmsh.
"""
struct GeoCode <: Geo
    points::UniqueVector{GeoPoint}
    lines::UniqueVector{GeoLine}
    curves::UniqueVector{GeoCurve}
    planes::UniqueVector{GeoPlane}
    surface_differences::UniqueVector{GeoSurfaceDifference}
    extrusions::UniqueVector{GeoSurfaceExtrusion}
    volume_differences::UniqueVector{GeoVolumeDifference}
    surface_groups::UniqueVector{GeoSurfaceGroup}
    volume_groups::UniqueVector{GeoVolumeGroup}
end

import Base: ==
function ==(g1::Geo, g2::Geo)
    if typeof(g1) != typeof(g2)
        return false
    end
    for n in fieldnames(typeof(g1))
        if getfield(g1, n) != getfield(g2, n)
            return false
        end
    end
    return true
end

GeoCode() = GeoCode([UniqueVector(eltype(t)[]) for t in fieldtypes(GeoCode)]...)

import Base: append!
function append!(g1::GeoCode, g2::GeoCode)
    for n in fieldnames(GeoCode)
        arr = getfield(g1, n)
        for e in getfield(g2, n)
            if !(e in arr)
                push!(arr, e)
            end
        end
    end
end

function compile(gc::GeoCode, p::GeoPoint)
    t = findfirst(isequal(p), gc.points)
    return "Point($t) = {$(p.x), $(p.y), $(p.z), $(p.mesh_size)};"
end

function compile(gc::GeoCode, l::GeoLine)
    lt = findfirst(isequal(l), gc.lines)
    pt1, pt2 = findfirst(isequal(l.p1), gc.points), findfirst(isequal(l.p2), gc.points)
    return "Line($lt) = {$pt1, $pt2};"
end

function vector_to_string(l::Vector)
    s = "$(tuple(l...))"[2:end-1]
    if length(l) == 1
        return s[1:end-1]
    else
        return s
    end
end

function compile(gc::GeoCode, c::GeoCurve)
    ct = findfirst(isequal(c), gc.curves)
    lts = [findfirst(isequal(l), gc.lines) for l in c.lines]
    return "Curve Loop($ct) = {" * vector_to_string(lts) * "};"
end

function compile(gc::GeoCode, p::GeoPlane)
    pt = findfirst(isequal(p), gc.planes)
    ct = findfirst(isequal(p.curve), gc.curves)
    return "Plane Surface($pt) = {$ct};"
end

function compile(gc::GeoCode, sd::GeoSurfaceDifference)
    pt1 = findfirst(isequal(sd.p1), gc.planes)
    pt2 = findfirst(isequal(sd.p2), gc.planes)
    ending = sd.keep_tool ? "" : " Delete;"
    return "BooleanDifference{Surface{$pt1}; Delete;}{Surface{$pt2};$ending}"
end

function compile(gc::GeoCode, se::GeoSurfaceExtrusion)
    pt = findfirst(isequal(se.plane), gc.planes)
    return "Extrude {$(se.x), $(se.y), $(se.z)} {Surface{$pt};}"
end

function compile(gc::GeoCode, vd::GeoVolumeDifference)
    vt1 = findfirst(isequal(vd.v1), gc.extrusions)
    vt2 = findfirst(isequal(vd.v2), gc.extrusions)
    ending = vd.keep_tool ? "" : " Delete;"
    return "BooleanDifference{Volume{$vt1}; Delete;}{Volume{$vt2};$ending}"
end

function compile(gc::GeoCode, group::GeoSurfaceGroup)
    ts = [findfirst(isequal(p), gc.planes) for p in group.planes]
    return "Physical Surface(\"$(group.name)\") = {" * vector_to_string(ts) * "};"
end

function compile(gc::GeoCode, group::GeoVolumeGroup)
    ts = [findfirst(isequal(p), gc.extrusions) for p in group.extrusions]
    return "Physical Volume(\"$(group.name)\") = {" * vector_to_string(ts) * "};"
end

function compile(gc::GeoCode)
    code = String["//points\n"]
    for p in gc.points
        push!(code, compile(gc, p) * "\n")
    end
    push!(code, "\n//lines\n")
    for l in gc.lines
        push!(code, compile(gc, l) * "\n")
    end
    push!(code, "\n//curves\n")
    for c in gc.curves
        push!(code, compile(gc, c) * "\n")
    end
    push!(code, "\n//planes\n")
    for p in gc.planes
        push!(code, compile(gc, p) * "\n")
    end
    push!(code, "\n//surface differences\n")
    for sd in gc.surface_differences
        push!(code, compile(gc, sd) * "\n")
    end
    push!(code, "\n//surface extrusions\n")
    for se in gc.extrusions
        push!(code, compile(gc, se) * "\n")
    end
    push!(code, "\n//volume differences\n")
    for vd in gc.volume_differences
        push!(code, compile(gc, vd) * "\n")
    end
    push!(code, "\n//surface groups\n")
    for group in gc.surface_groups
        push!(code, compile(gc, group) * "\n")
    end
    push!(code, "\n//volume groups\n")
    for group in gc.volume_groups
        push!(code, compile(gc, group) * "\n")
    end
    return code
end

export geo_write!
"""
    geo_write!(file_name::String, gc::GeoCode=GeoCode(); random_factor=1e-9,
        random_factor_3D=1e-9, characteristic_length_factor::Real=1, footer::String="")

Write the contents of a GeoCode to a file along with additional geo code in `footer`.
The `random_factor` and `random_factor_3D` arguments adjust how Gmsh uses randomness in its
2D and 3D meshing algorithms, respectively. `characteristic_length_factor` scales the
characteristic length of the mesh; smaller values yield denser meshes.
"""
function geo_write!(file_name::String, gc::GeoCode=GeoCode(); random_factor=1e-9,
    random_factor_3D=1e-9, characteristic_length_factor::Real=1, footer::String="")
    code_strings = ["Mesh.RandomFactor=$random_factor;\n";
                    "Mesh.RandomFactor3D=$random_factor_3D;\n";
                    "Mesh.CharacteristicLengthFactor=$characteristic_length_factor;\n";
                    "SetFactory(\"OpenCASCADE\");\n\n";
                    compile(gc);
                    "\n//additional\n";
                    footer]
    file = Base.open(file_name, write=true)
    for s in code_strings
        write(file, s)
    end
    close(file)
end

"""
    GeoCode(points::AbstractVector{<:NTuple{4, <:Real}})

Generate the GeoCode for a closed 3D curve. The 4th element of each point is a mesh size
that Gmsh uses to spatially modulate the mesh density.
"""
function GeoCode(points::AbstractVector{<:NTuple{4, <:Real}})
    gc = GeoCode()
    for p in points
        push!(gc.points, GeoPoint(p...))
    end
    n = length(points)
    for i in 1:n
        p1, p2 = gc.points[i], gc.points[mod1(i+1, n)]
        push!(gc.lines, GeoLine(p1, p2))
    end
    if n > 0
        geo_curve = GeoCurve(gc.lines)
        push!(gc.curves, geo_curve)
        push!(gc.planes, GeoPlane(geo_curve))
    end
    return gc
end
