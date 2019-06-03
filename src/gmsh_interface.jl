using DiscreteExteriorCalculus: Point, Simplex, TriangulatedComplex
using UniqueVectors: UniqueVector
using DotEnv

DotEnv.config(path=joinpath(@__DIR__, "../.env"))
include(ENV["gmshjlpath"])

# this code has been tested and works well on both Mac OSX Mojave and Ubuntu 18.10
function eventloop()
    while true
        gmsh.graphics.draw()
        gmsh.fltk.wait()
        sleep(.01)
    end
end

initialize!() = gmsh.initialize()
gui!() = schedule(Task(eventloop))
# TODO: what's the right way to end this Task and finalize the gmsh
# process without killing the main Julia process?
gmsh_open!(file_name::String) = gmsh.open(file_name)
# K is the dimension + 1 of the manifold (i.e. number of points per simplex)
mesh!(K::Int) = gmsh.model.mesh.generate(K-1)

# if tag is -1, create a new view and return the new tag
# if tag is greater than 0, replace the existing view with that tag
function add_field!(name::String, node_tags::AbstractVector{Int},
                   vector_field::AbstractVector{<:AbstractVector{<:Real}}, tag::Int=-1)
    view_tag = gmsh.view.add(name, tag)
    gmsh.view.addModelData(view_tag, 0, "", "NodeData", node_tags, vector_field)
    return view_tag
end

function add_field!(name::String, node_tags::AbstractVector{Int},
                   scalar_field::AbstractVector{<:Real}, tag::Int=-1)
    view_tag = gmsh.view.add(name, tag)
    gmsh.view.addModelData(view_tag, 0, "", "NodeData", node_tags, [[x] for x in scalar_field])
    return view_tag
end

function get_node_tags(N::Int)
    node_tags, coords = gmsh.model.mesh.getNodes()
    node_tags = Int.(node_tags)
    coords = reshape(Float64.(coords), 3, length(node_tags))[1:N, :]
    return UniqueVector(node_tags), coords
end

function get_points(N::Int, scale::Real=1.0)
    node_tags, coords = get_node_tags(N)
    coords *= scale
    points = [Point(coords[:, i]) for i in 1:size(coords,2)]
    return node_tags, points
end

# see gmsh.model.mesh.getElementProperties
# K is the number of points per simplex, i.e. the dimension + 1
# 2 => line, 3 => triangle, 4 => tetrahedra
const K_to_element_type = Dict(2 => 1, 3 => 2, 4 => 4)
function get_simplex_node_tags(K::Int, tag::Int=-1)
    simplex_tags, simplex_node_tags = gmsh.model.mesh.getElementsByType(K_to_element_type[K], tag)
    simplex_node_tags = reshape(Int.(simplex_node_tags), K, length(simplex_tags))
    return simplex_node_tags
end

function get_simplices(K::Int, node_tags::AbstractVector{Int},
    points::AbstractVector{Point{N}}, tag::Int=-1) where N
    simplex_node_tags = get_simplex_node_tags(K, tag)
    simplices = [Simplex([points[findfirst(isequal(nt), node_tags)]
        for nt in simplex_node_tags[:, i]]) for i in 1:size(simplex_node_tags, 2)]
    return simplices
end

function get_physical_group_tags()
    dim_tags = gmsh.model.getPhysicalGroups()
    Ks = [Int(dim_tags[1])+1 for dim_tags in dim_tags]
    names = [gmsh.model.getPhysicalName(dim_tag...) for dim_tag in dim_tags]
    entity_tags = [Int.(gmsh.model.getEntitiesForPhysicalGroup(dim_tag...)) for dim_tag in dim_tags]
    return names, Ks, entity_tags
end

function get_physical_groups(node_tags::AbstractVector{Int},
    points::AbstractVector{Point{N}}) where N
    names, Ks, entity_tags = get_physical_group_tags()
    groups = [vcat([get_simplices(K, node_tags, points, tag) for tag in tags]...)
        for (K, tags) in zip(Ks, entity_tags)]
    return Dict{String, Vector{Simplex{N, K}} where K}(n => g for (n, g) in zip(names, groups))
end

function get_triangulated_complex(N::Int, K::Int, scale::Real=1.0)
    node_tags, points = get_points(N, scale)
    simplices = get_simplices(K, node_tags, points)
    return node_tags, points, TriangulatedComplex(simplices)
end
