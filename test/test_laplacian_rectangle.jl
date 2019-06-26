using Test, DiscreteExteriorCalculus, DiscretePDEs
using UniqueVectors: UniqueVector
using LinearAlgebra: eigen
using AdmittanceModels: sparse_nullbasis
using Base.Iterators: product

r1, r2 = .5, .4
file_name = joinpath(@__DIR__, "rectangle.geo")
geo_write!(file_name, characteristic_length_factor=.2,
    footer="""
    Rectangle(1) = {0, 0, 0, $r1, $r2, 0};
    """)
initialize!()
gmsh_open!(file_name)
N, K = 2, 3
mesh!(K)

node_tags, points, tcomp = get_triangulated_complex(N, K)
@test typeof(tcomp) <: TriangulatedComplex{N, K}
comp = tcomp.complex
orient!(comp)
_, exterior = boundary_components_connected(comp)
m = Metric(N)
mesh = Mesh(tcomp, circumcenter(m))

laplacian = differential_operator(m, mesh, "Δ", 1, true)
constraint = zero_constraint(comp, exterior.cells[1], 1)
nullbasis = sparse_nullbasis(constraint)

vals, vects = eigen(collect(transpose(nullbasis) * laplacian * nullbasis))
inds = sortperm(vals)
vals, vects = vals[inds], vects[:, inds]

# solution is ψ(x,y) = sin(m*π*x/R1)*sin(n*π*y/R2) for integers m and n
rect_modes(m, n, R1, R2) = (m * π/R1)^2 + (n * π/R2)^2
n = 15
correct_modes = sort([rect_modes(i,j,r1,r2)
    for (i,j) in vcat(product(1:n, 1:n)...)])[1:n]

@testset "laplacian on a rectangle" begin
    @test isapprox(vals[1:n], correct_modes, rtol=1e-2)
end

if false # set to true to display plots with Gmsh
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    for i in 1:n
        v = nullbasis * vects[:,i]
        add_field!("Eigenvector", node_tags, v[ordering])
    end
    gui!()
end
rm(file_name)
