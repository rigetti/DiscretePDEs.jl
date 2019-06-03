using Test, DiscreteExteriorCalculus, DiscretePDEs
const DEC = DiscreteExteriorCalculus
const DPE = DiscretePDEs
using UniqueVectors: UniqueVector
using LinearAlgebra: issymmetric, eigen, norm
using AdmittanceModels: sparse_nullbasis, apply_transform
using Base.Iterators: product

r1, r2 = .5, .4
file_name = joinpath(@__DIR__, "rectangle.geo")
DPE.geo_write!(file_name, characteristic_length_factor=.2,
    footer="""
    Rectangle(1) = {0, 0, 0, $r1, $r2, 0};
    """)
DPE.initialize!()
DPE.gmsh_open!(file_name)
N, K = 2, 3
DPE.mesh!(K)

node_tags, points, tcomp = DPE.get_triangulated_complex(N, K)
@test typeof(tcomp) <: TriangulatedComplex{N, K}
comp = tcomp.complex
DEC.orient!(comp)
_, exterior = DEC.boundary_components_connected(comp)
m = Metric(N)
mesh = Mesh(tcomp, DEC.circumcenter(m))

laplacian = DEC.differential_operator(m, mesh, "Δ", 1, true)
constraint = DPE.zero_constraint(comp, exterior.cells[1], 1)
nullbasis = sparse_nullbasis(constraint)

vals, vects = eigen(collect(transpose(nullbasis) * laplacian * nullbasis))
inds = sortperm(vals)
vals, vects = vals[inds], vects[:, inds]

# solution is ψ(x,y) = sin(π*x/m)*sin(π*y/n) for integers m and n
rect_modes(m, n, R1, R2) = (m * π/R1)^2 + (n * π/R2)^2
n = 15
correct_modes = sort([rect_modes(i,j,r1,r2)
    for (i,j) in vcat(product(1:n, 1:n)...)])[1:n]

@testset "laplacian on a rectangle" begin
    @test isapprox(vals[1:n], correct_modes, rtol=1e-2)
end

if false
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    for i in 1:n
        v = nullbasis * vects[:,i]
        DPE.add_field!("Eigenvector", node_tags, v[ordering])
    end
    DPE.gui!()
end
rm(file_name)
