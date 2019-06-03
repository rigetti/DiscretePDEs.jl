using Test, DiscreteExteriorCalculus, DiscretePDEs
const DEC = DiscreteExteriorCalculus
const DPE = DiscretePDEs
using UniqueVectors: UniqueVector
using LinearAlgebra: issymmetric, eigen, norm
using AdmittanceModels: sparse_nullbasis, apply_transform
using Base.Iterators: product

file_name = joinpath(@__DIR__, "sphere.geo")
DPE.geo_write!(file_name, characteristic_length_factor=.3,
    footer="""
    Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
    """)
DPE.initialize!()
DPE.gmsh_open!(file_name)
N, K = 3, 3
DPE.mesh!(K)

node_tags, points, tcomp = DPE.get_triangulated_complex(N, K)
@test typeof(tcomp) <: TriangulatedComplex{N, K}
comp = tcomp.complex
DEC.orient!(comp)
m = Metric(N)
mesh = Mesh(tcomp, DEC.circumcenter(m))

laplacian = DEC.differential_operator(m, mesh, "Δ", 1, true)

vals, vects = eigen(collect(laplacian))
inds = sortperm(vals)
vals, vects = vals[inds], vects[:, inds]

# solutions are spherical harmonics (or combinations of them since they
# are degenerate)
correct_modes = vcat([[l*(l+1) for _ in -l:l] for l in 0:3]...)

@testset "laplacian on a sphere" begin
    @test isapprox(vals[1:length(correct_modes)], correct_modes, rtol=1e-2)
end

if false
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    for i in 1:length(correct_modes)
        v = vects[:,i]
        DPE.add_field!("Eigenvector", node_tags, v[ordering])
    end
    DPE.gui!()
end
rm(file_name)
