using Test, DiscreteExteriorCalculus, DiscretePDEs
using UniqueVectors: UniqueVector
using LinearAlgebra: eigen

file_name = joinpath(@__DIR__, "sphere.geo")
geo_write!(file_name, characteristic_length_factor=.3,
    footer="""
    Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
    """)
initialize!()
gmsh_open!(file_name)
N, K = 3, 3
mesh!(K)

node_tags, points, tcomp = get_triangulated_complex(N, K)
@test typeof(tcomp) <: TriangulatedComplex{N, K}
comp = tcomp.complex
orient!(comp)
m = Metric(N)
mesh = Mesh(tcomp, circumcenter(m))

laplacian = differential_operator(m, mesh, "Δ", 1, true)

vals, vects = eigen(collect(laplacian))
inds = sortperm(vals)
vals, vects = vals[inds], vects[:, inds]

# solutions are spherical harmonics (or combinations of them since they
# are degenerate)
correct_modes = vcat([[l*(l+1) for _ in -l:l] for l in 0:3]...)

@testset "laplacian on a sphere" begin
    @test isapprox(vals[1:length(correct_modes)], correct_modes, rtol=1e-2)
end

if false # set to true to display plots with Gmsh
    comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    for i in 1:length(correct_modes)
        v = vects[:,i]
        add_field!("Eigenvector", node_tags, v[ordering])
    end
    gui!()
end
rm(file_name)
