using Test, DiscreteExteriorCalculus, DiscretePDEs
const DEC = DiscreteExteriorCalculus
const DPE = DiscretePDEs
using UniqueVectors: UniqueVector
using LinearAlgebra: norm
using AdmittanceModels: lossless_modes_sparse, apply_transform, get_Y

# Modes of a rectangular box with constant μ and ϵ.
box_modes(m, n, l, a, b, c, μ, ϵ) = sqrt((m * π/a)^2 + (n * π/b)^2 +
    (l * π/c)^2)/(2π * sqrt(μ * ϵ))
# The lowest mode of our box is the 011 mode.
a, b, c = 10, 12, 14
correct_A(x, y, z) = sin((π/b) * y) * sin((π/c) * z) * [1, 0, 0]

file_name = joinpath(@__DIR__, "box.geo")
DPE.geo_write!(file_name, characteristic_length_factor=1,
    footer="""
    Box(1) = {0, 0, 0, $a, $b, $c};
    """)

DPE.initialize!()
DPE.gmsh_open!(file_name)
N, K = 3, 4
DPE.mesh!(K)
node_tags, points, tcomp = DPE.get_triangulated_complex(N, K)
group_dict = DPE.get_physical_groups(node_tags, points)
@test typeof(tcomp) <: TriangulatedComplex{N, K}
comp = tcomp.complex
DEC.orient!(comp)
boundary = DEC.boundary(comp)
m = Metric(N)
mesh = Mesh(tcomp, DEC.circumcenter(m))

μ⁻, ϵ = 2, 3
μ⁻_form = DPE.get_material(comp, μ⁻, 3)
Λ⁻_form = DPE.get_material(comp, 0, 2)
σ_form = DPE.get_material(comp, 0, 2)
ϵ_form = DPE.get_material(comp, ϵ, 2)

pso, null_basis = DPE.coulomb_pso(m, mesh, Vector{Cell{N}}[], boundary,
    μ⁻_form, Λ⁻_form, σ_form, ϵ_form)
constrained_pso = apply_transform(pso, null_basis)

density(mat) = count(!iszero, mat)/(size(mat, 1) * size(mat, 2))
@test density(null_basis) < .002
@test all(map(density, get_Y(constrained_pso)) .< .05)

λs, vs = lossless_modes_sparse(constrained_pso, maxiter=1e5)
freq = imag(λs[1])/(2π)
correct_freq = box_modes(0,1,1,a,b,c,1/μ⁻,ϵ)

comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
correct_vec_A = [correct_A(p.coords...) for p in comp_points]
correct_vec_A /= maximum(norm.(correct_vec_A))
vec_A = DEC.sharp(m, comp, null_basis * vs[:,1])
vec_A /= maximum(norm.(vec_A))
i = argmax(norm.(correct_vec_A))
vec_A *= sign(transpose(vec_A[i]) * correct_vec_A[i])

@testset "box modes" begin
    @test isapprox(freq, correct_freq, rtol=5e-3)
    @test maximum([norm(x - y) for (x, y) in zip(vec_A, correct_vec_A)]) < .25
end

if false
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    DPE.add_field!("Vector potential", node_tags, vec_A[ordering])
    DPE.add_field!("Correct vector potential", node_tags, correct_vec_A[ordering])
    DPE.gui!()
end
rm(file_name)
