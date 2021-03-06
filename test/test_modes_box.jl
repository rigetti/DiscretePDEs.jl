using Test, DiscreteExteriorCalculus, DiscretePDEs
using UniqueVectors: UniqueVector
using LinearAlgebra: norm
using AdmittanceModels: lossless_modes_dense, apply_transform, get_Y

# Modes of a rectangular box with constant μ and ϵ.
box_modes(m, n, l, a, b, c, μ, ϵ) = sqrt((m * π/a)^2 + (n * π/b)^2 +
    (l * π/c)^2)/(2π * sqrt(μ * ϵ))
# The lowest mode of our box is the 011 mode.
a, b, c = 10, 12, 14
correct_A(x, y, z) = sin((π/b) * y) * sin((π/c) * z) * [1, 0, 0]

file_name = joinpath(@__DIR__, "box.geo")
geo_write!(file_name, characteristic_length_factor=1,
    footer="""
    Box(1) = {0, 0, 0, $a, $b, $c};
    """)

initialize!()
gmsh_open!(file_name)
N, K = 3, 4
mesh!(K)
node_tags, points, tcomp = get_triangulated_complex(N, K)
group_dict = get_physical_groups(node_tags, points)
@test typeof(tcomp) <: TriangulatedComplex{N, K}
comp = tcomp.complex
orient!(comp)
m = Metric(N)
mesh = Mesh(tcomp, circumcenter(m))

μ⁻, ϵ = 2, 3
μ⁻_form = get_material(comp, μ⁻, 3)
Λ⁻_form = get_material(comp, 0, 2)
σ_form = get_material(comp, 0, 2)
ϵ_form = get_material(comp, ϵ, 2)

pso, null_basis = electrodynamics_pso(m, mesh, Vector{Cell{N}}[], boundary(comp),
    μ⁻_form, Λ⁻_form, σ_form, ϵ_form)
constrained_pso = apply_transform(pso, null_basis)

density(mat) = count(!iszero, mat)/(size(mat, 1) * size(mat, 2))
@test density(null_basis) < .002
@test all(map(density, get_Y(constrained_pso)) .< .05)

λs, vs = lossless_modes_dense(constrained_pso)
freq = imag(λs[1])/(2π)
correct_freq = box_modes(0,1,1,a,b,c,1/μ⁻,ϵ)

comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
correct_vec_A = [correct_A(p.coords...) for p in comp_points]
correct_vec_A /= maximum(norm.(correct_vec_A))
vec_A = sharp(m, comp, null_basis * vs[:,1])
vec_A /= maximum(norm.(vec_A))
i = argmax(norm.(correct_vec_A))
vec_A *= sign(transpose(vec_A[i]) * correct_vec_A[i])

@testset "box modes" begin
    @test isapprox(freq, correct_freq, rtol=5e-3)
    @test maximum([norm(x - y) for (x, y) in zip(vec_A, correct_vec_A)]) < .25
end

if false # set to true to display plots with Gmsh
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    add_field!("Vector potential", node_tags, vec_A[ordering])
    add_field!("Correct vector potential", node_tags, correct_vec_A[ordering])
    gui!()
end
rm(file_name)
