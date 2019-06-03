using Test, DiscreteExteriorCalculus, DiscretePDEs
const DEC = DiscreteExteriorCalculus
const DPE = DiscretePDEs
using UniqueVectors: UniqueVector
using LinearAlgebra: norm
using AdmittanceModels: lossless_modes_dense, apply_transform, get_Y

# Modes of a rectangular box with constant μ and spatially varying ϵ.
# Specifically, ϵ for y < y0 and ϵr * ϵ for y > y0. For ϵr close to 1, the
# lowest mode is still the 011 mode, but the field is now
# A = sin(α * y) * sin(π/c * z) * cos(ω * t) * x_unit for y < y0
# A = sin(β * y + γ) * sin(π/c * z) * cos(ω * t) * x_unit for y >= y0
# where α, β, γ are chosen so that A_x is continuous at y_0, goes to 0 at y=b,
# and the Helmholtz equation is satisfied with a single value of ω.
# These conditions give a quadratic for α which is solved in the code below.

a, b, c, y0 = 10, 12, 14, 8
μ, ϵ = 2, 3
ϵr = 1.5
# compute analytical mode
quadratic_coef = y0^2/(ϵr * (y0 - b)^2) - 1
linear_coef = -2 * y0 * π/(ϵr * (y0 - b)^2)
constant_coef = (π/c)^2 * (1/ϵr - 1) + π^2/(ϵr * (y0 - b)^2)
αs = DEC.solve_quadratic(quadratic_coef, linear_coef, constant_coef)
α = αs[argmin(abs.(αs .- π/b))]
β = (α * y0 - π)/(y0 - b)
γ = π - β * b
correct_freq = sqrt((π/c)^2 + α^2)/(2π * sqrt(μ * ϵ))
correct_A(x, y, z) = (y < y0 ? sin(α * y) : sin(β * y + γ)) * sin((π/c) * z) * [1, 0, 0]

# compute mode numerically
file_name = joinpath(@__DIR__, "box_dielectric.geo")
DPE.geo_write!(file_name, characteristic_length_factor=.8,
    footer="""
    Box(1) = {0, 0, 0, $a, $b, $c};
    Box(2) = {0, $y0, 0, $a, $(b-y0), $c};
    BooleanDifference{Volume{1}; Delete;}{Volume{2};}
    Physical Volume("1") = {1};
    Physical Volume("2") = {2};
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

μ⁻_form = DPE.get_material(comp, 1/μ, 3)
Λ⁻_form = DPE.get_material(comp, 0, 2)
σ_form = DPE.get_material(comp, 0, 2)
ϵ_form = DPE.get_material(comp, group_dict["1"], ϵ, 2) +
    DPE.get_material(comp, group_dict["2"], ϵr * ϵ, 2)

pso, null_basis = DPE.coulomb_pso(m, mesh, Vector{Cell{N}}[], boundary,
    μ⁻_form, Λ⁻_form, σ_form, ϵ_form)
constrained_pso = apply_transform(pso, null_basis)

λs, vs = lossless_modes_dense(constrained_pso)
freq = imag.(λs[1])/(2π)

comp_points = UniqueVector([c.points[1] for c in comp.cells[1]])
correct_vec_A = [correct_A(p.coords...) for p in comp_points]
correct_vec_A /= maximum(norm.(correct_vec_A))
vec_A = DEC.sharp(m, comp, null_basis * vs[:,1])
vec_A /= maximum(norm.(vec_A))
i = argmax(norm.(correct_vec_A))
vec_A *= sign(transpose(vec_A[i]) * correct_vec_A[i])

@testset "box modes with dielectric" begin
    @test isapprox(freq, correct_freq, rtol=2e-1)
    @test maximum([norm(x - y) for (x, y) in zip(vec_A, correct_vec_A)]) < .25
end

if false
    ordering = [findfirst(isequal(p), comp_points) for p in points]
    DPE.add_field!("Vector potential", node_tags, vec_A[ordering])
    DPE.add_field!("Correct vector potential", node_tags, correct_vec_A[ordering])
    DPE.gui!()
end
rm(file_name)
