using DiscreteExteriorCalculus: Metric, Simplex, Mesh
import DiscreteExteriorCalculus; const DEC = DiscreteExteriorCalculus
using SparseArrays: spdiagm
using AdmittanceModels: Blackbox, PSOModel, apply_transform, sparse_nullbasis,
    impedance_matrices, admittance_matrices

################################################################################
# Electrostatics
################################################################################

# Each source is a vector of cells with K = 1.
function electrostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary_points::AbstractVector{Cell{N}},
    ϵ::AbstractVector{<:Real}) where {N, K}
    k = 1
    pm = -DEC.hodge_square_sign(m, K, k)
    d₁, ★, d₀ = DEC.differential_operator_sequence(m, mesh, "d★d", k, true)
    Y = pm * d₁ * ★ * spdiagm(0 => ϵ) * d₀
    comp = mesh.primal.complex
    row_inds = [findfirst(isequal(s[1]), comp.cells[k]) for s in sources]
    n = length(sources)
    P = sparse(row_inds, 1:n, ones(n), length(comp.cells[k]), n)
    bbox = Blackbox([0.0], [Y], P, P, collect(1:n))
    conductors = filter(s -> length(s) > 1, sources)
    constraints = [zero_constraint(comp, boundary_points, k),
        [constant_constraint(comp, s, k) for s in conductors]...]
    return bbox, sparse_nullbasis(constraints)
end

electrostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary::CellComplex{N},
    ϵ::AbstractVector{<:Real}) where {N, K} =
    electrostatics_blackbox(m, mesh, sources, boundary.cells[1], ϵ)

get_charge_source(comp::CellComplex{N}, group::AbstractVector{Simplex{N, K}}) where {N, K} =
    DEC.submanifold(comp, group).cells[1]

################################################################################
# Electrodynamics without charge
################################################################################

# Each source is a vector of cells with K = 2. Note that each source
# must be consistently oriented.
function coulomb_pso(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary_points::AbstractVector{Cell{N}},
    boundary_edges::AbstractVector{Cell{N}},
    μ⁻::AbstractVector{<:Real}, Λ⁻::AbstractVector{<:Real},
    σ::AbstractVector{<:Real}, ϵ::AbstractVector{<:Real}) where {N, K}
    k = 2
    pm = DEC.hodge_square_sign(m, K, k)
    d₁, ★, d₀ = DEC.differential_operator_sequence(m, mesh, "d★d", k, true)
    K₀ = pm * d₁ * ★ * spdiagm(0 => μ⁻) * d₀
    d, ★ = DEC.differential_operator_sequence(m, mesh, "d★", k, true)
    K₁ = ★ * spdiagm(0 => Λ⁻)
    G = ★ * spdiagm(0 => σ)
    C = ★ * spdiagm(0 => ϵ)
    comp = mesh.primal.complex
    row_inds = vcat([[findfirst(isequal(c), comp.cells[k]) for c in s] for s in sources]...)
    col_inds = vcat([[i for _ in s] for (i, s) in enumerate(sources)]...)
    P = sparse(row_inds, col_inds, ones(length(row_inds)), length(comp.cells[k]), length(sources))
    pso = PSOModel(K₀ + K₁, G, C, P, P, collect(1:length(sources)))
    interior_inds = findall(i -> !(comp.cells[1][i] in boundary_points), 1:length(comp.cells[1]))
    d★ϵ_constraint = (d * C/maximum(abs.(ϵ)))[interior_inds, :] # only enforce d★ϵA = 0 in the interior
    constraint = vcat(zero_constraint(comp, boundary_edges, k), d★ϵ_constraint)
    return pso, sparse_nullbasis(constraint)
end

coulomb_pso(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary::CellComplex{N},
    μ⁻::AbstractVector{<:Real}, Λ⁻::AbstractVector{<:Real},
    σ::AbstractVector{<:Real}, ϵ::AbstractVector{<:Real}) where {N, K} =
    coulomb_pso(m, mesh, sources, boundary.cells[1], boundary.cells[2], μ⁻, Λ⁻, σ, ϵ)

get_current_source(comp::CellComplex{N}, group::AbstractVector{Simplex{N, K}}) where {N, K} =
    DEC.orient!(DEC.submanifold(comp, group).cells[2])

################################################################################
# Magnetostatics
################################################################################

# Each source is a vector of cells with K = 2. Note that each source
# must be consistently oriented.
function magnetostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary_points::AbstractVector{Cell{N}},
    boundary_edges::AbstractVector{Cell{N}},
    μ⁻::AbstractVector{<:Real}, Λ⁻::AbstractVector{<:Real}) where {N, K}
    pso, null_basis = coulomb_pso(m, mesh, sources, boundary_points,
        boundary_edges, μ⁻, Λ⁻, zeros(size(Λ⁻)), ones(size(Λ⁻)))
    bbox = Blackbox([0.0], [pso.K], pso.P, pso.Q, pso.ports)
    return bbox, null_basis
end

magnetostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary::CellComplex{N}, μ⁻::AbstractVector{<:Real},
    Λ⁻::AbstractVector{<:Real}) where {N, K} = magnetostatics_blackbox(m, mesh,
    sources, boundary.cells[1], boundary.cells[2], μ⁻, Λ⁻)

################################################################################
# Common
################################################################################

# solve for (φ, Q) or (A, I)
function solve_statics(bbox::Blackbox, null_basis::AbstractMatrix{<:Real},
    source_magnitudes::AbstractVector{<:Real})
    constrained_bbox = apply_transform(bbox, null_basis)
    field = null_basis * (collect(constrained_bbox.Y[1]) \ collect(constrained_bbox.P * source_magnitudes))
    source = bbox.Y[1] * field
    return field, source
end

# elastance, inductance
impedance_matrix(bbox::Blackbox, null_basis::AbstractMatrix{<:Real}) =
    impedance_matrices(apply_transform(bbox, null_basis))[1]

# capacitance, inv_inductance
admittance_matrix(bbox::Blackbox, null_basis::AbstractMatrix{<:Real}) =
    admittance_matrices(apply_transform(bbox, null_basis))[1]

# E = -dφ, estimate ♯E, a vector field
φ_to_vec_E(m::Metric{N}, mesh::Mesh{N}, φ::AbstractVector{<:Real}) where N =
    DEC.sharp(m, mesh.primal.complex, -DEC.differential_operator(m, mesh, "d", 1, true, φ))

# B = ★dA, estimate ♯B, a dual vector field
A_to_vec_B(m::Metric{N}, mesh::Mesh{N}, A::AbstractVector{<:Real}) where N =
    DEC.sharp(m, mesh.dual.complex, DEC.differential_operator(m, mesh, "★d", 2, true, A))

# Q = ★ρ so ρ = s★Q where s = ★★. Similarly, J = s★I.
function source_density(m::Metric{N}, mesh::Mesh{N, K},
    source::AbstractVector{<:Real}, k::Int) where {N, K}
    pm = DEC.hodge_square_sign(m, K, k)
    return pm * DEC.differential_operator(m, mesh, "★", K-k+1, false, source)
end

function get_material(comp::CellComplex{N}, group::AbstractVector{Simplex{N, K}},
    value::Real, k::Int) where {N, K}
    material = zeros(length(comp.cells[k]))
    for c in DEC.submanifold(comp, group).cells[k]
        material[findfirst(isequal(c), comp.cells[k])] = value
    end
    return material
end

get_material(comp::CellComplex, value::Real, k::Int) = ones(length(comp.cells[k])) * value

const speed_of_light = 299792458.0
const μ₀ = 4π/1e7
const ϵ₀ = 1/(μ₀ * speed_of_light^2)
