using DiscreteExteriorCalculus: Metric, Simplex, Cell, CellComplex, Mesh, sharp,
    subcomplex, orient!, differential_operator_sequence, differential_operator,
    zero_constraint, constant_constraint
using SparseArrays: sparse, spdiagm
using AdmittanceModels: Blackbox, PSOModel, apply_transform, sparse_nullbasis,
    impedance_matrices, admittance_matrices

export electrostatics_blackbox
"""
    electrostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
        sources::AbstractVector{<:AbstractVector{Cell{N}}},
        boundary_points::AbstractVector{Cell{N}},
        ϵ::AbstractVector{<:Real}) where {N, K}
    electrostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
        sources::AbstractVector{<:AbstractVector{Cell{N}}},
        boundary::CellComplex{N},
        ϵ::AbstractVector{<:Real}) where {N, K}

Each source is a vector of cells with `K == 1` representing a conductor. Discretize the
Poisson equation to find a Blackbox that maps charge to voltage as well as a sparse
constraint matrix that asserts that the voltage is constant on each source.
"""
function electrostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary_points::AbstractVector{Cell{N}},
    ϵ::AbstractVector{<:Real}) where {N, K}
    k = 1
    ★★, d₁, ★, d₀ = differential_operator_sequence(m, mesh, "★★d★d", k, true)
    Y = -★★ * d₁ * ★ * spdiagm(0 => ϵ) * d₀
    comp = mesh.primal.complex
    row_inds = [findfirst(isequal(s[1]), comp.cells[k]) for s in sources]
    n = length(sources)
    P = sparse(row_inds, 1:n, ones(n), length(comp.cells[k]), n)
    bbox = Blackbox([0.0], [Y], P, P, collect(1:n))
    conductors = filter(s -> length(s) > 1, sources)
    constraint = vcat(zero_constraint(comp, boundary_points, k),
        [constant_constraint(comp, s, k) for s in conductors]...)
    return bbox, sparse_nullbasis(constraint)
end

electrostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary::CellComplex{N},
    ϵ::AbstractVector{<:Real}) where {N, K} =
    electrostatics_blackbox(m, mesh, sources, boundary.cells[1], ϵ)

export get_charge_source
"""
    get_charge_source(comp::CellComplex{N},
        group::AbstractVector{Simplex{N, K}}) where {N, K}

Return a vector of cells representing the points in the given simplices.
"""
get_charge_source(comp::CellComplex{N}, group::AbstractVector{Simplex{N, K}}) where {N, K} =
    subcomplex(comp, group).cells[1]

export magnetostatics_blackbox
"""
    magnetostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
        sources::AbstractVector{<:AbstractVector{Cell{N}}},
        boundary_points::AbstractVector{Cell{N}},
        boundary_edges::AbstractVector{Cell{N}},
        μ⁻::AbstractVector{<:Real}, Λ⁻::AbstractVector{<:Real}) where {N, K}
    magnetostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
        sources::AbstractVector{<:AbstractVector{Cell{N}}},
        boundary::CellComplex{N}, μ⁻::AbstractVector{<:Real},
        Λ⁻::AbstractVector{<:Real}) where {N, K}

Each source is a vector of cells with `K == 2` representing a current line. In order to get
correct results, the current specified must have divergence 0. Discretize the Maxwell
equations in the Coulomb gauge to find a Blackbox mapping current to the line integral of
the vector potential (for a closed loop this is magnetic flux) as well as a sparse
constraint matrix that asserts that the vector potential is orthogonal to the boundary
and has divergence 0 in the interior.
"""
function magnetostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary_points::AbstractVector{Cell{N}},
    boundary_edges::AbstractVector{Cell{N}},
    μ⁻::AbstractVector{<:Real}, Λ⁻::AbstractVector{<:Real}) where {N, K}
    pso, null_basis = electrodynamics_pso(m, mesh, sources, boundary_points,
        boundary_edges, μ⁻, Λ⁻, zeros(size(Λ⁻)), ones(size(Λ⁻)))
    bbox = Blackbox([0.0], [pso.K], pso.P, pso.Q, pso.ports)
    return bbox, null_basis
end

magnetostatics_blackbox(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary::CellComplex{N}, μ⁻::AbstractVector{<:Real},
    Λ⁻::AbstractVector{<:Real}) where {N, K} = magnetostatics_blackbox(m, mesh,
    sources, boundary.cells[1], boundary.cells[2], μ⁻, Λ⁻)

export get_current_source
"""
    get_current_source(comp::CellComplex{N},
        group::AbstractVector{Simplex{N, K}}) where {N, K}

Return a vector of cells representing the lines in the given simplices.
"""
get_current_source(comp::CellComplex{N}, group::AbstractVector{Simplex{N, K}}) where {N, K} =
    orient!(subcomplex(comp, group).cells[2])

export electrodynamics_pso
"""
    electrodynamics_pso(m::Metric{N}, mesh::Mesh{N, K},
        sources::AbstractVector{<:AbstractVector{Cell{N}}},
        boundary_points::AbstractVector{Cell{N}},
        boundary_edges::AbstractVector{Cell{N}},
        μ⁻::AbstractVector{<:Real}, Λ⁻::AbstractVector{<:Real},
        σ::AbstractVector{<:Real}, ϵ::AbstractVector{<:Real}) where {N, K}
    electrodynamics_pso(m::Metric{N}, mesh::Mesh{N, K},
        sources::AbstractVector{<:AbstractVector{Cell{N}}},
        boundary::CellComplex{N},
        μ⁻::AbstractVector{<:Real}, Λ⁻::AbstractVector{<:Real},
        σ::AbstractVector{<:Real}, ϵ::AbstractVector{<:Real}) where {N, K}

Each source is a vector of cells with `K == 2` representing a current line. In order to get
correct results, the current specified must have divergence 0. Discretize the Maxwell
equations in the Coulomb gauge to find a Positive Second Order model mapping current to
voltage as well as a sparse constraint matrix that asserts that the vector potential is
orthogonal to the boundary and has divergence 0 in the interior.
"""
function electrodynamics_pso(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary_points::AbstractVector{Cell{N}},
    boundary_edges::AbstractVector{Cell{N}},
    μ⁻::AbstractVector{<:Real}, Λ⁻::AbstractVector{<:Real},
    σ::AbstractVector{<:Real}, ϵ::AbstractVector{<:Real}) where {N, K}
    k = 2
    ★★, d₁, ★, d₀ = differential_operator_sequence(m, mesh, "★★d★d", k, true)
    K₀ = ★★ * d₁ * ★ * spdiagm(0 => μ⁻) * d₀
    d, ★ = differential_operator_sequence(m, mesh, "d★", k, true)
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

electrodynamics_pso(m::Metric{N}, mesh::Mesh{N, K},
    sources::AbstractVector{<:AbstractVector{Cell{N}}},
    boundary::CellComplex{N},
    μ⁻::AbstractVector{<:Real}, Λ⁻::AbstractVector{<:Real},
    σ::AbstractVector{<:Real}, ϵ::AbstractVector{<:Real}) where {N, K} =
    electrodynamics_pso(m, mesh, sources, boundary.cells[1], boundary.cells[2], μ⁻, Λ⁻, σ, ϵ)

export solve_statics
"""
    solve_statics(bbox::Blackbox, null_basis::AbstractMatrix{<:Real},
        source_magnitudes::AbstractVector{<:Real})

Given a Blackbox and constraint matrix from `electrostatics_blackbox` or
`magnetostatics_blackbox`, return the scalar potential and charge or vector potential and
current, respectively.
"""
function solve_statics(bbox::Blackbox, null_basis::AbstractMatrix{<:Real},
    source_magnitudes::AbstractVector{<:Real})
    constrained_bbox = apply_transform(bbox, null_basis)
    field = null_basis * (collect(constrained_bbox.Y[1]) \ collect(constrained_bbox.P * source_magnitudes))
    source = bbox.Y[1] * field
    return field, source
end

export impedance_matrix
"""
    impedance_matrix(bbox::Blackbox, null_basis::AbstractMatrix{<:Real})

Given a Blackbox and constraint matrix from `electrostatics_blackbox` or
`magnetostatics_blackbox`, return the elastance matrix or inductance matrix, respectively.
"""
impedance_matrix(bbox::Blackbox, null_basis::AbstractMatrix{<:Real}) =
    impedance_matrices(apply_transform(bbox, null_basis))[1]

export admittance_matrix
"""
    admittance_matrix(bbox::Blackbox, null_basis::AbstractMatrix{<:Real})

Given a Blackbox and constraint matrix from `electrostatics_blackbox` or
`magnetostatics_blackbox`, return the capacitance matrix or inverse inductance matrix,
respectively.
"""
admittance_matrix(bbox::Blackbox, null_basis::AbstractMatrix{<:Real}) =
    admittance_matrices(apply_transform(bbox, null_basis))[1]

export electric_field
"""
    electric_field(m::Metric{N}, mesh::Mesh{N}, φ::AbstractVector{<:Real}) where N

Given the scalar potential from `solve_statics`, approximate the electric field at each
primal vertex.
"""
electric_field(m::Metric{N}, mesh::Mesh{N}, φ::AbstractVector{<:Real}) where N =
    sharp(m, mesh.primal.complex, -differential_operator(m, mesh, "d", 1, true, φ))

export magnetic_field
"""
    magnetic_field(m::Metric{N}, mesh::Mesh{N}, A::AbstractVector{<:Real}) where N

Given the scalar potential from `solve_statics`, approximate the magnetic field at each
dual vertex.
"""
magnetic_field(m::Metric{N}, mesh::Mesh{N}, A::AbstractVector{<:Real}) where N =
    sharp(m, mesh.dual.complex, differential_operator(m, mesh, "★d", 2, true, A))

export source_density
"""
    source_density(m::Metric{N}, mesh::Mesh{N, K},
        source::AbstractVector{<:Real}, k::Int) where {N, K}

Given the charge or current from `solve_statics`, find the charge density or current
density, respectively.
"""
function source_density(m::Metric{N}, mesh::Mesh{N, K},
    source::AbstractVector{<:Real}, k::Int) where {N, K}
    return differential_operator(m, mesh, "★★★", K-k+1, false, source)
end

export get_material
"""
    get_material(comp::CellComplex{N}, group::AbstractVector{Simplex{N, K}},
        value::Real, k::Int) where {N, K}

Return a vector representing the material distribution that is 0 everywhere except on the
simplices in `group` where it is `value`.
"""
function get_material(comp::CellComplex{N}, group::AbstractVector{Simplex{N, K}},
    value::Real, k::Int) where {N, K}
    material = zeros(length(comp.cells[k]))
    for c in subcomplex(comp, group).cells[k]
        material[findfirst(isequal(c), comp.cells[k])] = value
    end
    return material
end

"""
    get_material(comp::CellComplex, value::Real, k::Int)

Return a vector representing the material distribution that is `value` everywhere.
"""
get_material(comp::CellComplex, value::Real, k::Int) = ones(length(comp.cells[k])) * value

const speed_of_light = 299792458.0
const μ₀ = 4π/1e7
const ϵ₀ = 1/(μ₀ * speed_of_light^2)
