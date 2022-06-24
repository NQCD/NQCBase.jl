using Distances: Distances
using StaticArrays: SMatrix, SVector

mutable struct PeriodicCell{S,T<:AbstractFloat,L} <: AbstractMatrix{T}
    vectors::SMatrix{S,S,T,L}
    inverse::SMatrix{S,S,T,L}
    periodicity::SVector{S,Bool}
    function PeriodicCell(vectors::AbstractMatrix, inverse::AbstractMatrix, periodicity::AbstractVector)
        T = promote_type(eltype(vectors), eltype(inverse))
        S = size(vectors, 1)
        L = length(vectors)
        new{S,T,L}(vectors, inverse, periodicity)
    end
end

function PeriodicCell(vectors::AbstractMatrix, periodicity::AbstractVector)
    return PeriodicCell(vectors, inv(vectors), periodicity)
end

function PeriodicCell(vectors::AbstractMatrix)
    S = size(vectors, 1)
    periodicity = SVector{S}(true for _=1:S)
    return PeriodicCell(vectors, periodicity)
end

Base.getindex(cell::PeriodicCell, I::Vararg{Int,2}) = getindex(cell.vectors, I...)
Base.setindex!(cell::PeriodicCell, v, I::Vararg{Int,2}) = setindex!(cell.vectors, v, I...)
Base.eltype(::PeriodicCell{S,T}) where {S,T} = T
Base.size(cell::PeriodicCell) = size(cell.vectors)

Base.:*(cell::PeriodicCell, r::AbstractVector) = cell.vectors * r
Base.:\(cell::PeriodicCell, r::AbstractVector) = cell.inverse * r

function AtomsBase.bounding_box(cell::PeriodicCell{S}) where {S}
    return SVector{S}(vec * u"bohr" for vec in eachcol(cell.vectors))
end

function AtomsBase.boundary_conditions(cell::PeriodicCell{S}) where {S}
    return SVector{S}(
        bc ? AtomsBase.Periodic() : AtomsBase.DirichletZero() for bc in cell.periodicity
    )
end

AtomsBase.periodicity(cell::PeriodicCell) = cell.periodicity
AtomsBase.n_dimensions(::PeriodicCell{S}) where {S} = S

function set_periodicity!(cell::PeriodicCell, periodicity::AbstractVector{<:Bool})
    cell.periodicity = periodicity
end

function set_vectors!(cell::PeriodicCell, vectors::AbstractMatrix)
    cell.vectors = vectors
    cell.inverse = inv(cell.vectors)
end

function apply_cell_boundaries!(cell::PeriodicCell, R::AbstractMatrix)
    @views for r in eachcol(R)
        apply_cell_boundaries!(cell, r)
    end
end

function apply_cell_boundaries!(cell::PeriodicCell{S}, R::AbstractVector) where {S}
    fractional = cell \ R
    corrected = SVector{S}(
        cell.periodicity[j] ? mod(fractional[j], 1) : fractional[j] for j in eachindex(R)
    )
    R .= cell * corrected
end

function check_atoms_in_cell(cell::PeriodicCell, R::AbstractMatrix)::Bool
    @views for r in eachcol(R) # atoms
        fractional = cell \ r
        outside = @. (fractional > 1) | (fractional < 0)
        any(outside) && return false
    end
    return true
end

const periodic_distance = Distances.PeriodicEuclidean([1, 1, 1])
function evaluate_periodic_distance(cell::PeriodicCell, r1::AbstractVector, r2::AbstractVector)
    fractional1 = cell \ r1
    fractional2 = cell \ r2
    return Distances.evaluate(periodic_distance, fractional1, fractional2)
end

struct InfiniteCell{S} end
InfiniteCell() = InfiniteCell{1}()

apply_cell_boundaries!(::InfiniteCell, ::AbstractArray) = nothing
check_atoms_in_cell(::InfiniteCell, r) = true
AtomsBase.bounding_box(::InfiniteCell{S}) where {S} = AtomsBase.infinite_box(S)

function AtomsBase.boundary_conditions(::InfiniteCell{S}) where {S}
    return SVector{S}(AtomsBase.DirichletZero() for _ in 1:S)
end

AtomsBase.periodicity(::InfiniteCell{S}) where {S} = SVector{S}(false for _ in 1:S)
AtomsBase.n_dimensions(::InfiniteCell{S}) where {S} = S

