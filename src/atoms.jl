
"""
    Atoms{T<:AbstractFloat}

Basic atomic parameters: element symbols, numbers and masses

Masses are converted to atomic units.
Constructed using either element symbols or masses.

```jldoctest
julia> Atoms(:H)
Atoms{Float64}([:H], [1], [1837.4715941070515])

julia> Atoms([:H, :H, :H, :C])
Atoms{Float64}([:H, :H, :H, :C], [1, 1, 1, 6], [1837.4715941070515, 1837.4715941070515, 1837.4715941070515, 21894.713607956142])

julia> Atoms([100, 200])
Atoms{Float64}([:X, :X], [0, 0], [100.0, 200.0])
```
"""
struct Atoms{T<:AbstractFloat}
    types::Vector{Symbol}
    numbers::Vector{Int}
    masses::Vector{T}
    function Atoms{T}(types::AbstractVector) where {T<:AbstractFloat}
        numbers = [element.number for element in elements[types]]
        masses = [austrip(element.atomic_mass) for element in elements[types]]
        new{T}(types, numbers, masses)
    end
end
Atoms(types::AbstractVector) = Atoms{Float64}(types)
Atoms{T}(types...) where {T<:AbstractFloat} = Atoms{T}(collect(types))
Atoms(types...) = Atoms(collect(types))

Base.getindex(A::Atoms{T}, i) where {T} = Atoms{T}(A.types[i])

function Base.:(==)(a::Atoms, b::Atoms)
    return (a.types == b.types) && (a.numbers == b.numbers) && (a.masses ≈ b.masses)
end

function Base.show(io::IO, atoms::Atoms{T}) where {T}
    print(io, "Atoms{$T}($(atoms.types))")
end

AtomsBase.atomic_symbol(atoms::Atoms) = atoms.types
AtomsBase.atomic_mass(atoms::Atoms) = atoms.masses * u"me_au"
AtomsBase.atomic_number(atoms::Atoms) = atoms.numbers

AtomsBase.atomic_symbol(atoms::Atoms, i) = atoms.types[i]
AtomsBase.atomic_mass(atoms::Atoms, i) = atoms.masses[i] * u"me_au"
AtomsBase.atomic_number(atoms::Atoms, i) = atoms.numbers[i]

struct Particles{T<:AbstractFloat}
    masses::Vector{T}
    Particles{T}(masses::AbstractVector) where {T<:AbstractFloat} = new{T}(austrip.(masses))
end
Particles(masses::AbstractVector{T}) where {T<:AbstractFloat} = Particles{T}(masses)
Particles(masses::AbstractVector) = Particles{Float64}(masses)
Particles{T}(masses...) where {T<:AbstractFloat} = Particles{T}(collect(masses))
Particles(masses...) = Particles(collect(masses))
Atoms(masses::AbstractVector{<:Real}) = Particles(masses)

Base.getindex(A::Particles{T}, i) where {T} = Particles{T}(A.masses[i])

function Base.:(==)(a::Particles, b::Particles)
    return (a.masses ≈ b.masses)
end

AtomsBase.atomic_mass(atoms::Particles) = atoms.masses * u"me_au"
AtomsBase.atomic_mass(atoms::Particles, i) = atoms.masses[i] * u"me_au"
AtomsBase.atomic_symbol(::Particles) = throw(error("Particles do not have atomic symbols, use `Atoms` instead."))
AtomsBase.atomic_symbol(::Particles, _) = throw(error("Particles do not have atomic symbols, use `Atoms` instead."))
AtomsBase.atomic_number(::Particles) = throw(error("Particles do not have atomic numbers, use `Atoms` instead."))
AtomsBase.atomic_number(::Particles, _) = throw(error("Particles do not have atomic numbers, use `Atoms` instead."))

const AtomTypes = Union{Atoms, Particles}

Base.length(atoms::AtomTypes) = length(atoms.masses)
Base.range(atoms::AtomTypes) = range(1; length=length(atoms))
Base.IndexStyle(::AtomTypes) = IndexLinear()
masses(atoms::AtomTypes) = atoms.masses
