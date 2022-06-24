module NQCBase

using PeriodicTable
using Unitful, UnitfulAtomic
using Requires
using AtomsBase: AtomsBase

include("unit_conversions.jl")
include("atoms.jl")
export Atoms
export Particles
export masses

include("cells.jl")
export PeriodicCell
export InfiniteCell
export set_periodicity!
export set_vectors!

include("io/extxyz.jl")
export Cell

include("system.jl")
export System

function __init__()
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" @eval include("io/ase.jl")
end

end
