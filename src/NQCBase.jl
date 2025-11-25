module NQCBase

using PeriodicTable
using Unitful, UnitfulAtomic
using Requires

# Basic types for atomic structures
include("unit_conversions.jl")
include("atoms.jl")
include("cells.jl")

# Convenience Structure type

"""
    Structure{T}

Structures are storage types to keep atoms, positions, cells and further information in one place. 
Many functions in the NQCD packages use only parts of the structure information for efficiency, e.g. integration routines, but it can be convenient to have it all in one type. 

If you are developing in NQCD, please include multiple dispatch versions of your functions using Structure types where this would be convenient to the user. 
"""
struct Structure{T}
    atoms::NQCBase.Atoms # Atoms object
    positions::AbstractMatrix{T} # Positions matrix in atomic units, each column containing one atom's positions
    cell::AbstractCell # Unit cell object
    info::Dict{String, Any} # Other structure information, e.g. from an ExtXYZ header
end
Structure(atoms::NQCBase.Atoms, positions::AbstractMatrix, cell::AbstractCell) = Structure(atoms, positions, cell, Dict{String, Any}())

# I/O Interfaces
include("io/extxyz.jl")
include("atoms_base.jl")


export Cell
export System
export Trajectory
export Position, Velocity





end
