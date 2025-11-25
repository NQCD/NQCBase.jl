module NQCBaseMakieExt

using Makie
using PeriodicTable
using UnitfulAtomic

# Default atom colours
const default_atom_fills = Dict(
    [Symbol(el.symbol), parse(Makie.Colors.Colorant, el.cpk_hex) for el in elements]
)
# 
const default_atomic_radii = Dict() # ToDo: Fill these and decide what atomic radii to use for display


end