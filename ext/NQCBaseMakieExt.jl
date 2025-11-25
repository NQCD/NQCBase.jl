module NQCBaseMakieExt

using Makie
using PeriodicTable
using UnitfulAtomic
using Unitful
using NQCBase

# Default atom colours
const default_atom_fills = Dict(
    [Symbol(el.symbol), parse(Makie.Colors.Colorant, el.cpk_hex) for el in elements]
)
# Empirical atomic radii (pm) from Wikipedia data page:
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
const default_atomic_radii_pm = Dict{Symbol, Float64}(
    :H => 25,
    :He => 120,
    :Li => 145,
    :Be => 105,
    :B => 85,
    :C => 70,
    :N => 65,
    :O => 60,
    :F => 50,
    :Ne => 160,
    :Na => 180,
    :Mg => 150,
    :Al => 125,
    :Si => 110,
    :P => 100,
    :S => 100,
    :Cl => 100,
    :Ar => 71,
    :K => 220,
    :Ca => 180,
    :Sc => 160,
    :Ti => 140,
    :V => 135,
    :Cr => 140,
    :Mn => 140,
    :Fe => 140,
    :Co => 135,
    :Ni => 135,
    :Cu => 135,
    :Zn => 135,
    :Ga => 130,
    :Ge => 125,
    :As => 115,
    :Se => 115,
    :Br => 115,
    :Rb => 235,
    :Sr => 200,
    :Y => 180,
    :Zr => 155,
    :Nb => 145,
    :Mo => 145,
    :Tc => 135,
    :Ru => 130,
    :Rh => 135,
    :Pd => 140,
    :Ag => 160,
    :Cd => 155,
    :In => 155,
    :Sn => 145,
    :Sb => 145,
    :Te => 140,
    :I => 140,
    :Cs => 260,
    :Ba => 215,
    :La => 195,
    :Ce => 185,
    :Pr => 185,
    :Nd => 185,
    :Pm => 185,
    :Sm => 185,
    :Eu => 185,
    :Gd => 180,
    :Tb => 175,
    :Dy => 175,
    :Ho => 175,
    :Er => 175,
    :Tm => 175,
    :Yb => 175,
    :Lu => 175,
    :Hf => 155,
    :Ta => 145,
    :W => 135,
    :Re => 135,
    :Os => 130,
    :Ir => 135,
    :Pt => 135,
    :Au => 135,
    :Hg => 150,
    :Tl => 190,
    :Pb => 180,
    :Bi => 160,
    :Po => 190,
    :Ra => 215,
    :Ac => 195,
    :Th => 180,
    :Pa => 180,
    :U => 175,
    :Np => 175,
    :Pu => 175,
    :Am => 175,
    :Cm => 176,
 )

bettersphere = Makie.GeometryBasics.mesh(
     Makie.GeometryBasics.Tesselation(
         Makie.GeometryBasics.Sphere(Point3f(0.0,0.0,0.0), 1.0),
         128
     )
 )

@recipe AtomicStructure (structure) begin
    # Pass default atom fills and atomic radii to allow modification. 
    atom_fills = default_atom_fills
    atomic_radii_Å = Dict([el, ustrip(uconvert(u"Å", rad * u"pm")) for (el, rad) in default_atomic_radii_pm]) # Convert to Angstrom
    # Default edge width
    strokewidth = 1.5
    # Default higher-quality sphere model for raster exports. 
    marker = bettersphere
    Makie.mixin_generic_plot_attributes()...
end

# Plot in 3D by default
Makie.args_preferred_axis(::Type{<: AtomicStructure}) = Makie.Axis3

function Makie.plot!(
    structure::AtomicStructure{<:Tuple{NQCBase.Structure}},
)
    input_nodes = [:converted_1]
    output_nodes = [:positions, :markersizes, :interior_colors, :edge_colors, ]
    
    map!(structure.attributes, input_nodes, output_nodes) do nqcd_structure
        pos = ustrip.(auconvert.(u"Å", nqcd_structure.positions))
        positions = Point3f.(eachcol(positions))
        markersizes = [get()]
        return positions
    end
    
    
    
    return structure
end