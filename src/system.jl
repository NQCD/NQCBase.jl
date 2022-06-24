using AtomsBase: AtomsBase
using ExtXYZ: ExtXYZ

struct NoTemperature end
function (temperature::NoTemperature)(time=0)
    throw(error("""
    The temperature of the system is `NoTemperature` but you have attempted to access the temperature. \
    Be sure to set the temperature when you initialise the `System`.
    """
    ))
end

struct FixedTemperature{T}
    value::T
end
(temperature::FixedTemperature)(time=0) = austrip(temperature.value)

struct TimeDependentTemperature{T}
    temperature_function::T
end
function (temperature::TimeDependentTemperature)(time=0)
    return austrip(temperature.temperature_function(time))
end

struct ElectronPhononTemperature{E,P}
    electron_temperature::E
    phonon_temperature::P
end
function (temperature::ElectronPhononTemperature)(time=0)
    return austrip(temperature.phonon_temperature(time))
end

struct System{S,L<:Unitful.Length,V,A,C} <: AtomsBase.AbstractSystem{S}
    atoms::A
    cell::C
    position::Vector{SVector{S,L}}
    velocity::Vector{SVector{S,V}}
    data::Dict{Symbol,Any}
    function System(atoms, cell, position, velocity, data)
        S = AtomsBase.n_dimensions(cell)
        A = typeof(atoms)
        C = typeof(cell)
        L = eltype(first(position))
        V = eltype(first(velocity))
        new{S,L,V,A,C}(atoms, cell, position, velocity, data)
    end
end

function Base.show(io::IO, sys::System)
    print(io, "System")
    AtomsBase.show_system(io, sys)
end

Base.length(sys::System) = length(sys.atoms)

function System(filename)
    isxyzfile(filename) || throw(error("File must be a `.xyz`."))

    atoms = Atoms(filename)
    cell = Cell(filename)
    S = AtomsBase.n_dimensions(cell)
    dict = ExtXYZ.read_frame(filename)
    position = positions_from_extxyz_dict(dict)
    if haskey(dict, "vel")
        velocity = velocities_from_extxyz_dict(dict)
    else
        velocity = [zeros(SVector{S,typeof(1.0u"Å/fs")}) for _ in eachindex(position)]
    end

    data = Dict{Symbol,Any}()
    for key in keys(dict["info"])
        data[Symbol(key)] = dict["info"][key]
    end

    return System(atoms, cell, position, velocity, data)
end

function positions_from_extxyz_dict(dict::Dict{String,Any})
    return [SVector{length(col)}(col)*u"Å" for col in eachcol(dict["arrays"]["pos"])]
end

function velocities_from_extxyz_dict(dict::Dict{String,Any})
    return [SVector{length(col)}(col)*u"Å/fs" for col in eachcol(dict["arrays"]["vel"])]
end

AtomsBase.bounding_box(sys::System) = AtomsBase.bounding_box(sys.cell)
AtomsBase.boundary_conditions(sys::System) = AtomsBase.boundary_conditions(sys.cell)

AtomsBase.position(sys::System) = sys.position
AtomsBase.velocity(sys::System) = sys.velocity
AtomsBase.atomic_symbol(sys::System) = AtomsBase.atomic_symbol(sys.atoms)
AtomsBase.atomic_mass(sys::System) = AtomsBase.atomic_mass(sys.atoms)
AtomsBase.atomic_number(sys::System) = AtomsBase.atomic_number(sys.atoms)

AtomsBase.position(sys::System, i) = sys.position[i]
AtomsBase.velocity(sys::System, i) = sys.velocity[i]
AtomsBase.atomic_symbol(sys::System, i) = AtomsBase.atomic_symbol(sys.atoms, i)
AtomsBase.atomic_mass(sys::System, i) = AtomsBase.atomic_mass(sys.atoms, i)
AtomsBase.atomic_number(sys::System, i) = AtomsBase.atomic_number(sys.atoms, i)

function save(filename, sys::System)
    isxyzfile(filename) || throw(error("File must be a `.xyz`."))
    ExtXYZ.write_frame(filename, to_extxyz_dict(sys))
end

function to_extxyz_dict(sys::System)

    dict = Dict{String, Any}()
    dict["N_atoms"] = length(sys)

    dict["pbc"] = Array(AtomsBase.periodicity(sys))

    dict["cell"] = Array(permutedims(ustrip.(uconvert.(u"Å", reduce(hcat, AtomsBase.bounding_box(sys))))))

    dict["arrays"] = Dict{String,Any}()
    dict["arrays"]["species"] = String.(AtomsBase.atomic_symbol(sys))
    dict["arrays"]["pos"] = ustrip.(uconvert.(u"Å", reduce(hcat, AtomsBase.position(sys))))
    dict["arrays"]["vel"] = ustrip.(uconvert.(u"Å/fs", reduce(hcat, AtomsBase.velocity(sys))))

    dict["info"] = Dict{String,Any}()
    for key in keys(sys.data)
        dict["info"][String(key)] = sys.data[key]
    end

    return dict
end
