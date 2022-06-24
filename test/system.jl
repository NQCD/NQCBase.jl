using NQCBase
using Test
using Unitful
using AtomsBase
using StaticArrays
# NQCDynamics.System("test.xyz")

@test_throws ErrorException System("not_an_xyz_file.txt")


sys = System(joinpath(@__DIR__, "test.xyz"))
@test length(sys) == 8

@testset "AtomsBase" begin

    @test bounding_box(sys) isa SVector
    @test boundary_conditions(sys) isa SVector

    @test position(sys) isa Vector{<:SVector}
    @test velocity(sys) isa Vector{<:SVector}
    @test atomic_symbol(sys) isa Vector{Symbol}
    @test atomic_mass(sys) isa Vector{<:Unitful.Mass}
    @test atomic_number(sys) isa Vector{<:Integer}

    @test position(sys, 1) isa SVector
    @test velocity(sys, 1) isa SVector
    @test atomic_symbol(sys, 1) isa Symbol
    @test atomic_mass(sys, 1) isa Unitful.Mass
    @test atomic_number(sys, 1) isa Integer
end

NQCBase.save("output.xyz", sys)
NQCBase.load("output.xyz")
# atoms = Atoms("inputfile")
# cell = Cell("inputfile")
# model = MyModel()
# temperature = NoTemperature()

# sim = Simulation{Classical}(atoms, model; cell, temperature)

# configurations = DynamicalDistribution(velocity, position, size(sim)) * PureState(1)

# run_ensemble(sim, tspan, configurations; output)
