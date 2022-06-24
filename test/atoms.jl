using Test
using NQCBase

@testset "Atoms constructors" begin
    @test Atoms([:C, :H]) isa Atoms{Float64}
    @test Atoms(:C, :H) isa Atoms{Float64}
    @test Atoms{Float64}([:C, :H]) isa Atoms{Float64}
    @test Atoms{Float32}([:C, :H]) isa Atoms{Float32}
    @test Atoms{Float64}(:H) isa Atoms{Float64}
    @test Atoms{Float32}(:C, :H) isa Atoms{Float32}
end

@testset "Atoms" begin
    atoms = Atoms(:C, :H)
    @test length(atoms) == 2
    @test atoms.numbers == [6, 1]
    @test atoms.types == [:C, :H]
    @test atoms.masses isa Vector
    @test masses(atoms) == atoms.masses
    @test range(atoms) == 1:2
    @test atoms[1] isa Atoms{Float64}
    @test atoms[1:2] isa Atoms{Float64}
end

@testset "Particles constructors" begin
    @test Particles([1, 2, 3]) isa Particles{Float64}
    @test Particles([1.0f0, 2, 3]) isa Particles{Float32}
    @test Particles([1.0, 2, 3]) isa Particles{Float64}
    @test Particles(1.0, 2.0, 3.0) isa Particles{Float64}
    @test Particles(1.0f0, 2.0f0, 3.0f0) isa Particles{Float32}
    @test Particles{Float64}([1, 2, 3]) isa Particles{Float64}
    @test Particles{Float32}([1, 2, 3]) isa Particles{Float32}
    @test Particles{Float64}(1, 2, 3) isa Particles{Float64}
    @test Particles{Float32}(1, 2, 3) isa Particles{Float32}
    @test Atoms([1, 2, 3]) isa Particles{Float64}
end

@testset "Particles" begin
    atoms = Particles(1, 2)
    @test length(atoms) == 2
    @test atoms.masses == [1.0, 2.0]
    @test masses(atoms) == [1.0, 2.0]
    @test range(atoms) == 1:2
    @test atoms[1] isa Particles{Float64}
    @test atoms[1:2] isa Particles{Float64}
end
