
using Test
using NQCBase

atoms = Atoms([:H, :C, :O, :N])
cell = PeriodicCell(rand(3, 3) .* 10)
R = rand(3, 4) .* 10
structure = NQCBase.Structure(atoms, R, cell)

@testset "to/from_extxyz_dict (Atoms, Positions, Cell)" begin
    dict = NQCBase.to_extxyz_dict(atoms, R, cell) # Test conversion with atoms, R, cell and Structure methods. 
    @test dict["cell"] ≈ au_to_ang.(permutedims(cell.vectors, (2,1)))
    converted_structure = NQCBase.from_extxyz_dict(dict)
    @test converted_structure.cell.vectors ≈ cell.vectors
    @test converted_structure.cell.inverse ≈ cell.inverse
end
@testset "to/from_extxyz_dict (Structure)" begin
    dict = NQCBase.to_extxyz_dict(structure) # Test conversion with atoms, R, cell and Structure methods. 
    @test dict["cell"] ≈ au_to_ang.(permutedims(structure.cell.vectors, (2,1)))
    converted_structure = NQCBase.from_extxyz_dict(dict)
    @test converted_structure.cell.vectors ≈ structure.cell.vectors
    @test converted_structure.cell.inverse ≈ structure.cell.inverse
end

@testset "Single frame" begin
    file_buffer = "output.xyz"
    write_extxyz(file_buffer, atoms, R, cell) # Save a 1-structure file
    structure = read_extxyz(file_buffer) |> first # Load the 1-structure file
    @test structure.atoms == atoms
    @test structure.cell.vectors ≈ cell.vectors
    @test structure.cell.inverse ≈ cell.inverse
    @test structure.cell.periodicity ≈ cell.periodicity
    @test all(isapprox.(structure.positions, R, atol = 1e-8))
end

@testset "Multiple frames" begin
    fb = "output.xyz"
    atoms = Atoms([:H, :C, :O, :N])
    cell = PeriodicCell(rand(3, 3) .* 10)
    R = [rand(3, 4) .* 10 for _=1:100]
    write_extxyz(fb, atoms, R, cell)
    structures = read_extxyz(fb)
    for (i,structure) in enumerate(structures)
        @test structure.atoms == atoms
        @test structure.cell.vectors ≈ cell.vectors
        @test structure.cell.inverse ≈ cell.inverse
        @test structure.cell.periodicity ≈ cell.periodicity
        @test all(isapprox.(structure.positions, R[i], atol = 1e-8))
    end
end

rm("output.xyz")
