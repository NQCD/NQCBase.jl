using Test
using NQCBase
using AtomsBase
using Unitful, UnitfulAtomic

@testset "PeriodicCell" begin

    # 1D cells
    @test PeriodicCell(hcat(1)) isa PeriodicCell
    @test PeriodicCell(hcat(1.0)) isa PeriodicCell

    # 3D cells
    x = [1.0, 0.0, 0.0]
    y = [0.0, 1.0, 0.0]
    z = [0.0, 0.0, 1.0]
    a = PeriodicCell([x y z])
    @test a isa PeriodicCell
    @test a.inverse == a.vectors

    @testset "apply_cell_boundaries!" begin
        cell = PeriodicCell([1 0 0; 0 1 0; 0 0 1])
        R = rand(3, 4)
        A = copy(R)
        NQCBase.apply_cell_boundaries!(cell, A)
        @test R == A # Check unchanged when inside cell
        A += 2rand(3, 4) # Move atoms out of cell
        NQCBase.apply_cell_boundaries!(cell, A)
        @test all(0 .<= A .<= 1) # Check they're all back in
        A -= 2rand(3, 4) # Move atoms out of cell
        NQCBase.apply_cell_boundaries!(cell, A)
        @test all(0 .<= A .<= 1) # Check they're all back in
    end

    a = PeriodicCell(hcat(1))
    @test a.periodicity == [true]
    set_periodicity!(a, [false])
    @test a.periodicity == [false]

    @testset "evaluate_periodic_distance" begin
        cell = PeriodicCell([1 0 0; 0 1 0; 0 0 1])
        r1 = [0.1, 0.1, 0.1]
        r2 = [0.9, 0.9, 0.9]
        @test sqrt(3*0.2^2) ≈ NQCBase.evaluate_periodic_distance(cell, r1, r2)
        @test sqrt(3*0.2^2) ≈ NQCBase.evaluate_periodic_distance(cell, r2, r1)
        r1 = [0.4, 0.4, 0.4]
        r2 = [0.6, 0.6, 0.6]
        @test sqrt(3*0.2^2) ≈ NQCBase.evaluate_periodic_distance(cell, r1, r2)
        @test sqrt(3*0.2^2) ≈ NQCBase.evaluate_periodic_distance(cell, r2, r1)
    end

    @testset "check_atoms_in_cell" begin
        cell = PeriodicCell([1 0 0; 0 1 0; 0 0 1])
        R = rand(3, 10)
        @test NQCBase.check_atoms_in_cell(cell, R) # All atoms inside cell
        R[3, 5] += 1
        @test !NQCBase.check_atoms_in_cell(cell, R) # Some atoms not inside cell
        R[3, 5] -= 2
        @test !NQCBase.check_atoms_in_cell(cell, R) # Some atoms not inside cell
    end

    @testset "AtomsBase.jl" begin
        cell = PeriodicCell([1 0 0; 0 1 0; 0 0 1])
        box = bounding_box(cell)
        @test box[1] == [1, 0, 0]u"bohr"
        @test box[2] == [0, 1, 0]u"bohr"
        @test box[3] == [0, 0, 1]u"bohr"
        @test boundary_conditions(cell) == [Periodic(), Periodic(), Periodic()]
        @test periodicity(cell) == cell.periodicity
        @test n_dimensions(cell) == 3
    end

end

@testset "InfiniteCell" begin
    @test InfiniteCell() isa InfiniteCell
    @test InfiniteCell{3}() isa InfiniteCell
    cell = InfiniteCell()
    @test NQCBase.apply_cell_boundaries!(cell, rand(2,3)) === nothing
    @test NQCBase.check_atoms_in_cell(cell, rand(2,3)) == true
    @test bounding_box(cell) == [[Inf]]u"bohr"
    @test boundary_conditions(cell) == [DirichletZero()]
    @test periodicity(cell) == [false]
    @test n_dimensions(cell) == 1
end
