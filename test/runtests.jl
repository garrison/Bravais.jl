using Bravais
using Base.Test

lattice = BravaisLattice([4,6,8])
@test length(lattice) == 4*6*8
@test size(lattice) == (4*6*8,)
@test realspace(lattice, [1,4,3]) == [1,4,3]

lattice = HypercubicLattice([4,6,8]) # what is the difference btw this and regular BravaisLattice?

lattice = TriangularLattice([4,6])
@test length(lattice) == 4*6
@test realspace(lattice, [1,0]) == [1.0, 0]
@test realspace(lattice, [0,1]) == [0.5, sqrt(3)/2]
@test realspace(lattice, 2) == [1.0, 0]

lattice = HoneycombLattice([3,5])

lattice = HypercubicLattice([4,6])
for (i, site) in enumerate(lattice)
    @test lattice[i] == site
end
