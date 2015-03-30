using Bravais
using Base.Test

lattice = BravaisLattice([4,6,8])
@test length(lattice) == 4*6*8
@test size(lattice) == (4*6*8,)
# NOTE: assumes hypercubic by default, but does not define nearest neighbors
@test realspace(lattice, [1,4,3]) == [1,4,3]

lattice = TriangularLattice([4,6])
@test length(lattice) == 4*6
@test realspace(lattice, [1,0]) == [1.0, 0]
@test realspace(lattice, [0,1]) == [0.5, sqrt(3)/2]
@test realspace(lattice, 7) == [1.0, 0]

@test length(BravaisLattice([4,6])) == 4*6
@test length(HypercubicLattice([4,6,8,10])) == 4*6*8*10
@test length(TriangularLattice([4,6])) == 4*6
@test length(HoneycombLattice([3,5])) == 3*5*2
@test length(KagomeLattice([8,6])) == 8*6*3

@test_throws BoundsError checkbounds(HypercubicLattice([4,6]), 0)
@test_throws BoundsError checkbounds(HypercubicLattice([4,6]), 4*6+1)

lattices = (
    BravaisLattice([4,6]),
    HypercubicLattice([8]),
    HypercubicLattice([4,6]),
    HypercubicLattice([4,2,6]),
    TriangularLattice([4,6]),
    HoneycombLattice([3,5]),
    KagomeLattice([8,6]),
)

for lattice in lattices
    @test size(lattice) == (length(lattice),)

    last = false
    for (i, site) in enumerate(lattice)
        checkbounds(lattice, i)
        @test lattice[i] == site
        @test site in lattice
        @test findfirst(lattice, site) == i
        if i == length(lattice)
            last = true
        end
    end
    @test last == true


end
