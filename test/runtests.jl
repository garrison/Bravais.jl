using Bravais
using Test

import StaticArrays
using LinearAlgebra

debug = false

let lattice = BravaisLattice{3}([4,6,8])
    @test length(lattice) == 4*6*8
    @test size(lattice) == (4*6*8,)

    # NOTE: assumes hypercubic by default, but BravaisLattice does not
    # define any concept of nearest neighbors
    @test realspace(lattice, [1,4,3]) == [1,4,3]
end

let lattice = TriangularLattice([4,6])
    @test length(lattice) == 4*6
    @test realspace(lattice, [1,0]) == [1.0, 0]
    @test realspace(lattice, [0,1]) == [0.5, sqrt(3)/2]
    @test realspace(lattice, 7) == [1.0, 0]
end

@test length(BravaisLattice(StaticArrays.@SVector[4,6])) == 4*6
@test length(HypercubicLattice{4}([4,6,8,10])) == 4*6*8*10
@test length(TriangularLattice([4,6])) == 4*6
@test length(HoneycombLattice([3,5])) == 3*5*2
@test length(KagomeLattice([8,6])) == 8*6*3

@test_throws BoundsError checkbounds(SquareLattice([4,6]), 0)
@test_throws BoundsError checkbounds(SquareLattice([4,6]), 4*6+1)

lattices = Any[
    "Basic BravaisLattice",
    @inferred(BravaisLattice{2}([4,6])),
    @inferred(BravaisLattice(StaticArrays.@SVector[4,6])),

    "1D and effectively 1D",
    @inferred(ChainLattice([8])),
    @inferred(CubicLattice([1, 1, 8])),
    @inferred(CubicLattice([1, 8, 1])),

    "1D Twisted",
    @inferred(ChainLattice([6], Diagonal([6]), [0//1])),
    @inferred(ChainLattice([6], Diagonal([6]), [1//2])),
    @inferred(ChainLattice([6], Diagonal([6]), [1//5])),

    "2D",
    @inferred(SquareLattice([4,6])),
    @inferred(SquareLattice([6,4])),

    "2D twisted",
    @inferred(SquareLattice([4,6], Diagonal([4,6]), [1//7, 3//13])),
    @inferred(SquareLattice([6,4], Diagonal([6,4]), [1//7, 3//13])),

    "2D twisted cylinder",
    @inferred(SquareLattice([6,4], Diagonal([6,0]), [1//7, 0//1])),

    "2D twisted helix",
    @inferred(SquareLattice([6,4], [6 0; 1 4], [1//7, 3//13])),

    "3D like",
    @inferred(CubicLattice([4,2,6])),
    @inferred(HypercubicLattice{4}([4,2,6,1])),
    @inferred(HypercubicLattice(StaticArrays.@MVector[4,2,6,1])),

    "3D open twisted helix",
    @inferred(CubicLattice([4,6,8], [4 0 0; 0 6 0; 1 0 8], [1//7, 0//1, 3//13])),

    "Other lattices",
    @inferred(TriangularLattice([4,6])),
    @inferred(HoneycombLattice([3,5])),
    @inferred(KagomeLattice([8,6])),
]

for lattice in lattices
    if isa(lattice, AbstractString)
        debug && println(lattice)
        continue
    end

    len = @inferred length(lattice)
    @test size(lattice) == (len,)
    @test eachindex(lattice) == 1:len

    @test_throws BoundsError lattice[0]
    @test_throws BoundsError lattice[len+1]

    @static if VERSION >= v"1.2.0-DEV"
        # Relies on https://github.com/JuliaLang/julia/pull/27516
        @test @inferred(Nothing, findfirst(isequal(maxcoords(lattice)), lattice)) === nothing
    else
        @test findfirst(isequal(maxcoords(lattice)), lattice) === nothing
    end
    @inferred(in(@inferred(getindex(lattice, 1)),lattice))

    last = false
    for (i, site) in enumerate(lattice)
        checkbounds(lattice, i)
        @test lattice[i] == site
        @test site in lattice
        @test findfirst(isequal(site), lattice) == i
        @test last == false
        if i == length(lattice)
            last = true
        end

        @inferred wraparound(lattice, site) # FIXME: actually test something here

        @inferred realspace(lattice, site) # FIXME: actually test something here
        realspace(lattice, site, ones(Int, ndimensions(lattice))) # FIXME: actually test something here
    end
    @test last == true

    @inferred dimensions(lattice)
    d = @inferred ndimensions(lattice)
    η = @inferred twist(lattice)
    M = @inferred repeater(lattice)

    @test ishelical(lattice) == !isdiag(repeater(lattice))

    @test isbravais(lattice) == (lattice === bravais(lattice))

    # Check that translating multiple times, in total by an $\mathbf(A)_i$
    # vector, brings us back to where we started (and picking up a phase due
    # to $η_i$).
    sz = length(lattice)
    for i in 1:d
        if M[i,i] == 0
            @test_throws ArgumentError LatticeTranslationCache(lattice, i)
            continue
        end

        # Initialize our array of each index, along with a beginning
        # phase of zero
        indices = Tuple{Int, Rational{Int}}[(z, 0//1) for z in 1:sz]

        # Translate as many times as we need to in each dimension
        for j in 1:d
            if M[i,j] != 0
                ltrc = LatticeTranslationCache(lattice, j)
                for z in 1:M[i,j]
                    new_indices = Tuple{Int, Rational{Int}}[]
                    sizehint!(new_indices, sz)
                    for (s, old_η) in indices
                        newidx, wrap_η = @inferred translateη(lattice, s, j)
                        @test (newidx, wrap_η) == @inferred translateη(ltrc, s)
                        push!(new_indices, (newidx, wrap_η + old_η))
                    end
                    indices = new_indices
                end
            end
        end

        # assert that everything is as expected (i.e., overall we got the
        # identity element while picking up some phase)
        for z in 1:sz
            @test indices[z] == (z, η[i])
        end
    end

    if !isbravais(lattice)
        continue
    end

    @test transpose(primvecs(lattice)) * recivecs(lattice) ≈ Matrix(2π * I, ndimensions(lattice), ndimensions(lattice))

    # Check the momenta across the boundary conditions
    n_k_idx = @inferred nmomenta(lattice)
    @test 1:n_k_idx == @inferred eachmomentumindex(lattice)
    for site in lattice
        for k_idx in 1:n_k_idx
            for i in 1:d
                if M[i,i] != 0
                    site2 = site + vec(M[i,:])
                    exp1 = exp(im * (kdotr(lattice, k_idx, site) + 2π * η[i]))
                    exp2 = exp(im * @inferred(kdotr(lattice, k_idx, site2)))
                    idx, η_wrap = @inferred wraparoundη(lattice, site2)
                    exp3 = exp(im * (kdotr(lattice, k_idx, idx) + 2π * η_wrap))
                    @test exp1 ≈ exp2
                    @test exp1 ≈ exp3
                    for charge in (1, 3)
                        k_total = @inferred momentum(lattice, k_idx, charge)
                        exp1 = exp(im * (kdotr(lattice, k_total, site) + 2π * η[i] * charge))
                        exp2 = exp(im * kdotr(lattice, k_total, site2))
                        exp3 = exp(im * (kdotr(lattice, k_total, idx) + 2π * η_wrap * charge))
                        @test exp1 ≈ exp2
                        @test exp1 ≈ exp3
                    end
                    k = @inferred momentumspace(lattice, k_idx)
                    exp1 = exp(im * (dot(k, realspace(lattice, site)) + 2π * η[i]))
                    exp2 = exp(im * dot(k, realspace(lattice, site2)))
                    site3, wrap = @inferred wraparound_site(lattice, site2)
                    exp3 = exp(im * (dot(k, realspace(lattice, site3)) + 2π * dot(wrap, twist(lattice))))
                    @test exp1 ≈ exp2
                    @test exp1 ≈ exp3
                end
            end
        end
    end
end

let lattice = KagomeLattice([2,3])
    # try a valid wraparound
    @inferred wraparound(lattice, [4,0,0])
    # invalid basis index
    @test_throws ArgumentError wraparound(lattice, [4,0,3])
    # attempt to "translate" in direction of the basis index
    @test_throws ArgumentError translateη(lattice, 1, 3)
end

# Invalid wraparound for OBC
let lattice = ChainLattice([8], Diagonal([0]))
    @inferred wraparound(lattice, [0])
    @inferred wraparound(lattice, [7])
    @test_throws ArgumentError wraparound(lattice, [19])
    @test_throws ArgumentError wraparound(lattice, [-1])
end

function test_neighbors(lattice, pairs, neigh=Val{1})
    mypairs = Set(pairs)
    neighbors(lattice, neigh) do i, j, η
        pop!(mypairs, (i, j))
        # FIXME: check that η is correct
    end
    @test isempty(mypairs)
end

pairs = [
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5),
    (5, 6),
    (6, 1),
]
test_neighbors(ChainLattice([6], Diagonal([6]), [1//5]), pairs)
test_neighbors(SquareLattice([6, 1]), pairs)
test_neighbors(SquareLattice([1, 6]), pairs)
test_neighbors(CubicLattice([6, 1, 1]), pairs)
test_neighbors(HypercubicLattice{4}([1, 1, 6, 1]), pairs)

pairs = [
    (1, 3),
    (2, 4),
    (3, 5),
    (4, 6),
    (5, 7),
    (6, 8),
    (7, 1),
    (8, 2),
    (1, 2),
    (3, 4),
    (5, 6),
    (7, 8),
    ]
# XXX:
#test_neighbors(SquareLattice([4, 2]), pairs)
#test_neighbors(HypercubicLattice{4}([4, 1, 2, 1]), pairs)
# make sure open boundary conditions in the 2 direction is equivalent
# (unless i decide to remove this feature that doesn't double count these
# links)
test_neighbors(SquareLattice([4, 2], Diagonal([4, 0])), pairs)

# 4x3 Open
pairs = [
    (1, 4),
    (2, 5),
    (3, 6),
    (4, 7),
    (5, 8),
    (6, 9),
    (7, 10),
    (8, 11),
    (9, 12),
    (10, 1),
    (11, 2),
    (12, 3),
    (1, 2),
    (2, 3),
    (4, 5),
    (5, 6),
    (7, 8),
    (8, 9),
    (10, 11),
    (11, 12),
]
test_neighbors(SquareLattice([4, 3], Diagonal([4, 0])), pairs)

# 4x4
pairs = [
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 1),
    (5, 6),
    (6, 7),
    (7, 8),
    (8, 5),
    (9, 10),
    (10, 11),
    (11, 12),
    (12, 9),
    (13, 14),
    (14, 15),
    (15, 16),
    (16, 13),
    (1, 5),
    (2, 6),
    (3, 7),
    (4, 8),
    (5, 9),
    (6, 10),
    (7, 11),
    (8, 12),
    (9, 13),
    (10, 14),
    (11, 15),
    (12, 16),
    (13, 1),
    (14, 2),
    (15, 3),
    (16, 4),
]
test_neighbors(SquareLattice([4, 4]), pairs)
test_neighbors(HypercubicLattice{4}([1, 4, 1, 4]), pairs)
test_neighbors(HypercubicLattice{4}([4, 1, 1, 4]), pairs)

#=
TEST(HypercubicLattice, NearestNeighbors_2x4) (
    std::set<std::pair<size_t, size_t> > pairs = (
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7)
    )
    test_neighbors(HypercubicLattice((2, 4)), pairs)
    test_neighbors(HypercubicLattice((1, 2, 4, 1)), pairs)
)

template <class LatticeType>
static void test_next_nearest_neighbors(const LatticeType &lattice, const std::set<std::pair<size_t, size_t> > &pairs)
(
    auto mypairs = pairs
    lattice.next_nearest_neighbors([&mypairs](size_t i, size_t j, boost::rational<int> eta) (
        const size_t erased = mypairs.erase((i, j))
        if (!erased)
            std::cerr << i << ' ' << j << '\n'
        ASSERT_EQ(erased, 1)

        // fixme: check that eta is correct
        (void) eta
#if 0
        std::cerr << eta << std::endl
#endif
    ))
    ASSERT_EQ(mypairs.size(), 0)
)

TEST(HypercubicLattice, NextNearestNeighbors_OBC) (
    std::set<std::pair<size_t, size_t> > pairs = (
        (0, 2),
        (1, 3),
        (2, 4),
        (3, 5)
    )
    Eigen::Matrix<int, 1, 1> repeater
    repeater(0, 0) = 0
    test_next_nearest_neighbors(HypercubicLattice((6), (0), repeater), pairs)
)

TEST(HypercubicLattice, NextNearestNeighbors_PBC) (
    std::set<std::pair<size_t, size_t> > pairs = (
        (0, 2),
        (1, 3),
        (2, 4),
        (3, 5),
        (4, 0),
        (5, 1)
    )
    test_next_nearest_neighbors(HypercubicLattice((6)), pairs)
)
=#

function test_neighbor_sublattices(lattice, allowed, neigh=Val{1})
    @inferred isbipartite(lattice)
    @inferred istripartite(lattice)

    neighbors(lattice, neigh) do i, j, wrap
        ind1 = @inferred sublattice_index(lattice, i)
        ind2 = sublattice_index(lattice, j)
        @test ind1 in allowed
        @test ind2 in allowed
        @test ind1 != ind2
    end
end

@test isbipartite(SquareLattice([4,4]))
@test !isbipartite(SquareLattice([3,3]))
@test isbipartite(SquareLattice([3,3], Diagonal([0,0])))
@test !isbipartite(SquareLattice([4,3]))
@test !isbipartite(SquareLattice([4,1])) # if two dimensions are specified, we don't consider this bipartite.
@test !istripartite(SquareLattice([4,4]))

let lattice = SquareLattice([4, 3], Diagonal([4,0]))
    test_neighbor_sublattices(lattice, [0,1])
    @test sublattice_index(lattice, 1) == 0
    @test sublattice_index(lattice, 2) == 1
    @test sublattice_index(lattice, 3) == 0
    @test sublattice_index(lattice, 4) == 1
    @test sublattice_index(lattice, 5) == 0
    @test sublattice_index(lattice, 6) == 1
    @test sublattice_index(lattice, 7) == 0
    @test sublattice_index(lattice, 8) == 1
    @test sublattice_index(lattice, 9) == 0
    @test sublattice_index(lattice, 10) == 1
    @test sublattice_index(lattice, 11) == 0
    @test sublattice_index(lattice, 12) == 1
end

let lattice = SquareLattice([4, 3])
    @test_throws ArgumentError sublattice_index(lattice, 1)
end

let lattice = SquareLattice([4, 2])
    test_neighbor_sublattices(lattice, [0,1])
    @test sublattice_index(lattice, 1) == 0
    @test sublattice_index(lattice, 2) == 1
    @test sublattice_index(lattice, 3) == 1
    @test sublattice_index(lattice, 4) == 0
    @test sublattice_index(lattice, 5) == 0
    @test sublattice_index(lattice, 6) == 1
    @test sublattice_index(lattice, 7) == 1
    @test sublattice_index(lattice, 8) == 0
end

let lattice = ChainLattice([6])
    test_neighbor_sublattices(lattice, [0,1])
    @test sublattice_index(lattice, 1) == 0
    @test sublattice_index(lattice, 2) == 1
    @test sublattice_index(lattice, 3) == 0
    @test sublattice_index(lattice, 4) == 1
    @test sublattice_index(lattice, 5) == 0
    @test sublattice_index(lattice, 6) == 1
end

let lattice = SquareLattice([6,1])
    @test_throws ArgumentError sublattice_index(lattice, 1)
end

let lattice = TriangularLattice([3,3])
    @test !isbipartite(lattice)
    @test istripartite(lattice)
end
let lattice = TriangularLattice([4,3])
    @test !isbipartite(lattice)
    @test !istripartite(lattice)
    @test_throws ArgumentError sublattice_index(lattice, 1)
end
let lattice = TriangularLattice([4,3], Diagonal([0,3]))
    @test !isbipartite(lattice)
    @test istripartite(lattice)
end
let lattice = TriangularLattice([4,3], Diagonal([0,0]))
    @test !isbipartite(lattice)
    @test istripartite(lattice)
end

let lattice = TriangularLattice([3,3])
    test_neighbor_sublattices(lattice, [0,1,2])
    @test sublattice_index(lattice, 1) == 0
    @test sublattice_index(lattice, 2) == 1
    @test sublattice_index(lattice, 3) == 2
    @test sublattice_index(lattice, 4) == 2
    @test sublattice_index(lattice, 5) == 0
    @test sublattice_index(lattice, 6) == 1
    @test sublattice_index(lattice, 7) == 1
    @test sublattice_index(lattice, 8) == 2
    @test sublattice_index(lattice, 9) == 0
end

let lattice = HoneycombLattice([5,5])
    @test isbipartite(lattice)
    @test !istripartite(lattice)
    test_neighbor_sublattices(lattice, [0,1])
    for (i, site) in enumerate(lattice)
        @test sublattice_index(lattice, i) == sublattice_index(lattice, site) == site[end] == ((i ⊻ 1) & 1)
    end
end

let lattice = KagomeLattice([4,3])
    @test !isbipartite(lattice)
    @test istripartite(lattice)
    test_neighbor_sublattices(lattice, [0,1,2])
    for (i, site) in enumerate(lattice)
        @test sublattice_index(lattice, i) == sublattice_index(lattice, site) == site[end]
    end
end

function test_neighbor_distances(lattice, neigh=Val{1}, expected_squared=1)
    neighbors(lattice, neigh) do i, j, wrap
        diff = realspace(lattice, i) - realspace(lattice, j) - realspace(bravais(lattice), wrap .* dimensions(lattice))
        dist_squared = dot(diff, diff)
        @test expected_squared ≈ dist_squared atol=1e-8
    end
end

test_neighbor_distances(ChainLattice([8]))
test_neighbor_distances(ChainLattice([8], Diagonal([0])))
test_neighbor_distances(ChainLattice([8]), Val{1}, 1)
test_neighbor_distances(ChainLattice([8], Diagonal([0])), Val{1}, 1)
test_neighbor_distances(ChainLattice([8]), Val{2}, 4)
test_neighbor_distances(ChainLattice([8], Diagonal([0])), Val{2}, 4)
test_neighbor_distances(ChainLattice([8]), Val{5}, 25)
test_neighbor_distances(ChainLattice([8], Diagonal([0])), Val{5}, 25)
test_neighbor_distances(SquareLattice([4,6]))
test_neighbor_distances(SquareLattice([4,6], Diagonal([0,0])))
test_neighbor_distances(SquareLattice([4,6]), Val{2}, 2)
test_neighbor_distances(SquareLattice([4,6], Diagonal([0,0])), Val{2}, 2)
test_neighbor_distances(TriangularLattice([4,6]))
test_neighbor_distances(TriangularLattice([4,6], Diagonal([0,0])))
test_neighbor_distances(TriangularLattice([4,6]), Val{2}, 3)
test_neighbor_distances(TriangularLattice([4,6], Diagonal([0,0])), Val{2}, 3)
test_neighbor_distances(TriangularLattice([4,6]), Val{3}, 4)
test_neighbor_distances(TriangularLattice([4,6], Diagonal([0,0])), Val{3}, 4)
test_neighbor_distances(HoneycombLattice([4,6]))
test_neighbor_distances(HoneycombLattice([4,6], Diagonal([0,0])))
test_neighbor_distances(KagomeLattice([4,6]))
test_neighbor_distances(KagomeLattice([4,6], Diagonal([0,0])))
