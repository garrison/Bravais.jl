using Bravais
using Base.Test

using Compat

lattice = BravaisLattice([4,6,8])
@test length(lattice) == 4*6*8
@test size(lattice) == (4*6*8,)

# NOTE: assumes hypercubic by default, but BravaisLattice does not
# define any concept of nearest neighbors
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

lattices = Any[
    # Basic BravaisLattice
    BravaisLattice([4,6]),

    # 1D and effectively 1D
    HypercubicLattice([8]),
    HypercubicLattice([1, 1, 8]),
    HypercubicLattice([1, 8, 1]),

    # 1D Twisted
    HypercubicLattice([6], diagm([6]), [0//1]),
    HypercubicLattice([6], diagm([6]), [1//2]),
    HypercubicLattice([6], diagm([6]), [1//5]),

    # 2D
    HypercubicLattice([4,6]),
    HypercubicLattice([6,4]),

    # 2D twisted
    HypercubicLattice([4,6], diagm([4,6]), [1//7, 3//13]),
    HypercubicLattice([6,4], diagm([6,4]), [1//7, 3//13]),

    # 2D twisted cylinder
    HypercubicLattice([6,4], diagm([6,0]), [1//7, 0//1]),

    # 2D twisted helix
    HypercubicLattice([6,4], [6 0; 1 4], [1//7, 3//13]),

    # 3D like
    HypercubicLattice([4,2,6]),
    HypercubicLattice([4,2,6,1]),

    # 3D open twisted helix
    HypercubicLattice([4,6,8], [4 0 0; 0 6 0; 1 0 8], [1//7, 0//1, 3//13]),

    # Other lattices
    TriangularLattice([4,6]),
    HoneycombLattice([3,5]),
    KagomeLattice([8,6]),
]

for lattice in lattices
    len = length(lattice)
    @test size(lattice) == (len,)

    @test_throws BoundsError lattice[0]
    @test_throws BoundsError lattice[len+1]

    @test findfirst(lattice, maxcoords(lattice)) == 0

    last = false
    for (i, site) in enumerate(lattice)
        checkbounds(lattice, i)
        @test lattice[i] == site
        @test site in lattice
        @test findfirst(lattice, site) == i
        @test last == false
        if i == length(lattice)
            last = true
        end

        wraparound(lattice, site) # FIXME: actually test something here

        realspace(lattice, site) # FIXME: actually test something here
        realspace(lattice, site, ones(Int, ndimensions(lattice))) # FIXME: actually test something here
    end
    @test last == true

    d = ndimensions(lattice)
    η = twist(lattice)
    M = repeater(lattice)

    @test ishelical(lattice) == !isdiag(repeater(lattice))

    if !isbravais(lattice)
        # FIXME!!
        continue
    end

    # Check the momenta across the boundary conditions
    n_k_idx = nmomenta(lattice)
    for site in lattice
        for k_idx in 1:n_k_idx
            for i in 1:d
                if M[i,i] != 0
                    site2 = site + vec(M[i,:])
                    exp1 = exp(im * (kdotr(lattice, k_idx, site) + 2π * η[i]))
                    exp2 = exp(im * kdotr(lattice, k_idx, site2))
                    idx, η_wrap = wraparoundη(lattice, site2)
                    exp3 = exp(im * (kdotr(lattice, k_idx, idx) + 2π * η_wrap))
                    @test_approx_eq exp1 exp2
                    @test_approx_eq exp1 exp3
                end
            end
        end
    end

    # Check that translating multiple times, in total by an $\mathbf(A)_i$
    # vector, brings us back to where we started (and picking up a phase due
    # to $η_i$).
    sz = length(lattice)
    for i in 1:d
        if M[i,i] == 0
            continue
        end

        # Initialize our array of each index, along with a beginning
        # phase of zero
        indices = @compat Tuple{Int, Rational{Int}}[(z, 0//1) for z in 1:sz]

        # Translate as many times as we need to in each dimension
        for j in 1:d
            for z in 1:M[i,j]
                new_indices = @compat Tuple{Int, Rational{Int}}[]
                sizehint!(new_indices, sz)
                for (s, old_η) in indices
                    newidx, wrap_η = translateη(lattice, s, j)
                    push!(new_indices, (newidx, wrap_η + old_η))
                end
                indices = new_indices
            end
        end

        # assert that everything is as expected (i.e., overall we got the
        # identity element while picking up some phase)
        for z in 1:sz
            @test indices[z] == (z, η[i])
        end
    end

end

function test_nearest_neighbors(lattice, pairs)
    mypairs = Set(pairs)
    neighbors(lattice, :nearest) do i, j, η
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
test_nearest_neighbors(HypercubicLattice([6], diagm([6]), [1//5]), pairs)
test_nearest_neighbors(HypercubicLattice([6, 1]), pairs)
test_nearest_neighbors(HypercubicLattice([1, 6]), pairs)
test_nearest_neighbors(HypercubicLattice([6, 1, 1]), pairs)
test_nearest_neighbors(HypercubicLattice([1, 1, 6, 1]), pairs)

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
#test_nearest_neighbors(HypercubicLattice([4, 2]), pairs)
#test_nearest_neighbors(HypercubicLattice([4, 1, 2, 1]), pairs)
# make sure open boundary conditions in the 2 direction is equivalent
# (unless i decide to remove this feature that doesn't double count these
# links)
test_nearest_neighbors(HypercubicLattice([4, 2], diagm([4, 0])), pairs)

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
test_nearest_neighbors(HypercubicLattice([4, 3], diagm([4, 0])), pairs)

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
test_nearest_neighbors(HypercubicLattice([4, 4]), pairs)
test_nearest_neighbors(HypercubicLattice([1, 4, 1, 4]), pairs)
test_nearest_neighbors(HypercubicLattice([4, 1, 1, 4]), pairs)

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
    test_nearest_neighbors(HypercubicLattice((2, 4)), pairs)
    test_nearest_neighbors(HypercubicLattice((1, 2, 4, 1)), pairs)
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

TEST(HypercubicLattice, SublatticeIndex_4_3) (
    const HypercubicLattice lattice((4, 3))
    ASSERT_EQ(lattice.sublattice_index(0), 0)
    ASSERT_EQ(lattice.sublattice_index(1), 1)
    ASSERT_EQ(lattice.sublattice_index(2), 0)
    ASSERT_EQ(lattice.sublattice_index(3), 1)
    ASSERT_EQ(lattice.sublattice_index(4), 0)
    ASSERT_EQ(lattice.sublattice_index(5), 1)
    ASSERT_EQ(lattice.sublattice_index(6), 0)
    ASSERT_EQ(lattice.sublattice_index(7), 1)
    ASSERT_EQ(lattice.sublattice_index(8), 0)
    ASSERT_EQ(lattice.sublattice_index(9), 1)
    ASSERT_EQ(lattice.sublattice_index(10), 0)
    ASSERT_EQ(lattice.sublattice_index(11), 1)
)

TEST(HypercubicLattice, SublatticeIndex_4_2) (
    const HypercubicLattice lattice((4, 2))
    ASSERT_EQ(lattice.sublattice_index(0), 0)
    ASSERT_EQ(lattice.sublattice_index(1), 1)
    ASSERT_EQ(lattice.sublattice_index(2), 1)
    ASSERT_EQ(lattice.sublattice_index(3), 0)
    ASSERT_EQ(lattice.sublattice_index(4), 0)
    ASSERT_EQ(lattice.sublattice_index(5), 1)
    ASSERT_EQ(lattice.sublattice_index(6), 1)
    ASSERT_EQ(lattice.sublattice_index(7), 0)
)

TEST(HypercubicLattice, SublatticeIndex_8_1) (
    const HypercubicLattice lattice((8, 1))
    ASSERT_EQ(lattice.sublattice_index(0), 0)
    ASSERT_EQ(lattice.sublattice_index(1), 1)
    ASSERT_EQ(lattice.sublattice_index(2), 0)
    ASSERT_EQ(lattice.sublattice_index(3), 1)
    ASSERT_EQ(lattice.sublattice_index(4), 0)
    ASSERT_EQ(lattice.sublattice_index(5), 1)
    ASSERT_EQ(lattice.sublattice_index(6), 0)
    ASSERT_EQ(lattice.sublattice_index(7), 1)
)
=#
