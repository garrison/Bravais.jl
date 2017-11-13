# Bravais.jl -- Copyright 2015 James R. Garrison and contributors.
# Provided under the MIT License.

# NOTE: Most symbols used below are defined in the accompanying
# documentation.

# FIXME: need other constructors for e.g. OBC

# FIXME: have a way of getting the neighbors of a given site (instead
# of iterating through all of them).  for this, we may actually want
# to double count (or not).

# FIXME: in lattice w/ basis, be sure to error out anytime the "basis
# index" is out of range

__precompile__()

module Bravais

using StaticArrays

include("delegate.jl")

# Julia currently lacks row-major versions of sub2ind and ind2sub.  We
# want to use these so the ordering for returned tuples is correct;
# see discussion at <https://github.com/JuliaLang/julia/pull/10337#issuecomment-78206799>.
#
# Also, once julia easily supports it we can speed up ind2sub by
# precomputing some "magic numbers" for given lattice dimensions (see
# https://github.com/JuliaLang/julia/issues/8188#issuecomment-56763806).
rowmajor_ind2sub(dims, index) = reverse(ind2sub(reverse(dims), index))

abstract type AbstractSiteNetwork <: AbstractVector{Vector{Int}} end
abstract type AbstractLattice{D} <: AbstractSiteNetwork end

abstract type AbstractBravaisLattice{D} <: AbstractLattice{D} end
abstract type AbstractLatticeWithBasis{D,Dp1} <: AbstractLattice{D} end

abstract type WrappedBravaisLattice{D} <: AbstractBravaisLattice{D} end
abstract type WrappedLatticeWithBasis{D,Dp1} <: AbstractLatticeWithBasis{D,Dp1} end

# NOTE: In the following, D is the number of dimensions, Dsq == D^2,
# and Dp1 == D+1.

struct BravaisLattice{D,Dsq} <: AbstractBravaisLattice{D}
    N_tot::Int  # total number of sites
    N::SVector{D,Int}  # lattice extent in each dimension
    M::SMatrix{D,D,Int,Dsq}  # "M" matrix (will be diagonal for non-helical boundary)
    η::SVector{D,Rational{Int}}  # twist in each direction (rational multiple of 2π)
    a::SMatrix{D,D,Float64,Dsq}  # primitive lattice vectors
    b::SMatrix{D,D,Float64,Dsq}  # reciprocal lattice vectors
    momenta::Vector{SVector{D,Rational{Int}}}
    strides::SVector{D,Int}

    function BravaisLattice{D,Dsq}(N::AbstractVector{Int},
                                   M::AbstractMatrix{Int}=Diagonal(N), # assumes pbc
                                   η::AbstractVector{Rational{Int}}=zeros(SVector{D,Rational{Int}}),
                                   a::AbstractMatrix{Float64}=eye(SMatrix{D,D})) where {D,Dsq}
        @assert Dsq == D * D

        # check N
        length(N) == D || throw(ArgumentError(""))
        all(d_i > 0 for d_i in N) || throw(ArgumentError(""))
        N_tot = prod(N)

        # check M
        size(M) == (D, D) || throw(ArgumentError(""))
        istril(M) || throw(ArgumentError(""))
        for i in 1:D
            M[i,i] == 0 || M[i,i] == N[i] || throw(ArgumentError(""))
            if M[i,i] != N[i]
                all(s == 0 for s in M[i,:]) || throw(ArgumentError(""))
                all(s == 0 for s in M[:,i]) || throw(ArgumentError(""))
            end
        end

        # check η
        length(η) == D || throw(ArgumentError(""))
        all(0 <= η_i < 1 for η_i in η) || throw(ArgumentError(""))
        all(M[i,i] != 0 || η[i] == 0 for i in 1:D) || throw(ArgumentError(""))

        # check a and generate b
        size(a) == (D, D) || throw(ArgumentError(""))
        b = 2pi * transpose(inv(a))

        # generate allowed momenta
        momenta_range = [s != 0 ? s : 1 for s in diag(M)]
        nmomenta = prod(momenta_range)
        momenta = sizehint!(SVector{D,Rational{Int}}[], nmomenta)
        momentum = zeros(MVector{D,Rational{Int}})
        for idx in 1:nmomenta
            # FIXME: getting idx should be easier, once
            # multidimensional iteration is supported
            ñ = [rowmajor_ind2sub(tuple(momenta_range...), idx)...] .- 1
            for i in 1:D
                if M[i,i] == 0
                    momentum[i] = 0
                else
                    tmp = sum(Rational{Int}[M[i,j] * momentum[j] for j in 1:i-1])
                    momentum[i] = (ñ[i] + η[i] - tmp) // M[i,i]
                end
            end
            push!(momenta, momentum)
        end

        # calculate strides
        strides = zeros(Int, D)
        s = 1
        for i in D:-1:1
            strides[i] = s
            s *= N[i]
        end

        new(N_tot, copy(N), copy(M), copy(η), copy(a), b, momenta, strides)
    end
end

Base.@pure square(x) = x * x

BravaisLattice{D}(args...) where {D} =
    BravaisLattice{D,square(D)}(args...)

BravaisLattice(N::StaticVector{D,Int}, args...) where {D} =
    BravaisLattice{D}(N, args...)

struct LatticeWithBasis{D,Dsq,Dp1} <: AbstractLatticeWithBasis{D,Dp1}
    N_tot::Int
    maxcoords::SVector{Dp1,Int}
    bravaislattice::BravaisLattice{D,Dsq}
    basis::Vector{SVector{D,Float64}}
    strides::SVector{Dp1,Int}

    function LatticeWithBasis{D,Dsq,Dp1}(N::AbstractVector{Int},
                                         M::AbstractMatrix{Int}=Diagonal(N), # assumes pbc
                                         η::AbstractVector{Rational{Int}}=zeros(SVector{D,Rational{Int}}),
                                         a::AbstractMatrix{Float64}=eye(SMatrix{D,D}),
                                         basis::Vector{SVector{D,Float64}}=[zeros(SVector{D,Float64})]) where {D,Dsq,Dp1}
        @assert Dp1 == D + 1
        bravaislattice = BravaisLattice{D,Dsq}(N, M, η, a)

        # check basis
        nbasis = length(basis)
        nbasis > 0 || throw(ArgumentError(""))

        # determine N_tot and maxcoords, now that we know the basis size
        N_tot = prod(N) * nbasis
        maxcoords = vcat(N, [nbasis])

        # calculate strides
        strides = zeros(Int, length(maxcoords))
        strides[end] = 1
        s = nbasis
        for i in D:-1:1
            strides[i] = s
            s *= N[i]
        end

        new(N_tot, maxcoords, bravaislattice, copy(basis), strides)
    end
end

isbravais(::AbstractBravaisLattice) = true
isbravais(::AbstractLatticeWithBasis) = false

@doc doc"Returns the underlying Bravais lattice" ->
bravais(lattice::AbstractBravaisLattice) = lattice
bravais(lattice::LatticeWithBasis) = lattice.bravaislattice
bravais(lattice::WrappedLatticeWithBasis) = bravais(lattice.lattice)

const LatticeImplUnion{D,Dsq} = Union{BravaisLattice{D,Dsq}, LatticeWithBasis{D,Dsq}}
const WrappedLatticeUnion{D} = Union{WrappedBravaisLattice{D}, WrappedLatticeWithBasis{D}}

_strides(lattice::LatticeImplUnion) = lattice.strides
_strides(lattice::WrappedLatticeUnion) = _strides(lattice.lattice)

maxcoords(lattice::BravaisLattice) = lattice.N
maxcoords(lattice::LatticeWithBasis) = lattice.maxcoords
maxcoords(lattice::WrappedLatticeUnion) = maxcoords(lattice.lattice)

Base.length(lattice::LatticeImplUnion) = lattice.N_tot

Base.eachindex(lattice::AbstractLattice) = Base.OneTo(length(lattice))

Base.size(lattice::AbstractLattice) = (length(lattice),)

function Base.getindex(lattice::Union{AbstractBravaisLattice{Dprime},AbstractLatticeWithBasis{D,Dprime} where D}, index::Integer) where {Dprime}
    checkbounds(lattice, index)
    # Alternatively: return [rowmajor_ind2sub(tuple(maxcoords(lattice)...), index)...] - 1
    strides = _strides(lattice)
    @assert length(strides) == Dprime
    rv = MVector{Dprime,Int}()
    r = index - 1
    for i in 1:Dprime-1
        d, r = divrem(r, strides[i])
        rv[i] = d
    end
    rv[end] = r
    return SVector(rv)
end

function Base.in(site::AbstractVector{Int}, lattice::AbstractLattice)
    mc = maxcoords(lattice)
    length(mc) == length(site) || return false
    for i in 1:length(site)
        0 <= site[i] < mc[i] || return false
    end
    return true
end

Base.findfirst(lattice::AbstractLattice, site::AbstractVector{Int}) = site ∉ lattice ? 0 : dot(site, _strides(lattice)) + 1

function Base.start(lattice::Union{AbstractBravaisLattice{Dprime},AbstractLatticeWithBasis{D,Dprime} where D}) where {Dprime}
    zeros(SVector{Dprime,Int})
end

function Base.done(lattice::AbstractLattice, site::AbstractVector{Int})
    mc = maxcoords(lattice)
    @assert length(site) == length(mc)
    for i in 1:length(mc)
        site[i] < mc[i] || return true
    end
    return false
end

function Base.next(lattice::Union{AbstractBravaisLattice{Dprime},AbstractLatticeWithBasis{D,Dprime} where D}, site::SVector{Dprime,Int})::NTuple{2,SVector{Dprime,Int}} where {Dprime}
    newsite = MVector{Dprime,Int}(site)
    mc = maxcoords(lattice)
    # XXX FIXME @assert something # otherwise BoundsError!
    @assert length(mc) == Dprime >= 1
    for i in Dprime:-1:2
        newsite[i] += 1
        if newsite[i] == mc[i]
            newsite[i] = 0
        else
            return site, SVector(newsite)
        end
    end
    newsite[1] += 1
    return site, SVector(newsite)
end

# FIXME: can we return readonly array views/proxies of some of the
# following things?  It seems to be difficult if not impossible at the
# moment -- even base library functions like eigfact return mutable
# arrays.

dimensions(lattice::BravaisLattice) = lattice.N  # FIXME: rename this `extent`?
dimensions(lattice::LatticeWithBasis) = dimensions(bravais(lattice))
dimensions(lattice::WrappedLatticeUnion) = dimensions(lattice.lattice)

ndimensions(::AbstractLattice{D}) where {D} = D

twist(lattice::BravaisLattice) = lattice.η
twist(lattice::LatticeWithBasis) = bravais(lattice).η
twist(lattice::WrappedLatticeUnion) = twist(lattice.lattice)
repeater(lattice::BravaisLattice) = lattice.M
repeater(lattice::LatticeWithBasis) = bravais(lattice).M
repeater(lattice::WrappedLatticeUnion) = repeater(lattice.lattice)

ishelical(lattice::AbstractLattice) = !isdiag(repeater(lattice))

# We intentionally define many of the following functions for Bravais
# lattices only.  If one wants to query them for a lattice with a
# basis, call e.g. primvecs(bravais(lattice)).

@doc doc"Returns the primitive vectors of the direct lattice" ->
primvecs(lattice::BravaisLattice) = lattice.a

@doc doc"Returns the primitive vectors of the reciprocal lattice" ->
recivecs(lattice::BravaisLattice) = lattice.b

nmomenta(lattice::BravaisLattice) = length(lattice.momenta)
nmomenta(lattice::WrappedBravaisLattice) = nmomenta(lattice.lattice)

eachmomentumindex(lattice::AbstractBravaisLattice) = 1:nmomenta(lattice)

# FIXME: need a way to iterate momenta

momentum(lattice::BravaisLattice, idx::Integer) = lattice.momenta[idx]
function momentum(lattice::BravaisLattice{D}, idx::Integer, charge::Integer)::SVector{D,Rational{Int}} where {D}
    # "total momentum", really.  note that this may return things greater than one.
    x1 = momentum(lattice, idx)
    if charge == 1
        return x1
    end
    offsets = zeros(MVector{D,Rational{Int}})
    for i in 1:D
        if lattice.M[i,i] != 0
            offsets[i] = lattice.η[i] * (charge - 1)
            for j in 1:i-1
                offsets[i] -= lattice.M[i,j] * offsets[j]
            end
            offsets[i] = offsets[i] // lattice.M[i,i]
        end
    end
    return x1 + offsets
end

kdotr(lattice::BravaisLattice, ksite::AbstractVector{Rational{Int}}, site::AbstractVector{Int}) = 2pi * dot(ksite, site)
kdotr(lattice::BravaisLattice, kidx::Integer, site::AbstractVector{Int}) = kdotr(lattice, momentum(lattice, kidx), site)
kdotr(lattice::BravaisLattice, k, ridx::Integer) = kdotr(lattice, k, lattice[ridx])

momentumspace(lattice::BravaisLattice, k::AbstractVector{Rational{Int}}) = lattice.b * k
momentumspace(lattice::BravaisLattice, kidx::Integer) = momentumspace(lattice, momentum(lattice, kidx))
# FIXME: make it possible to get momentum-space points translated to the first brillouin zone (momentumspace_bz)

realspace(lattice::BravaisLattice, site::AbstractVector{Int}) = lattice.a * site
# NOTE: the following function does not error out if wrap makes no sense due to OBC.
realspace(lattice::BravaisLattice, site::AbstractVector{Int}, wrap::AbstractVector{Int}) = lattice.a * (site + transpose(lattice.M) * wrap)

function realspace(lattice::LatticeWithBasis, site::AbstractVector{Int}, args...)
    length(site) == length(lattice.maxcoords) || throw(ArgumentError(""))
    return realspace(bravais(lattice), site[1:end-1], args...) + lattice.basis[site[end]+1]
end

realspace(lattice::LatticeImplUnion, ridx::Integer, args...) = realspace(lattice, lattice[ridx], args...)

function wraparound_site!(lattice::LatticeImplUnion{D}, site::AbstractVector{Int}) where {D}
    mc = maxcoords(lattice)
    length(site) == length(mc) || throw(ArgumentError(""))

    # For lattice w/ basis, make sure the last index is in range, as
    # we cannot wrap it around.
    if isa(lattice, LatticeWithBasis)
        if !(0 <= site[end] < mc[end])
            throw(ArgumentError("Site has invalid basis index."))
        end
    end

    wrap = zeros(Int, D)

    N = bravais(lattice).N
    M = bravais(lattice).M

    # FIXME: use a @generated function to unroll
    for i in D:-1:1
        if !(0 <= site[i] < N[i])
            if M[i,i] != 0
                # periodic/twisted BC's in this direction
                wrap[i] = fld(site[i], N[i])
                for j in 1:i
                    site[j] -= wrap[i] * M[i,j]
                end
            else
                # OBC in this direction
                #
                # XXX: This does not provide any exception guarantee.
                # If we wanted, we could easily undo what we had done
                # before throwing.  (Additionally, any mutating
                # functions that call this may also wish to undo their
                # side effects before propagating the exception.)
                throw(ArgumentError("This lattice has open boundary conditions in the direction we are trying to wraparound!"))
            end
        end
    end

    return site, wrap
end

wraparound_site!(lattice::WrappedLatticeUnion, site::AbstractVector{Int}) = wraparound_site!(lattice.lattice, site)

function wraparound_site(lattice::Union{AbstractBravaisLattice{Dprime},AbstractLatticeWithBasis{D,Dprime} where D}, site::AbstractVector{Int}) where {Dprime}
    site_, wrap = wraparound_site!(lattice, MVector{Dprime,Int}(site))
    return SVector(site_), wrap
end

wraparound_site(lattice::AbstractLattice, index::Integer) = wraparound_site(lattice, lattice[index])

function wraparound(lattice::AbstractLattice, site_or_index::Union{AbstractVector{Int}, Integer})
    site, wrap = wraparound_site(lattice, site_or_index)
    idx = findfirst(lattice, site)
    return idx, wrap
end

function wraparoundη(lattice::AbstractLattice, site_or_index::Union{AbstractVector{Int}, Integer})
    idx, wrap = wraparound(lattice, site_or_index)
    η = dot(wrap, twist(lattice))
    return idx, η
end

function translate_site!(lattice::AbstractLattice, site::AbstractVector{Int}, direction::Integer)
    if !isbravais(lattice) && direction > ndimensions(lattice)
        throw(ArgumentError("Cannot translate in the 'direction' of the basis index."))
    end
    site[direction] += 1
    return wraparound_site!(lattice, site)
end

function translate_site(lattice::Union{AbstractBravaisLattice{Dprime},AbstractLatticeWithBasis{D,Dprime} where D}, site::AbstractVector{Int}, direction::Integer) where {Dprime}
    site_, wrap = translate_site!(lattice, MVector{Dprime,Int}(site), direction)
    return SVector(site_), wrap
end

translate_site(lattice::AbstractLattice, index::Integer, direction::Integer) = translate_site(lattice, lattice[index], direction)

function translate(lattice::AbstractLattice, site_or_index::Union{AbstractVector{Int}, Integer}, direction::Integer)
    site, wrap = translate_site(lattice, site_or_index, direction)
    idx = findfirst(lattice, site)
    return idx, wrap
end

function translateη(lattice::AbstractLattice, site_or_index::Union{AbstractVector{Int}, Integer}, direction::Integer)
    idx, wrap = translate(lattice, site_or_index, direction)
    η = dot(wrap, twist(lattice))
    return idx, η
end

function neighbors(f, lattice::AbstractLattice, neigh=Val{1})
    for ridx in 1:length(lattice)
        siteneighbors(f, lattice, ridx, neigh)
    end
end

function neighborsη(f, lattice::AbstractLattice, neigh=Val{1})
    neighbors(lattice, neigh) do idx1, idx2, wrap
        η = dot(wrap, twist(lattice))
        f(idx1, idx2, η)
    end
end

#= Begin common code that will help with specific Bravais lattice implementations =#

sublattice_index(lattice::AbstractLattice, ridx::Int) = sublattice_index(lattice, lattice[ridx])

siteneighbors(f, lattice::AbstractLattice, site_or_index::Union{AbstractVector{Int}, Integer}) = siteneighbors(f, lattice, site_or_index, 1)

function _siteneighbors2d(f, lattice::AbstractLattice{2}, ridx::Integer, offsets)
    M = repeater(lattice)
    mc = maxcoords(lattice)

    site = lattice[ridx]

    for offset in offsets
        newsite = site + offset
        g(i) = M[i,i] == 0 && !(0 <= newsite[i] < mc[i])
        if g(1) || g(2)
            # In directions w/ OBC, do not provide the neighbors!
            #
            # NOTE: In e.g. a 3d lattice, similar logic will fail
            # if there's one OBC direction and helical BCs in the
            # two periodic directions.
            continue
        end
        newidx, wrap = wraparound(lattice, newsite)
        f(ridx, newidx, wrap)
    end
    nothing
end

#= Begin specific Bravais lattice implementations =#

function _hypercubic_sublattice_index(site::AbstractVector{Int})
    parity = 0
    for x in site
        parity = parity ⊻ x
    end
    return parity & 1
end

struct HypercubicLattice{D,Dsq} <: WrappedBravaisLattice{D}
    lattice::BravaisLattice{D,Dsq}
    bipartite::Bool

    function HypercubicLattice{D,Dsq}(N::AbstractVector{Int},
                                      M::AbstractMatrix{Int}=Diagonal(N), # assumes pbc
                                      η::AbstractVector{Rational{Int}}=zeros(SVector{D,Rational{Int}})) where {D,Dsq}
        bravaislattice = BravaisLattice{D,Dsq}(N, M, η)
        bipartite = true
        for i in 1:D
            if M[i,i] != 0
                # Attempt to go across the boundary in the `i`
                # direction, and test if the site is on the same
                # sublattice after being wrapped around.
                site1 = zeros(Int, D)
                site1[i] = N[i]
                site2, = wraparound_site(bravaislattice, site1)
                if _hypercubic_sublattice_index(site1) != _hypercubic_sublattice_index(site2)
                    bipartite = false
                    break
                end
            end
        end
        new(bravaislattice, bipartite)
    end
end

HypercubicLattice{D}(args...) where {D} =
    HypercubicLattice{D,square(D)}(args...)

HypercubicLattice(N::StaticVector{D,Int}, args...) where {D} =
    HypercubicLattice{D}(N, args...)

isbipartite(lattice::HypercubicLattice) = lattice.bipartite
istripartite(::HypercubicLattice) = false

function sublattice_index(lattice::HypercubicLattice, site::AbstractVector{Int})
    isbipartite(lattice) || throw(ArgumentError("Hypercubic lattice must be bipartite for it to have sublattice indices."))
    return _hypercubic_sublattice_index(site)
end

const ChainLattice = HypercubicLattice{1,1}
const SquareLattice = HypercubicLattice{2,4}
const CubicLattice = HypercubicLattice{3,9}

# FIXME: do the wrapping etc in a common function for all lattice
# types.  then the siteneighbors function (or `siteneighbordsimpl`)
# can be specified even more simply.

# BASICALLY: give two things: single count and double count offsets.
# also, for each basis index.  then again, this generic code is not
# going to be quite as optimized as e.g. the hypercubic lattice.  but
# that is okay, i can just have special code for that.

# XXX FIXME: we want to be able (in the end) to choose different
# primitive vectors so we can have a weird helical lattice.  how are
# we going to support this??

function siteneighbors(f, lattice::ChainLattice, ridx::Integer, ::Type{Val{1}})
    M = lattice.lattice.M
    mc = maxcoords(lattice)
    site = lattice[ridx]
    newsite = @SVector [site[1] + 1]
    if M[1,1] <= 1 && newsite[1] >= mc[1]
        # Do not provide neighbors across an open boundary.
        return
    end
    newidx, wrap = wraparound(lattice, newsite)
    f(ridx, newidx, wrap)
    nothing
end

function siteneighbors(f, lattice::ChainLattice, ridx::Integer, ::Type{Val{N}}) where {N}
    M = lattice.lattice.M
    mc = maxcoords(lattice)
    site = lattice[ridx]
    newsite = @SVector [site[1] + N]
    if M[1,1] <= 1 && newsite[1] >= mc[1]
        # Do not provide neighbors across an open boundary.
        return
    end
    newidx, wrap = wraparound(lattice, newsite)
    f(ridx, newidx, wrap)
    nothing
end

function siteneighbors(f, lattice::HypercubicLattice{D}, ridx::Integer, ::Type{Val{1}}) where {D}
    M = lattice.lattice.M
    mc = maxcoords(lattice)

    site = lattice[ridx]

    for i in 1:D
        newsite = MVector(site)
        newsite[i] += 1
        if M[i,i] <= 1 && newsite[i] >= mc[i]
            # In directions w/ OBC or length one, do not provide the neighbors!
            continue
        end
        newidx, wrap = wraparound(lattice, newsite)
        f(ridx, newidx, wrap)
    end
end

function siteneighbors(f, lattice::SquareLattice, ridx::Integer, ::Type{Val{2}})
    offsets = @SVector [@SVector([1,1]), @SVector([-1,1])]
    _siteneighbors2d(f, lattice, ridx, offsets)
end

function _triangular_sublattice_index(site::AbstractVector{Int})
    @assert length(site) == 2
    # FIXME: mod is slow
    return mod(site[2] - site[1], 3)
end

struct TriangularLattice <: WrappedBravaisLattice{2}
    lattice::BravaisLattice{2,4}
    tripartite::Bool

    function TriangularLattice(N::AbstractVector{Int},
                               M::AbstractMatrix{Int}=Diagonal(N), # assumes pbc
                               η::AbstractVector{Rational{Int}}=zeros(SVector{2,Rational{Int}}))
        bravaislattice = BravaisLattice{2,4}(N, M, η, [1.0 0; 0.5 sqrt(3)/2]')
        tripartite = true
        for i in 1:2
            if M[i,i] != 0
                # Attempt to go across the boundary in the `i`
                # direction, and test if the site is on the same
                # sublattice after being wrapped around.
                site1 = [0, 0]
                site1[i] = N[i]
                site2, = wraparound_site(bravaislattice, site1)
                if _triangular_sublattice_index(site1) != _triangular_sublattice_index(site2)
                    tripartite = false
                    break
                end
            end
        end
        new(bravaislattice, tripartite)
    end
end

isbipartite(::TriangularLattice) = false # although "technically" it might be if it's just a chain
istripartite(lattice::TriangularLattice) = lattice.tripartite

function sublattice_index(lattice::TriangularLattice, site::AbstractVector{Int})
    istripartite(lattice) || throw(ArgumentError("Triangular lattice must be tripartite for it to have sublattice indices."))
    return _triangular_sublattice_index(site)
end

# FIXME: many of these neighbor functions can use special handling when the lattice height or width is 1 in a direction.  or we could just forbid this.

function siteneighbors(f, lattice::TriangularLattice, ridx::Integer, ::Type{Val{1}})
    offsets = @SVector [@SVector([1,0]), @SVector([0,1]), @SVector([-1,1])]
    _siteneighbors2d(f, lattice, ridx, offsets)
end

function siteneighbors(f, lattice::TriangularLattice, ridx::Integer, ::Type{Val{2}})
    offsets = @SVector [@SVector([2,-1]), @SVector([1,1]), @SVector([-1,2])]
    _siteneighbors2d(f, lattice, ridx, offsets)
end

function siteneighbors(f, lattice::TriangularLattice, ridx::Integer, ::Type{Val{3}})
    offsets = @SVector [@SVector([2,0]), @SVector([0,2]), @SVector([-2,2])]
    _siteneighbors2d(f, lattice, ridx, offsets)
end

function rhombi(lattice::TriangularLattice)
    @SVector [
        @SVector([@SVector([0,0]), @SVector([1,-1]), @SVector([1,0]), @SVector([0,1])]),
        @SVector([@SVector([0,0]), @SVector([0,1]), @SVector([-1,1]), @SVector([-1,0])]),
        @SVector([@SVector([0,0]), @SVector([-1,0]), @SVector([0,-1]), @SVector([1,-1])]),
    ]
end

@delegate WrappedBravaisLattice.lattice [ Base.length, ndimensions, primvecs, recivecs, momentum, kdotr, momentumspace, realspace ]

#= End specific Bravais lattice implementations =#

#= Begin specific lattice w/ basis implementations =#

struct HoneycombLattice <: WrappedLatticeWithBasis{2,3}
    lattice::LatticeWithBasis{2,4,3}

    function HoneycombLattice(N::AbstractVector{Int},
                              M::AbstractMatrix{Int}=Diagonal(N), # assumes pbc
                              η::AbstractVector{Rational{Int}}=zeros(SVector{2,Rational{Int}}))
        a = @SMatrix([1.5 sqrt(3)/2; 0 sqrt(3)])'
        basis = [@SVector([0.0, 0]), @SVector([1.0, 0])]
        return new(LatticeWithBasis{2,4,3}(N, M, η, a, basis))
    end
end

isbipartite(::HoneycombLattice) = true
istripartite(::HoneycombLattice) = false

function sublattice_index(::HoneycombLattice, site::AbstractVector{Int})
    retval = typeassert(site[end], Int)
    @assert 0 <= retval < 2 # for bipartite lattice
    return retval
end

function siteneighbors(f, lattice::HoneycombLattice, ridx::Integer, ::Type{Val{1}})
    site = lattice[ridx]
    if site[end] == 0
        offsets = @SVector [@SVector([0, 0, 1]), @SVector([-1, 0, 1]), @SVector([-1, 1, 1])]
        _siteneighbors2d(f, lattice, ridx, offsets)
    else
        @assert site[end] == 1
        # Commented out, otherwise we are double counting!
        #offsets = @SVector [@SVector([0, 0, -1]), @SVector([1, 0, -1]), @SVector([1, -1, -1])]
        #_siteneighbors2d(f, lattice, ridx, offsets)
    end
end

struct KagomeLattice <: WrappedLatticeWithBasis{2,3}
    lattice::LatticeWithBasis{2,4,3}

    function KagomeLattice(N::AbstractVector{Int},
                           M::AbstractMatrix{Int}=Diagonal(N), # assumes pbc
                           η::AbstractVector{Rational{Int}}=zeros(SVector{2,Rational{Int}}))
        a = @SMatrix([2 0; 1 sqrt(3)])'
        basis = [@SVector([0.0, 0]), @SVector([0.5, √3/2]), @SVector([1.0, 0])]
        return new(LatticeWithBasis{2,4,3}(N, M, η, a, basis))
    end
end

isbipartite(::KagomeLattice) = false
istripartite(::KagomeLattice) = true

function sublattice_index(::KagomeLattice, site::AbstractVector{Int})
    retval = typeassert(site[end], Int)
    @assert 0 <= retval < 3 # for tripartite lattice
    return retval
end

function siteneighbors(f, lattice::KagomeLattice, ridx::Integer, ::Type{Val{1}})
    site = lattice[ridx]
    if site[end] == 0
        offsets = @SVector [@SVector([0, 0, 1]), @SVector([0, -1, 1])]
    elseif site[end] == 1
        offsets = @SVector [@SVector([0, 0, 1]), @SVector([-1, 1, 1])]
    else
        @assert site[end] == 2
        offsets = @SVector [@SVector([0, 0, -2]), @SVector([1, 0, -2])]
    end
    _siteneighbors2d(f, lattice, ridx, offsets)
end

@delegate WrappedLatticeWithBasis.lattice [ Base.length, realspace ]

#= End specific lattice w/ basis implementations =#

#= Begin cache objects =#

@doc doc"""
Precomputes the results of the translateη function on each site.

ltrc = LatticeTranslationCache(lattice, direction)
translateη(ltrc, site)
""" ->
struct LatticeTranslationCache{LatticeType<:AbstractLattice}
    lattice::LatticeType
    direction::Int
    cache::Vector{Tuple{Int,Rational{Int}}}

    function LatticeTranslationCache(lattice::LatticeType, direction::Integer) where {LatticeType<:AbstractLattice}
        cache = Tuple{Int,Rational{Int}}[]
        sizehint!(cache, length(lattice))
        for i in 1:length(lattice)
            push!(cache, translateη(lattice, i, direction))
        end
        new{LatticeType}(lattice, Int(direction), cache)
    end
end

translateη(ltrc::LatticeTranslationCache, j::Integer) = ltrc.cache[j]
translateη(ltrc::LatticeTranslationCache, site::AbstractVector{Int}) = translateη(ltrc, findfirst(ltrc.lattice, site))

#= End cache objects =#

export
    AbstractSiteNetwork,
    AbstractLattice,
    AbstractBravaisLattice,
    AbstractLatticeWithBasis,
    BravaisLattice,
    HypercubicLattice,
    ChainLattice,
    SquareLattice,
    CubicLattice,
    TriangularLattice,
    LatticeWithBasis,
    HoneycombLattice,
    KagomeLattice,
    bravais,
    isbravais,
    maxcoords,
    ndimensions,
    dimensions,
    twist,
    repeater,
    ishelical,
    primvecs,
    recivecs,
    nmomenta,
    eachmomentumindex,
    momentum,
    kdotr,
    realspace,
    wraparound_site!,
    wraparound_site,
    wraparound,
    wraparoundη,
    translate_site!,
    translate_site,
    translate,
    translateη,
    momentumspace,
    siteneighbors,
    neighbors,
    neighborsη,
    isbipartite,
    istripartite,
    sublattice_index,
    LatticeTranslationCache

end # module
