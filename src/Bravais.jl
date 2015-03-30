# Bravais.jl -- Copyright 2015 James R. Garrison and contributors.
# Provided under the MIT License.

# NOTE: Most symbols used below are defined in the accompanying
# documentation.

# FIXME: perhaps change to use tuples, not arrays.  they are ordered,
# but they are also immutable (which might be bad for some things, but
# we can always convert to array and back.)  Or we could just use
# statically-sized arrays eventually (which would require templating
# on d), which seems more appropriate.  Then anyone who wants to order
# on them must convert to a tuple, but this is fine..

# FIXME: need other constructors for e.g. OBC

# TODO: neighbors(f, lattice, :first)

# FIXME: have a way of getting the neighbors of a given site (instead
# of iterating through all of them).  for this, we may actually want
# to double count (or not).

# FIXME: a wraparound function that returns the number of wraps
# instead of the phase picked up by the wraps.  Same with translate.

module Bravais

VERSION < v"0.4-" && using Docile

# This macro was originally written by JMW for DataStructures.jl.
# It has been modified for use here.
macro delegate(source, targets)
    funcnames = map(x -> esc(x), targets.args)
    fdefs = Any[]
    typename = esc(source.args[1])
    fieldname = esc(Expr(:quote, source.args[2].args[1]))
    for funcname in funcnames
        push!(fdefs, quote
                         ($funcname)(a::($typename), args...) =
                             ($funcname)(a.($fieldname), args...)
                     end)
    end
    return Expr(:block, fdefs...)
end

# Elementwise floating division is not (yet) implemented in Julia; see
# https://github.com/JuliaLang/julia/issues/9113
elementwise_fld(a, b) = [fld(x, y) for (x, y) in zip(a, b)]

# Julia currently lacks row-major versions of sub2ind and ind2sub.  We
# want to use these so the ordering for returned tuples is correct;
# see discussion at <https://github.com/JuliaLang/julia/pull/10337#issuecomment-78206799>.
rowmajor_sub2ind(dims, indices...) = sub2ind(reverse(dims), reverse(indices)...)
rowmajor_ind2sub(dims, index) = reverse(ind2sub(reverse(dims), index))

abstract AbstractLattice <: AbstractVector{Vector{Int}}

abstract AbstractBravaisLattice <: AbstractLattice
abstract AbstractLatticeWithBasis <: AbstractLattice

abstract WrappedBravaisLattice <: AbstractBravaisLattice
abstract WrappedLatticeWithBasis <: AbstractLatticeWithBasis

immutable BravaisLattice <: AbstractBravaisLattice
    N_tot::Int  # total number of sites
    N::Vector{Int}  # lattice extent in each dimension
    M::Matrix{Int}  # "M" matrix (will be diagonal for non-helical boundary)
    η::Vector{Rational{Int}}  # twist in each direction (rational multiple of 2π)
    a::Matrix{Float64}  # primitive lattice vectors
    b::Matrix{Float64}  # reciprocal lattice vectors
    momenta::Matrix{Rational{Int}}  # FIXME: rename x

    function BravaisLattice(N::Vector{Int},
                            M::Matrix{Int}=diagm(N), # assumes pbc
                            η::Vector{Rational{Int}}=zeros(Rational{Int}, length(N)),
                            a::Matrix{Float64}=eye(length(N)))

        # check N
        @assert all(d_i -> d_i>0, N)
        d = length(N)
        N_tot = prod(N)

        # check M
        @assert size(M) == (d, d)
        @assert M == tril(M)
        for i in 1:d
            @assert M[i,i] == 0 || M[i,i] == N[i]
            if M[i,i] != N[i]
                @assert all(s->(s == 0), M[i,:])
                @assert all(s->(s == 0), M[:,i])
            end
        end

        # check η
        @assert length(η) == d
        @assert all(η_i -> 0<=η_i<1, η)
        @assert all([M[i,i] != 0 || η[i] == 0 for i in 1:d])

        # check a and generate b
        @assert size(a) == (d, d)
        b = 2pi * transpose(inv(a))

        # generate allowed momenta
        momenta_range = [s != 0 ? s : 1 for s in diag(M)]
        nmomenta = prod(momenta_range)
        momenta = Array(Rational{Int}, d, nmomenta)
        for idx in 1:nmomenta
            # FIXME: getting idx should be easier, once
            # multidimensional iteration is supported
            ñ = [ind2sub(tuple(momenta_range...), idx)...] - 1
            for i in 1:d
                if M[i,i] == 0
                    momenta[i,idx] = 0
                else
                    tmp = sum([M[i,j] * momenta[j,idx] for j in 1:i-1])
                    momenta[i,idx] = (ñ[i] + η[i] - tmp) // M[i,i]
                end
            end
        end

        new(N_tot, copy(N), copy(M), copy(η), copy(a), b, momenta)
    end
end

immutable LatticeWithBasis <: AbstractLatticeWithBasis
    N_tot::Int
    maxcoords::Vector{Int}
    bravaislattice::BravaisLattice
    basis::Matrix{Float64}

    function LatticeWithBasis(N::Vector{Int},
                              M::Matrix{Int}=diagm(N), # assumes pbc
                              η::Vector{Rational{Int}}=zeros(Rational{Int}, length(N)),
                              a::Matrix{Float64}=eye(length(N)),
                              basis::Matrix{Float64}=zeros(length(N), 1))
        bravaislattice = BravaisLattice(N, M, η, a)

        # check basis
        nbasis = size(basis)[2]
        @assert nbasis > 0
        @assert size(basis)[1] == length(N)

        # determine N_tot and maxcoords, now that we know the basis size
        N_tot = prod(N) * nbasis
        maxcoords = vcat(N, [nbasis])

        new(N_tot, maxcoords, bravaislattice, copy(basis))
    end
end

isbravais(lattice::AbstractBravaisLattice) = true
isbravais(lattice::AbstractLatticeWithBasis) = false

@doc doc"Returns the underlying Bravais lattice" ->
bravais(lattice::AbstractBravaisLattice) = lattice
bravais(lattice::LatticeWithBasis) = lattice.bravaislattice
bravais(lattice::WrappedLatticeWithBasis) = bravais(lattice.lattice)

# This is not exported (at the moment).  It is used to simplify some
# code below.
maxcoords(lattice::BravaisLattice) = lattice.N
maxcoords(lattice::LatticeWithBasis) = lattice.maxcoords
maxcoords(lattice::WrappedBravaisLattice) = maxcoords(lattice.lattice)
maxcoords(lattice::WrappedLatticeWithBasis) = maxcoords(lattice.lattice)

Base.length(lattice::BravaisLattice) = lattice.N_tot
Base.length(lattice::LatticeWithBasis) = lattice.N_tot

Base.size(lattice::AbstractLattice) = (length(lattice),)

function Base.getindex(lattice::AbstractLattice, index::Integer)
    checkbounds(lattice, index)
    return [rowmajor_ind2sub(tuple(maxcoords(lattice)...), index)...] - 1
end

Base.findfirst(lattice::AbstractLattice, site::Vector{Int}) = rowmajor_sub2ind(tuple(maxcoords(lattice)...), (site + 1)...)

Base.in(site::Vector{Int}, lattice::AbstractLattice) = length(maxcoords(lattice)) == length(site) && all(0 .<= site .< maxcoords(lattice))

Base.start(lattice::AbstractLattice) = zeros(Int, length(maxcoords(lattice)))
Base.done(lattice::AbstractLattice, site::Vector{Int}) = !all(site .< maxcoords(lattice))
function Base.next(lattice::AbstractLattice, site::Vector{Int})
    newsite = copy(site)
    mc = maxcoords(lattice)
    # XXX FIXME @assert something # otherwise BoundsError!
    @assert length(mc) >= 1
    for i in length(mc):-1:2
        newsite[i] += 1
        if newsite[i] == mc[i]
            newsite[i] = 0
        else
            return site, newsite
        end
    end
    newsite[1] += 1
    return site, newsite
end

# FIXME: can we return readonly views of some of the following things?

# FIXME: define (some of) these for LatticeWithBasis.
dimensions(lattice::BravaisLattice) = lattice.N  # FIXME: rename this `extent`?
ndimensions(lattice::BravaisLattice) = length(lattice.N)
twist(lattice::BravaisLattice) = lattice.η
repeater(lattice::BravaisLattice) = lattice.M

# We intentionally define many of the following functions for Bravais
# lattices only.  If one wants to query them for a lattice with a
# basis, call e.g. primvecs(bravais(lattice)).

@doc doc"Returns the primitive vectors of the Bravais lattice" ->
primvecs(lattice::BravaisLattice) = lattice.a

@doc doc"Returns the primitive vectors of the reciprocal lattice" ->
recivecs(lattice::BravaisLattice) = lattice.b

momentum(lattice::BravaisLattice, idx) = lattice.momenta[:, idx]
function momentum(lattice::BravaisLattice, idx, charge::Int)
    # "total momentum", really.  note that this may return things greater than one.
    rv = momentum(lattice, idx)
    if charge == 1
        return rv
    end
    for i in 1:length(lattice.N)
        if lattice.M[i,i] != 0
            rv[i] += lattice.η[i] * (charge - 1) // lattice.M[i,i]
        end
    end
    return rv
end

kdotr(lattice::BravaisLattice, kidx, r::Vector{Int}) = 2pi * dot(momentum(lattice, kidx), r)
kdotr(lattice::BravaisLattice, kidx, ridx::Integer) = kdotr(kidx, lattice[ridx])

momentumspace(lattice::BravaisLattice, k::Vector{Int}) = lattice.b * k
momentumspace(lattice::BravaisLattice, kidx::Integer) = momentumspace(lattice, momentum(lattice, kidx)) 
# FIXME: make it possible to get points translated to the first brillouin zone

realspace(lattice::BravaisLattice, r::Vector{Int}) = lattice.a * r
realspace(lattice::BravaisLattice, ridx::Integer) = realspace(lattice, lattice[ridx])
# FIXME: realspace with M offset

# XXX FIXME: test this!!
function realspace(lattice::LatticeWithBasis, r::Vector{Int})
    @assert length(r) == length(lattice.maxcoords)
    return realspace(bravais(lattice), r[1:end-1]) + lattice.basis[:, r[end]+1]
end
realspace(lattice::LatticeWithBasis, ridx::Integer) = realspace(lattice, lattice[ridx])

function wraparound!(lattice::BravaisLattice, site::Vector{Int})
    d = ndimensions(lattice)
    @assert length(site) == d

    wraps = zeros(Int, d)

    for i in d:-1:1
        if !(0 <= site < lattice.N[i])
            if lattice.M[i,i] != 0
                # periodic/twisted BC's in this direction
                wraps[i] = fld(site[j], lattice.N[i])
                for j in 1:i
                    site[j] -= wraps[i] * lattice.M[i,j]
                end
            else
                # OBC in this direction
                #
                # XXX: This does not provide any exception guarantee.
                # If we wanted, we could easily undo what we had done
                # before throwing.
                throw(ArgumentError("This lattice has open boundary conditions in the direction we are trying to wraparound!"))
            end
        end
    end

    # figure out the phase picked up
    η = wraps .* lattice.η

    return site, η
end

wraparound(lattice::BravaisLattice, site::Vector{Int}) = wraparound!(lattice, copy(site))
wraparound(lattice::BravaisLattice, index::Integer) = wraparound(lattice, lattice[index])

# FIXME: wraparound for lattice w/ basis

function translate!(lattice::BravaisLattice, site::Vector{Int}, direction::Integer)
    # FIXME: possibly refuse to translate if M[i,i] == 0
    site[direction] += 1
    return wraparound!(lattice, site)
end
translate(lattice::BravaisLattice, site::Vector{Int}, direction::Integer) = translate!(lattice, copy(site), direction)
# FIXME: make the following one return (newridx, η)
translate(lattice::BravaisLattice, ridx::Integer, direction::Integer) = translate(lattice, lattice[ridx], direction)

#= End BravaisLattice =#

#= Begin specific Bravais lattice implementations =#

immutable HypercubicLattice <: WrappedBravaisLattice
    lattice::BravaisLattice

    function HypercubicLattice(N::Vector{Int},
                               M::Matrix{Int}=diagm(N), # assumes pbc
                               η::Vector{Rational{Int}}=zeros(Rational{Int}, length(N)))
        new(BravaisLattice(N, M, η))
    end
end

function nearestneighbors(f, lattice::HypercubicLattice)
    # NOTE: iterates each bond once.
    for site in lattice
        for i in 1:ndimensions(lattice)
            newsite = copy(site)
            newsite[i] += 1
            # FIXME: instead of `wrap`, why not eta?  well, "wrap"
            # allows me to do something special for 2-leg ladder.
            #
            # FIXME: in directions w/ OBC, do not provide the neighbors!
            #
            # FIXME: also, this elementwise_fld does not consider
            # helical boundaries (correct?)
            wrap = elementwise_fld(newsite, dimensions(lattice))
            newsite -= wrap .* dimensions(lattice)
            f(findfirst(lattice, site), findfirst(lattice, newsite), wrap)
        end
    end
end

immutable TriangularLattice <: WrappedBravaisLattice
    lattice::BravaisLattice

    function TriangularLattice(N::Vector{Int},
                               M::Matrix{Int}=diagm(N), # assumes pbc
                               η::Vector{Rational{Int}}=zeros(Rational{Int}, length(N)))
        @assert length(N) == 2
        new(BravaisLattice(N, M, η, [1.0 0; 0.5 sqrt(3)/2]'))
    end
end

function nearestneighbors(f, lattice::TriangularLattice)
    # FIXME: make it work like the hypercubic thing above
    for site in lattice
        for offset in ([1,0], [0,1], [-1,1])
            newsite = site + offset
            wrap = elementwise_fld(newsite, dimensions(lattice))
            newsite -= wrap .* dimensions(lattice)
            f(findfirst(lattice, site), findfirst(lattice, newsite), wrap)
        end
    end
end

@delegate WrappedBravaisLattice.lattice [ Base.length, dimensions, ndimensions, twist, repeater, primvecs, recivecs, momentum, kdotr, momentumspace, realspace, wraparound!, wraparound, translate!, translate ]

#= End specific Bravais lattice implementations =#

#= Begin specific lattice w/ basis implementations =#

immutable HoneycombLattice <: WrappedLatticeWithBasis
    lattice::LatticeWithBasis

    function HoneycombLattice(N::Vector{Int},
                              M::Matrix{Int}=diagm(N), # assumes pbc
                              η::Vector{Rational{Int}}=zeros(Rational{Int}, length(N)))
        @assert length(N) == 2
        a = [1.5 sqrt(3)/2; 0 sqrt(3)]'
        basis = [0 0; 1.0 0]'
        return new(LatticeWithBasis(N, M, η, a, basis))
    end
end

immutable KagomeLattice <: WrappedLatticeWithBasis
    lattice::LatticeWithBasis

    function KagomeLattice(N::Vector{Int},
                           M::Matrix{Int}=diagm(N), # assumes pbc
                           η::Vector{Rational{Int}}=zeros(Rational{Int}, length(N)))
        @assert length(N) == 2
        a = [2 0; 1 sqrt(3)]'
        basis = [0 0; 0.5 sqrt(3)/2; 1.0 0]'
        return new(LatticeWithBasis(N, M, η, a, basis))
    end
end

@delegate WrappedLatticeWithBasis.lattice [ Base.length, realspace ]

#= End specific lattice w/ basis implementations =#

export AbstractLattice, BravaisLattice, HypercubicLattice, TriangularLattice, LatticeWithBasis, HoneycombLattice, KagomeLattice, dimensions, ndimensions, momentum, kdotr, realspace, momentumspace, nearestneighbors

end # module
