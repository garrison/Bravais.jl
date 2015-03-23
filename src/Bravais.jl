# Bravais.jl -- Copyright 2015 James R. Garrison and contributors.
# Provided under the MIT License.

# NOTE: Most symbols used below are defined in the accompanying
# documentation.

module Bravais

VERSION < v"0.4-" && using Docile

# This macro was originally written by JMW for DataStructures.jl.
# It has been modified for use here.
macro delegate(sources, targets)
    funcnames = map(x -> esc(x), targets.args)
    fdefs = Any[]
    for source in sources.args
        typename = esc(source.args[1])
        fieldname = esc(Expr(:quote, source.args[2].args[1]))
        for funcname in funcnames
            push!(fdefs, quote
                             ($funcname)(a::($typename), args...) =
                                 ($funcname)(a.($fieldname), args...)
                         end)
        end
    end
    return Expr(:block, fdefs...)
end

abstract AbstractLattice
abstract AbstractBravaisLattice <: AbstractLattice
abstract AbstractLatticeWithBasis <: AbstractLattice

@doc doc"Returns the underlying Bravais lattice" ->
bravais(lattice::AbstractBravaisLattice) = lattice

# Elementwise floating division is not (yet) implemented in Julia; see
# https://github.com/JuliaLang/julia/issues/9113
elementwise_fld(a, b) = [fld(x, y) for (x, y) in zip(a, b)]

#= Begin BravaisLattice =#

immutable type BravaisLattice <: AbstractLattice
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

# FIXME: need other constructors for e.g. OBC

function Base.getindex(lattice::BravaisLattice, index)
    if !(0 < index <= lattice.N_tot)
        throw(BoundsError())
    end
    return [ind2sub(tuple(lattice.N...), index)...] - 1
end


#function Base.findfirst(lattice::BravaisLattice, site)
#
#end
#Base.find()

# XXX FIXME: possibly Base.checkbounds()

Base.length(lattice::BravaisLattice) = lattice.N_tot
Base.size(lattice::BravaisLattice) = (lattice.N_tot,)

Base.start(lattice::BravaisLattice) = zeros(Int, length(lattice.N))
Base.done(lattice::BravaisLattice, state) = !all(state .< lattice.N)
function Base.next(lattice::BravaisLattice, state)
    newstate = copy(state)
    # XXX FIXME @assert something # otherwise BoundsError!
    @assert length(lattice.N) >= 1
    for i in 1:length(lattice.N)-1
        newstate[i] += 1
        if newstate[i] == lattice.N[i]
            newstate[i] = 0
        else
            return state, newstate
        end
    end
    newstate[length(lattice.N)] += 1
    return state, newstate
end

# FIXME: methods to access the "members" (N, η, M)
dimensions(lattice::BravaisLattice) = lattice.N
ndimensions(lattice::BravaisLattice) = length(lattice.N)

# FIXME: can we return readonly views of these?
primvecs(lattice::BravaisLattice) = lattice.a
recivecs(lattice::BravaisLattice) = lattice.b

# FIXME: implement in()  (how? which way?)

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

realspace(lattice::BravaisLattice, r::Vector{Int}) = lattice.a * r
realspace(lattice::BravaisLattice, ridx::Integer) = realspace(lattice, lattice[ridx])
# FIXME: realspace with M offset

momentumspace(lattice::BravaisLattice, k::Vector{Int}) = lattice.b * k
momentumspace(lattice::BravaisLattice, kidx::Integer) = momentumspace(lattice, momentum(lattice, kidx)) 
# FIXME: make it possible to get points translated to the first brillouin zone

#= End BravaisLattice =#

#= Begin specific Bravais lattice implementations =#

immutable type HypercubicLattice <: AbstractBravaisLattice
    lattice::BravaisLattice

    function HypercubicLattice(N::Vector{Int},
                               M::Matrix{Int}=diagm(N), # assumes pbc
                               η::Vector{Rational{Int}}=zeros(Rational{Int}, length(N)))
        new(BravaisLattice(N, M, η))
    end
end

function nearestneighbors(f, lattice::HypercubicLattice)
    for site in lattice
        for i in 1:ndimensions(lattice)
            newsite = copy(site)
            newsite[i] += 1
            wrap = elementwise_fld(newsite, dimensions(lattice))
            newsite -= wrap .* dimensions(lattice)
            f(findfirst(lattice, site), findfirst(lattice, newsite), wrap)
        end
    end
end

immutable type TriangularLattice <: AbstractBravaisLattice
    lattice::BravaisLattice

    function TriangularLattice(N::Vector{Int},
                               M::Matrix{Int}=diagm(N), # assumes pbc
                               η::Vector{Rational{Int}}=zeros(Rational{Int}, length(N)))
        @assert length(N) == 2
        new(BravaisLattice(N, M, η, [1.0 0; 0.5 sqrt(3)/2]'))
    end
end

function nearestneighbors(f, lattice::TriangularLattice)
    for site in lattice
        for offset in ([1,0], [0,1], [-1,1])
            newsite = site + offset
            wrap = elementwise_fld(newsite, dimensions(lattice))
            newsite -= wrap .* dimensions(lattice)
            f(findfirst(lattice, site), findfirst(lattice, newsite), wrap)
        end
    end
end

@delegate [ HypercubicLattice.lattice, TriangularLattice.lattice ] [ Base.getindex, Base.length, Base.size, Base.start, Base.done, Base.next, dimensions, ndimensions, momentum, kdotr, realspace, momentumspace ]

#= End specific Bravais lattice implementations =#

#= Begin LatticeWithBasis =#

immutable type LatticeWithBasis <: AbstractLatticeWithBasis
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

bravais(lattice::LatticeWithBasis) = lattice.bravaislattice

function Base.getindex(lattice::LatticeWithBasis, index)
    @assert 0 < index <= lattice.N_tot
    return [ind2sub(tuple(lattice.maxcoords...), index)...] - 1
end

Base.length(lattice::LatticeWithBasis) = lattice.N_tot
Base.size(lattice::LatticeWithBasis) = (lattice.N_tot,)

# XXX FIXME: test this!!
function realspace(lattice::LatticeWithBasis, r::Vector{Int})
    @assert length(r) == length(lattice.maxcoords)
    return realspace(bravais(lattice), r[1:end-1]) + lattice.basis[:, r[end]+1]
end
realspace(lattice::LatticeWithBasis, ridx::Integer) = realspace(lattice, lattice[ridx])
# FIXME: realspace with M offset

#= End LatticeWithBasis =#

#= Begin specific lattice w/ basis implementations =#

immutable type HoneycombLattice <: AbstractLatticeWithBasis
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

immutable type KagomeLattice <: AbstractLatticeWithBasis
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

@delegate [ HoneycombLattice.lattice, KagomeLattice.lattice ] [ Base.getindex, Base.length, Base.size, realspace, bravais ]

#= End specific lattice w/ basis implementations =#

export AbstractLattice, BravaisLattice, HypercubicLattice, TriangularLattice, LatticeWithBasis, HoneycombLattice, KagomeLattice, dimensions, ndimensions, momentum, kdotr, realspace, momentumspace, bravais, nearestneighbors

end # module
