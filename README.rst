=======
Bravais
=======

.. image:: https://travis-ci.org/garrison/Bravais.jl.svg?branch=master
    :target: https://travis-ci.org/garrison/Bravais.jl

.. image:: https://coveralls.io/repos/garrison/Bravais.jl/badge.svg
    :target: https://coveralls.io/r/garrison/Bravais.jl

A `Julia <http://julialang.org/>`_ package for working with lattices in `condensed matter physics <http://en.wikipedia.org/wiki/Condensed_matter_physics>`_.

.. NOTE:: This package is a work in progress.

Lattices are incredibly important in condensed matter physics.  Atoms in a `solid <http://en.wikipedia.org/wiki/Solid>`_ are typically arranged in a lattice, and we can even trap atoms in `optical lattices <http://en.wikipedia.org/wiki/Optical_lattice>`_.  One especially nice thing for computer simulation is that problems on a finite lattice often have finite-dimensional Hilbert spaces.  And for people who are more inclined toward thinking in the continuum: remember that to regularize a quantum field theory, it needs to be able to be defined on a lattice!

Documentation at http://bravaisjl.readthedocs.org/

Goals/Features
==============

We want to support a variety of different lattice types in an arbitrary number of dimensions, including:

- Bravais lattices (e.g. hypercubic, triangular, face-centered cubic, body-centered cubic) as well as lattices with a basis (e.g. honeycomb, kagome, diamond, pyrochlore).
- Arbitrary boundary conditions (open, periodic, antiperiodic, or [more generally] twisted), independently in each dimension.
- Helical boundaries, where translating the length of the lattice in one dimension also leads to an offset in another dimension(s).
- Other, user-configurable lattices, provided the primitive vectors and basis vectors are given.

For a given lattice, we want to be able to easily query/calculate the following:

- Provide a systematic enumeration of the sites of the lattice, allowing mapping between indices and site labels (aka coordinates in terms of the primitive vectors).
- Ability to query the actual real-space coordinates of each lattice site.
- Ability to perform translation operations on a site, getting both the new site and what wrap(s) around the boundary were necessary.  This is useful for computational methods that are able to take advantage of the discrete translational symmetry of the lattice.
- Ability to determine (and enumerate) the allowed momenta of a system, based on the boundary conditions.
- Mechanism for querying the nearest neighbors of a site for common lattices.
- Given a momentum :math:`\mathbf{k}_i` and site :math:`\mathbf{r}_j`, calculate dot products :math:`\mathbf{k}_i \cdot \mathbf{r}_j`.
- Routines for plotting a lattice, its reciprocal lattice, and the first `Brillouin Zone (BZ) <http://en.wikipedia.org/wiki/Brillouin_zone>`_.

Getting started
===============

To install, test, and import::

  Pkg.clone("https://github.com/garrison/Bravais.jl.git")
  Pkg.test("Bravais")
  using Bravais

Example lattices
----------------

Lattice construction is fairly simple.  Some examples below.

A 1D chain with open boundary conditions::

  ChainLattice([12], Diagonal([0]))

A 1D chain with periodic boundary conditions::

  ChainLattice([12])

A 2D square lattice on a "cylinder"::

  SquareLattice([12,4], Diagonal([12,0]))

A 2D triangular (aka hexagonal) lattice with PBC::

  TriangularLattice([4,6])

Honeycomb with PBC::

  HoneycombLattice([3,5])

4D hypercube::

  HypercubicLattice{4}([8,8,8,8])

Documentation
=============

See http://bravaisjl.readthedocs.org/
