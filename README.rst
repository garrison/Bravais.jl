=======
Bravais
=======

A package for working with lattices in `condensed matter physics <http://en.wikipedia.org/wiki/Condensed_matter_physics>`_.

Lattices are incredibly important.  Atoms in a `solid <http://en.wikipedia.org/wiki/Solid>`_ are typically arranged in a lattice, and we can even construct `optical lattices <http://en.wikipedia.org/wiki/Optical_lattice>`_.  In addition, problems in condensed matter physics on a finite lattice often have finite-size Hilbert spaces.  And to regularize a quantum field theory, it needs to be defined on a lattice!

Goals/Features
==============

We want to support a variety of different lattice types in an arbitrary number of dimensions, including:

- Bravais lattices (e.g. hypercubic, triangular) as well as lattices with a basis (e.g. honeycomb, kagome, pyrochlore).
- Arbitrary boundary conditions (open, periodic, antiperiodic, or [more generally] twisted), independently in each dimension.
- Helical boundaries, where translating the length of the lattice in one dimension also leads to an offset in another dimension(s).
- Other, user-configurable lattices, provided the primitive vectors and basis vectors are given.

For a given lattice, we want to be able to easily query/calculate the following:

- Provide a systematic enumeration of the sites of the lattice, allowing mapping between indices and site labels.
- Ability to query coordinates a lattice site in real space.
- Allows one to perform translation operations on a site, getting both the new site and what wraps were done.  This is useful for computational methods that are able to take advantage of the discrete translational symmetry of the lattice.
- Allows one to determine (and enumerate) the allowed momenta of a system, based on the boundary conditions.
- Query the nearest neighbors of a site for common lattices
- Given a momentum :math:`\mathbf{k}_i` and site :math:`\mathbf{r}_j`, calculate dot products :math:`\mathbf{k}_i \cdot \mathbf{r}_j`.
- Plotting a lattice, its reciprocal lattice, and the first `Brillouin Zone (BZ) <http://en.wikipedia.org/wiki/Brillouin_zone>`_.

Getting started
===============

To install and import::

  Pkg.clone("https://github.com/garrison/Bravais.jl.git")
  Pkg.test("Bravais")
  using Bravais

Example lattices
----------------

A 1D chain with open boundary conditions::

  HypercubicLattice([12], diagm([0]))

A 1D chain with periodic boundary conditions::

  HypercubicLattice([12])

A 2D "cylinder"::

  HypercubicLattice([12, 4], diagm([12, 0]))

A 2D triangular (aka hexagonal) lattice::

  TriangularLattice([4,6])

Honeycomb::

  HoneycombLattice([3,5])

Details
=======

A Bravais lattice in :math:`d` dimensions is defined to be all points :math:`\mathbf{R} = \sum_{i=1}^d n_i \mathbf{a}_i` where each :math:`n_i \in \mathbb{Z}`.  The vectors :math:`\mathbf{a}_i` are called the "primitive vectors" of the lattice, and there is one for each dimension in which the lattice extends.

We cannot address (or store information about) an infinite number of points on a computer, so instead we label a finite number of points on the lattice.  Then we may wish to use "periodic boundary conditions" in one or more directions to emulate an infinite number of points.  (This paragraph needs more work.)

Given positive integers :math:`N_1, \ldots, N_d` along with the primitive vectors :math:`\mathbf{a}_i`, we define a finite lattice with :math:`N_\mathrm{tot}=\prod_{i=1}^d N_i` sites, each site given by :math:`\mathbf{r}_\alpha = \sum_{i=1}^d n_i \mathbf{a}_i` where :math:`n_i \in \mathbb{Z}_{N_i}` and :math:`\alpha = 1, \ldots, N_\mathrm{tot}`.

At this point, the finite lattice only represents a subset of the original lattice sites.  But we would like to in many cases (fully open boundary conditions being the notable exception) repeat our "finite" lattice, tiling it throughout the lattice, such that every site of our original infinite lattice is represented by some site on our finite lattice.  (A diagram may help explain.)

We define vectors :math:`\mathbf{A}_1, \ldots, \mathbf{A}_d`, each of which is some linear combination of the primitive vectors, such that each site on the infinite lattice can be *uniquely* written as :math:`\mathbf{R} = \mathbf{r}_\alpha + \sum_{i=1}^d \tilde{N}_i \mathbf{A}_i`, where :math:`\tilde{N}_i \in \mathbb{Z}`.  There will often (particularly for :math:`d>1`) be infinitely many different ways of choosing these, some of which will result in "helical" (right?) boundary conditions.  (May help to give a few examples here, and diagrams!  May also be good to mention above that in a finite lattice, the actual choice of the primitive vectors is no longer arbitrary.)

The reciprocal lattice
----------------------

The reciprocal lattice is defined as all wave vectors :math:`\mathbf{K}` satisfying :math:`e^{i\mathbf{K}\cdot\mathbf{R}}=1` for all points :math:`\mathbf{R}` in the infinite Bravais lattice.  In other words, we require :math:`\mathbf{K} \cdot \mathbf{R} = 2\pi M` for some :math:`M \in \mathbb{Z}`.  This can best be achieved by choosing the reciprocal lattice's basis vectors :math:`\mathbf{b}_j` such that :math:`\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi \delta_{ij}`.  Then each :math:`\mathbf{K}` can be written as :math:`\mathbf{K} = \sum_{j=1}^d m_j \mathbf{b}_j` with :math:`m_j \in \mathbb{Z}`.  This gives :math:`\mathbf{K} \cdot \mathbf{R} = \sum_{i=1}^d\sum_{j=1}^d n_i m_j \, \mathbf{a}_i \cdot \mathbf{b}_j = 2\pi \sum_{i=1}^d m_i n_i`, which will always satisfy the original condition :math:`e^{i\mathbf{K}\cdot\mathbf{R}}=1`.  For more details, see Ashcroft and Mermin pages 86-87.

Allowed momenta
---------------

We can choose either periodic/twisted or open boundary conditions (PBC and OBC, respectively) independently in each direction :math:`\mathbf{A}_i`.  We will begin by studying the case in which the boundary conditions in each direction are periodic/twisted.  With periodic/twisted boundary conditions on a finite lattice, only certain momenta are possible in the system.  After exploring the case of fully periodic/twisted boundary conditions, we will extend our reasoning to include the (somewhat simpler) case in which one or all dimensions have open boundary conditions.

Recall that `Bloch's theorem <http://en.wikipedia.org/wiki/Bloch_wave>`_ says the eigenstates of a Hamiltonian can be chosen such that each :math:`\psi` is associated with a wave vector :math:`\mathbf{k}` such that :math:`\psi(\mathbf{r} + \mathbf{R}) = e^{i\mathbf{k} \cdot \mathbf{R}}\psi(\mathbf{r})` for every :math:`\mathbf{R}` in the lattice.  (See e.g. Ashcroft and Mermin, page 134.)  Our goal in the following is to determine, given some boundary conditions on a finite lattice, what wave vectors :math:`\mathbf{k}` are allowed.

There will be as many allowed momenta as there are points on the finite lattice ASSUMING PBC IN EACH DIRECTION.  Typically allowed momenta are given by points within the first Brillouin Zone.  We want to uniquely label them, but for simplicity we will label them systematically without the requirement that they be in the *first* Brillouin Zone.

We define the values :math:`\theta_i` such that

.. math::
   \psi(\mathbf{r} + \mathbf{A}_i) = e^{i\theta_i}\psi(\mathbf{r})

for all :math:`i`.  We can combine our knowledge that :math:`\mathbf{A}_i` is in the lattice with Bloch's theorem to give :math:`e^{i\mathbf{k} \cdot \mathbf{A}_i}\psi(\mathbf{r}) = e^{i\theta_i}\psi(\mathbf{r})`, or equivalently :math:`e^{i\left[ \mathbf{k} \cdot \mathbf{A}_i - \theta_i \right]} = 1`, for all :math:`i`.

We know that the :math:`\mathbf{A}_i`'s must be linear combinations of the primitive vectors, so we can write them as :math:`\mathbf{A}_i = \sum_{j=1}^d M_{ij} \mathbf{a}_j`, where each :math:`M_{ij}` is an integer.  (For periodic/twisted boundary conditions, our diagonal elements must be :math:`M_{ii} = N_i`, the lattice extent in each direction.  We will see later that for any dimension :math:`i` in which we have open boundary conditions, we instead have :math:`M_{ii} = 0`.)  We will also write our wave vector in terms of fractions of the reciprocal lattice's basis vectors: :math:`\mathbf{k} = \sum_{h=1}^d x_h \mathbf{b}_h`.  Then,

.. math::
   \mathbf{k} \cdot \mathbf{A}_i &= \sum_{h=1}^d \sum_{j=1}^d x_h M_{ij} \mathbf{b}_h \cdot \mathbf{a}_j \\
   &= 2\pi \sum_{j=1}^d M_{ij} x_j

With this, our requirement becomes

.. math::
   \left[ -\frac{\theta_i}{2\pi} + \sum_{j=1}^d M_{ij} x_j \right] = \tilde{n}_i

for all :math:`i`, where each :math:`\tilde{n}_i` is some nonnegative integer less than :math:`N_i`.  We define :math:`\eta_i = \theta_i/2\pi` to give

.. math::
   \left[ -\eta_i + \sum_{j=1}^d M_{ij} x_j \right] = \tilde{n}_i ,

which can also be written as a matrix equation, :math:`Mx = \tilde{n} + \eta`.

Let us assume, for vast simplification, that :math:`M_{ij}` is lower triangular (i.e. only the values for which :math:`i \ge j` are allowed to be nonzero).  (This is not a significant restriction, and in many cases the matrix will actually be diagonal.)  We can then solve the above equation iteratively for each :math:`i` beginning with :math:`i=0`.  Rewriting it with this assumption gives:

.. math::
   \sum_{j=1}^{i} M_{ij} x_j = \tilde{n}_i + \eta_i

We then solve for :math:`x_i` to give

.. math::
   x_i = \frac{1}{M_{ii}} \left[ \tilde{n}_i + \eta_i - \sum_{j=1}^{i-1} M_{ij} x_j \right]

which holds for any dimension in which there are periodic/twisted boundary conditions.

Now we briefly consider the case of open boundary conditions.  For any direction :math:`i` in which there is open boundary conditions, set :math:`M_{ij}=M_{ji}=0\ \forall j` (i.e. the corresponding row and column of the matrix :math:`M` must be zero) and :math:`\eta_i=0`.  Then :math:`x_i=0` (zero momentum) is the only unique solution (is it?) in that direction, as we expect.  For the directions in which there are periodic boundary conditions (or, more generally, twisted boundary conditions), the allowed momenta are must be determined, as we now explain.

Number of allowed moment: product over all dimensions with periodic/twisted BC's (FIXME)

For a lattice with a basis, the allowed momenta are given entirely by the underlying Bravais lattice.

Just like the lattice sites themselves, the `Bravais` package provides enumeration of the allowed momenta in a system.

Allowed total momenta
---------------------

The above considers the allowed momenta of the single particle problem.  If we have multiple particles, we may wish to determine the possible *total momenta*.  They are given as follows, where :math:`c` is the "charge" (i.e. particle count).

.. math::
   x_i^\prime = x_i + (c-1) \frac{\eta_i}{M_{ii}}

For OBC, the denominator blows up, but it should be obvious that :math:`x_i^\prime = 0`.

Generic lattice code
--------------------

OK, so what do we need to determine a lattice?  :math:`\mathbf{a}_i`, :math:`\mathbf{b}_i`, :math:`N_i`, :math:`\eta_i`, and the lower triangular matrix :math:`M_{ij}`.  Note for the diagonal elements that :math:`M_{ii} = N_i` (for periodic or twisted boundary conditions) or :math:`M_{ii} = 0` (for open boundary conditions).  We also rely on the user implementing the lattice type to specify the concept of "nearest neighbors", as what is meant by the :math:`n`'th nearest neighbors depends on the details of the lattice spacing in each direction.

Here's a table for our variables and what symbols are used in the code

+------------------------+------------------------+---------------------------------+----------------------------------+
| Symbol                 | Internal variable name |                                 | Description                      |
+========================+========================+=================================+==================================+
| :math:`N_i`            | ``N[i]``               | ``dimensions(lattice)[i]``      | lattice extent in each direction |
+------------------------+------------------------+---------------------------------+----------------------------------+
| :math:`d`              | ``d``                  | ``length(dimensions(lattice))`` | number of dimensions             |
|                        |                        | or ``ndimensions(lattice)``     |                                  |
+------------------------+------------------------+---------------------------------+----------------------------------+
| :math:`N_\mathrm{tot}` | ``N_tot``              | ``length(lattice)``             | total number of sites            |
+------------------------+------------------------+---------------------------------+----------------------------------+

And we are going to want to be able to talk about realizations of these lattice points in real space, so the following things matter.

+----------------------+------------------------+------------------------------+--------------------------------------+
| Symbol               | Internal variable name |                              | Description                          |
+======================+========================+==============================+======================================+
| :math:`\mathbf{a}_i` | ``a[:,i]``             | ``primvecs(lattice)[:,i]``   | primitive vectors                    |
+----------------------+------------------------+------------------------------+--------------------------------------+
| :math:`\mathbf{b}_i` | ``b[:,i]``             | ``recivecs(lattice)[:,i]``   | reciprocal lattice primitive vectors |
+----------------------+------------------------+------------------------------+--------------------------------------+

FIXME: also something here for the points of the different bravais sites.

As soon as we want to start talking about allowed momenta, the following two things matter as well.

+----------------+-----------------------+
| Symbol         | Variable name         |
+================+=======================+
| :math:`\eta_i` | ``eta[i]``            |
+----------------+-----------------------+
| :math:`M_{ij}` | ``repeater(i, j)``    |
+----------------+-----------------------+

Our basic ``BravaisLattice`` type contains all of these things.

We have a ``wraparound()``  (and ``wraparound!``) function, which takes a site that may or may not be on the actual finite lattice, and returns its lattice index along with the phase that it picks up.  So for instance given the site :math:`\mathbf{r}_\alpha + \mathbf{A}_i`, it returns the site index :math:`\alpha` of :math:`\mathbf{r}_\alpha` along with the phase :math:`\eta_i` picked up when [un]wrapping the boundary conditions.  As above, the phase :math:`\eta_i` returned is defined by

.. math::
   \psi(\mathbf{r} + \mathbf{A}_i) = e^{2\pi i\eta_i}\psi(\mathbf{r})

There is also a ``translation_operators()`` method, which returns a "translation operator" (really a vector meant for mapping) for each dimension in which :math:`M_{ii}` is nonzero (i.e. for each direction that is not OBC).  So, for instance, ``translation_operators()[i][alpha]`` returns the new site index :math:`\beta` (along with any phase picked up :math:`\eta`) of the site :math:`\mathbf{r}_\alpha + \mathbf{a}_i` such that

.. math::
   \psi(\mathbf{r}_\alpha + \mathbf{a}_i) = e^{2\pi i\eta}\psi(\mathbf{r}_\beta).

Wrapping condition in second quantization
-----------------------------------------

We wish to generalize the above wrapping equation to second quantization.  Note that :math:`\psi(\mathbf{r}) = \langle \mathbf{r} \vert \psi \rangle = \langle 0 \vert c_\mathbf{r} \vert \psi \rangle`.  Using this, we get

.. math::
   \psi(\mathbf{r} + \mathbf{A}_i) = \langle 0 \vert c_{\mathbf{r} + \mathbf{A}_i} \vert \psi \rangle

.. math::
   \psi(\mathbf{r} + \mathbf{A}_i) = e^{i\theta_i} \langle 0 \vert c_{\mathbf{r}} \vert \psi \rangle

Together, these imply

.. math::
   c_{\mathbf{r} + \mathbf{A}_i} &= e^{i\theta_i} c_{\mathbf{r}} \\
   c_{\mathbf{r} + \mathbf{A}_i}^\dagger &= e^{-i\theta_i} c_{\mathbf{r}}^\dagger

As a result of this,

.. math::
   T_i^L \vert \psi \rangle = e^{-i\theta_i N_c} \vert \psi \rangle

when working in second quantization.  (Explain this.)  where :math:`N_c` is the "charge" (poorly chosen name, which should be updated.)

nearest_neighbors() functions
-----------------------------

Returns (via a callback) :math:`i`, :math:`j`, and :math:`\eta`, such that the relevant hopping term would be :math:`e^{2\pi\eta}c_i^\dagger c_j`.

Specific lattice implementations
--------------------------------

Hypercubic
~~~~~~~~~~

- works in any dimension
- does not double count bonds on a two-leg ladder (fixme: do we really want this?)
- when considering nearest neighbors, do we really want it to be this general?  oh well, we can have subclasses that specialize it, since next-nearest neighbors will mean something different depending on dimension.
