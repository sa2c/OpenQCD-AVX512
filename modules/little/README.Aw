
********************************************************************************

                   Definition of the little Dirac operator

********************************************************************************

The little Dirac operator was introduced in

  M. Luescher: "Local coherence and deflation of the low quark modes
                in lattice QCD", JHEP 0707 (2007) 081

Here the data structures used to represent the operator are described.


Definition
----------

The deflation subspace is spanned by the fields (*b).sd[1],..,(*b).sd[Ns] on
the blocks b of the DFL_BLOCKS grid. When the subspace is created, the basis
fields are orthonormalized on each block. The restriction of the Wilson-Dirac
operator Dw+i*mu*gamma_5 to the deflation subspace is referred to as the
little Dirac operator. It is completely specified by the matrix elements

  A_{n,k;m,l}=(v_{n,k},(Dw+i*mu*gamma_5)*v_{m,l})

where v_{n,0},v_{n,1},..,v_{n,Ns-1} are the basis vectors on the block with
index n.


Matrix arrays
-------------

The DFL_BLOCKS grid consists of the sublattices of the even and odd blocks
(see dfl/dfl_geometry.c). In each local lattice, there are nb blocks, half of
which are even and half odd. The number nbb of blocks on the exterior boundary
of the local lattice also divides into equal numbers of even and odd blocks.

The matrix A_{n,k;m,l} decomposes into four parts Aee, Aoo, Aoe and Aeo in the
obvious way. Each of these parts may be stored in the form of one-dimensional
arrays of complex Ns x Ns matrices. Explicitly

  Aee[n][Ns*k+l] = (v_{n,k},(Dw+i*mu*gamma_5)*v_{n,l}),

  Aoo[n][Ns*k+l] = (v_{m,k},(Dw+i*mu*gamma_5)*v_{m,l}), m=n+nb/2,

  Aoe[8*n+ifc][Ns*k+l] = (v_{m,k},(Dw+i*mu*gamma_5)*v_{inn[m][ifc],l}),

  Aeo(8*n+ifc][Ns*k+l] = (v_{inn[m][ifc],k},(Dw+i*mu*gamma_5)*v_{m,l}),

where n=0,..,nb/2-1 labels the even blocks, m=nb/2,..,nb-1 the odd blocks,
ifc=0,..,7 the 8 coordinate directions -0,+0,..,-3,+3, while inn[m][ifc] is
the index of the block in direction ifc of the block with index m.

In the case of the double-precision operator, the length of the arrays Aoe and
Aeo is 4*nb+nbb/2 rather than 4*nb. The additional nbb/2 elements at end of
the arrays are used at intermediate stages of the computations as buffers for
the matrices on the odd exterior boundary points of the block lattice. These
are stored in the order of the boundary points (see dfl/dfl_geometry.c and
README.Aw_com).


Even-odd flag
-------------

The even-odd flag can be set and unset by calling set_tm_parms() (see
flags/lat_parms.c). Initially the flag is not set.

The programs for the Dirac operator and thus those that construct the little
Dirac operator apply the twisted mass term i*mu*gamma_5 on the even sites of
the lattice only if the flag is set. The associated deflation projectors are
suitable for the solution of the Dirac operator with such a twisted-mass term.


Data structure
--------------

The single- and double-precision arrays representing the little Dirac operator
are collected in the structures Aw_t and Aw_dble_t (see include/little.h). The
elements of these structures are

 Ns,nb
 Aee[nb/2][Ns*Ns]
 Aoo[nb/2][Ns*Ns]
 Aoe[4*nb][Ns*Ns]
 Aeo[4*nb][Ns*Ns]

As already mentioned, the length of the last two arrays is 4*nb+nbb/2
rather than 4*nb in the case of the double precision operator.


Even-odd preconditioned operator
--------------------------------

The even-odd preconditioned little operator Ahat acts on fields supported on
the even blocks. It is related to the little operator A through

  Ahat=1-Aee^(-1)*Aeo*Aoo^(-1)*Aoe

The preconditioned operator may be represented by Aw_t and Aw_dble_t
structures containing the matrix arrays

  Aee^(-1), Aoo^(-1), Aee^(-1)*Aeo and Aoo^(-1)*Aoe

instead of Aee, Aoo, Aoe and Aeo.

In the case of the double-precision preconditioned operator, the array of the
even-even terms is of length nb/2+nbb/2 instead of nb/2. The additional nbb/2
elements at end of the array are used as buffers for the matrices on the even
interior boundary points of the block lattice (see README.Aw_com).


Little-little Dirac operator
----------------------------

The deflation subspace is constructed by projecting Ns global spinor fields to
the blocks of the DFL_BLOCKS grid. These global fields are linear combinations
of the basis fields v_{n,k} and span a subspace of dimension Ns within the
deflation subspace (which has dimension nb*Ns).

The even-odd preconditioned little Dirac operator may be deflated using the
restriction of the global modes to the even blocks as the deflation modes.
These fields (which are also contained in the deflation subspace) are referred
to as the little modes, and the restriction of the even-odd preconditioned
little Dirac operator to the space spanned by them as the little-little Dirac
operator. Its action is completely specified by its matrix elements in the
space of the little modes, i.e. by a complex Ns x Ns matrix (the program
set_ltl_modes() orthonormalizes the little modes before the little-little
Dirac operator is calculated).

The single-precision little modes md_k (k=0,..,Ns-1) and Awhat*md_k are stored
in the first and second half of the first Ns fields returned by vflds(). The
double-precision fields are stored in the same way in the Ns fields returned
by vdflds().


Boundary conditions
-------------------

In the case of boundary conditions of type 0,1 and 2, the hopping terms Aeo
and Aoe that go across the boundaries of the lattice at global time 0 and
NPROC0*L0-1 are equal to zero.

The programs in this directory obtain the little Dirac operator always in the
same way and thus effectively as if periodic boundary conditions were imposed
in the time direction. Since the quark fields vanish at time 0 when boundary
conditions of type 0,1 or 2 are chosen, the calculation gives the correct
result (i.e. vanishing hopping terms across the lattice boundaries) also in
these cases.
