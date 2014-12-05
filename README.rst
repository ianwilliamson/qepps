=====================================================
QEPPS: Quadratic eigenvalue problem parameter sweeper
=====================================================

:Author:       Ian Williamson <ian.williamson@utexas.edu>
:Organization: Microelectronics Research Center, The University of Texas at Austin    


Background
----------
Comsol offers a GUI in which problems can be modeled, meshed, solved, and visualized. This is the way in which most people use the software, however Comsol also supplies a Matlab API exposing most of the features of the GUI so that sequences of operations may be scripted. We have successfully utilized this API in the past to perform advanced parameter sweeps and to automate complex solver sequences that would be extremely tedious in the GUI. The API also exposes Comsol’s internal linear algebra data structures such as the stiffness matrix, mass matrix, force vector, and solution vector. This means that the underlying linear algebra problem could be solved entirely in Matlab, though no advantage is typically gained by doing this, especially for large problems. 


Motivation
----------
Modal studies in electromagnetics are quadratic eigenvalue problems. This means that they can be represented as::

   (λ^2*E + λ*D + K)*U = 0

where E, D, and K are matrices, U is the eigenvector, and λ is the eigenvalue. Physically, U corresponds to the electric or magnetic field distribution over the discretized domain, with each element corresponding to one of the field components at a location within the 2D or 3D mesh. The eigenvalue, λ, corresponds to either the mode's guided effective index (in waveguides) or to the bloch wave vector in photonic crystals and other periodic geometries.

In the context of electromagnetics/optics, we are often interested in sweeping frequency to obtain broadband dispersion of the structure’s mode(s). This is useful for photonic band gap engineering, understanding signal attenuation characterization, as well as other studies.

It would be highly impractical to generate and upload the full E, D, and K matrices for every frequency value of interest, so these matrices will be componentized. Comsol allows the underlying weak form (i.e. the differential equations that are being solved) to be modified. Typically this feature is used to specify some custom physics that isn’t prepackaged by Comsol, but in this case, it will be used to factor out any expressions involving frequency from the various terms in the governing equations. This includes dispersive materials.


Building
--------
QEPPS has been developed using TACC resources. Accordingly, most of its dependencies can be satisfied by loading the prepackaged TACC modules. These include:

- PETSc 3.5 (complex)
- SLEPc 3.5 (complex)
- MUMPS 4.10 (complex)
- libgrvy 0.32
- LUA 5.2

PETSC, SLEPc, MUMPS, and libgrvy can all be loaded on TACC resources with the following command::

   module load petsc/3.5-complex slepc/3.5-complex mumps/4.10.0-complex grvy/0.32.0

LUA 5.2 is included in the source of QEPPS under src/lua and the QEPPS Makefile is already configured to build and link against this LUA.

Once the dependencies are satisfied, simply cd to the top directory and run::

   make


Usage
-----
The polynomial eigenvalue problem (PEP) context from PETSc/SLEPc 3.5 is used to solve the linear algebra problem. There are a variety of approaches for solving quadratic eigenvalue problems, including direct linearization as well as iterative approaches. PETSc also supports a variety of external solvers such as MUMPS. The solver type in PETSc can be selected and configured at run time via the options database.

The LUA scripting language is used to handle the runtime configuration for QEPPS. This allows the user to easily, and programmatically modify the parameter sweep values and can easily expand the classes of problems to be solved. QEPPS assembles the quadratic eigenvalue problem matrices, E, D, and K in the following forms::

   E = E0*e0(f) + E1*e1(f) + E2*e2(f) + ...
   D = D0*d0(f) + D1*d1(f) + D2*d2(f) + ...
   K = K0*k0(f) + K1*k1(f) + K2*k2(f) + ...

where Ei, Di, and Ki are component matrices and ei(f), di(f), and ki(f) are arbitrary functions of the sweep parameter. After parsing the input script, the LUA interpreter handles all math function evaluations of ei(f), di(f), and ki(f).
