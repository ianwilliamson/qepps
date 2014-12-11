=====================================================
QEPPS: Quadratic eigenvalue problem parameter sweeper
=====================================================

:Author:       Ian Williamson <ian.williamson@utexas.edu>
:Organization: Microelectronics Research Center, The University of Texas at Austin    


Background
----------

Comsol offers a GUI in which problems can be modeled, meshed, solved, and visualized. This is the way in which most people use the software, however Comsol also supplies a Matlab API exposing most of the features of the GUI so that sequences of operations may be scripted. This API has been successfully used in our group to perform advanced parameter sweeps and to automate complex solver sequences that would be extremely tedious in the GUI. The API also exposes Comsol’s internal linear algebra data structures such as the stiffness matrix, mass matrix, force vector, and solution vector. This means that the underlying linear algebra problem could be solved entirely in Matlab, though no advantage is typically gained by doing this, especially for large problems. QEPPS has been developed to solve a subset of the problems that we encounter in computational nanophotonics.


Motivation
----------
Modal studies in electromagnetics are quadratic eigenvalue problems. This means that they can be represented as

.. math::

  (  \lambda^2 \textbf{E} + \lambda \textbf{D} + \textbf{K}  ) \textbf{u} = 0

where **E**, **D**, and **K** are matrices, **u** is the eigenvector, and λ is the eigenvalue. Physically, **u** corresponds to the electric or magnetic field distribution over the discretized domain, with each element corresponding to one of the field components at a location within the 2D or 3D mesh. The eigenvalue, λ, corresponds to either the mode's guided effective index (in waveguides) or to the bloch wave vector in photonic crystals and other periodic geometries.

In the context of electromagnetics/optics, we are often interested in sweeping frequency to obtain broadband dispersion of the structure’s mode(s). This is useful for photonic band gap engineering, understanding signal attenuation, and many other studies.


Building
--------
QEPPS has been developed using TACC resources. Accordingly, most of the dependencies can be satisfied by loading the prepackaged TACC modules. The full list of dependencies is:

- PETSc 3.5 (complex)
- SLEPc 3.5 (complex)
- MUMPS 4.10 (complex)
- libgrvy 0.32
- LUA 5.2

The appropriate versions of PETSC, SLEPc, MUMPS, and libgrvy can all be added to the user env on at TACC with the following command::

   module load petsc/3.5-complex slepc/3.5-complex mumps/4.10.0-complex grvy/0.32.0

LUA 5.2 is included in the source of QEPPS under src/lua and the QEPPS makefile is already configured to build and link against LUA in this location.

After the dependencies have been satisfied, all that is needed to build QEPPS is the command::

   make


Test problems
-------------
The configuration LUA scripts and data files for several test problems are provided under the tests/ subdirectory. These can be run in their current form, without modification on a single TACC stampede dev node. Launcher bash scripts are also included for running each problem. Currently two problems are provided and both are relatively small; the entire parameter sweep for each should complete in less than a minute.

These can be used to validate the results that are obtained after modifying QEPPS or trying different solver options.


Usage
-----
It is highly likely that the end user will want to solve their own problems. As demonstrated by the provided test problems, a LUA script file along with command line arguments to the PETSc options database control all runtime configuration of QEPPS. This approach affords the user maximal flexibility in modifying the parameter sweep values and changing the problem configuration.

QEPPS assembles the quadratic eigenvalue problem matrices, **E**, **D**, and **K** in the following way

.. math::
  \textbf{E} = \textbf{E0} e_0(f) + \textbf{E1} e_1(f) + \textbf{E2} e_2(f) + \ldots \\

  \textbf{D} = \textbf{D0} d_0(f) + \textbf{D1} d_1(f) + \textbf{D2} d_2(f) + \ldots \\

  \textbf{K} = \textbf{K0} k_0(f) + \textbf{K1} k_1(f) + \textbf{K2} k_2(f) + \ldots

Ei, Di, and Ki are component matrices and ei(f), di(f), and ki(f) are scaling functions of the sweep parameter. The scaling functions are specified in the LUA configuration script and their evaluation is handled at run time by the embedded LUA engine. The locations of the data files are also specified in the LUA script. Additionally, various options for controling QEPPS behavior are also specified in the LUA script. Please see the example problems under the tests/ subdirectory for detailed explanations and examples.

For detailed documentation on the PETSc and SLEPc command line arguments and options, as well as the MUMPS solver, please reference the respective user manuals at

- http://www.mcs.anl.gov/petsc/petsc-3.5/docs/manual.pdf
- http://www.grycap.upv.es/slepc/documentation/slepc.pdf
- http://mumps.enseeiht.fr/doc/userguide_4.10.0.pdf