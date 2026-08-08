# Welcome to xfemm

Welcome to the xfemm project. xfemm is a software project intended to
create a direct interface to a high quality magnetics finite element code
based on FEMM. The objective of xfemm is to create a cross-platform
command line magnetics finite element solver written in standard C++, a
set of magnetics problem definition and post-processing functions in
native Matlab/Octave code, and a mex interface to the solvers

## NOTE TO USERS

If you use xfemm, particularly for industrial work, but also academic, 
it will be greatly appreciated if you could write an email stating this 
and how it has supported your work. This is a low-cost way to ensure 
further development and maintenance will continue! Contact the authors 
on the discussion forum, or you will find an email address in the source 
files.

If you wish to cite xfemm in your work, please use the following:

Crozier, R, Mueller, M., "A New MATLAB and Octave Interface to a 
Popular Magnetics Finite Element Code", Proceedings of the 22nd 
International Conference on Electric Machines (ICEM 2016), September 
2016.

We would also suggest you cite the original FEMM program.

## Installation and Setup

There are two ways to make use of the xfemm project. One is to use it as
a collections of standalone programs run from the command line. The other
is to use it through the Matlab/Octave programming language.

The standalone programs are fmesher and fsolver. There is also a library
of post-processing functions called fpproc, but no standalone program
interface is provided to this at this time, you will have to create your
own (note that the Matlab/Octave interface does, however, provide full
access to fpproc). 

## Compiling Standalone Binary Programs

Released versions of xfemm come with pre-built binaries. But if you want
to compile xfemm on your platform, you can do so quite easily with cmake
and your compiler of choice. Run cmake on the CMakeLists.txt in the cfemm
directory to create the build system, and then build the project. On Linux
this would be done as

    cd <install dir>/xfemm/cfemm
    cmake .
    make

the binary files are found in the xfemm/cfemm/bin directory


### Side-note: Compiling with an external triangle

Xfemm ships with triangle 1.6 by Jonathan Shewchuk of the Carnegie Mellon
University, which was released in 2005.  There is a newer (unofficial) version
available, that also incorporates the aCute mesher developed at the University
of Florida, Gainesville.  This version has the main goal "to turn Triangle into
a re-usable library and the introduction of a simplified C API."

You can get the new triangle code here:
https://github.com/wo80/Triangle

If you install the external triangle to a standard location (e.g. /usr/local),
xfemm will automatically pick it up. You can notice this by the line "Found
triangle <version>" when you run cmake.  If you install the external triangle
in a different location, you'll need to tell cmake where to find it:

    cd <install dir>/xfemm/cfemm
    cmake . -DCMAKE_PREFIX_PATH=<triangle install dir>
    make

## Optional PETSc linear solver backend

The linear systems that the solvers build can optionally be solved with
PETSc instead of the built-in conjugate-gradient solver.  This gives access
to the PETSc solver/preconditioner catalogue (CG, GMRES, BiCGStab, ... ;
Jacobi, SOR, ILU, Hypre AMG, and direct solves via LU/MUMPS/SuperLU) and to
better convergence diagnostics.

A PETSc build with complex scalars is required (e.g. the Debian/Ubuntu
package `libpetsc-complex-dev`).  Enable it at configure time:

    cd <install dir>/xfemm/cfemm
    cmake . -DXFEMM_USE_PETSc=ON
    make

The PETSc backend is then selected at run time with the
`XFEMM_SOLVER_BACKEND` environment variable:

    XFEMM_SOLVER_BACKEND=petsc fsolver <problem>

The legacy solver remains the default, so existing scripts are unaffected.
The PETSc solver/PC choices can be tuned with the standard `PETSC_OPTIONS`
environment variable, e.g.:

    PETSC_OPTIONS="-ksp_type cg -pc_type hypre -pc_hypre_type boomeramg" \
        XFEMM_SOLVER_BACKEND=petsc fsolver <problem>
    PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps" \
        XFEMM_SOLVER_BACKEND=petsc fsolver <problem>

Known limitation: harmonic problems that use the Newton iteration
(`ACSolver = 1` in the problem file) are not yet supported by the PETSc
backend; use the legacy backend for those.

## Compiling Matlab Interface

Detailed instructions for compiling the Matlab inteface can be found in
the README file provided in the mfemm directory. The process is fairly well 
automated.
