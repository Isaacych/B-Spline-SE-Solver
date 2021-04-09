# B-Spline-SE-Solver
This suite of codes solve the Schroedinger equation using a spectral representation method with the B-spline functions.
It is currently coded to solve multichannel bounded systems. Algorithm for calculating the continnum state is not yet implemented.
Users can defined their own potential in mod_PEC.f90.
This code is based on the B-spline codes from Jacob Williams. Additional info can be found in mod_bspline.f90 and in https://github.com/jacobwilliams/bspline-fortran.

How to compile:
This code requires Lapack and BLAS libraries. Please replace the path for the libraries for your machine in the sample compile.sh file, if necessary.

How to use:
This code takes the input file bse.inp. The user can select the order of BSpline, define the grid for radius, define the masses of the system, and select the type of calculation.
So far, it is limited to solve bounded system. But user can choose to solve the Schroedinger equation using banded matrix or full matrix.
Wave function (eigenvector) of the system can be printed out by choosing the number of wave functions and which wave functions to be printed.
