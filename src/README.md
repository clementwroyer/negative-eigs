# Source files directory

Source code for the project - Approximate negative eigenvalues of Hessian matrices from CUTEst.

# Contents

## Code files

### FindBestOrder.m

MATLAB procedure to find the best sampling order among a certain set of orderings.

## findiffapprox.m 

MATLAB function to compute a finite-difference approximation of a Hessian matrix.

## HESSIANSNCFD.mat

MAT file containing all the Hessians (exact versions+finite-difference approximations).

### HessPbmsNCFinDiff

List of Eigenvalues of the matrices.

### initHessCUTEST.m

MATLAB testing script used to generate our pool of matrices - produces HESSIANSNCFD.mat and HessPbmsNCFinDiff.
Requires CUTEst/MATLAB interface.

### ListPbmsNC

Text file - List of CUTEst problems with negative curvature used in our experiments

### NES.m

Main MATLAB routine (Negative Eigenvalue Search) to find negative eigenvalues of a matrix through submatrices.

### SetOrder.m

Matlab function computing orderings for the coefficients in our matrix approximation.

### testfindiff.m

MATLAB testing script for the findiffapprox function.

### testsNES.m

MATLAB testing script for the NES function - used to produce our results.

