# C_D_RBF_PU
A Compact Radial Basis Function Partition of Unity Method

    By: Sara Arefian and Davoud Mirzaei, September 2022

------------

This file contains the Matlab code for the rational RBF-PU method of
 
"S. Arefian, D. Mirzaei, A Compact Radial Basis Function Partition of Unity Method, Comput. Math. Appl. 2022."

------------

The user just needs to run 'StartRun.m' to see all the results in the paper, and even more examples. 

------------

-Initially, this code provides the results for 2D examples on the unit square but all functions except PolyMat and PointsInSquare work in all dimensions. 

-To generalize PolyMat for other dimensions you just need to switch over different cases for the MultiIndex vector. Or, you may use an integer partitioning algorithm to produce this vector in arbitrary dimensions for a given polynomial order. 

-Of course, similar functions can be developed by the user to produce scattered points in other dimensions. 

- As is pointed out in the paper, the polyharmoic spline kernels (PHS) are implemented in this work. However, other RBFs can be simply replaced. You just need to ignore the scaling commands in the 'LagrangeMat.m' file and implement your favorite RBF in 'Frbf.m'. 

------------
