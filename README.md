Sampling Uniformly From the Set of Positive Definite Matrices With Trace Constraint
==================

This is the MATLAB implementation of the sampling algorithm referenced in:

Martin Mittelbach, Bho Matthiesen, and Eduard Jorswieck, "[Sampling Uniformly From the Set of Positive Definite Matrices With Trace Constraint](http://dx.doi.org/10.1109/TSP.2012.2186447)," IEEE Transactions on Signal Processing, vol. 60, no. 5, pp. 2167-2179, May 2012.


## Abstract of Article

We derive a parameterization of positive definite matrices using the Cholesky decomposition in combination with hyperspherical coordinates. Based on the parameterization we develop a simple and efficient method to randomly generate positive definite matrices with constant or bounded trace according to a uniform distribution. Further, we present an efficient implementation using the inversion method and either rejection sampling or transforming a beta distribution. The matrix parameterization might be of independent interest, whereas the random sampling algorithm finds applications in Monte Carlo simulations, testing of algorithms, and performance studies. With the help of an abstract example we describe how the sampling method can be used to approximate the optimum in a difficult, e.g., nonconvex, optimization problem for which no solution or efficient global optimization algorithm is known. In this paper we consider real as well as complex matrices.


## Acknowledgements

This research was supported in part by the Deutsche Forschungsgemeinschaft (DFG) under grant JO 801/3-1.


## License and Referencing

This code is licensed under the GPLv3 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
