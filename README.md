# Projects
FYS3150 Projects Fall 2019
 
This repository contains the source files for all the projects in FYS3150 in the fall of 2019.
 
 ## Project 1
This project solves the Poisson equation with the use of a tridiagonal matrix. There is a general and a special algorithm, which is compared to the speed of a LU-Decomposition. To execute algorithms; run fys3150_project1_algorithms.py and follow the user interface in terminal. To see relative error i special algorithm; execute fys3150_project1_relative_error.py

## Project 2
This project solves the Schroedinger's equation for one and two electrons in the harmonic oscillator potential well and finds the eigenvalues of the systems. The algorithm used is the so called Jacobi's rotational algorithm. 

## Project 3
This report shows how the Monte Carlo integration algorithm is a superior algorithm compared to the Gauss Guadrature integration method in terms of calculating a multi-dimensional integral. The superiority is both in terms of accuracy and calculation speed. The report also discuss how both methods can be sped up and improved in terms of expressing the integrand in a different basis (spherical coordinates) and/or using a suitable Probability Density Function or implementing parallelization. The observations in this report are useful to determine which method to use when facing other integration problems later on.

## Project 5
This report demonstrates how the Variational Monte Carlo can be used with trail wave functions to find an upper bound on the ground state of the three dimensional system with two-electron in the Harmonic Oscillator potential. It shows that introducing the Jastrow factor for electron correlations to the Hydrogen wave function makes the energy drop significantly and a relative error of just 0.04898 is found. The Variational Monte Carlo method is proven to be superior to that of Jacobis rotational method with the Töplitz matrix (described in Project 2 here on GitHub), but only when choosing the right wave function.

The time consuming part of the Variational Monte Carlo method was the minimization in terms of the variational parameters. This minimization was done in a slow and brute force manner and with several parameters; this was the bottle neck. After the variational parameters was found, the expectation values was calculated in a matter of seconds.

Furthermore, the systems are shown to partly coincide with the Virial Theorem. The Hydrogen-like wave function seems to have almost perfect validity without the repulsive interactions of the electrons. The introduction of the Jastrow factor makes the Virial Theorem invalid for low frequencies in every case, but valid for high frequencies. The over all trend is them; Virial Theorem hold for high HO frequencies.
