#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#define   NULL_PTR   (void *) 0
#define   ZERO       1.0E-10
#define   MAXIT      10
#define   EPS        3.0e-14
#define   PI         3.14159265359

using namespace std;

double int_spherical_function(double, double, double, double, double, double);
void gauleg(double, double, double *, double *, int);
void gauss_laguerre(double *, double *, int, double);
double gammln(double);

double brute_force_MC(double *);
double exp_MC(double,double,double,double,double,double);





