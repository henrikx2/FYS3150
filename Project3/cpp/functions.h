#include <iostream>
#include <new>
#include <cmath>
#include <cstring>

#define   ZERO       1.0E-7
#define   MAXIT      10
#define   EPS        3.0e-14
#define   PI         acos(-1.)
#define   EXACT      5*PI*PI/(16.0*16.0)

using namespace std;

double int_spherical_function(double, double, double, double, double, double);
void gauleg(double, double, double *, double *, int);
void gauss_laguerre(double *, double *, int, double);
double gammln(double);

double brute_force_MC(double *);
double exp_MC(double,double,double,double,double,double);





