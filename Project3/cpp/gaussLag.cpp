//   This is a simple program which tests Gaussian quadrature using
//   Legendre and Laguerre polynomials
//   It integrates the simple function x* exp(-x) for the interval
//   x \in [0,infty). The exact result is 1. For Legendre based quadrature a
//   tangent mapping is also used.

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include "functions.h"
using namespace std;
using namespace std::chrono;

double int_spherical_function(double, double, double, double, double, double);

//     Here we define various functions called by the main program

//   Main function begins here
int main()
{
     int n;
     cout << "Read in the number of integration points" << endl;
     cin >> n;
//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method

     double *x_t = new double [n];
     double *w_t = new double [n];
     double *x_p = new double [n];
     double *w_p = new double [n];
     double *x_r = new double [n+1];
     double *w_r = new double [n+1];

//   Start timing
     high_resolution_clock::time_point time1 = high_resolution_clock::now();

//   set up the mesh points and weights
     gauleg(0, PI, x_t, w_t, n);
     gauleg(0, 2*PI, x_p, w_p, n);

//   set up the mesh points and weights and the power of x^alf
     double alf = 0.0;
     double alpha = 2.0; //For determinant
     gauss_laguerre(x_r, w_r, n, alf);

//   Initializing the brute force gauleg sum
     double int_gauss = 0.;
//   Run six for-loops to evaluate every integral
     int i,j,k,l,m,o;
     for (i = 1;i<n+1;i++){
     for (j = 1;j<n+1;j++){
     for (k = 0;k<n;k++){
     for (l = 0;l<n;l++){
     for (m = 0;m<n;m++){
     for (o = 0;o<n;o++){
     int_gauss += w_r[i]*w_r[j]*w_t[k]*w_t[l]*w_p[m]*w_p[o]*int_spherical_function(x_r[i],x_r[j],x_t[k],x_t[l],x_p[m],x_p[o]);
     }}}}}}

//    final time
      high_resolution_clock::time_point time2 = high_resolution_clock::now();
      double time_used = duration_cast<nanoseconds>(time2-time1).count();

//    final output
      cout << setiosflags(ios::showpoint | ios::uppercase);
      cout << "Gaussian-Laguerre Quadrature = "<< setw(20) << setprecision(15)  << int_gauss << endl;
      cout << "Error = "<< setw(20) << setprecision(15)  << abs((5*PI*PI/(16.0*16.0))-int_gauss) << endl;
      cout << "Time = "<< setw(20) << setprecision(15)  << time_used/(1.0E+9) << " s" << endl;
      delete [] x_r;
      delete [] w_r;
      delete [] x_t;
      delete [] w_t;
      delete [] x_p;
      delete [] w_p;
      return 0;
}  // end of main program

//  this function defines the function to integrate
double int_spherical_function(double r1, double r2, double t1, double t2, double p1, double p2)
{
   
   double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
   double lent = sqrt(r1*r1+r2*r2-2.0*r1*r2*cosb);
   double f = exp(-3*(r1+r2))*r1*r1*r2*r2*sin(t1)*sin(t2)/lent;
   if(lent > ZERO)
   return f;
   else 
   return 0;
   /*
   double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
   double lent = sqrt(r1*r1+r2*r2-2.0*r1*r2*cosb);
   double f = sin(t1)*sin(t2)/lent;
   if(lent > ZERO)
   return f;
   else 
   return 0;*/
} // end of function to evaluate
