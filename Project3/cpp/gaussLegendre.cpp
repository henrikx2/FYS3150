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
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10
#define PI 3.14159265359
using namespace std;
using namespace std::chrono;

//     Here we define various functions called by the main program

double int_function(double, double, double, double, double, double);
void gauleg(double, double, double *, double *, int);

//   Main function begins here
int main()
{
     int n;
     double a, b;
     cout << "Read in the number of integration points" << endl;
     cin >> n;
     cout << "Read in integration limits" << endl;
     cin >> a >> b;
//   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method
     double *x = new double [n];
     double *w = new double [n];
     
//   Start timing

     high_resolution_clock::time_point t1 = high_resolution_clock::now();

//   set up the mesh points and weights
     gauleg(a, b,x,w, n);
//   Initializing the brute force gauleg sum
     double int_gauss = 0.;
//   Run six for-loops to evaluate every integral
     int h,i,j,k,l,m;
     for (h = 0;h<n;h++){
     for (i = 0;i<n;i++){
     for (j = 0;j<n;j++){
     for (k = 0;k<n;k++){
     for (l = 0;l<n;l++){
     for (m = 0;m<n;m++){
     int_gauss += w[h]*w[i]*w[j]*w[k]*w[l]*w[m]*int_function(x[h],x[i],x[j],x[k],x[l],x[m]);
     }}}}}}
     
//    final time
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      double time_used = duration_cast<nanoseconds>(t2-t1).count();
      
//    final output
      cout  << setiosflags(ios::showpoint | ios::uppercase);
      cout << "Gaussian-Legendre Quadrature = "<< setw(20) << setprecision(15)  << int_gauss << endl;
      cout << "Error = "<< setw(20) << setprecision(15)  << abs((5*PI*PI/(16.0*16.0))-int_gauss) << endl;
      cout << "Time = "<< setw(20) << setprecision(15)  << time_used/(1.0E+9) << " s" << endl;
      delete [] x;
      delete [] w;
      return 0;
}  // end of main program

//  this function defines the function to integrate
double int_function(double x1, double x2, double y1, double y2, double z1, double z2)
{
  double alf = 2.;
  // Express integrand in terms of six variables
  double expo = exp(-2*alf*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2)));
  double diff = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  // Checking if diff > 0 by approximation
  if(diff > ZERO){return expo/diff;}
  else return 0;
} // end of function to evaluate

void gauleg(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()