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
void gauss_laguerre(double *, double *, int, double);
double gammln(double);

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
     double *x_r = new double [n];
     double *w_r = new double [n];

//   Start timing

     high_resolution_clock::time_point t1 = high_resolution_clock::now();

//   set up the mesh points and weights
     gauleg(0, PI, x_t, w_t, n/2);
     gauleg(0, 2*PI, x_p, w_p, n/2);

//   set up the mesh points and weights and the power of x^alf
     double alf = 2.0;
     double alpha = 2.0;
     gauss_laguerre(x_r, w_r, n, alf);

//   Initializing the brute force gauleg sum
     double int_gauss = 0.;
//   Run six for-loops to evaluate every integral
     int h,i,j,k,l,m;
     for (h = 0;h<n;h++){
     for (i = 0;i<n;i++){
     for (j = 0;j<n/2;j++){
     for (k = 0;k<n/2;k++){
     for (l = 0;l<n/2;l++){
     for (m = 0;m<n/2;m++){
     int_gauss += w_r[h]*w_r[i]*w_t[j]*w_t[k]*w_p[l]*w_p[m]*int_function(x_r[h],x_r[i],x_t[j],x_t[k],x_p[l],x_p[m]);
     }}}}}}
     
     int_gauss /= pow(2*alpha,4)*2*alpha;

//    final time
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      double time_used = duration_cast<nanoseconds>(t2-t1).count();

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
double int_function(double r1, double r2, double t1, double t2, double p1, double p2){
   double cosb = cos(t1)*cos(t2) + sin(t1)*sin(t2)*cos(p1-p2);
   double f = exp(-3*(r1+r2))*r1*r1*r2*r2*sin(t1)*sin(t2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
   if(abs(r1*r1+r2*r2-2*r1*r2*cosb) > ZERO)
   return f;
   else 
   return 0;
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

void gauss_laguerre(double *x, double *w, int n, double alf)
{
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}
// end function gaulag

double gammln( double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

// end function gammln
//#undef EPS
//#undef MAXIT
