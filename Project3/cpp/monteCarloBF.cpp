#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "functions.h"
using namespace std;
using namespace std::chrono;

//     Main function begins here     
int main()
{
     int n;
     double a, b;
     cout << "Read in the number of Monte-Carlo samples" << endl;
     cin >> n;
     cout << "Read in integration limits" << endl;
     cin >> a >> b;
     double x[6], y, fx; 
     double int_mc = 0.;  double variance = 0.; double sigma;
     double sum_sigma= 0. ; long idum=-1 ;  
     double length = b-a; // we fix the max size of the box
     double jacobidet = pow((length),6);
     double inverse_period = 1./RAND_MAX;
     double exact = 5*PI*PI/(16.0*16.0);

//   Start timing

     high_resolution_clock::time_point t1 = high_resolution_clock::now();

//   evaluate the integral with importance sampling    
     for (int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
       for (int j = 0; j< 6; j++) {
           //x[j]=-length+2*length*ran0(&idum);
           x[j] = double(rand())*inverse_period*(b-a)+a;
       }
       fx = brute_force_MC(x); 
       int_mc += fx;
       sum_sigma += fx*fx;
     }
     int_mc = int_mc/((double) n );
     sum_sigma = sum_sigma/((double) n );
     variance=sum_sigma-int_mc*int_mc;

     sigma = jacobidet*sqrt(variance/((double) n ));

//    final time
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      double time_used = duration_cast<nanoseconds>(t2-t1).count();

//   final output 
      cout << setiosflags(ios::showpoint | ios::uppercase);
      cout << "Monte Carlo Brute Force = " << setw(10) << setprecision(8) << jacobidet*int_mc << endl;
      cout << "Error = "<< setw(20) << setprecision(15)  << abs(exact-jacobidet*int_mc) << endl;
      cout << "Standard deviation = " << setw(10) << setprecision(8) << sigma << endl;
      cout << "Variance = " << setw(10) << setprecision(8) << sigma*sigma << endl;    
      cout << "Time = "<< setw(20) << setprecision(15)  << time_used/(1.0E+9) << " s" << endl;

     return 0;

//   write to files

}  // end of main program 
