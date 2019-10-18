// importance sampling with gaussian deviates
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
     cout << "Read in the number of Monte-Carlo samples" << endl;
     cin >> n;
     double fx; 
     double int_mc = 0.;  double variance = 0.; double sigma = 0.;
     double sum_sigma= 0. ; long idum=-1 ;
     double alpha = 2.0;  
     double jacobidet = 4*pow(acos(-1.),4.)/pow(2.0*alpha,2);
     double inverse_period = 1./RAND_MAX;
     double *f = new double [n];

//   Start timing

     high_resolution_clock::time_point t1 = high_resolution_clock::now();

//   evaluate the integral with importance sampling    
     for (int i = 1; i <= n; i++){

     double r1 = -0.25*log((1-double(rand())*inverse_period));
     double r2 = -0.25*log((1-double(rand())*inverse_period));
     double t1 = double(rand())*inverse_period*PI;
     double t2 = double(rand())*inverse_period*PI;
	double p1 = double(rand())*inverse_period*2*PI;
	double p2 = double(rand())*inverse_period*2*PI;
       
       fx = exp_MC(r1,r2,t1,t2,p1,p2); 
       int_mc += fx;
       f[i] = fx;
     }
     int_mc = int_mc/((double) n );
     for (int i = 1; i <= n; i++){
          variance = (f[i]-int_mc)*(f[i]-int_mc);
     }
     variance *= (jacobidet/n);
     
     sigma = sqrt(variance/((double) n ));

//    final time
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      double time_used = duration_cast<nanoseconds>(t2-t1).count();

//   final output 
      cout << setiosflags(ios::showpoint | ios::uppercase);
      cout << "Monte carlo result= " << setw(10) << setprecision(8) << jacobidet*int_mc << endl;
      cout << "Error = "<< setw(20) << setprecision(15)  << abs((5*PI*PI/(16.0*16.0))-jacobidet*int_mc) << endl;
      cout << "Standard deviation = " << setw(10) << setprecision(8) << sigma << endl;
      cout << "Variance = " << setw(10) << setprecision(8) << variance << endl;   
      cout << "Time = "<< setw(20) << setprecision(15)  << time_used/(1.0E+9) << " s" << endl;
     return 0;
}  // end of main program 
