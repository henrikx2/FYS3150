#include <iostream>
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <armadillo>
#include "vmcsolver.h"
#include <omp.h>

using namespace std;
using namespace arma;

#define   ZERO       1.0E-10

vec int_vec(double min, double max,double nSteps);
vector<int> readvalues(string file);

void best_steplength();
void plot_minima(int trail);
void minimize_alpha_beta(VMCSolver *solver,int trail, double omega);
void minimize_beta(VMCSolver *solver,double& variance_min, double tol, double int_tol);
void plot_stability(int trail);
void calculate_optimal(int trail,vec omegas);
void plot_virial(int trail,vec omegas);
void test_known_wave_func();
