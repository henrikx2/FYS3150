//The code is a modification of the example by Morten Hjorth-Jensen found at "http://compphysics.github.io/ComputationalPhysics/doc/pub/vmc/html/vmc.html".

#include "vmcsolver.h"
#include "lib.h"
#include "cppFunctions.h"
#include <iostream>
#include <chrono>
using namespace std;
using namespace chrono;

int main(){

    int input;
    cout << "1 - plot_minima(1st wave function) (approx. 8 min)" << endl;
    cout << "2 - plot_minima(2nd wave function) (approx. 12 min)" << endl;
    cout << "3 - minimize_alpha_beta(1st wave function) (approx. 2 min)" << endl;
    cout << "4 - minimize_alpha_beta(2nd wave function) (approx 30 min)" << endl;
    cout << "5 - plot_stability(1st wave function) (approx. 20 min)" << endl;
    cout << "6 - plot_stability(2nd wave function) (approx. 30 min)" << endl;
    cout << "7 - calculate_optimal(1st wave function)" << endl;
    cout << "8 - calculate_optimal(2nd wave function)" << endl;
    cout << "9 - plot_virial(1st wave function)" << endl;
    cout << "10 - plot_virial(2nd wave function)" << endl;
    cout << "11 - test_known_wave_function()" << endl;
    cout << "Enter number of calculation to perform from list above (1-11):" << endl;
    cin >> input;

    if(input == 1){plot_minima(1);}
    if(input == 2){plot_minima(2);}
    if(input == 3){minimize_alpha_beta(1);}
    if(input == 4){minimize_alpha_beta(2);}
    if(input == 5){plot_stability(1);}
    if(input == 6){plot_stability(2);}
    if(input == 7){
        vec omegas = vec("0.01 0.5 1.0");
        calculate_optimal(1,omegas);
    }
    if(input == 8){
        vec omegas = vec("0.01 0.5 1.0");
        calculate_optimal(2,omegas);
    }
    if(input == 9){
        vec omegas_virial = int_vec(0.01,1,100);
        plot_virial(1,omegas_virial);
    }
    if(input == 10){
        vec omegas_virial = int_vec(0.01,1,100);
        plot_virial(2,omegas_virial);
    }
    if(input == 11){
        test_known_wave_func();
    }
    
    // Run single calculation with timing
    /*
    // Start timing
    high_resolution_clock::time_point time1 = high_resolution_clock::now();
    
    VMCSolver *solver = new  VMCSolver();
    solver->trail = 2;
    solver->alpha = 0.99474245;
    solver->beta = 0.28324775;
    solver->nCycles = 100000000;
    solver->runMonteCarloIntegration();

    // Stops the clock
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(time2-time1);
    double runtime = time_span.count();
    cout << "RUNTIME: " << to_string(runtime) << " s" << endl;
    */
    return 0;

} // end of main
