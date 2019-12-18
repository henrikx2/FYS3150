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
    int all = 0;
    cout << "0 - Run all calculations (approx. XXX min)" << endl;
    cout << "1 - plot_minima(1st wave function) (approx. 8 min)" << endl;
    cout << "2 - plot_minima(2nd wave function) (approx. 12 min)" << endl;
    cout << "3 - minimize_alpha_beta(1st wave function) (approx. 2 min)" << endl;
    cout << "4 - minimize_alpha_beta(2nd wave function) (approx 2 min)" << endl;
    cout << "5 - plot_stability(1st wave function) (approx. 2 min)" << endl;
    cout << "6 - plot_stability(2nd wave function) (approx. 3 min)" << endl;
    cout << "7 - calculate_optimal(1st wave function)" << endl;
    cout << "8 - calculate_optimal(2nd wave function)" << endl;
    cout << "9 - plot_virial(1st wave function) (approx. 18 min)" << endl;
    cout << "10 - plot_virial(2nd wave function) (approx. 20 min)" << endl;
    cout << "11 - test_known_wave_function()" << endl;
    cout << "12 - best_steplength() (approx. min)" << endl;
    cout << "Enter number of calculation to perform from list above (0-11):" << endl;
    cin >> input;
    
    if(input == 0){
        input += 1;
        all += 1;
    }
    if(input == 1){
        plot_minima(1);
        if(all == 1){input += 1;}
    }
    if(input == 2){
        plot_minima(2);
        if(all == 1){input += 3;}
    }
    if(input == 3){
        VMCSolver *solver = new VMCSolver();
        // Start timing
        high_resolution_clock::time_point time1 = high_resolution_clock::now();

        vec omegas = vec("0.01 0.02 0.5 1.0 2.0");
        for(int i = 0; i < omegas.n_elem; i++){
        minimize_alpha_beta(solver,1,omegas[i]);
        }

        // Stops the clock
        high_resolution_clock::time_point time2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double> >(time2-time1);
        double total_runtime = time_span.count();
        total_runtime = total_runtime/60; // Make runtime in minutes
        cout << "TOTAL RUNTIME: " << to_string(total_runtime) << " min" << endl << endl;
    }
    if(input == 4){
        VMCSolver *solver = new VMCSolver();
        // Start timing
        high_resolution_clock::time_point time1 = high_resolution_clock::now();

        vec omegas = vec("0.01 0.02 0.5 1.0 2.0");
        for(int i = 0; i < omegas.n_elem; i++){
        minimize_alpha_beta(solver,2,omegas[i]);
        }

        // Stops the clock
        high_resolution_clock::time_point time2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double> >(time2-time1);
        double total_runtime = time_span.count();
        total_runtime = total_runtime/60; // Make runtime in minutes
        cout << "TOTAL RUNTIME: " << to_string(total_runtime) << " min" << endl << endl;
    }
    if(input == 5){
        plot_stability(1);
        if(all == 1){input += 1;}    
    }
    if(input == 6){
        plot_stability(2);
        if(all == 1){input += 1;}
    }
    if(input == 7){
        vec omegas = vec("0.01 0.02 0.5 1.0 2.0");
        calculate_optimal(1,omegas);
        if(all == 1){input += 1;}
    }
    if(input == 8){
        vec omegas = vec("0.01 0.02 0.5 1.0 2.0");
        calculate_optimal(2,omegas);
        if(all == 1){input += 1;}
    }
    if(input == 9){
        //vec omegas_virial = int_vec(0.01,1.11,11);
        vec omegas_virial = vec("0.01 0.02 0.5 1.0 2.0");
        plot_virial(1,omegas_virial);
        if(all == 1){input += 1;}
    }
    if(input == 10){
        //vec omegas_virial = int_vec(0.01,1.11,11);
        vec omegas_virial = vec("0.01 0.02 0.5 1.0 2.0");
        plot_virial(2,omegas_virial);
        if(all == 1){input += 1;}
    }
    if(input == 11){
        test_known_wave_func();
        if(all == 1){input += 1;}
    }
    if(input == 12){
        best_steplength();
    }
    if(all == 1){
        cout << "All calculations run!" << endl;
    }
    return 0;

} // end of main
