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
    cout << "1 - plot_minima(1st wave function) (approx. 1 min)" << endl;
    cout << "2 - plot_minima(2nd wave function) (approx. 2 min)" << endl;
    cout << "3 - minimize_alpha_beta(1st wave function) (approx. 1 min)" << endl;
    cout << "4 - minimize_alpha_beta(2nd wave function) (approx 1 min)" << endl;
    cout << "5 - plot_stability(1st wave function) (approx. 2 min)" << endl;
    cout << "6 - plot_stability(2nd wave function) (approx. 3 min)" << endl;
    cout << "7 - calculate_optimal(1st wave function)" << endl;
    cout << "8 - calculate_optimal(2nd wave function)" << endl;
    cout << "9 - plot_virial(1st wave function)" << endl;
    cout << "10 - plot_virial(2nd wave function)" << endl;
    cout << "11 - test_known_wave_function()" << endl;
    cout << "12 - best_steplength() (approx. 27 min)" << endl;
    cout << "Enter number of calculation to perform from list above (0-12):" << endl;
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
        if(all == 1){input += 1;}
    }
    if(input == 3){
        minimize_alpha_beta(1);
        if(all == 1){input += 1;}
    }
    if(input == 4){
        minimize_alpha_beta(2);
        if(all == 1){input += 1;}
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
        vec omegas_virial = int_vec(0.01,1.11,100);
        plot_virial(1,omegas_virial);
        if(all == 1){input += 1;}
    }
    if(input == 10){
        vec omegas_virial = int_vec(0.01,1.11,100);
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
