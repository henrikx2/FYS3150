#include "cppFunctions.h"
#include "vmcsolver.h"
#include <chrono>

using namespace std;
using namespace arma;
using namespace chrono;

// Similar to linspace in python
vec int_vec(double min,double max,double nSteps){
    double h = static_cast<double>((max-min)/nSteps);
    vec Linspace = vec(nSteps,fill::zeros);
    Linspace[0] = min;
    for (int i = 1; i < nSteps; i++){
        Linspace[i] = Linspace(i-1) + h;
    }
    return Linspace;
} // end of int_vec

// Read a Coloumn Vector from a file
vector<int> readvalues(string file){
    string line;
    vector<int> N_values;
    ifstream values;
    values.open(file,ios::in);
    while (getline(values,line)){
        if (atoi(line.c_str()) != 0){
            N_values.push_back(atoi(line.c_str()));
        }
        else{}
    }
    values.close();
    return N_values;
} // end of readvalues

// Calculate energies and variances as function of alpha
void plot_minima(int trail){
    cout << "Calculating: plot_minima(trail = " << to_string(trail) <<")" << endl;
    VMCSolver *solver = new VMCSolver();
    solver->print = "no print";
    solver->trail = trail;
    if (trail == 1)
    {
        // Initial alphas
        double nSteps = 2000.0;
        double alf_max = 2.1;
        double alf_min = 0.1;

        vec alpha_initial = int_vec(alf_min,alf_max,nSteps);
        vec energy = vec(nSteps,fill::zeros);
        vec energyVariance = vec(nSteps,fill::zeros);

        // Start timing
        high_resolution_clock::time_point time1 = high_resolution_clock::now();

        for(int i = 0; i < alpha_initial.n_elem; i++){
            solver->alpha = alpha_initial[i];
            solver->runMonteCarloIntegration();
            energy[i] = solver->energy;
            energyVariance[i] = solver->energyVariance;
        }
        uvec find_index;
        int index;
        find_index = find(energyVariance == min(energyVariance));
        index = find_index[0];
        cout << "Lowest alpha: " << alpha_initial[index] << endl;

        // Stops the clock
        high_resolution_clock::time_point time2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double> >(time2-time1);
        double runtime = time_span.count();
        runtime = runtime/60; // Make runtime in minutes
        cout << "RUNTIME: " << to_string(runtime) << " min" << endl;

        // Write values to file
        string alphas_file = "../data/Results_5_plot_minima_trail_" + to_string(solver->trail) + ".txt" ;
        ofstream output_values;
        output_values.open(alphas_file,ios::out);
        output_values << alpha_initial << endl;
        output_values << energy << endl;
        output_values << energyVariance << endl;
        output_values.close();
        cout << "Data saved to: " << alphas_file << endl;
    }
    
    if (trail == 2)
    {
        // Initial betas
        double nSteps = 2000.0;
        double beta_max = 1.8;
        double beta_min = 0.1;
        solver->alpha = 0.99451448; //Using optimal alpha for 2nd wave function

        vec beta_initial = int_vec(beta_min,beta_max,nSteps);
        vec energy = vec(nSteps,fill::zeros);
        vec energyVariance = vec(nSteps,fill::zeros);

        // Start timing
        high_resolution_clock::time_point time1 = high_resolution_clock::now();

        for(int i = 0; i < beta_initial.n_elem; i++){
            solver->beta = beta_initial[i];
            solver->runMonteCarloIntegration();
            energy[i] = solver->energy;
            energyVariance[i] = solver->energyVariance;
        }
        uvec find_index;
        int index;
        find_index = find(energyVariance == min(energyVariance));
        index = find_index[0];
        cout << "Lowest beta: " << beta_initial[index] << endl;

        // Stops the clock
        high_resolution_clock::time_point time2 = high_resolution_clock::now();
        duration<double> time_span = duration_cast<duration<double> >(time2-time1);
        double runtime = time_span.count();
        runtime = runtime/60; // Make runtime in minutes
        cout << "RUNTIME: " << to_string(runtime) << " min" << endl;

        // Write values to file
        string betas_file = "../data/Results_5c_plot_minima_trail_" + to_string(solver->trail) + ".txt" ;
        ofstream output_values;
        output_values.open(betas_file,ios::out);
        output_values << beta_initial << endl;
        output_values << energy << endl;
        output_values << energyVariance << endl;
        output_values.close();
        cout << "Data saved to: " << betas_file << endl;
    }
} // end of plot_minima

// Minimise energy and variance, and find best optimized alpha (and beta)
void minimize_alpha_beta(int trail){
    cout << "Calculating: minimize_alpha_beta(trail = " << trail << ")" << endl;
    VMCSolver *solver = new VMCSolver();
    solver->print = "no print";
    solver->trail = 1;
    
    double nSteps = 100.0;
    double alf_max;
    double alf_min;

    double tol = 1E-4;
    double int_tol = 1E-5;
    double variance_min = 1E10;
    double variance_old;

    // Initial alphas
    vec alpha_initial = int_vec(alf_min,alf_max,nSteps);
    vec energy = vec(nSteps,fill::zeros);
    vec energyVariance = vec(nSteps,fill::zeros);

    // Start timing
    high_resolution_clock::time_point time1 = high_resolution_clock::now();

    cout << "Minimizing " << solver->trail << ". wave function with ALPHA" << endl << endl;
    while(variance_min >= tol){
        variance_min = 1E10;
        alf_max = 1.1;
        alf_min = 0.5;
        alpha_initial = int_vec(alf_min,alf_max,nSteps);
        uvec find_index;
        int index;
        int runN = 1;
        while(variance_min >= tol){
            runN += 1;
            for(int i = 0; i < alpha_initial.n_elem; i++){
                solver->alpha = alpha_initial[i];
                solver->runMonteCarloIntegration();
                energy[i] = solver->energy;
                energyVariance[i] = solver->energyVariance;
            }
            find_index = find(energyVariance == min(energyVariance));
            index = find_index[0];
            solver->alpha = alpha_initial[index];
            double interval = abs(alf_max-alf_min);
            if(interval <= int_tol){
                variance_min = min(energyVariance);
                variance_old = variance_min;
                cout << "Variance MINIMIZED with ALPHA!" << endl << endl;
                cout << "Alpha = " << setprecision(8) << solver->alpha << endl;
                if(solver->trail == 2){
                    cout << "Beta = " << setprecision(8) << solver->beta << endl;
                }
                cout << "Energy: " << energy[index] << endl;
                cout << "Energy Variance: " << setprecision(16) << energyVariance[index] << endl << endl;
                break;
            }
            else{
                // Set new alpha interval
                variance_min = min(energyVariance);
                alf_min = alpha_initial[index-1];
                alf_max = alpha_initial[index+1];
                alpha_initial = int_vec(alf_min,alf_max,nSteps);
            }
        }
        if(variance_min <= tol){
            cout << "Variance MINIMIZED to tolerance: " << tol << endl << endl;
            cout << "MC-cycles: " << solver->nCycles << endl;
            cout << "Alpha = " << setprecision(8) << solver->alpha << endl;
            cout << "Beta = " << setprecision(8) << solver->beta << endl;
            cout << "Energy: " << energy[index] << endl;
            cout << "Energy Variance: " << setprecision(16) << energyVariance[index] << endl << endl;
            break;
        }
        if(trail == 1){break;}
        if(trail == 2){
            solver->trail = trail;
            cout << "Minimizing 2. wave function with BETA" << endl << endl;
            minimize_beta(solver, variance_min, tol, int_tol);
            if(variance_min <= tol || abs(variance_min-variance_old) <= 1E-6){
                cout << "Change in variance small (" << abs(variance_min-variance_old) << "), assuming minima." << endl;
                break;
            }
            cout << "Minimizing " << solver->trail << ". wave function with ALPHA" << endl << endl;

        }
    }
    // Stops the clock
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(time2-time1);
    double runtime = time_span.count();
    runtime = runtime/60; // Make runtime in minutes
    cout << "RUNTIME: " << to_string(runtime) << " min" << endl;

} // end of minimize_alpha_beta

// Minimize beta
void minimize_beta(VMCSolver *solver,double& variance_min, double tol, double int_tol){

    // Initial betas
    double nSteps = 100.0;
    double beta_max = 0.9;
    double beta_min = 0.1;

    vec beta_initial = int_vec(beta_min,beta_max,nSteps);
    vec energy = vec(nSteps,fill::zeros);
    vec energyVariance = vec(nSteps,fill::zeros);

    variance_min = 1E10;
    uvec find_index;
    int index;
    int runN = 1;
    while(variance_min > tol){
        runN += 1;
        for(int i = 0; i < beta_initial.n_elem; i++){
            solver->beta = beta_initial[i];
            solver->runMonteCarloIntegration();
            energy[i] = solver->energy;
            energyVariance[i] = solver->energyVariance;
        }
        find_index = find(energyVariance == min(energyVariance));
        index = find_index[0];
        solver->beta = beta_initial[index];
        double interval = abs(beta_max-beta_min);
        if(interval <= int_tol){
            variance_min = min(energyVariance);
            cout << "Variance MINIMIZED with BETA!" << endl << endl;
            cout << "Alpha = " << setprecision(8) << solver->alpha << endl;
            cout << "Beta = " << setprecision(8) << solver->beta << endl;
            cout << "Energy: " << energy[index] << endl;
            cout << "Energy Variance: " << setprecision(16) << energyVariance[index] << endl << endl;
            break;
        }
        else{
            // Set new beta interval
            variance_min = min(energyVariance);
            beta_min = beta_initial[index-1];
            beta_max = beta_initial[index+1];
            beta_initial = int_vec(beta_min,beta_max,nSteps);
        }
        
    }
    if(variance_min <= tol){
        cout << "Variance MINIMIZED to tolerance: " << tol << endl << endl;
        cout << "MC-cycles: " << solver->nCycles << endl;
        cout << "Alpha = " << setprecision(8) << solver->alpha << endl;
        cout << "Beta = " << setprecision(8) << solver->beta << endl;
        cout << "Energy: " << energy[index] << endl;
        cout << "Energy Variance: " << setprecision(16) << energyVariance[index] << endl << endl;
    }
} // end of minimize_beta

// Calculate energy for optimized alpha as function of MC-cycles
void plot_stability(int trail){
    cout << "Calculating: plot_stability()" << endl;

    VMCSolver *solver = new VMCSolver();
    solver->print = "no print";
    solver->trail = trail;

    // Set optimal variational parameters
    if(trail == 1){
        solver->alpha = 0.87909613;
    }
    if(trail == 2){
        solver->alpha = 0.99451448;
        solver->beta = 0.28375046;
    }
    
    int MC_min = 10;
    int MC_max = 1000000;
    int steps = 10000;

    vec Nvalues = int_vec(MC_min,MC_max,steps);
    vec energy = vec(Nvalues.n_elem,fill::zeros);
    vec energyVariance = vec(Nvalues.n_elem,fill::zeros);

    // Start timing
    high_resolution_clock::time_point time1 = high_resolution_clock::now();
    
    for(int i = 0; i < steps; i++){
        solver->nCycles = Nvalues[i];
        solver->runMonteCarloIntegration();
        energy[i] = solver->energy;
        energyVariance[i] = solver->energyVariance;

        if(i % (steps/10) == 0 && i != 0){
            cout << "n = " << solver->nCycles << endl;
        }
    }

    // Stops the clock
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(time2-time1);
    double runtime = time_span.count();
    runtime = runtime/60; // Make runtime in minutes
    cout << "RUNTIME: " << to_string(runtime) << " min" << endl;

    // Write values to file
    string file = "../data/Results_5_plot_stability_trail_"+to_string(trail)+"_MC_"+to_string(MC_max)+".txt" ;
    ofstream output_values;
    output_values.open(file,ios::out);
    output_values << Nvalues << endl;
    output_values << energy << endl;
    output_values << energyVariance << endl;
    output_values.close();
} // end of plot_stability

// Calculate optimal values
void calculate_optimal(int trail,vec omegas){
    cout << "Calculating: calculate_optimal(trail = " << to_string(trail) << ")" << endl;
    VMCSolver *solver = new VMCSolver();
    solver->nCycles = 1000000;
    solver->print = "print";
    solver->trail = trail;

    // Set optimal variational parameters
    if(trail == 1){
        solver->alpha = 0.87909613;
    }
    if(trail == 2){
        solver->alpha = 0.99451448;
        solver->beta = 0.28375046;
    }
    // Start timing
    high_resolution_clock::time_point time1 = high_resolution_clock::now();

    for(int i = 0; i < omegas.n_elem; i++){
        solver->omega = omegas[i];
        cout << "OMEGA = " << to_string(omegas[i]) << " (Omega_r = "<< to_string(omegas[i]/2) << ")" << endl << endl;
        solver->runMonteCarloIntegration();
    }
    // Stops the clock
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(time2-time1);
    double runtime = time_span.count();
    cout << "RUNTIME: " << to_string(runtime) << " s" << endl;

} //end of calculate_optimal

// Plot values to test Virial Theorem
void plot_virial(int trail, vec omegas){
    cout << "Calculating: plot_virial(trail = " << trail << ")" << endl;
    VMCSolver *solver = new VMCSolver();
    solver->print = "no print";
    solver->trail = trail;

    // Set optimal variational parameters
    if(trail == 1){
        solver->alpha = 0.87909613;
    }
    if(trail == 2){
        solver->alpha = 0.99451448;
        solver->beta = 0.28375046;
    }

    vec Epot_HO = vec(omegas.n_elem,fill::zeros);
    vec Epot_electrons = vec(omegas.n_elem,fill::zeros);
    vec Ekin = vec(omegas.n_elem,fill::zeros);

    for(int i = 0; i < omegas.n_elem; i++){
        solver->omega = omegas[i];
        solver->runMonteCarloIntegration();
        Epot_HO[i] = solver->Ep_HO;
        Epot_electrons[i] = solver->Ep_electron;
        Ekin[i] = solver->Ek;
    }

    // Write values to file
    string file = "../data/Results_5_plot_virial_trail_"+to_string(trail)+"_MC_"+to_string(solver->nCycles)+".txt" ;
    ofstream output_values;
    output_values.open(file,ios::out);
    output_values << omegas << endl;
    output_values << Ekin << endl;
    output_values << Epot_HO << endl;
    output_values << Epot_electrons << endl;
    output_values.close();

} // end plot_virial

// Test a wave function with known energies
void test_known_wave_func(){    
    cout << "Calculating: test_known_wave_func()" << endl;
    VMCSolver *solver = new VMCSolver();
    solver->print = "print";
    solver->trail = 0;
    solver->nParticles = 1;
    solver->alpha = 1.0;
    solver->beta = 1.0;
    solver->omega = 1.0;
    solver->nCycles = 1000000;

    solver->runMonteCarloIntegration();
} // end test_known_wave_func
