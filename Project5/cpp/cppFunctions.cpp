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

// Find best steplength
void best_steplength(){
    cout << "Calculating: best_steplength()" << endl;
    VMCSolver *solver = new VMCSolver;
    solver->nCycles = 100000;
    solver->omega = 1.0;
    solver->trail = 1;
    solver->print = "no print";
    solver->optimize_steplength = 0;

    int size = 50;
    vec omegas = int_vec(0.01,1.1,size);
    vec alphas = int_vec(0.1,1.6,size);
    mat steps_best = mat(size,size,fill::zeros);
    vec steps_best_2D = vec(size,fill::zeros);

    double diff_new;
    double diff_old;

    // Start timing
        high_resolution_clock::time_point time1 = high_resolution_clock::now();
    
    cout << "3D:" << endl;
    for(int o = 0; o < omegas.n_elem; o++){
        cout << o+1 << "/" << size << endl;
        
        solver->omega = omegas[o];
        
        for(int i = 0; i < alphas.n_elem; i++){
            
            solver->alpha = alphas[i];
            
            diff_old = 1;
            for(double step = 0.1; step < 100; step += 0.4){
                
                solver->stepLength = step;
                solver->runMonteCarloIntegration();
                        
                diff_new = abs(0.5 - solver->AC);
                if(diff_new < diff_old){
                    diff_old = diff_new;
                    steps_best(o,i) = step;
                }
                if(abs(diff_new-diff_old) > 0.1){break;}
            }
        }
    }
    cout << "2D:" << endl;
    for(int i = 0; i < alphas.n_elem; i++){
        cout << i+1 << "/" << size << endl;
        solver->alpha = alphas[i];
        
        diff_old = 1;
        
        for(double step = 0.1; step < 10; step += 0.1){
            solver->stepLength = step;
            solver->runMonteCarloIntegration();
                        
            diff_new = abs(0.5 - solver->AC);
            if(diff_new < diff_old){
                diff_old = diff_new;
                steps_best_2D[i] = step;
            }
        }
    }

    // Stops the clock
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(time2-time1);
    double runtime = time_span.count();
    runtime = runtime/60; // Make runtime in minutes
    cout << "RUNTIME: " << to_string(runtime) << " min" << endl;

    // Write 2D-values to file
    string file2D = "../data/Results_5_best_steplength_MC_"+to_string(solver->nCycles)+"_2D.txt" ;
    ofstream output_values_2D;
    output_values_2D.open(file2D,ios::out);
    output_values_2D << alphas << endl;
    output_values_2D << steps_best_2D << endl;
    output_values_2D.close();
    cout << "Data saved to: " << file2D << endl;
    /*
    // Write 3D-values to file
    string file3D = "../data/Results_5_best_steplength_MC_"+to_string(solver->nCycles)+"_3D.txt" ;
    ofstream output_values;
    output_values.open(file3D,ios::out);
    output_values << omegas << endl;
    output_values << alphas << endl;
    output_values << steps_best << endl;
    output_values.close();
    cout << "Data saved to: " << file3D << endl;*/
} // end of best_step

// Calculate energies and variances as function of alpha
void plot_minima(int trail){
    cout << "Calculating: plot_minima(trail = " << to_string(trail) <<")" << endl;
    VMCSolver *solver = new VMCSolver();
    solver->print = "no print";
    solver->trail = trail;
    solver->nCycles = 100000;
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
    cout << "---------------------------------------------------------------------" << endl;
    cout << "Calculating: minimize_alpha_beta(trail = " << trail << ")" << endl << endl;
    VMCSolver *solver = new VMCSolver();
    solver->print = " no print";
    solver->trail = 1;
    solver->omega = 1.0;
    solver->nCycles = 100000;

    double nSteps = 50.0;
    if(trail == 1){nSteps = 200;}
    double alf_max;
    double alf_min;

    double tol = 1E-3;
    double int_tol = 1E-5;
    double variance_min = 1E10;
    double variance_old;

    // Initial alphas
    vec alpha_initial;
    vec energy_vec = vec(nSteps,fill::zeros);
    vec energyVariance_vec = vec(nSteps,fill::zeros);

    // Start timing
    high_resolution_clock::time_point time1 = high_resolution_clock::now();

    //cout << "Minimizing " << solver->trail << ". wave function with ALPHA" << endl << endl;
    cout << "Minimizing..." << endl;
    while(variance_min >= tol){
        variance_min = 1E10;
        alf_max = 1.1;
        alf_min = 0.7;
        if(trail == 1){
            nSteps = 200;
            alf_min = 0.1;
        }
        alpha_initial = int_vec(alf_min,alf_max,nSteps);
        uvec find_index;
        int index;
        int runN = 1;
        while(variance_min >= tol){
            runN += 1;
            for(int i = 0; i < alpha_initial.n_elem; i++){
                solver->alpha = alpha_initial[i];
                solver->runMonteCarloIntegration();
                energy_vec[i] = solver->energy;
                energyVariance_vec[i] = solver->energyVariance;
            }
            find_index = find(energyVariance_vec == min(energyVariance_vec));
            index = find_index[0];
            solver->alpha = alpha_initial[index];
            double interval = abs(alf_max-alf_min);
            if(interval <= int_tol){
                variance_min = min(energyVariance_vec);
                variance_old = variance_min;
                cout << "." << endl;
                /*
                cout << "Variance MINIMIZED with ALPHA!" << endl << endl;
                cout << "Alpha = " << setprecision(8) << solver->alpha << endl;
                if(solver->trail == 2){
                    cout << "Beta = " << setprecision(8) << solver->beta << endl;
                }
                cout << "Energy: " << energy_vec[index] << endl;
                cout << "Energy Variance: " << setprecision(16) << energyVariance_vec[index] << endl << endl;
                */
                break;
            }
            else{
                // Set new alpha interval
                variance_min = min(energyVariance_vec);
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
            cout << "Energy: " << energy_vec[index] << endl;
            cout << "Energy Variance: " << setprecision(16) << energyVariance_vec[index] << endl << endl;
            break;
        }
        if(trail == 1){
            cout << "Energy MINIMIZED with variance: " << endl << endl;
            cout << "MC-cycles: " << solver->nCycles << endl;
            cout << "Alpha = " << setprecision(8) << solver->alpha << endl;
            cout << "Energy: " << energy_vec[index] << endl;
            cout << "Energy Variance: " << setprecision(16) << energyVariance_vec[index] << endl << endl;
            break;
        }
        if(trail == 2){
            solver->trail = trail;
            //cout << "Minimizing 2. wave function with BETA" << endl << endl;
            minimize_beta(solver, variance_min, tol, int_tol);
            if(variance_min <= tol){
                break;
            }
            if(abs(variance_min-variance_old) <= 1E-6){
                cout << "Change in variance small (" << abs(variance_min-variance_old) << "), assuming minima." << endl << endl;
                cout << "Variance MINIMIZED to tolerance: " << tol << endl << endl;
                cout << "MC-cycles: " << solver->nCycles << endl;
                cout << "Alpha = " << setprecision(8) << solver->alpha << endl;
                cout << "Beta = " << setprecision(8) << solver->beta << endl;
                cout << "Energy: " << energy_vec[index] << endl;
                cout << "Energy Variance: " << setprecision(16) << energyVariance_vec[index] << endl << endl;
                break;
                }
            //cout << "Minimizing " << solver->trail << ". wave function with ALPHA" << endl << endl;

        }
    }
    // Stops the clock
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(time2-time1);
    double runtime = time_span.count();
    runtime = runtime/60; // Make runtime in minutes
    cout << "RUNTIME: " << to_string(runtime) << " min" << endl << endl;

} // end of minimize_alpha_beta

// Minimize beta
void minimize_beta(VMCSolver *solver,double& variance_min, double tol, double int_tol){

    // Initial betas
    double nSteps = 50.0;
    double beta_max = 0.5;
    double beta_min = 0.01;

    vec beta_initial = int_vec(beta_min,beta_max,nSteps);
    vec energy_beta_vec = vec(nSteps,fill::zeros);
    vec energyVariance_beta_vec = vec(nSteps,fill::zeros);

    variance_min = 1E10;
    uvec find_index;
    int index;
    int runN = 1;
    while(variance_min > tol){
        runN += 1;
        for(int i = 0; i < beta_initial.n_elem; i++){
            solver->beta = beta_initial[i];
            solver->runMonteCarloIntegration();
            energy_beta_vec[i] = solver->energy;
            energyVariance_beta_vec[i] = solver->energyVariance;
        }
        find_index = find(energyVariance_beta_vec == min(energyVariance_beta_vec));
        index = find_index[0];
        solver->beta = beta_initial[index];
        double interval = abs(beta_max-beta_min);
        if(interval <= int_tol){
            variance_min = min(energyVariance_beta_vec);
            cout << "." << endl;
            /*
            cout << "Variance MINIMIZED with BETA!" << endl << endl;
            cout << "Alpha = " << setprecision(8) << solver->alpha << endl;
            cout << "Beta = " << setprecision(8) << solver->beta << endl;
            cout << "Energy: " << energy_beta_vec[index] << endl;
            cout << "Energy Variance: " << setprecision(16) << energyVariance_beta_vec[index] << endl << endl;
            */
            break;
        }
        else{
            // Set new beta interval
            variance_min = min(energyVariance_beta_vec);
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
        cout << "Energy: " << energy_beta_vec[index] << endl;
        cout << "Energy Variance: " << setprecision(16) << energyVariance_beta_vec[index] << endl << endl;
    }
} // end of minimize_beta

// Calculate energy for optimized alpha as function of MC-cycles
void plot_stability(int trail){
    cout << "Calculating: plot_stability(trail = " << trail << ")" << endl;

    VMCSolver *solver = new VMCSolver();
    solver->print = "no print";
    solver->trail = trail;

    // Set optimal variational parameters
    if(trail == 1){
        solver->alpha = 0.85538944;
    }
    if(trail == 2){
        solver->alpha = 0.98197771;
        solver->beta = 0.304;
    }

    int MC_max = 100000;
    int h = 10;
    int nSteps = MC_max/h;
    vec Nvalues = vec(nSteps,fill::zeros);
    Nvalues[0] = 0.001;
    for (int i = 1; i < nSteps; i++){
        Nvalues[i] = Nvalues(i-1) + h;
    }    

    vec energy = vec(Nvalues.n_elem,fill::zeros);
    vec energyVariance = vec(Nvalues.n_elem,fill::zeros);

    // Start timing
    high_resolution_clock::time_point time1 = high_resolution_clock::now();
    
    for(int i = 0; i < nSteps; i++){
        solver->nCycles = Nvalues[i];
        solver->runMonteCarloIntegration();
        energy[i] = solver->energy;
        energyVariance[i] = solver->energyVariance;

        if(i % (nSteps/10) == 0 && i != 0){
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
    cout << "Data saved to: " << file << endl;

} // end of plot_stability

// Calculate optimal values
void calculate_optimal(int trail,vec omegas){
    cout << "Calculating: calculate_optimal(trail = " << to_string(trail) << ")" << endl;
    VMCSolver *solver = new VMCSolver();
    solver->nCycles = 100000;
    solver->print = "print";
    solver->trail = trail;

    // Start timing
    high_resolution_clock::time_point time1 = high_resolution_clock::now();

    // Set optimal variational parameters
    if(trail == 1){
        solver->alpha = 0.85538944;
    }
    if(trail == 2){
        solver->alpha = 0.98197771;
        solver->beta = 0.304;
    }

    for(int i = 0; i < omegas.n_elem; i++){
        solver->omega = omegas[i];
        cout << "OPTIMAL VALUES, OMEGA = " << to_string(omegas[i]) << " (Omega_r = "<< to_string(omegas[i]/2) << ")" << endl << endl;
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
    solver->nCycles = 100000;

    vec Epot_HO = vec(omegas.n_elem,fill::zeros);
    vec Epot_electrons = vec(omegas.n_elem,fill::zeros);
    vec Ekin = vec(omegas.n_elem,fill::zeros);

    // Start timing
    high_resolution_clock::time_point time1 = high_resolution_clock::now();

     // Set optimal variational parameters
    if(trail == 1){
        solver->alpha = 0.85538944;
    }
    if(trail == 2){
        solver->alpha = 0.98197771;
        solver->beta = 0.304;
    }

    for(int i = 0; i < omegas.n_elem; i++){
        solver->omega = omegas[i];
        solver->runMonteCarloIntegration();
        Epot_HO[i] = solver->Ep_HO;
        Epot_electrons[i] = solver->Ep_electron;
        Ekin[i] = solver->Ek;
    }

    // Stops the clock
    high_resolution_clock::time_point time2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(time2-time1);
    double total_runtime = time_span.count();
    total_runtime = total_runtime; //make time in minutes
    cout << "TOTAL RUNTIME: " << to_string(total_runtime) << " s" << endl;

    // Write values to file
    string file = "../data/Results_5_plot_virial_trail_"+to_string(trail)+"_MC_"+to_string(solver->nCycles)+".txt" ;
    ofstream output_values;
    output_values.open(file,ios::out);
    output_values << omegas << endl;
    output_values << Ekin << endl;
    output_values << Epot_HO << endl;
    output_values << Epot_electrons << endl;
    output_values.close();
    cout << "Data saved to: " << file << endl;

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
    solver->nCycles = 100000;

    solver->runMonteCarloIntegration();
} // end test_known_wave_func
