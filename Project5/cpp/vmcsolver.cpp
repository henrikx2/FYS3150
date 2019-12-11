//The code is a modification of the example by Morten Hjorth-Jensen found at "http://compphysics.github.io/ComputationalPhysics/doc/pub/vmc/html/vmc.html".

#include <iostream>
#include <fstream>
#include <iomanip>
#include "vmcsolver.h"
#include "lib.h"

//Initialize solver
VMCSolver::VMCSolver() :
    
    nDimensions(3),
    nParticles(2),
    idum(-1),

    stepLength(1),
    alpha(0.873424),
    beta(0.5),
    omega(1.0),
    nCycles(1000000),
    trail(1),
    energy(0.0),
    energySquared(0.0),
    energyVariance(0.0),
    print("print")
{
}

//Different functions here

void VMCSolver::runMonteCarloIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    double waveFunctionOld = 0;
    double waveFunctionNew = 0;
    double energySum = 0;
    double energySquaredSum = 0;
    double distanceSum = 0;
    double deltaE;

    AC = 0;
    distance = 0;
    Ek = 0;
    Ep_HO = 0;
    Ep_electron = 0;
    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {
        // Store the current value of the wave function
        waveFunctionOld = waveFunction(rOld);
        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }
            // Recalculate the value of the wave function
            waveFunctionNew = waveFunction(rNew);
            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
                }
                AC += 1;
            }
            else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }
            // update energies
            deltaE = localEnergy(rNew);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }

    double div = nCycles*nParticles;
    energy = energySum/div;
    energySquared = energySquaredSum/div;
    energyVariance = (energySquared-energy*energy);
    AC /= div;
    distance /= div;
    Ek /= div;
    Ep_HO /= div;
    Ep_electron /= div;

    // Print values
    if(print == "print"){
        cout << "MC-cycles: " << nCycles << endl;
        cout << "Alpha = " << alpha << endl;
        if(trail == 2){cout << "Beta = " << beta << endl;}
        cout << "Total Energy: " << energy << endl;
        cout << "Energy Variance: " << energyVariance << endl;
        cout << "Kinetic Energy: " << Ek << endl;
        cout << "Potential Energy: " << Ep_HO+Ep_electron << endl;
        cout << "Mean distance r12: " << distance << endl;
        cout << "Accepted Configs: " << AC << endl << endl;
    }
}

double VMCSolver::localEnergy(const mat &r)
{
    double kineticEnergy = 0;
    double potentialEnergy = 0;
    double electronelectron = 0;
    double r12 = 0;

    // Potential and kinetic energy
    if(trail == 1 || trail == 2){
        kineticEnergy += 3*alpha*omega;
        double rParticles;
        for(int i = 0; i < nParticles; i++) {
            rParticles = 0;
            for(int j = 0; j < nDimensions; j++){
                rParticles += r(i,j)*r(i,j);
            }
            potentialEnergy += 0.5*omega*omega*rParticles;
            kineticEnergy -= 0.5*alpha*alpha*omega*omega*rParticles;
        }
        // Contribution from electron-electron potential
        for(int k = 0; k < nDimensions; k++) {
            r12 += (r(0,k) - r(1,k)) * (r(0,k) - r(1,k));
        }
        electronelectron += 1 / sqrt(r12);
        
    }
    // If we have second trail wavefunction add contributions
    if(trail == 2){
        double denom = 1+beta*sqrt(r12);
        kineticEnergy += 1/(2*denom*denom)*(alpha*omega*sqrt(r12)-1/(2*denom*denom)-2/sqrt(r12)+2*beta/denom);
    }
    // Test known wavefunction
    if(trail == 0){
        kineticEnergy = 0;
        potentialEnergy = 0;
        double rho = r(0,0)*r(0,0)+r(0,1)*r(0,1)+r(0,2)*r(0,2);
        kineticEnergy += -1/sqrt(rho)-alpha/2*(alpha-2/sqrt(rho));
    }
    
    distance += r12;
    Ek += kineticEnergy;
    Ep_HO += potentialEnergy;
    Ep_electron += electronelectron;
    return kineticEnergy + potentialEnergy + electronelectron;
}

double VMCSolver::waveFunction(const mat &r)
{
    double argument = 0;
    double psi_1 = 0;
    if(trail == 1 || trail == 2){
        double rSingleParticle;
        for(int i = 0; i < nParticles; i++) {
            rSingleParticle = 0;
            for(int j = 0; j < nDimensions; j++) {
                rSingleParticle += r(i,j) * r(i,j);
            }
            argument += rSingleParticle;
        }
        psi_1 = exp(-0.5*alpha*omega*argument);
    }
    if(trail == 2){
        argument = 0;
        double r12 = 0;
        for(int j = 0; j < nDimensions; j++) {
            r12 += (r(0,j) - r(1,j)) * (r(0,j) - r(1,j));
        }
        argument += sqrt(r12)/(1+beta*sqrt(r12));
        psi_1 *= exp(0.5*argument);
    }
    if(trail == 0){
        double rParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rParticle += r(0,j) * r(0,j);
        }
        psi_1 = alpha*rParticle*exp(-alpha*rParticle);
    }
    return psi_1;
}