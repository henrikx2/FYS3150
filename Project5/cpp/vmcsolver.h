//The code is a modification of the example by Morten Hjorth-Jensen found at "http://compphysics.github.io/ComputationalPhysics/doc/pub/vmc/html/vmc.html".

#ifndef VMCSOLVER_H
#define VMCSOLVER_H
#include <iomanip>
#include <armadillo>
#include <omp.h>
using namespace arma;
using namespace std;
class VMCSolver
{
public:
    VMCSolver();
    void runMonteCarloIntegration();
    double stepLength;
    double alpha;
    double beta;
    double omega;
    int trail;
    int nCycles;
    int nParticles;

    double Ek;
    double Ep_HO;
    double Ep_electron;
    double energy;
    double energySquared;
    double energyVariance;
    double AC;
    double distance;

    string print;

private:
    double waveFunction(const mat &r);
    double localEnergy(const mat &r);
    int nDimensions;
    long idum;
    mat rOld;
    mat rNew;
};
#endif // VMCSOLVER_H
