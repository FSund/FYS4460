#ifndef CSTATISTICSSAMPLER_H
#define CSTATISTICSSAMPLER_H
#include <armadillo>
#include "CState.h"
#include "CAtom.h"
#include <cstdio>

using namespace std;
using namespace arma;

class CStatisticsSampler
{
friend class MainApplication;
private:
    const CState* state;
public:
    CStatisticsSampler(const CState &state);
    void sample(const CState &state, bool MDunits, const double t, bool save=1);

    double kineticEnergy(bool MDunits);
    double potentialEnergy(bool MDunits);
    double temperature(bool MDunits);
    double pressure(bool MDunits);
    double diffusion(bool MDunits);

    void initialize_pairCorrelation(int nBins_, bool MDunits);
    imat pairCorrelation();
    mat pairCorrelation_manual(string filename, bool MDunits, int nBins_);

    void cylinder_flow(bool &save);
    void reset_flow();

    void print(bool MDunits);

    FILE *velocityFile;
    FILE *temperatureFile;
    FILE *pressureFile;
    FILE *energyFile;
    FILE *diffusionFile;
    FILE *pairCorrelationFile;
    FILE *flowFile;
    bool outputFolderExists;

    int nAtoms;
    int nMatrixAtoms;
    int nMovingAtoms;
    int nBoxes;
    vec3 systemDim;
    double volume;
    double rho;
    mat initPos;
    double K;
    double U;
    double T;
    double P;
    double rsquared;

    int nBins;
    double L_bins;

    double tCurr;
    double tPrev;
    mat currPosMoving;
    mat prevPosMoving;

    int nFlowBins;
    double rMaxFlow;
    double drFlow;
//    vec rBins;
    vec vxAverage;

    vector<const CAtom*> movingAtoms;
};

#endif // CSTATISTICSSAMPLER_H
