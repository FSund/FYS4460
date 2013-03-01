#include "CStatisticsSampler.h"

CStatisticsSampler::CStatisticsSampler(const CState &state_):
    state(&state_),
    nAtoms(state->getnAtoms()),
    nBoxes(state->getnBoxes()),
    systemDim(state->getSize()),
    volume(prod(systemDim)),
    rho(nAtoms/volume)
{
    initPos = mat(3, nAtoms);
    for (int i = 0; i < nAtoms; i++)
    {
        initPos.col(i) = state->getAtom(i).getPosition();
    }

    temperatureFile = fopen("./output/temperature.dat","w");
    pressureFile = fopen("./output/pressure.dat","w");
    energyFile = fopen("./output/energy.dat","w");
    velocityFile = fopen("./output/velocity.dat","w");
    diffusionFile = fopen("./output/diffusion.dat","w");
    pairCorrelationFile = fopen("./output/pairCorrelation.dat", "w");
}

void CStatisticsSampler::sample(
        const CState &state_,
        const bool MDunits,
        const double t)
{
    state = &state_;
    K = kineticEnergy(MDunits);
    U = potentialEnergy(MDunits);
    T = temperature(MDunits);
    P = pressure(MDunits);
    rsquared = diffusion(MDunits);

    fprintf(energyFile, "%f %f %f \n", t, K, U);
    fprintf(temperatureFile, "%f %f \n", t, T);
    fprintf(pressureFile, "%f %f \n", t, P);
    fprintf(diffusionFile, "%f %f \n", t, rsquared);

    fflush(0); // flush files (write now instead of waiting for program to finish)
}

double CStatisticsSampler::kineticEnergy(bool MDunits)
{
    vec3 v;
    K = 0.0;
    for (int i = 0; i < nAtoms; i++)
    {
        v = state->getAtom(i).getVelocity();
        K += v(0)*v(0) + v(1)*v(1) + v(2)*v(2); // MD units
    }

    K *= 0.5; // still MD units

    if (MDunits) return K; // MD units
    else return K*E0; // SI units
}

double CStatisticsSampler::potentialEnergy(bool MDunits)
{
    U = 0.0;
    for (int i = 0; i < nAtoms; i++)
    {
        U += state->getAtom(i).getPotEn();
    }
    U *= 4; // MD units, factor 4 from Lennard-Jones potential

    if (MDunits) return U; // MD ubits
    else return U*E0; // SI units
}

double CStatisticsSampler::temperature(bool MDunits)
{
    T = 2.0*K/(3.0*nAtoms); // MD units

    if (MDunits) return T; // MD units
    else return T*T0; // SI units
}

double CStatisticsSampler::pressure(bool MDunits)
{
    double sum = 0.0;
    for (int i = 0; i < nAtoms; i++)
    {
        sum += state->getAtom(i).getPressure(); // MD units
    }

    P = rho*T + sum/(3.0*volume); // MD units  

    if (MDunits) return P; // MD units
    else return P*(F0/(L0*L0)); // SI units
}

double CStatisticsSampler::diffusion(bool MDunits)
{
    vec3 dr;
    rsquared = 0.0;
    for (int i = 0; i < nAtoms; i++)
    {
        dr = state->getAtom(i).getPosition() - initPos.col(i);
        dr += state->getAtom(i).getBoundaryCrossings()%systemDim;;
        rsquared += dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2);
    }
    rsquared /= nAtoms;

    if (MDunits) return rsquared;
    else return rsquared*L0*L0;
}

void CStatisticsSampler::initialize_pairCorrelation(int nBins_, bool MDunits)
{
    nBins = nBins_;
    L_bins = norm(systemDim/2.0, 2)/nBins;
    vec bins(nBins);

    for (int i = 0; i < nBins; i++)
        bins(i) = i*L_bins;

    if (!MDunits) bins *= L0;
    for (int i = 0; i < nBins; i++)
        fprintf(pairCorrelationFile, "%f ", bins(i));
    fprintf(pairCorrelationFile, " \n");
}

imat CStatisticsSampler::pairCorrelation()
{
    double r;
    vec3 rvec;
    ivec count = zeros<ivec>(nBins,1);
    vec3 systemDim2 = systemDim/2.0;

    L_bins = norm(systemDim2, 2)/nBins;
    for (int i = 0; i < nAtoms; i++)
    {
        for (int j = i + 1; j < nAtoms; j++)
        {
            rvec = state->getAtom(i).getPosition() - state->getAtom(j).getPosition();
            for (int k = 0; k < 3; k++)
                if (abs(rvec(k)) >= systemDim2(k))
                    rvec(k) -= int(rvec(k)/systemDim2(k))*systemDim(k);
            r = sqrt(rvec(0)*rvec(0) + rvec(1)*rvec(1) + rvec(2)*rvec(2));
            count(int(r/L_bins)) += 2;
        }
    }

    for (int i = 0; i < nBins; i++)
        fprintf(pairCorrelationFile, "%d ", count(i));
    fprintf(pairCorrelationFile, " \n ");

    return count;
}

mat CStatisticsSampler::pairCorrelation_manual(string filename, bool MDunits, int nBins_)
{
    double r, L_bins_;
    vec3 rvec;
    ivec count = zeros<ivec>(nBins_,1);
    vec3 systemDim2 = systemDim/2.0;

    L_bins_ = norm(systemDim2, 2)/nBins_;
    for (int i = 0; i < nAtoms; i++)
    {
        for (int j = i + 1; j < nAtoms; j++)
        {
            rvec = state->getAtom(i).getPosition() - state->getAtom(j).getPosition();
            for (int k = 0; k < 3; k++)
                if (abs(rvec(k)) >= systemDim2(k))
                    rvec(k) -= int(rvec(k)/systemDim2(k))*systemDim(k);
            r = sqrt(rvec(0)*rvec(0) + rvec(1)*rvec(1) + rvec(2)*rvec(2));
            count(int(r/L_bins_)) += 2;
        }
    }
    mat g(nBins_, 2);
    for (int i = 0; i < nBins_; i++)
        g(i, 0) = i*L_bins_;
    if (!MDunits)
        g *= L0;
    g.col(1) = conv_to<mat>::from(count);

    FILE* manualPairCorrelationFile = fopen(filename.c_str(), "w");

    for (int i = 0; i < nBins_; i++)
        fprintf(manualPairCorrelationFile, "%f ", g(i,0));
    for (int i = 0; i < nBins_; i++)
        fprintf(manualPairCorrelationFile, "%d ", count(i));
    fprintf(manualPairCorrelationFile, " \n");

    return g;
}

void CStatisticsSampler::print(bool MDunits)
{
    if (MDunits)
    {
        cout << "Kinetic energy   = " << K << endl;
        cout << "Potential energy = " << U << endl;
        cout << "Temperature      = " << T << endl;
        cout << "Pressure         = " << P << endl;
    }
    else
    {
        cout << "Kinetic energy   = " << K*E0 << endl;
        cout << "Potential energy = " << U*E0 << endl;
        cout << "Temperature      = " << T*T0 << endl;
        cout << "Pressure         = " << P*(F0/(L0*L0)) << endl;
    }
};

