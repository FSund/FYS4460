#include "mainapplication.h"

MainApplication::MainApplication()
{
}

void MainApplication::runApplication(int argc, char *argv[])
{
//    cout << "MainApplication::runApplication()" << endl;

    double dt, T, Tbath, tau, L;
    int nSteps, thermostat, N, calculateStatistics, saveStates;

    if (argc != 11)
    {
        cout << "Usage: dt nSteps T(MD) thermostat#(0,1,2) Tbath(MD) tau L(SI) N calculateStatistics saveStates" << endl;
        if (argc == 1)
        {
            dt = 0.005;
            nSteps = 500;
            T = 1.0;
            thermostat = 1;
            Tbath = 3.0;
            tau = 15;
            L = 5.260;
            N = 12;
            calculateStatistics = 1;
            saveStates = 1;

            cout << "Using default of: dt = " << dt << ", nSteps = " << nSteps << ", T = " << T
                 << ", thermostat# = " << thermostat << ", Tbath = " << Tbath
                 << ", tau = " << tau << ", L = " << L
                 << ", N = " << N << ", statistics? = " << calculateStatistics
                 << ", save states? = " << saveStates << endl << endl;
        } else
            exit(1);
    }
    else if (argc == 11)
    {
        dt = atof(argv[1]);
        nSteps = atoi(argv[2]);
        T = atof(argv[3]);
        thermostat = atoi(argv[4]);
        Tbath = atof(argv[5]);
        tau = atof(argv[6]);
        L = atof(argv[7]);
        N = atoi(argv[8]);
        calculateStatistics = atoi(argv[9]);
        saveStates = atoi(argv[10]);
    }

    double tt = 1.0/10;
    long seed = -1;
    cout << "dt                 :    " << dt << endl;
    cout << "Tbath              :    " << Tbath << endl;

    CState state;
    state = initialize(T, L, N, &seed);
    CStatisticsSampler sampler(state);
    sampler.initialize_pairCorrelation(200, 1);

    ostringstream filename;
//    ostringstream velocityFilename;
//    ofstream fout;
//    fout.open("./output/output5.5.txt");

//    double dt = 0.005; // in MD units
//    bool calculateStatistics = 1;
//    double Tbath = 1.0; // MD units
//    double tau = 10.0;

//    cout << "Initial conditions (r, v, a) (MD units):" << endl;
//    cout << state.getAtom(1).getPosition().t();
//    cout << state.getAtom(1).getVelocity().t();
//    cout << state.getAtom(1).getForce().t();
//    cout << endl;

//    state.average();

//    state.load("./output/states/final_T5.5_N1000_dt0.005.xyz");

    for (int i = 0; i < nSteps; i++)
    {
        cout << "n = " << i << endl;
        filename.str(string());
        filename << "./output/states/state." << setfill('0') << setw(4) << i << ".xyz";

        if (saveStates) state.save(filename.str(), 1, 1);

        if (calculateStatistics) sampler.sample(state, 1, dt*i);
        //if (i>250) sampler.pairCorrelation();
        cout << "P = " << sampler.P << ", T = " << sampler.T << endl;

        if (thermostat == 1) state.berendsen(Tbath, sampler.T, tt);
        else if (thermostat == 2) state.andersen(Tbath, sampler.T, tt, &seed);

        state.move(dt, calculateStatistics);

        //cout << "Intermediate conditions (r, v, a) (MD units):" << endl;
        //cout << state.getAtom(95).getPosition().t();
        //cout << state.getAtom(95).getVelocity().t();
        //cout << state.getAtom(95).getForce().t();
    }
//    state.save("./output/states/final_T5.5_N2000_dt0.005.xyz", 0, 0);

    sampler.sample(state, 1, nSteps*dt);
    sampler.pairCorrelation_manual("./output/pairCorrelation_final", 1, 200);

//    cout << endl << "Final conditions (r, v, a) (MD units): " << endl;
//    cout << state.getAtom(1).getPosition().t();
//    cout << state.getAtom(1).getVelocity().t();
//    cout << state.getAtom(1).getForce().t() << endl;
//    cout << "Boundary crossings: " << endl;
//    cout << state.getAtom(1).getBoundaryCrossings().t() << endl;
//    cout << "System dimensions (MD units): " << endl;
//    cout << state.size.t() << endl;

    //state.average();
}

CState MainApplication::initialize(double T_, double L_, int N_, long* seed)
{
//    cout << "MainApplication::Initialize()" << endl;

    ivec3 N;
    vec3 L;

    N << N_ << N_ << N_;
    //L << 5.260 << 5.260 << 5.260;
    L << L_ << L_ << L_;
    L /= L0; // converting to MD units
    double T = T_;

    //double T = 100;           // Kelvin
    //T = T/T0;                 // converting to MD units
    //double T = 5.5; // MD units;

    double sigma_v = sqrt(T); // MD units
    double interactionLength = 3.0; // MD units
    //long seed = -1; // CHANGE THIS

    cout << "T_initial (MD)     :    " << T << endl;
    cout << "sigma_v (MD)       :    " << sigma_v << endl;

    CState state("Ar", "fcc", N, L, interactionLength);
    state.randnVelocity(0.0, sigma_v, seed);
    //state.randuVelocity(0.0, sigma_v, &seed);

    return state;
}
