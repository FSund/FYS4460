#include "mainapplication.h"

MainApplication::MainApplication()
{
}

void MainApplication::runApplication(int argc, char *argv[])
{
    double dt, T, Tbath, tau, L;
    int nSteps, thermostat, N, calculateStatistics, saveStates;

    if (argc != 11)
    {
        cout << endl << "! Usage: 'dt nSteps T(MD) thermostat#(0,1,2) Tbath(MD) tau L(SI) N calculateStatistics saveStates'" << endl;
        if (argc == 1)
        {
            dt = 0.005;
            nSteps = 200;
            T = 1.7; // to get 0.851
            thermostat = 1;
            Tbath = 0.851;
            tau = 15;
            L = 5.720;
            N = 20;
            calculateStatistics = 1;
            saveStates = 0;

            cout << endl << "! Using default of: dt = " << dt << ", nSteps = " << nSteps << ", T = " << T
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

//    state.load("./output/states/state.1000.xyz");
    for (int i = 0; i < nSteps; i++)
    {
        if (i%50==0) cout << "n = " << i << " of " << nSteps << endl;

        filename.str(string());
        filename << "./output/states/state." << setfill('0') << setw(4) << i << ".xyz";

        if (saveStates) state.save(filename.str(), 0, 1, 1, 0);

        if (calculateStatistics) sampler.sample(state, 1, dt*i);
//        if (i>250) sampler.pairCorrelation();

        if (thermostat == 1 && calculateStatistics) state.berendsen(Tbath, sampler.T, tt);
        else if (thermostat == 2) state.andersen(Tbath, tt, &seed);

        state.move(dt, calculateStatistics);
    }
    sampler.sample(state, 1, nSteps*dt);
//    filename.str(string());
//    filename << "./output/states/state." << setfill('0') << setw(4) << nSteps << ".xyz";
//    state.save(filename.str(), 0, 0);

    //    sampler.pairCorrelation_manual("./output/pairCorrelation_final.dat", 1, 200);
}

void MainApplication::thermalized_system()
{
    cout << "MainApplication::porous_system" << endl;

    double T = 1.8;
    double dt = 0.01;
    int nSteps = 5000;
    double Tbath = 0.851;
    int calculateStatistics = 1;
    long seed = -1;
    double tt = 1.0/15;

    double L_ = 5.720;
    int N_ = 20;

    ostringstream filename;

    vec3 L;
    ivec3 N;
    L << L_ << L_ << L_;
    N << N_ << N_ << N_;
    double interactionlength = 3.0;

    CState state = initialize(T, L_, N_, &seed);

    CStatisticsSampler sampler(state);
    sampler.initialize_pairCorrelation(200, 1);
    sampler.sample(state, 1, 0, 1);
    cout << "Initial temp = " << sampler.T << endl;

    for (int i = 0; i < nSteps; i++)
    {
        if (i%25==0) cout << "n = " << i << " of " << nSteps << endl;

        filename.str(string());
        filename << "./output/states/state." << setfill('0') << setw(5) << i << ".xyz";

//        sampler.sample(state, 1, dt*i, 0);
//        cout << "T (before thermostat) = " << sampler.T << endl;
        if (i < 3000) state.berendsen(Tbath, sampler.T, tt);
        sampler.sample(state, 1, dt*i, 1);
        cout << "T (after thermostat)  = " << sampler.T << endl;

        if (i%100 == 0 || i == 0) state.save(filename.str(), 0, 0, 0, 0);
        state.move(dt, calculateStatistics);
    }
    // save the final state
    filename.str(string());
    filename << "./output/states/state." << setfill('0') << setw(5) << nSteps << ".xyz";
    state.save(filename.str(), 0, 0, 0, 0);
}

void MainApplication::porous_system3()
{
    cout << "MainApplication::porous_system3" << endl;

    double dt = 0.02;
    int nSteps = 5000;
    double Tbath = 1.5;
    bool calculateStatistics = 1;
    long seed = -1;
    double tt = 1.0/15;

    double L_ = 5.720/L0;
    int N_ = 20;

    ostringstream filename;

    vec3 L;
    ivec3 N;
    L << L_ << L_ << L_;
    N << N_ << N_ << N_;
    double interactionlength = 3.0;

    CState state("state.05000.xyz", N, L, interactionlength); // OK

    state.generate_cylindrical_pore(20.0);
//    state.generate_spherical_pores(&seed, 20, 20.0, 30.0);
//    state.FILIP_pores();
    state.remove_half_the_atoms();
    state.saveMatrix("./output/states/matrix.xyz");

    CStatisticsSampler sampler(state);
    sampler.initialize_pairCorrelation(200, 1);
    sampler.sample(state, 1, 0, 1);
    cout << "Initial temp = " << sampler.T << endl;

    for (int i = 0; i < nSteps; i++)
    {
        if (i%25==0) cout << "n = " << i << " of " << nSteps << endl;

        sampler.sample(state, 1, dt*i, 1);
        cout << "T (before thermostat) = " << sampler.T << endl;
//        if (i < 3000) state.berendsen(Tbath, sampler.T, tt);

//        if (i%100 == 0 || i == 0)
        {
            filename.str(string());
            filename << "./output/states/state." << setfill('0') << setw(5) << i << ".xyz";
            state.save(filename.str(), 0, 0, 0, 0);
        }

        state.move(dt, calculateStatistics);
    }
    // save the final state
    filename.str(string());
    filename << "./output/states/state." << setfill('0') << setw(5) << nSteps << ".xyz";
    state.save(filename.str(), 0, 0, 0, 0);
}

void MainApplication::test()
{
    cout << "MainApplication::porous_system3" << endl;

    double dt = 0.01;
    int nSteps = 1;
    double Tbath = 1.5;
    int calculateStatistics = 1;
    long seed = -1;
    double tt = 1.0/15;

    double L_ = 5.720/L0;
    int N_ = 20;

    ostringstream filename;

    vec3 L;
    ivec3 N;
    L << L_ << L_ << L_;
    N << N_ << N_ << N_;
    double interactionlength = 3.0;

    CState state("state.05000.xyz", N, L, interactionlength);

//    state.generate_cylindrical_pore(20.0);
    state.generate_spherical_pores(&seed, 20, 20.0, 30.0);
//    state.FILIP_pores();
    state.remove_half_the_atoms();
    state.saveMatrix("./output/states/matrix.xyz");

    CStatisticsSampler sampler(state);
    sampler.initialize_pairCorrelation(200, 1);
    sampler.sample(state, 1, 0, 1);
    cout << "Initial temp = " << sampler.T << endl;

    for (int i = 0; i < nSteps; i++)
    {
        if (i%25==0) cout << "n = " << i << " of " << nSteps << endl;

        filename.str(string());
        filename << "./output/states/state." << setfill('0') << setw(5) << i << ".xyz";

        sampler.sample(state, 1, dt*i, 0);
        cout << "T (before thermostat) = " << sampler.T << endl;
        if (i < 3000) state.berendsen(Tbath, sampler.T, tt);

//        sampler.sample(state, 1, dt*i, 1);
//        cout << "T (after thermostat)  = " << sampler.T << endl;

        if (i%100 == 0 || i == 0) state.save(filename.str(), 0, 0, 0, 0);
        state.move(dt, calculateStatistics);
    }
    // save the final state
    filename.str(string());
    filename << "./output/states/state." << setfill('0') << setw(5) << nSteps << ".xyz";
    state.save(filename.str(), 0, 0, 0, 0);
}

CState MainApplication::initialize(double T_, double L_, int N_, long* seed)
{
    ivec3 N;
    vec3 L;

    N << N_ << N_ << N_;
    L << L_ << L_ << L_;
    L /= L0; // converting to MD units
    double T = T_;
    double sigma_v = sqrt(T); // MD units
    double interactionLength = 3.0; // MD units

    cout << "T_initial (MD)     :    " << T << endl;
    cout << "sigma_v (MD)       :    " << sigma_v << endl;

    CState state("Ar", "fcc", N, L, interactionLength);
    state.randnVelocity(0.0, sigma_v, seed);
    //state.randuVelocity(0.0, sigma_v, seed);

    return state;
}
