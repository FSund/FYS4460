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

void MainApplication::porous_system_d()
{
//    cout << "MainApplication::porous_system" << endl;

////    double T_, Tbath, tau;
////    int thermostat, saveStates;
//    double L_, dt;
//    int N_, nSteps, calculateStatistics;

//    dt = 0.005;
//    nSteps = 200;
////    thermostat = 1;
////    Tbath = 0.851;
////    tau = 15;
//    calculateStatistics = 1;
////    saveStates = 0;
////    long seed = -1;

//    L_ = 5.720/L0;
//    N_ = 20;

//    ////
//    // cout << L_*N_ << endl;
//    ////

//    vec3 L;
//    ivec3 N;
//    L << L_ << L_ << L_;
//    N << N_ << N_ << N_;

//    double interactionlength = 3.0;

//    CState state;
//    state.load("N20_at_T0.845.xyz");

////    CStatisticsSampler sampler(state);
////    sampler.initialize_pairCorrelation(200, 1);
//    ostringstream filename;

//    double r = 20.0/L0; // Angstrom, radius of sphere in center
//    cout << "r = " << r << endl;

//    vec3 center = (N%L)/2.0;
//    vec3 rvec;
//    double dr;

//    ////
//    cout << center << endl;
//    ////

////    for (int i = 0; i < state.getnAtoms(); i++)
////    {
////        rvec = state.getAtom(i).getPosition() - center;
////        // cout << rvec.t();
////        dr = norm(rvec, 2);
////        if (dr >= r)
////        {
////            // cout << "i = " << i << ", size = " << state.atoms.size() << endl;
////            state.atoms[i]->setMatrixAtom();
////        }
////    }

//    state = CState(state, N, L, interactionlength);

//    for (int i = 0; i < nSteps; i++)
//    {
//        filename.str(string());
//        filename << "./output/states/state." << setfill('0') << setw(4) << i << ".xyz";

//        state.save(filename.str(), 1, 1, 1, 1);
//        state.move(dt, calculateStatistics);
//    }
}

void MainApplication::porous_system()
{
//    cout << "MainApplication::porous_system" << endl;

////    double dt, T_, Tbath, tau, L_;
////    int nSteps, thermostat, N_, calculateStatistics, saveStates;

//    double dt = 0.005;
//    int nSteps = 1000;
////    int thermostat = 1;
//    double Tbath = 1.05;
////    double tau = 15;
//    int calculateStatistics = 1;
////    int saveStates = 0;
//    long seed = -1;
//    double tt = 1.0/10;

////    double T_ = 1.0;
//    double L_ = 5.720/L0;
//    int N_ = 20;

//    vec3 L;
//    ivec3 N;
//    L << L_ << L_ << L_;
//    N << N_ << N_ << N_;

//    double interactionlength = 3.0;

//    CState state;
//    state.load("N20_at_T0.845.xyz");
//    state = CState(state, N, L, interactionlength);

//    ostringstream filename;

//    generate_pores(state, &seed);

//    ////
////    for (int i = 0; i < state.nAtoms; i++)
////        cout << state.elements[i];
//    ////

//    CStatisticsSampler sampler(state);
//    sampler.initialize_pairCorrelation(200, 1);

//    for (int i = 0; i < nSteps; i++)
//    {
//        if (i%25==0) cout << "n = " << i << " of " << nSteps << endl;

//        filename.str(string());
//        filename << "./output/states/state." << setfill('0') << setw(4) << i << ".xyz";

//        sampler.sample(state, 1, dt*i);
//        if (i < 300) state.berendsen(Tbath, sampler.T, tt);

//        state.save(filename.str(), 0, 0, 0, 0);
//        state.move(dt, calculateStatistics);
//    }
}

void MainApplication::porous_system2()
{
    cout << "MainApplication::porous_system2" << endl;

//    double dt, T_, Tbath, tau, L_;
//    int nSteps, thermostat, N_, calculateStatistics, saveStates;

    double dt = 0.005;
    int nSteps = 1000;
//    int thermostat = 1;
    double Tbath = 1.5;
//    double tau = 15;
    int calculateStatistics = 1;
//    int saveStates = 0;
    long seed = -1;
    double tt = 1.0/10;

//    double T_ = 1.0;
    double L_ = 5.720/L0;
    int N_ = 20;

    ostringstream filename;

    vec3 L;
    ivec3 N;
    L << L_ << L_ << L_;
    N << N_ << N_ << N_;
    double interactionlength = 3.0;

    CState state("N20_at_T0.845.xyz", N, L, interactionlength);

    state.generate_spherical_pores(&seed, 20, 20.0, 30.0);
    state.remove_half_the_atoms();
    state.saveMatrix("./output/states/matrix.xyz");

    CStatisticsSampler sampler(state);
//    sampler.initialize_pairCorrelation(200, 1);

    for (int i = 0; i < nSteps; i++)
    {
        if (i%25==0) cout << "n = " << i << " of " << nSteps << endl;

        filename.str(string());
        filename << "./output/states/state." << setfill('0') << setw(4) << i << ".xyz";

        sampler.sample(state, 1, dt*i);
        cout << "T = " << sampler.T << ", K = " << sampler.K << endl;
        if (i < 300) state.berendsen(Tbath, sampler.T, tt);

        state.save(filename.str(), 0, 0, 0, 0);
        state.move(dt, calculateStatistics);
    }
}

void MainApplication::porous_system3()
{
    cout << "MainApplication::porous_system3" << endl;

    double dt = 0.005;
    int nSteps = 1000;
    double Tbath = 1.5;
    int calculateStatistics = 1;
    long seed = -1;
    double tt = 1.0/10;

    double L_ = 5.720/L0;
    int N_ = 20;

    ostringstream filename;

    vec3 L;
    ivec3 N;
    L << L_ << L_ << L_;
    N << N_ << N_ << N_;
    double interactionlength = 3.0;

    CState state("N20_at_T0.845.xyz", N, L, interactionlength);

    state.generate_cylindrical_pore(20.0);
//    state.generate_spherical_pores(&seed, 50, 20, 30);
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
        filename << "./output/states/state." << setfill('0') << setw(4) << i << ".xyz";

        sampler.sample(state, 1, dt*i, 1);
        cout << "T (before thermostat) = " << sampler.T << endl;
        if (i < 100) state.berendsen(Tbath, sampler.T, tt);
//        sampler.sample(state, 1, dt*i, 0);
//        cout << "T (after thermostat)  = " << sampler.T << endl;

        state.save(filename.str(), 0, 0, 0, 0);
        state.move(dt, calculateStatistics);

//        vec3 vel;
//        for (int i = 0; i < state.getnAtoms(); i++)
//        {
//            if (state.getAtom(i).matrixAtom)
//            {
//                vel = state.getAtom(i).getVelocity();
//                if (norm(vel, 2) != 0)
//                    cout << vel.t();
//            }
//        }
    }
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
