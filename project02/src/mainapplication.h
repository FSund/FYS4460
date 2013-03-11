#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include <cmath>
#include "CAtom.h"
#include "CState.h"
#include "CStatisticsSampler.h"

using namespace std;
using namespace arma;

class MainApplication
{
public:
    MainApplication();
    void runApplication(int argc, char *argv[]);
    void porous_system_d();
    void porous_system();
    void porous_system2();
private:
    CState initialize(double T_, double L_, int N_, long *seed);

    void test();
    double U(vec3 r0, vec3 r1);
};

#endif // MAINAPPLICATION_H
