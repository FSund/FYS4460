#ifndef CSTATE_H
#define CSTATE_H

//class CAtom;
//class CBox;

#include <armadillo>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "CAtom.h"
#include "CBox.h"
#include "lib.h"
#include "inlines.h"
//#include <omp.h>

//using namespace std;
//using namespace arma;

class CState
{
friend class mainapplication;
public:
    CState();
    CState(string element, string structure, ivec3 N, vec3 L, double interactionLength);
//    CState(const CState &state, const ivec3 &N, const vec3 &L, const double &interactionLength);
    CState(const string filename, ivec3 &N, vec3 &L, double &interactionLength);

    void makeAtoms(ivec3 N, vec3 L, mat r, double sites_per_cell);
    void makeBoxes(const vec3 size, const double interactionLength);
    void fillBoxes();
    void randnVelocity(double mean, double sigma, long *seed);
    void randuVelocity(double mean, double vmax, long *seed);
    void generate_spherical_pores(long* seed, int N_pores, double min_pore_radius, double max_pore_radius);

    void remove_half_the_atoms();

    void save(string filename, bool saveSpeed, bool saveForces, bool indexing, bool markMatrixAtoms);
    void saveMatrix(string filename);
    void saveVelocity(string filename, bool MDUnits);
    void load(string filename);

    void move(const double dt, const bool statistics);
    void newForces();
    void newForcesAndStatistics();
    void berendsen(const double Tbath, const double T, const double tt);
    void andersen(const double Tbath, const double tt, long *seed);

    void average();

    int getnAtoms() const;
    int getnBoxes() const;
    vec3 getSize() const;
    const CAtom &getAtom(int i) const;

    // unused
//    CBox getBox(int i) const;
//    vector<CAtom*> getAtoms() const;
//    vector<CBox*> getBoxes() const;

    friend class MainApplication;
protected:
    vector<CAtom*> atoms;
//    vector<string> elements;
    vector<CBox*> boxes;
    imat boxIndexes;
    //string element;
    string structure;
    int nAtoms;
    int nBoxes;
    int atoms_per_cell;
    int nSites;
    int sites_per_cell;
    ivec3 N;
    ivec3 NBoxes;
    vec3 L;
    vec3 size;
    vec3 boxDimensions;
    double kinetic;
    double potential;

//    struct data
//    {
//        vector<CAtom*> atoms;
//    };
//    data matrix;

//    vector<CAtom*> moving_atoms;
};

// inline stuff for better performance
inline const CAtom &CState::getAtom(int i) const
{
    return *atoms[i];
}

#endif // CSTATE_H
