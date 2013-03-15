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
    ~CState();

    void makeAtoms(ivec3 N, vec3 L, mat r, double sites_per_cell, string atomType);
    void makeBoxes(const vec3 size, const double interactionLength);
    void fillBoxes();
    void randnVelocity(double mean, double sigma, long *seed);
    void randuVelocity(double mean, double vmax, long *seed);

    void generate_spherical_pores(long* seed, int N_pores, double min_pore_radius, double max_pore_radius);
    void generate_cylindrical_pore(const double radius);
    void FILIP_pores();
    void remove_half_the_atoms();
//    void decrease_density_by_factor(double &factor);

    void save(string filename, bool saveSpeed=0, bool saveForces=0, bool indexing=0, bool markMatrixAtoms=0);
    void saveMatrix(string filename);
    void saveVelocity(string filename, bool MDUnits);

    void move(const double dt, const bool statistics);
    void newForces();
    void newForcesAndStatistics();
    void gravity();
    void berendsen(const double Tbath, const double T, const double tt);
    void andersen(const double Tbath, const double tt, long *seed);

    void average();

    int getnAtoms() const;
    int getnMovingAtoms() const;
    int getnMatrixAtoms() const;
    int getnBoxes() const;
    vec3 getSize() const;
    const CAtom &getAtom(const int &i) const;
    const CAtom *getAtomPtr(const int &i) const;

    // unused
//    CBox getBox(int i) const;
//    vector<CAtom*> getAtoms() const;
//    vector<CBox*> getBoxes() const;

    friend class MainApplication;
protected:
    int load(string filename);
    vector<CAtom*> atoms;
    vector<CAtom*> movingAtoms;
//    vector<string> elements;
    vector<CBox*> boxes;
    imat boxIndexes;
    //string element;
//    string structure;
    int nAtoms;
    int nMovingAtoms;
    int nMatrixAtoms;
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
inline const CAtom &CState::getAtom(const int &i) const
{
    return *atoms[i];
}

inline const CAtom *CState::getAtomPtr(const int &i) const
{
    return atoms[i];
}

#endif // CSTATE_H
