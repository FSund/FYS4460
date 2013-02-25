#ifndef CBOX_H
#define CBOX_H
#include "CAtom.h"
#include "linkedList.h"
#include "inlines.h"
#include <vector>
#include <armadillo>

//#include "CState.h"
//class CState;

using namespace std;
using namespace arma;

class CBox
{
public:
    CBox();
    CBox(const int nBoxes,
         const vec3 position,
         const vec3 size,
         const ivec3 boxIndex,
         const ivec3 NBoxes);
    CBox(const int nBoxes,
         const vec3 position,
         const vec3 size,
         const ivec3 boxIndex,
         const ivec3 NBoxes,
         CAtom *&firstAtomptr);

    //void addAtomByRef(CAtom* ptratom);
    void addAtom(CAtom *&atomptr);
    void findNeighbours(const vec3 &systemSize);
    void purgeAtoms(linkedList<CAtom*> &boxlessAtoms);
    void calculateForces(const vector<CBox*> boxes);
    void calculateForcesAndStatistics(const vector<CBox*> boxes);
    void resetForcesBool();

    linkedList<CAtom*> readAtomList() const;
    void getFirstAtom(linkedList<CAtom*> &atomList);
    CAtom* readFirstAtomPtr() const;

    void printAtomAddresses();
    int length() const;

    vec3 forceFromBox(const vec3 r0) const;
    void forceFromBox(const vec3 r0, vec3 &force, double &potSum, double &pressureSum) const;
protected:
    int nBoxes;
    int nNeighbours;
    int boxNumber;
    vec3 position;
    vec3 size;
    ivec3 index;
    ivec3 NBoxes;
    ivec neighbourNumbers;
    imat neighbourIndices;
    //vector<ivec3> neighbourIndices;
    mat displacementVectors;
    bool empty;
    bool forcesAreCalculated;
    linkedList<CAtom*> atomPtrList;

    int nAtoms;
};

#endif // CBOX_H
