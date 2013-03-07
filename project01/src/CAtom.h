#ifndef CATOM_H
#define CATOM_H
#include <armadillo>

//const double SIGMA = 3.405; // Angstrom
//const double MASS = 39.948; // amu
//const double EPSILON = 0.010318; // eV

const double L0 = 3.405;    // Angstrom
const double t0 = 2156.9;   // fs
const double F0 = .30303;   // eV
const double E0 = 0.01038;  // eV
const double T0 = 119.74;   // K
const double pi = atan(1)*4;

using namespace std;
using namespace arma;

class CAtom
{
// default member is private -> 'r' and 'v' will be private
//    vec3 r; // position
//    vec3 v; // velocity
public:
    CAtom();
    CAtom(const vec3 &position, const vec3 &velocity);

    vec3   getPosition() const;
    vec3   getVelocity() const;
    vec3      getForce() const;
    double    getPotEn() const;
    double getPressure() const;
    ivec3 getBoundaryCrossings() const;

//    vec3 getNewPosition() const;
    vec3 getNewVelocity() const;
    vec3    getNewForce() const;

    void setPosition(const vec3 &newPosition);
    void setVelocity(const vec3 &newVelocity);
//    void    setForce(const vec3 &newForce);

    void setNewPosition(const vec3 &newNewPosition);
    void setNewVelocity(const vec3 &newNewVelocity);
    void    setNewForce(const vec3 &newNewForce);
//    void    setNewPotEn(const double &newNewPotEn);

    void addToNewForce(const vec3 &addForce);
//    void addToNewPotEn(const double &addPot);
    void addToBoundaryCrossings(const ivec3 addBoundaryCrossings);
    void addToNewStatistics(const double addPot, const double addPressure);

    void getData(vec3 &Position, vec3 &Velocity, vec3 &Force);
    void setData(const vec3 &newPosition, const vec3 &newVelocity, const vec3 &newForce);
    void setNewData(const vec3 &newNewPosition, const vec3 &newNewVelocity, const vec3 &newNewForce);

    void forward();
    void resetStatistics();
    ~CAtom();

    const vec3 &getNewPosition();
protected:
    vec3 position;
    vec3 velocity;
    vec3 force;
    double potEn;
    double pressure;

//    vec3 oldposition;
//    vec3 oldvelocity;
//    vec3 oldforce;

    vec3 newposition;
    vec3 newvelocity;
    vec3 newforce;
    double newPotEn;
    double newPressure;

    ivec3 boundaryCrossings;
};

// inline stuff for better performance
inline const vec3 &CAtom::getNewPosition()
{
    return newposition;
}

inline void CAtom::addToNewForce(const vec3 &addForce)
{
    newforce(0) += addForce(0);
    newforce(1) += addForce(1);
    newforce(2) += addForce(2);
}

inline void CAtom::addToNewStatistics(const double addPot, const double addPressure)
{
    newPotEn    += addPot;
    newPressure += addPressure;
}

#endif // CATOM_H
