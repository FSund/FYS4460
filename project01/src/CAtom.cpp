#include "CAtom.h"

CAtom::CAtom():
    position(zeros<vec>(3,1)),
    velocity(zeros<vec>(3,1)),
    force(zeros<vec>(3,1)),
    potEn(0.0),
    pressure(0.0),
    boundaryCrossings(zeros<ivec>(3,1))
{
    // cout << "CAtom: Using the 'default' constructor that initializes to zeros<vec>(3)" << endl;
     // default constructor
}

CAtom::CAtom(const vec3 &position, const vec3 &velocity):
    position(position),
    velocity(velocity),
    force(zeros<vec>(3)),
    potEn(0.0),
    pressure(0.0),
    boundaryCrossings(zeros<ivec>(3,1))
{
    // efficient constructor that assigns 'r' and 'v' meaningful values
    // ('position' and 'velocity'). This doesn't do implicit conversions
}

//CAtom::CAtom(vec3 position = zeros(3), vec3 velocity = zeros(3)):r(position),v(velocity)
//{
//    // smart and efficient constructor that assigns 'r' and 'v' meaningful values
//    // ('position' and 'velocity'). If one or both of the arguments are missing,
//    // it assigns zeros(3) to the data fields with missing arguments. This might
//    // be a bit inefficient if constructing a lot of "atoms" ???
//}

vec3 CAtom::getPosition() const
{
    return position;
    //return this->r;
    //return (*this).r;
}

vec3 CAtom::getVelocity() const
{
    return velocity;
    //return this->v;
    //return (*this).v;
}

vec3 CAtom::getForce() const
{
    return force;
}

double CAtom::getPotEn() const
{
    return potEn;
}

double CAtom::getPressure() const
{
    return pressure;
}

ivec3 CAtom::getBoundaryCrossings() const
{
    return boundaryCrossings;
}

/* Inlined in header file*/
//vec3 CAtom::getNewPosition() const
//{
//    return newposition;
//}

vec3 CAtom::getNewVelocity() const
{
    return newvelocity;
}

vec3 CAtom::getNewForce() const
{
    return newforce;
}

void CAtom::getData(vec3 &Position, vec3 &Velocity, vec3 &Force)
{
    Position = position;
    Velocity = velocity;
    Force = force;
}

void CAtom::setPosition(const vec3 &newPosition)
{
    position = newPosition;
}

void CAtom::setVelocity(const vec3 &newVelocity)
{
    velocity = newVelocity;
}

/* Unused */
//void CAtom::setForce(const vec3 &newForce)
//{
//    force = newForce;
//}

/* Inlined in header file*/
//void CAtom::addToNewForce(const vec3 &addForce)
//{
//    //newforce += addForce;
//    newforce(0) += addForce(0);
//    newforce(1) += addForce(1);
//    newforce(2) += addForce(2);
//}

/* Inlined in header file*/
//void CAtom::addToNewPotEn(const double &addPot)
//{
//    newPotEn += addPot;
//}

void CAtom::setNewPosition(const vec3 &newNewPosition)
{
    newposition = newNewPosition;
}

void CAtom::setNewVelocity(const vec3 &newNewVelocity)
{
    newvelocity = newNewVelocity;
}

void CAtom::setNewForce(const vec3 &newNewForce)
{
    newforce = newNewForce;
}

/* Unused*/
//void CAtom::setNewPotEn(const double &newNewPotEn)
//{
//    newPotEn = newNewPotEn;
//}

void CAtom::setData(const vec3 &newPosition, const vec3 &newVelocity, const vec3 &newForce)
{
    position = newPosition;
    velocity = newVelocity;
    force    = newForce;
}

void CAtom::setNewData(const vec3 &newNewPosition, const vec3 &newNewVelocity, const vec3 &newNewForce)
{
    newposition = newNewPosition;
    newvelocity = newNewVelocity;
    newforce    = newNewForce;
}

void CAtom::addToBoundaryCrossings(const ivec3 addBoundaryCrossings)
{
    boundaryCrossings += addBoundaryCrossings;
}

/* Inlined in header file*/
//void CAtom::addToNewStatistics(const double addPot, const double addPressure)
//{
//    newPotEn    += addPot;
//    newPressure += addPressure;
//}

void CAtom::forward()
{
    //cout << "CAtom::forward()" << endl;
//    oldposition = position;
//    oldvelocity = velocity;
//    oldforce    = force;
    position    = newposition;
    velocity    = newvelocity;
    force       = newforce;
    potEn       = newPotEn;
    pressure    = newPressure;
    newposition = zeros<vec>(3,1);
    newvelocity = zeros<vec>(3,1);
    newforce    = zeros<vec>(3,1);
    newPotEn    = 0.0;
    newPressure = 0.0;
}

void CAtom::resetStatistics()
{
    newPotEn = 0.0;
    newPressure = 0.0;
}

CAtom::~CAtom()
{
    // default destructor
}
