#include "CAtom.h"

CAtom::CAtom():
    position(zeros<vec>(3,1)),
    velocity(zeros<vec>(3,1)),
    force(zeros<vec>(3,1)),
    potEn(0.0),
    pressure(0.0),
    boundaryCrossings(zeros<ivec>(3,1)),
    matrixAtom(0)
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
    boundaryCrossings(zeros<ivec>(3,1)),
    matrixAtom(0)
{
    // efficient constructor that assigns 'r' and 'v' meaningful values
    // ('position' and 'velocity'). This doesn't do implicit conversions
}

CAtom::CAtom(const vec3 &position, const vec3 &velocity, string atomType_):
    position(position),
    velocity(velocity),
    force(zeros<vec>(3)),
    potEn(0.0),
    pressure(0.0),
    boundaryCrossings(zeros<ivec>(3,1)),
    matrixAtom(0),
    atomType(atomType_)
{
    setAtomType(atomType_);
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

/* Inlined in header file !!!*/

//vec3 CAtom::getNewPosition() const
//{
//    return newposition;
//}

/* Unused */
//vec3 CAtom::getNewVelocity() const
//{
//    return newvelocity;
//}

vec3 CAtom::getNewForce() const
{
    return newforce;
}

void CAtom::getData(vec3 &position_, vec3 &velocity_, vec3 &force_)
{
    position_ = position;
    velocity_ = velocity;
    force_ = force;
}

void CAtom::setPosition(const vec3 &position_)
{
    if (matrixAtom)
        cout << "! ny setPosition for matrixAtom !" << endl;
    position = position_;
}

void CAtom::setVelocity(const vec3 &velocity_)
{
//    if (matrixAtom && norm(velocity_, 2) != 0)
//    {
//        cout << "! setVelocity for matrixAtom !" << endl;
//        cout << velocity_.t();
//    }
    velocity = velocity_;
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

void CAtom::setNewPosition(const vec3 &newPosition_)
{
//    if (matrixAtom)
//        cout << "! setNewPosition() for matrixAtom !" << endl;
    newposition = newPosition_;
}

void CAtom::setNewVelocity(const vec3 &newVelocity_)
{
//    if (matrixAtom)
//    {
//        cout << "! setNewVelocity() for matrixAtom !" << endl;
//        cout << newVelocity_.t();
//    }
//    if (!matrixAtom)
        newvelocity = newVelocity_;
}

void CAtom::setNewForce(const vec3 &newForce_)
{
    newforce = newForce_;
}

/* Unused*/
//void CAtom::setNewPotEn(const double &newNewPotEn)
//{
//    newPotEn = newNewPotEn;
//}

//void CAtom::setData(const vec3 &newPosition, const vec3 &newVelocity, const vec3 &newForce)
//{
//    position = newPosition;
//    velocity = newVelocity;
//    force    = newForce;
//}

//void CAtom::setNewData(const vec3 &newNewPosition, const vec3 &newNewVelocity, const vec3 &newNewForce)
//{
//    newposition = newNewPosition;
//    newvelocity = newNewVelocity;
//    newforce    = newNewForce;
//}

//void CAtom::setMatrixAtom()
//{
//    matrixAtom = 1;
//    velocity.zeros();
//    newvelocity.zeros();
//}

int CAtom::setAtomType(const string &atomType_)
{
    if (atomType_ == "Ar")
    {
        atomType = atomType_;
        matrixAtom = 0;

        return 0;
    }
    else if (atomType_ == "Ar_m")
    {
        atomType = atomType_;
        matrixAtom = 1;
        velocity.zeros();
        newvelocity.zeros();
        newposition = position;

        return 1;
    }
    else // default
    {
        cout << endl << "! Unknown atom type, exiting !" << endl << endl;
        exit(1);
    }
}

string CAtom::getAtomType() const
{
    return atomType;
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

void CAtom::checkForceCutoff()
{
    if (!matrixAtom)
    {
        for (int i = 0; i < 3; i++)
        {
            if (newforce(i) > forcemax)
            {
                cout << "! force cutoff MAX ! for " << atomType;
                cout << newforce(i) << " " << endl;
                newforce(i) = forcemax;
            }
            else
            if (newforce(i) < -forcemax)
            {
                cout << "! force cutoff MIN ! for " << atomType;
                cout << newforce(i) << " " << endl;
                newforce(i) = -forcemax;
            }
        }
    }
}

void CAtom::forward()
{
    if (!matrixAtom)
    {
        position    = newposition;
        velocity    = newvelocity;
        newposition.zeros();// = zeros<vec>(3,1);
        newvelocity.zeros();// = zeros<vec>(3,1);
    }
    force = newforce;
    newforce.zeros();//    = zeros<vec>(3,1);

    potEn       = newPotEn;
    pressure    = newPressure;
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
