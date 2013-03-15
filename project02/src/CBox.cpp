#include "CBox.h"

CBox::CBox():empty(1), nAtoms(0)
{
    //atomPtrList = linkedList<CAtom*>(); // may be redundant
}

CBox::CBox(const int nBoxes,
           const vec3 position,
           const vec3 size,
           const ivec3 boxIndex,
           const ivec3 NBoxes) :
    nBoxes(nBoxes),
    boxNumber(calculate_box_number(boxIndex, NBoxes)),
    position(position),
    size(size),
    index(boxIndex),
    NBoxes(NBoxes),
    empty(1),
    nAtoms(0)
{
    //atomPtrList = linkedList<CAtom*>(); // may be redundant
    //cout << "Creating empty box" << endl;
}

CBox::CBox(const int nBoxes,
           const vec3 position,
           const vec3 size,
           const ivec3 boxIndex,
           const ivec3 NBoxes,
           CAtom* &firstAtomptr) :
    nBoxes(nBoxes),
    boxNumber(calculate_box_number(boxIndex, NBoxes)),
    position(position),
    size(size),
    index(boxIndex),
    NBoxes(NBoxes),
    empty(0),
    nAtoms(1)
{
    atomPtrList = linkedList<CAtom*>(firstAtomptr); // there might be a better way to do this!!
    //cout << "Creating box with first atomptr in it" << endl;
}

void CBox::addAtom(CAtom* &atomptr)
{
    if (empty)
    {
        atomPtrList.item = atomptr;
        empty = 0;
    }
    else
    {
        atomPtrList.insertFirstItem(atomptr);
    }
    nAtoms++;
}

void CBox::findNeighbours(const vec3 &systemSize)
{
    int i = 0;
    bool alreadyInList;
    ivec3 boxIndex, dIndex;
    vec3 displacement;

    for (int di = -1; di <= 1; di++)
    {
        for (int dj = -1; dj <= 1; dj++)
        {
            for (int dk = -1; dk <= 1; dk++)
            {
                dIndex << di << dj << dk;
                boxIndex = index + dIndex;

                // making sure the index is within the box-space we have,
                // ie. applying periodic boundary conditions for the boxes
                for (int dim = 0; dim < 3; dim++)
                {
                    //add(dim) = int(floor(double(boxIndex(dim))/double(NBoxes(dim))))*NBoxes(dim); // PBC
                    boxIndex(dim) -= int(floor(double(boxIndex(dim))/double(NBoxes(dim)))) *NBoxes(dim); // PBC
                }

                // if we have the same box as we're in, don't add it to the lists
                if (boxIndex(0) == index(0) &&
                        boxIndex(1) == index(1) &&
                        boxIndex(2) == index(2))
                    continue;

                // check to see if we already have this box/index in the list
                // (only applies for systems with less than 27 boxes, meaning
                // smaller than 6*6*6*(3 sigma))
                alreadyInList = 0;
                for (int j = 0; j < i; j++)
                {
                    if (conv_to<int>::from(sum(boxIndex == neighbourIndices.col(j))) == 3) // if the compiler gives you trouble here, use the line below instead
                    //if (sum(boxIndex == neighbourIndices.col(j)) == 3)
                    {
                        alreadyInList = 1;
                        continue; // break out of this for loop
                    }
                    else if (!alreadyInList) alreadyInList = 0;
                }

                if (alreadyInList) continue;

                //neighbourIndices.push_back(boxIndex);
                neighbourIndices = join_rows(neighbourIndices, boxIndex);

                displacement(0) = (-(boxIndex(0) < index(0)+di) + (boxIndex(0) > index(0)+di))*systemSize(0);
                displacement(1) = (-(boxIndex(1) < index(1)+dj) + (boxIndex(1) > index(1)+dj))*systemSize(1);
                displacement(2) = (-(boxIndex(2) < index(2)+dk) + (boxIndex(2) > index(2)+dk))*systemSize(2);

                displacementVectors = join_rows(displacementVectors, displacement);
                i++; // keeping track of the number of neighbours we have added
                // this is most important when we have small systems (with less
                // than 6x6x6 unit cells for the Argon system), when we have less
                // than 26 neighbours
            }
        }
    }
    nNeighbours = i;

    // filling the "neighbourNumbers" 1D (armadillo) vector as well
    // (haven't decided which of the lists to use yet)
    neighbourNumbers = zeros<ivec>(nNeighbours);
    for (int i = 0; i < nNeighbours; i++)
        neighbourNumbers(i) = calculate_box_number(neighbourIndices.col(i), NBoxes);
}

void CBox::purgeAtoms(linkedList<CAtom*> &boxlessAtoms)
{
    if (empty)
    {
        // cout << "! Atomlist/box number " << boxNumber << " empty, exiting purge !" << endl;
        return;
    }

    linkedList<CAtom*>* runner;
    linkedList<CAtom*> newList(atomPtrList);
    runner = &newList; // making runner operate on the elements of newList

    vec3 atomPosition, maxBoundaries;
    CAtom* atomptr;
    maxBoundaries = position + size;

    // running through the list
    while (runner->readNext() != 0)
    {
        atomptr = runner->next->item; // the item/atom of next
        atomPosition = atomptr->getNewPosition();
        // checking if the atom is outside the box
        if (atomPosition(0) <  position(0) ||
            atomPosition(1) <  position(1) ||
            atomPosition(2) <  position(2) ||
            atomPosition(0) >= maxBoundaries(0) ||
            atomPosition(1) >= maxBoundaries(1) ||
            atomPosition(2) >= maxBoundaries(2) )
        {
            // putting the atom in the list of "boxless" atoms
            boxlessAtoms.insertFirstItem(atomptr);

            // removing the atom from the new list
            runner->dropNextItem(); // dropping the next item/atom from the list
            nAtoms--; // reducing the counter of number of atoms in the box

            // don't want to go advance in the list since we dropped the next item
            // must check the _new_ next item first
            continue;
        }

        runner = runner->next;
    }

    // checking the first item of the list manually, because the while-loop
    // above won't check this (because we have a stupid list)
    atomptr = newList(); // the first item/atom
    atomPosition = atomptr->getNewPosition();
    // checking if the atom is outside the box
    if (atomPosition(0) <  position(0) ||
        atomPosition(1) <  position(1) ||
        atomPosition(2) <  position(2) ||
        atomPosition(0) >= maxBoundaries(0) ||
        atomPosition(1) >= maxBoundaries(1) ||
        atomPosition(2) >= maxBoundaries(2) )
    {
        // putting the atom in the list of "boxless" atoms
        boxlessAtoms.insertFirstItem(atomptr);

        // removing the atom from the new list
        if (nAtoms == 1)
        {
            cout << "! Dropping the last atom in a box !" << endl;
            empty = 1;
            newList.item = 0;
            newList.next = 0;
            nAtoms = 0;
        }
        else
        {
            newList.dropFirstItem();
            nAtoms--;
        }
    }

    // replacing the list of this box with the new one we have made
    atomPtrList = newList;
    if (!empty) nAtoms = atomPtrList.length();
}

void CBox::calculateForces(const vector<CBox*> &boxes)
{
    if (empty)
    {
        // cout << "! Atomlist/box number " << boxNumber << " empty, exiting calculateForce !" << endl;
        return;
    }

    vec3 force, newPosition, rvec, forceComponent;
    CBox* neighbourBoxPtr;
    const linkedList<CAtom*>* atomList;
    atomList = &atomPtrList;

    // finding the forces from all neighbouring boxes
    while (atomList != 0)
    {
        newPosition = atomList->item->getNewPosition();
        force.zeros();

        for (int i = 0; i < nNeighbours; i++)
        {
            neighbourBoxPtr = boxes[neighbourNumbers(i)];

            // if we have already calculated the forces for this box (and its
            // neighbours) we just skip the force-calculation
            if (neighbourBoxPtr->forcesAreCalculated) continue;

            // else we calculate the force from the neighbouring box on the
            // current atom, and add the force from this atom to the forces on
            // the atoms in the neighbouring box

            rvec(0) = newPosition(0) + displacementVectors(0,i);
            rvec(1) = newPosition(1) + displacementVectors(1,i);
            rvec(2) = newPosition(2) + displacementVectors(2,i);

            forceComponent = neighbourBoxPtr->forceFromBox(rvec, atomList->item->matrixAtom);
            force(0) += forceComponent(0);
            force(1) += forceComponent(1);
            force(2) += forceComponent(2);
        }
        // adding the force from the neighbouring boxes on this atom, to the atom
        atomList->item->addToNewForce(force);

        // advancing to the next atom
        atomList = atomList->readNext();
    }

    //////////////////////////////////////////////////////
    // finding the forces from atoms in this box itself //
    const linkedList<CAtom*>* atomList2;
    vec3 newPosition2, drvec;
    double dr2, dr6, LJ;

    atomList = &atomPtrList;
    atomList2 = &atomPtrList;
    // loop over all atoms in this box
    while (atomList != 0)
    {
        newPosition = atomList->item->getNewPosition();
        force.zeros();

        // loop over all atoms with "index" higher than the current atom, making
        // use of Newton's third law for the rest of the forces
        atomList2 = atomList->readNext();
        while (atomList2 != 0)
        {
            if (atomList->item->matrixAtom && atomList2->item->matrixAtom)
            {
                atomList2 = atomList2->readNext();
                continue;
            }
            newPosition2 = atomList2->item->getNewPosition();
            drvec(0) = newPosition(0) - newPosition2(0);
            drvec(1) = newPosition(1) - newPosition2(1);
            drvec(2) = newPosition(2) - newPosition2(2);

            dr2 = drvec(0)*drvec(0) + drvec(1)*drvec(1) + drvec(2)*drvec(2);
            dr6 = dr2*dr2*dr2;
            LJ = 24.0*(2.0 - dr6)/(dr6*dr6*dr2);

            forceComponent(0) = LJ*drvec(0);
            forceComponent(1) = LJ*drvec(1);
            forceComponent(2) = LJ*drvec(2);

            force(0) += forceComponent(0);
            force(1) += forceComponent(1);
            force(2) += forceComponent(2);

            ////
//            if (!isfinite(force(0))) cout << "! infinite force, atoms in box = " << force.t();
            ////

            // using Newton's third law for the opposite force
            atomList2->item->addToNewForce(-forceComponent);

            // advancing to the next atom in the list
            atomList2 = atomList2->readNext();
        }
        // adding the force from all atoms in the box on the atom, to the atom
        atomList->item->addToNewForce(force);

        // advancing to the next atom in the list
        atomList = atomList->readNext();
    }

    // flagging this box, so that we know we can skip the forces from this box
    // on all neighbouring boxes. This is possible by using Newton's third law
    forcesAreCalculated = 1;
}

void CBox::calculateForcesAndStatistics(const vector<CBox*> &boxes)
{
    if (empty)
    {
        // cout << "! Atomlist/box number " << boxNumber << " empty, exiting calculateForce !" << endl;
        return;
    }

    vec3 force, newPosition, rvec, forceComp;
    CBox* neighbourBoxPtr;
    const linkedList<CAtom*>* atomList;
    double potComp, potSum, pressureComp, pressureSum;
    atomList = &atomPtrList;
    bool matrixAtom;

    // finding the forces from all neighbouring boxes
    while (atomList != 0)
    {
        newPosition = atomList->item->getNewPosition();
//        matrixAtom = atomList->item->matrixAtom;

        force.zeros();
        potSum = 0.0;
        pressureSum = 0.0;

        for (int i = 0; i < nNeighbours; i++)
        {
            neighbourBoxPtr = boxes[neighbourNumbers(i)];

            // if we have already calculated the forces for this box (and its
            // neighbours) we just skip the force-calculation
            if (neighbourBoxPtr->forcesAreCalculated) continue;

            // else we calculate the force from the neighbouring box on the
            // current atom, and add the force from this atom to the forces on
            // the atoms in the neighbouring box

            rvec(0) = newPosition(0) + displacementVectors(0,i);
            rvec(1) = newPosition(1) + displacementVectors(1,i);
            rvec(2) = newPosition(2) + displacementVectors(2,i);

            neighbourBoxPtr->forceFromBox(rvec, forceComp, potComp, pressureComp);

            force(0) += forceComp(0);
            force(1) += forceComp(1);
            force(2) += forceComp(2);

            potSum += potComp;
            pressureSum += pressureComp;
        }
        // adding the force from the neighbouring boxes on this atom, to the atom
        atomList->item->addToNewForce(force);
        atomList->item->addToNewStatistics(potSum, pressureSum);

        // advancing to the next atom
        atomList = atomList->readNext();
    }

    //////////////////////////////////////////////////////
    // finding the forces from atoms in this box itself //
    const linkedList<CAtom*>* atomList2;
    vec3 newPosition2, drvec;
    double dr2, dr6, LJ;

    atomList = &atomPtrList;
    atomList2 = &atomPtrList;
    // loop over all atoms in this box
    while (atomList != 0)
    {
        newPosition = atomList->item->getNewPosition();
        force.zeros();
        potSum = 0.0;
        pressureSum = 0.0;

        // loop over all atoms with "index" higher than the current atom, making
        // use of Newton's third law for the rest of the forces
        atomList2 = atomList->readNext();
        while (atomList2 != 0)
        {
            newPosition2 = atomList2->item->getNewPosition();
            drvec(0) = newPosition(0) - newPosition2(0);
            drvec(1) = newPosition(1) - newPosition2(1);
            drvec(2) = newPosition(2) - newPosition2(2);

            dr2 = drvec(0)*drvec(0) + drvec(1)*drvec(1) + drvec(2)*drvec(2);
            dr6 = dr2*dr2*dr2;
            LJ = 24.0*(2.0 - dr6)/(dr6*dr6*dr2);

            //dr12_inv = 1.0/(dr6*dr6);
            //LJ = 24*(2.0 - dr6)*dr12_inv/dr2;
            //potEn = (1.0 - dr6)*dr12_inv;

            forceComp(0) = LJ*drvec(0);
            forceComp(1) = LJ*drvec(1);
            forceComp(2) = LJ*drvec(2);
            potComp = (1.0 - dr6)/(dr6*dr6); // adding factor 4 later
            pressureComp = forceComp(0)*drvec(0) + forceComp(1)*drvec(1) + forceComp(2)*drvec(2);

            force(0) += forceComp(0);
            force(1) += forceComp(1);
            force(2) += forceComp(2);
            potSum += potComp;
            pressureSum += pressureComp;

            ////
//            if (!isfinite(force(0))) cout << "! infinite force, atoms in box = " << force.t();
            ////

            // using Newton's third law for the opposite force
            atomList2->item->addToNewForce(-forceComp);
            atomList2->item->addToNewStatistics(potComp, pressureComp);

            // advancing to the next atom in the list
            atomList2 = atomList2->readNext();
        }
        // adding the force from all atoms in the box on the atom, to the atom
        atomList->item->addToNewForce(force);
        atomList->item->addToNewStatistics(potSum, pressureSum);

        // advancing to the next atom in the list
        atomList = atomList->readNext();
    }

    // flagging this box, so that we know we can skip the forces from this box
    // on all neighbouring boxes. This is possible by using Newton's third law
    forcesAreCalculated = 1;
}

void CBox::resetForcesBool()
{
    forcesAreCalculated = 0;
}

void CBox::flush()
{
    empty = 1;
    atomPtrList.next = 0;
    atomPtrList.item = 0;
    nAtoms = 0;
}

vec3 CBox::forceFromBox(const vec3 &r0) const
{
    if (empty)
    {
        //cout << "! Atomlist/box number " << boxNumber << " empty, exiting forceFromBox and returning zerovec3 !" << endl;
        return zeros<vec>(3,1);
    }

    const linkedList<CAtom*>* atomList;
    atomList = &atomPtrList;
    vec3 r, drvec, forceComponent, force;
    double dr2, dr6, LJ;

    force.zeros();
    while (atomList != 0)
    {
        r = atomList->item->getNewPosition();
        drvec = r0 - r;

        dr2 = drvec(0)*drvec(0) + drvec(1)*drvec(1) + drvec(2)*drvec(2);
        dr6 = dr2*dr2*dr2;
        LJ = 24.0*(2.0 - dr6)/(dr6*dr6*dr2);

        forceComponent(0) = LJ*drvec(0);
        forceComponent(1) = LJ*drvec(1);
        forceComponent(2) = LJ*drvec(2);

        force(0) += forceComponent(0);
        force(1) += forceComponent(1);
        force(2) += forceComponent(2);

        atomList->item->addToNewForce(-forceComponent);

        atomList = atomList->readNext();
    }

    return force;
}

vec3 CBox::forceFromBox(const vec3 &r0, const bool &matrixAtom) const
{
    if (empty)
    {
        //cout << "! Atomlist/box number " << boxNumber << " empty, exiting forceFromBox and returning zerovec3 !" << endl;
        return zeros<vec>(3,1);
    }

    const linkedList<CAtom*>* atomList;
    atomList = &atomPtrList;
    vec3 r, drvec, forceComponent, force;
    double dr2, dr6, LJ;

    force.zeros();
    while (atomList != 0)
    {
        if (matrixAtom && atomList->item->matrixAtom)
        {
            atomList = atomList->readNext();
            continue;
        }

        r = atomList->item->getNewPosition();
        drvec = r0 - r;

        dr2 = drvec(0)*drvec(0) + drvec(1)*drvec(1) + drvec(2)*drvec(2);
        dr6 = dr2*dr2*dr2;
        LJ = 24.0*(2.0 - dr6)/(dr6*dr6*dr2);

        forceComponent(0) = LJ*drvec(0);
        forceComponent(1) = LJ*drvec(1);
        forceComponent(2) = LJ*drvec(2);

        force(0) += forceComponent(0);
        force(1) += forceComponent(1);
        force(2) += forceComponent(2);

        atomList->item->addToNewForce(-forceComponent);

        atomList = atomList->readNext();
    }

    return force;
}

void CBox::forceFromBox(
        const vec3 &r0,
        vec3 &force,
        double &potSum,
        double &pressureSum) const
{
    if (empty)
    {
        //cout << "! Atomlist/box number " << boxNumber << " empty, exiting forceFromBox and returning zerovec3 !" << endl;
        force.zeros();
        potSum = 0.0;
        pressureSum = 0.0;
        return;
    }

    const linkedList<CAtom*>* atomList;
    atomList = &atomPtrList;
    vec3 r, drvec, forceComp;
    double dr2, dr6, dr12_inv, LJ, potComp, pressureComp;

    force.zeros();
    potSum = 0.0;
    pressureSum = 0.0;
    while (atomList != 0)
    {
        r = atomList->item->getNewPosition();

        drvec = r0 - r;
        dr2 = drvec(0)*drvec(0) + drvec(1)*drvec(1) + drvec(2)*drvec(2);

        dr6 = dr2*dr2*dr2;
        dr12_inv = 1.0/(dr6*dr6);
        LJ = 24.0*(2.0 - dr6)*dr12_inv/dr2;

        forceComp(0) = LJ*drvec(0);
        forceComp(1) = LJ*drvec(1);
        forceComp(2) = LJ*drvec(2);
        potComp = dr12_inv - 1.0/dr6; // adding factor 4 later
        pressureComp = forceComp(0)*drvec(0) + forceComp(1)*drvec(1) + forceComp(2)*drvec(2);

        force(0) += forceComp(0);
        force(1) += forceComp(1);
        force(2) += forceComp(2);
        potSum += potComp;
        pressureSum += pressureComp;

        atomList->item->addToNewForce(-forceComp);
        atomList->item->addToNewStatistics(potComp, pressureComp);

        atomList = atomList->readNext();
    }
}

//void CBox::forceFromBox(
//        const vec3 &r0,
//        vec3 &force,
//        double &potSum,
//        double &pressureSum,
//        const bool &matrixAtom) const
//{
//    if (empty)
//    {
//        //cout << "! Atomlist/box number " << boxNumber << " empty, exiting forceFromBox and returning zerovec3 !" << endl;
//        force.zeros();
//        potSum = 0.0;
//        pressureSum = 0.0;
//        return;
//    }

//    const linkedList<CAtom*>* atomList;
//    atomList = &atomPtrList;
//    vec3 r, drvec, forceComp;
//    double dr2, dr6, dr12_inv, LJ, potComp, pressureComp;

//    force.zeros();
//    potSum = 0.0;
//    pressureSum = 0.0;
//    while (atomList != 0)
//    {
//        if (matrixAtom && atomList->item->matrixAtom)
//        {
//            atomList = atomList->readNext();
//            continue;
//        }

//        r = atomList->item->getNewPosition();

//        drvec = r0 - r;
//        dr2 = drvec(0)*drvec(0) + drvec(1)*drvec(1) + drvec(2)*drvec(2);

//        dr6 = dr2*dr2*dr2;
//        dr12_inv = 1.0/(dr6*dr6);
//        LJ = 24.0*(2.0 - dr6)*dr12_inv/dr2;

//        forceComp(0) = LJ*drvec(0);
//        forceComp(1) = LJ*drvec(1);
//        forceComp(2) = LJ*drvec(2);
//        potComp = dr12_inv - 1.0/dr6; // adding factor 4 later
//        pressureComp = forceComp(0)*drvec(0) + forceComp(1)*drvec(1) + forceComp(2)*drvec(2);

//        force(0) += forceComp(0);
//        force(1) += forceComp(1);
//        force(2) += forceComp(2);
//        potSum += potComp;
//        pressureSum += pressureComp;

//        atomList->item->addToNewForce(-forceComp);
//        atomList->item->addToNewStatistics(potComp, pressureComp);

//        atomList = atomList->readNext();
//    }
//}

linkedList<CAtom*> CBox::readAtomList() const
{
    return atomPtrList;
}

void CBox::getFirstAtom(linkedList<CAtom*> &atomList)
{
    atomList = atomPtrList;
}

CAtom *CBox::readFirstAtomPtr() const
{
    return atomPtrList(); // reading item field ( = CAtom*) of atomPtrList;
}

void CBox::printAtomAddresses()
{
    linkedList<CAtom*> atom(atomPtrList);

    while (atom.readNext() != 0)
    {
        cout << "address of atom in box    = " << atom.item << endl;
        atom = *atom.readNext();
    }
}

int CBox::length() const
{
    return atomPtrList.length();
}
