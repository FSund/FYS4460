#include "CState.h"

CState::CState()
{
} // default constructor

CState::CState(string element, string structure, ivec3 N, vec3 L, double interactionLength):
    structure(structure),
    N(N),
    L(L)
{
    int atoms_per_site, sites_per_cell;
    mat r; // positions;

    if (structure == "fcc")
    {
        sites_per_cell = 4;
        atoms_per_site = 1;
        //CAtom atom;
        nSites = prod(N)*sites_per_cell;
        nAtoms = nSites*atoms_per_site;
        for (int i = 0; i < nAtoms; i++)
        {
            elements.push_back(element);
        }
        size = N%L;

        r << 0.0 << 0.0 << 0.0 << endr
          << 0.5 << 0.5 << 0.0 << endr
          << 0.0 << 0.5 << 0.5 << endr
          << 0.5 << 0.0 << 0.5 << endr;
        r.col(0) *= L(0);
        r.col(1) *= L(1);
        r.col(2) *= L(2);

        makeAtoms(N, L, r, sites_per_cell);
        makeBoxes(size, interactionLength);
        fillBoxes();

    }
    else
    {
        cout << "Unknown structure" << endl;
        exit(1);
    }
}

void CState::makeAtoms(ivec3 N, vec3 L, mat r, double sites_per_cell)
{
    vec3 cellPos;

    atoms = vector<CAtom*>();
    atoms.reserve(nAtoms); // don't want resize, since we're using "push_back" to add item
    for (int ix = 0; ix < N(0); ix++)
    {
        for (int iy = 0; iy < N(1); iy++)
        {
            for (int iz = 0; iz < N(2); iz++)
            {
                for (int i = 0; i < sites_per_cell; i++)
                {
                    CAtom* atom = new CAtom;
                    cellPos << L(0)*ix << L(1)*iy << L(2)*iz;
                    atom->setPosition(cellPos + r.row(i).t());
                    atoms.push_back(atom);
                }
            }
        }
    }
}

void CState::makeBoxes(const vec3 systemSize, const double interactionLength)
{
    ivec3 boxIndex;
    vec3 boxPos;

    NBoxes = conv_to<ivec>::from(floor(systemSize/interactionLength));
    boxDimensions = systemSize/NBoxes;
    nBoxes = prod(NBoxes);

    boxes = vector<CBox*>();
    boxes.reserve(nBoxes);
    boxes.resize(nBoxes); // using resize since we aren't using push_back to add items
    for (int ix = 0; ix < NBoxes(0); ix++)
    {
        for (int iy = 0; iy < NBoxes(1); iy++)
        {
            for (int iz = 0; iz < NBoxes(2); iz++)
            {
                boxIndex << ix << iy << iz;
                boxPos = boxIndex%boxDimensions;
                CBox* box = new CBox(nBoxes, boxPos, boxDimensions, boxIndex, NBoxes);
                boxes[calculate_box_number(boxIndex, NBoxes)] = box;
            }
        }
    }

    // finding all the neighbouring boxes of each box, taking into account the
    // periodic boundary conditions
    for (int i = 0; i < nBoxes; i++)
        boxes[i]->findNeighbours(systemSize);

    cout << "box dimensions (MD): " << boxDimensions.t()
         << "               (SI): " << boxDimensions.t()*L0
         << "system size    (MD): " << systemSize.t()
         << "               (SI): " << systemSize.t()*L0
         << "number of boxes    :    " << nBoxes << endl << endl;
}

void CState::fillBoxes()
{
    vec3 atomPosition;
    ivec3 boxIndex;
    CAtom* atomptr;

    for (int i = 0; i < nAtoms; i++)
    {
        atomptr = atoms[i];
        atomPosition = atomptr->getPosition();
        for (int j = 0; j < 3; j++)
        {
            boxIndex(j) = int(floor(atomPosition(j)/boxDimensions(j)));
        }
        boxes[calculate_box_number(boxIndex, NBoxes)]->addAtom(atomptr);
    }
}

void CState::randnVelocity(double mean, double sigma_v, long *seed)
{
    vec2 randomNormal, randomUniform;
    mat velocities(nAtoms, 3);
    rowvec3 momentum;
    double R;

    for (int i = 0; i < nAtoms; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            // Box-Muller transform
            randomUniform << ran2(seed) << ran2(seed);
            R = sqrt(-2*log(randomUniform(0)));
            randomNormal(0) = R*cos(2*pi*randomUniform(1));
            velocities(i, j) = randomNormal(0)*sigma_v + mean;

            // unused random number... should probably use this one for something
            // randomNormal(1) = R*sin(2*pi*randomUniform(1))*sigma_v + mean;
        }
    }

    momentum = sum(velocities)/nAtoms; // finding the linear momentum of the system
    for (int i = 0; i < nAtoms; i++)
    {
        velocities.row(i) -= momentum; // removing initial linear momentum from the system
    }

    // sending the generated velocities to the atoms
    for (int i = 0; i < nAtoms; i++)
    {
        atoms[i]->setVelocity(velocities.row(i).t());
    }
}

void CState::randuVelocity(double mean, double vmax, long *seed)
{
    mat velocities(nAtoms, 3);
    rowvec momentum(nAtoms);

    // generating random, uniform velocities
    for (int i = 0 ; i < nAtoms; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            velocities(i, j) = (ran2(seed)-0.5)*vmax + mean;
        }
    }
    // removing any linear momentum from the system
    momentum = sum(velocities)/nAtoms;
    for(int i = 0; i < nAtoms; i++)
    {
        velocities.row(i) -= momentum;
    }
    // sending the velocities to the atoms
    for (int i = 0; i < nAtoms; i++)
    {
        atoms[i]->setVelocity(velocities.row(i).t());
    }
}

void CState::save(string filename, bool forces, bool indexing)
{
    // // check and fix the extension of the filename
    // if (filename.substr(filename.find_last_of(".")) != ".xyz")
    // {
    //     filename.append(".xyz");
    // }

    double force;

    // open file stream
    ofstream ofile;
    ofile.open(filename.c_str());

    ofile << nAtoms << endl;
    ofile << "Comment" << endl;

    //ofile.setf(ios::right); // I have no idea what this might do
    for (int i = 0; i < nAtoms; i++) {
        ofile << elements[i];
        //ofile << scientific; // uncomment if you want "0.000000e+00" formatting
        ofile << setw(16) << setprecision(8);
        for (int j = 0; j < 3; j++)
        {
            ofile << setw(18) << setprecision(8) << (atoms[i]->getPosition()(j))*L0;
        }
        for (int j = 0; j < 3; j++)
        {
            ofile << setw(16) << setprecision(8) << (atoms[i]->getVelocity()(j))*(L0/t0);
        }
        if (forces)
        {
            for (int j = 0; j < 3; j++)
            {
                ofile << setw(16) << setprecision(8) << (atoms[i]->getForce()(j))*F0;
            }
            force = norm(atoms[i]->getForce(), 2);
            ofile << setw(16) << setprecision(8) << force*F0;
        }
        if (indexing)
        {
            ofile << setw(16) << setprecision(8) << i+1;
        }
        ofile << endl;
    }
    ofile.close();

    //cout << "Exiting CState::save" << endl;
}

void CState::saveVelocity(string filename, bool MDUnits)
{
    vec3 velocity;

    // open file stream
    ofstream ofile;
    ofile.open(filename.c_str());

    //ofile.setf(ios::right); // I have no idea what this might do
    for (int i = 0; i < nAtoms; i++) {
        if (MDUnits)
            velocity = (atoms[i]->getVelocity());
        else
            velocity = (atoms[i]->getVelocity())*(L0/t0);

        //ofile << scientific; // uncomment if you want "0.000000e+00" formatting
        ofile << setw(16) << setprecision(8);
        for (int j = 0; j < 3; j++)
        {
            ofile << setw(16) << setprecision(8) << velocity(j);
        }
        ofile << setw(16) << setprecision(8) << norm(velocity, 2);
        ofile << endl;
    }
    ofile.close();
}

void CState::load(string filename)
{
    cout << "CState::load" << endl;

    ifstream ofile;
    string element, comment;
    int counter = 0;
    vec3 position, velocity;

    ofile.open(filename.c_str());
    if (ofile.is_open())
    {
        cout << "Managed to open the file!" << endl;
        if (ofile.good())
        {
            ofile >> nAtoms;
            ofile >> comment;
        }
        cout << "nAtoms = " << nAtoms << endl << "comment = " << comment << endl;
        while (counter < nAtoms*7)
        {
            if (counter%7 == 0 && ofile.good())
            {
                ofile >> element;
                elements.push_back(element);
                // cout << "i = " << counter << ", " << element << endl;
                counter++;
            }
            else if (ofile.good())
            {
                ofile >> position(0);
                ofile >> position(1);
                ofile >> position(2);
                ofile >> velocity(0);
                ofile >> velocity(1);
                ofile >> velocity(2);
                counter += 6;
                // cout << counter/7-1 << endl;
                atoms[counter/7-1]->setPosition(position/L0);
                atoms[counter/7-1]->setVelocity(velocity/(L0/t0));
            }
        }
    }
    else cout << "Unable to open file!" << endl;
}

void CState::move(const double dt, const bool statistics)
{
    vec3 position, velocity, force;
    vec3 newposition, newvelocity;
    vec3 vdt2;
    mat vHalf(3, nAtoms);
    ivec3 iBoundaryCrossings;
    vec3 fBoundaryCrossings;

    // first we find the new positions for all the atoms, and v(t + dt/2)
    for (int i = 0; i < nAtoms; i++)
    {
        // get data from the atom
        atoms[i]->getData(position, velocity, force);
        // first half of Verlet algorithm
        vdt2 = velocity + (force/2.0)*dt; // find v(t + dt/2)
        vHalf.col(i) = vdt2;                // store this for later
        newposition = position + vdt2*dt;   // new position

        // check if the new position is outside the periodic boundaries, count
        // how many borders it has crossed, and translate it inside the box
        for (int j = 0; j < 3; j++)
        {
            fBoundaryCrossings(j) = floor(newposition(j)/size(j));
            iBoundaryCrossings(j) = int(fBoundaryCrossings(j));
            // we could just use the integer version, but this gives overflow
            // much faster than the float version
        }
        if (iBoundaryCrossings(0) != 0 || iBoundaryCrossings(1) != 0 || iBoundaryCrossings(2) != 0)
        {
            //newposition -= iBoundaryCrossings%size;
            newposition -= fBoundaryCrossings%size; // see comment above
            atoms[i]->addToBoundaryCrossings(iBoundaryCrossings);
        }

        // sending the new position to the atom
        atoms[i]->setNewPosition(newposition);
    }

    // looping through all boxes, checking which (if any) atoms have moved
    // outside the boundaries of their box
    linkedList<CAtom*> boxlessAtoms; // linked list of atoms that have moved outside their boxes
    for (int i = 0; i < nBoxes; i++)
    {
        // this method checks if any of the atoms in a box has moved outside the
        // box, if so removes it from the box, and adds the offending atoms to
        // the linked list we supply
        boxes[i]->purgeAtoms(boxlessAtoms);
    }

    vec3 atomPosition;
    ivec3 boxIndex;
    CAtom* atomptr;
    // traversing the linked list of atoms that have moved outside their boxes,
    // which we have removed from the box (above), and adding them to their
    // correct boxes
    int imax = boxlessAtoms.length();
    // since I have a stupid linked list that can't be empty, I have to avoid
    // accessing the item generated when initializing a list without any
    // arguments. Therefore I avoid the _last_ item in the list (new items are
    // pushed to the top when adding to the list)
    for (int i = 0; i < imax-1; i++)
    {
        atomptr = boxlessAtoms(); // the current atom/item in the list
        atomPosition = atomptr->getNewPosition();

        // finding the index of the box this atom belongs in
        for (int j = 0; j < 3; j++)
            boxIndex(j) = int(floor(atomPosition(j)/boxDimensions(j)));

        // adding the atom to the box it belongs in
        boxes[calculate_box_number(boxIndex, NBoxes)]->addAtom(atomptr);

        // moving one step forward in the list
        boxlessAtoms = *boxlessAtoms.readNext();
    }

    // then we find the new forces for all the atoms
    // this needs to be done _after_ finding all the new positions, since we
    // have to find the forces on each atom (in it's new position) from all the
    // other atoms (in their new positions)
    if (statistics) newForcesAndStatistics();
    else newForces();

    // then we find the new velocities for all the atoms, and move them forwards in time
    for (int i = 0; i < nAtoms; i++)
    {
        newvelocity = vHalf.col(i) + (atoms[i]->getNewForce()/2.0)*dt;
        atoms[i]->setNewVelocity(newvelocity);
        // actually moving the atom forward one timestep, and zeroing out the "new" stuff
        atoms[i]->forward();
    }

    //cout << "Exiting CState::move" << endl;
}

void CState::newForces()
{
    vec3 zerovec3 = zeros<vec>(3);

    // first we zero out all the new forces (from last time)
    for (int i = 0; i < nAtoms; i ++)
        atoms[i]->setNewForce(zerovec3);
    for (int i = 0; i < nBoxes; i++)
        boxes[i]->resetForcesBool();

    // then we find the forces between all atoms in the box of each atom on the
    // atom itself, and the force from all boxes on all atoms, remembering
    // Newton's third law
    for (int i = 0; i < nBoxes; i++)
        boxes[i]->calculateForces(boxes);
}

void CState::newForcesAndStatistics()
{
    vec3 zerovec3 = zeros<vec>(3);

    // first we zero out all the new forces (from last time)
    for (int i = 0; i < nAtoms; i ++)
    {
        atoms[i]->resetStatistics();
        atoms[i]->setNewForce(zerovec3);
    }
    for (int i = 0; i < nBoxes; i++)
        boxes[i]->resetForcesBool();

    // then we find the forces between all atoms in the box of each atom on the
    // atom itself, and the force from all boxes on all atoms, remembering
    // Newton's third law
    for (int i = 0; i < nBoxes; i++)
        boxes[i]->calculateForcesAndStatistics(boxes);
}

void CState::berendsen(const double Tbath, const double T, const double tt)
{
    double gamma = sqrt(1 + tt*(Tbath/T - 1));
    vec3 velocity;

    if(!isfinite(gamma)) {
        cout << endl;
        cout << "! You're trying to use the berendsen thermostat without having sampled ";
        cout << "the temperature with the statistics sampler. This will give ";
        cout << "infinite velocities...! Exiting now." << endl << endl;
        exit(1);
    }

    for (int i = 0; i < nAtoms; i++)
    {
        velocity = atoms[i]->getVelocity();
        velocity *= gamma;
        atoms[i]->setVelocity(velocity);
    }
}

void CState::andersen(const double Tbath, const double tt, long* seed)
{
    double sigma = sqrt(Tbath), R;
    vec3 velocity;
    double randU[2];
    double randN[2];

    for (int i = 0; i < nAtoms; i++)
    {
        if (ran2(seed) < tt)
        {
            for (int j = 0; j < 3; j++)
            {
                randU[0] = ran2(seed);
                randU[1] = ran2(seed);
                R = sqrt(-2*log(randU[0]));
                randN[0] = R*cos(2.0*pi*randU[1]);
                // randN[1] = R*sin(2.0*pi*randU[1]);
                velocity(j) = randN[0]*sigma;
            }
            atoms[i]->setVelocity(velocity);
        }
    }
}

void CState::average()
{
    double vAverage = 0, vSigma = 0;
    double xAverage = 0, xSigma = 0;
    vec3 velocity;

    for (int i = 0; i < nAtoms; i++)
    {
        velocity = atoms[i]->getVelocity();
        vAverage += norm(velocity, 2);
        xAverage += velocity(0);
    }
    vAverage /= nAtoms;
    xAverage /= nAtoms;
    for (int i = 0; i < nAtoms; i++)
    {
        velocity = atoms[i]->getVelocity();
        vSigma += pow(norm(velocity, 2) - vAverage, 2.0);
        xSigma += pow(velocity(0) - xAverage, 2.0);
    }
    vSigma = sqrt(vSigma/nAtoms);
    xSigma = sqrt(xSigma/nAtoms);

    cout << "average velocity (r)       = " << vAverage << endl;
    cout << "standard deviation of v(r) = " << vSigma << endl;
    cout << "average v_x                = " << xAverage << endl;
    cout << "standard deviation of v_x  = " << xSigma << endl << endl;
}

int CState::getnAtoms() const
{
    return nAtoms;
}

int CState::getnBoxes() const
{
    return nBoxes;
}

vec3 CState::getSize() const
{
    return size;
}

CAtom CState::getAtom(int i) const
{
    return *atoms[i];
}

vector<CAtom*> CState::getAtoms() const
{
    return atoms;
}

vector<CBox*> CState::getBoxes() const
{
    return boxes;
}

CBox CState::getBox(int i) const
{
    return *boxes[i];
}
