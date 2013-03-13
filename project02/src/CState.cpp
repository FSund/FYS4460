#include "CState.h"

CState::CState()
{
} // default constructor

CState::CState(string element, string structure, ivec3 N, vec3 L, double interactionLength):
//    structure(structure),
    nMovingAtoms(0),
    nMatrixAtoms(0),
    N(N),
    L(L)
{
    int atoms_per_site, sites_per_cell;
    mat r; // positions;

    if (structure == "fcc")
    {
        sites_per_cell = 4;
        atoms_per_site = 1;

        nSites = prod(N)*sites_per_cell;
        nAtoms = nSites*atoms_per_site;
        nMovingAtoms = nAtoms;
        size = N%L;

//        for (int i = 0; i < nAtoms; i++)
//        {
//            elements.push_back(element);
//        }

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

//CState::CState(const CState &state, const ivec3 &N, const vec3 &L, const double &interactionLength):
//    atoms(state.atoms),
////    elements(state.elements),
//    nAtoms(state.atoms.size()),
//    N(N),
//    L(L),
//    size(N%L)
//{
//    makeBoxes(size, interactionLength);
//    fillBoxes();
//}

CState::CState(const string filename, ivec3 &N, vec3 &L, double &interactionLength):
    nMovingAtoms(0),
    nMatrixAtoms(0),
    N(N),
    L(L),
    size(N%L)
{
    if (load(filename)) // takes care of nAtoms, atoms, atomType and matrixAtom
        exit(1);
    makeBoxes(size, interactionLength);
    fillBoxes();
}

CState::~CState()
{
    for (int i = 0; i < nAtoms; i++)
        delete atoms[i];
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
    cout << "CState::fillBoxes()" << endl;

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

    cout << "Exiting CState::fillBoxes()" << endl << endl;
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

void CState::generate_spherical_pores(long *seed, int N_pores, double min_pore_radius, double max_pore_radius)
{
    double pore_radius = min_pore_radius/L0; // converting from Angstrom to MD units
    double pore_radius_dr = (max_pore_radius - min_pore_radius)/L0;
    mat pore_centers(N_pores, 3);
    vec pore_radi(N_pores);

    for (int i = 0; i < N_pores; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            pore_centers(i, j) = ran2(seed)*size(j);
        }
        pore_radi(i) = pore_radius + ran2(seed)*pore_radius_dr;
    }

    // PBC
    int ii = 0;
    mat displacementMat(27, 3);
    vec3 displacement;
    for (int i = -1; i <= 1; i++)
    {
        for (int j = -1; j <= 1; j++)
        {
            for (int k = -1; k <= 1; k++)
            {
                displacement << i*size(0) << j*size(1) << k*size(2);
                displacementMat.row(ii) = displacement;
                ii++;
            }
        }
    }

    vec3 rvec;
    double dr;
    string atomType;
    for (int poreNr = 0; poreNr < N_pores; poreNr++)
    {
        for (int i = 0; i < nAtoms; i++)
        {
            for (int j = 0; j < 27; j++)
            {
                rvec = atoms[i]->getPosition().t() - pore_centers.row(poreNr)
                        + displacementMat.row(j);
                dr = norm(rvec, 2);
                if (dr <= pore_radi(poreNr) && !atoms[i]->matrixAtom)
                {
                    atomType = atoms[i]->getAtomType();
                    atomType.append("_m");
                    atoms[i]->setAtomType(atomType);
                    nMatrixAtoms++;
                    nMovingAtoms--;
                }
            }
        }
    }
}

void CState::generate_cylindrical_pore(const double radius)
{
//    double pore_radius = min_pore_radius/L0; // converting from Angstrom to MD units
//    double pore_radius_dr = (max_pore_radius - min_pore_radius)/L0;

    int N_pores = 1;
//    mat pore_centers(N_pores, 3);

//    vec3 pore_radi;
//    pore_radi = size/12;
    double pore_dr = radius/L0;
    vec3 pore_thing;
    pore_thing << 0 << size(1)/2.0 << size(2)/2.0;

    vec3 rvec;
    string atomType;
    double dr;
    for (int poreNr = 0; poreNr < N_pores; poreNr++)
    {
        for (int i = 0; i < nAtoms; i++)
        {
            rvec = atoms[i]->getPosition() - pore_thing;
            dr = norm(rvec.rows(1,2),2);
            if (dr > pore_dr && !atoms[i]->matrixAtom)
            {
                atomType = atoms[i]->getAtomType();
                atomType.append("_m");
                atoms[i]->setAtomType(atomType);
                nMatrixAtoms++;
                nMovingAtoms--;
            }
        }
    }
}

void CState::FILIP_pores()
{
//    int N_pores = 1;
//    nMatrixAtoms???
//    double pore_dr = size(0)/6;
//    vec3 pore_thing;
//    pore_thing << 0 << size(1)/2.0 << size(2)/2.0;
//    vec3 rvec;
//    string atomType;
//    double dr;
//    for (int poreNr = 0; poreNr < N_pores; poreNr++)
//    {
//        for (int i = 0; i < nAtoms; i++)
//        {
//            rvec = atoms[i]->getPosition() - pore_thing;
//            dr = norm(rvec.rows(1,2),2);
//            rvec = atoms[i]->getPosition();
//            if (dr > pore_dr && !atoms[i]->matrixAtom && rvec(0)<(2.0*size(0)/3.0) && rvec(0)>(size(0)/3.0))
//            {
//                atomType = atoms[i]->getAtomType();
//                atomType.append("_m");
//                atoms[i]->setAtomType(atomType);
//            }
//        }
//    }
}

void CState::remove_half_the_atoms()
{
    cout << "CState::remove_half_the_atoms" << endl;

    vector<CAtom*> newAtoms;
    long idum = -1;
    for (int i = 0; i < nAtoms; i++)
    {
        if (!atoms[i]->matrixAtom && ran2(&idum) < 0.5)
            newAtoms.push_back(atoms[i]);
        else if (atoms[i]->matrixAtom)
            newAtoms.push_back(atoms[i]);
        else
            delete atoms[i];
    }

    atoms = newAtoms;
    nAtoms = atoms.size();
    nMovingAtoms = nAtoms - nMatrixAtoms;
    for (int i = 0; i < nBoxes; i++)
    {
        boxes[i]->flush();
    }
    fillBoxes();

    cout << "New nAtoms = " << nAtoms << endl;
    cout << "Exiting CState::remove_half_the_atoms" << endl << endl;
}

void CState::decrease_density_by_factor(double &factor)
{
    if (factor > 1.0 || factor < 0.0)
    {
        cout << "! Invalid factor to decrease density by. I won't do it !" << endl;
        return;
    }

    vector<CAtom*> newAtoms;
    long idum = -1;
    for (int i = 0; i < nAtoms; i++)
    {
        if (!atoms[i]->matrixAtom && ran2(&idum) < factor)
            newAtoms.push_back(atoms[i]); // keeping some "fluid" atoms
        else if (atoms[i]->matrixAtom)
            newAtoms.push_back(atoms[i]); // keping all matrix atoms
        else
            delete atoms[i]; // deleting all other atoms
    }

    atoms = newAtoms;
    nAtoms = atoms.size();
    nMovingAtoms = nAtoms - nMatrixAtoms;
    for (int i = 0; i < nBoxes; i++)
    {
        boxes[i]->flush();
    }
    fillBoxes();
}

void CState::save(string filename, bool saveSpeed, bool saveForces, bool indexing, bool markMatrixAtoms)
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

    if (!ofile)
    {
        cout << endl;
        cout << "! It seems like the '.xyz'-file couldn't be created. ";
        cout << "No .xyz-files will be saved. Please fix the folder or filename and rerun if you ";
        cout << "want to save these files." << endl << endl;
        return;
    }

    ofile << nAtoms << endl;
    ofile << "Comment" << endl;

    for (int i = 0; i < nAtoms; i++) {
        ofile << left << setfill(' ') << setw(16) << atoms[i]->getAtomType();
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
        if (saveSpeed)
        {
            ofile << setw(16) << setprecision(8) << norm(atoms[i]->getVelocity(), 2)*(L0/t0);
        }
        if (saveForces)
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
        if (markMatrixAtoms)
        {
            ofile << setw(16) << setprecision(8) << atoms[i]->matrixAtom;
        }

        // print "pressure" on each atom
        ofile << setw(16) << setprecision(8) << atoms[i]->getPressure();

        ofile << endl;
    }
    ofile.close();

    //cout << "Exiting CState::save" << endl;
}

void CState::saveMatrix(string filename)
{
    // open file stream
    ofstream ofile;
    ofile.open(filename.c_str());

    if (!ofile)
    {
        cout << endl;
        cout << "! It seems like the '.xyz'-file couldn't be created. ";
        cout << "No .xyz-files will be saved. Please fix the folder or filename and rerun if you ";
        cout << "want to save these files." << endl << endl;
        return;
    }

    // counting number of matrix atoms
    int nMatrixAtoms = 0;
    for (int i = 0; i < nAtoms; i++)
        if (atoms[i]->matrixAtom)
            nMatrixAtoms++;

    ofile << nMatrixAtoms << endl;
    ofile << "MatrixAtoms" << endl;

    for (int i = 0; i < nAtoms; i++) {
        if (atoms[i]->matrixAtom)
        {
            ofile << left << setfill(' ') << setw(16) << atoms[i]->getAtomType();
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
            ofile << endl;
        }
    }
    ofile.close();
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

int CState::load(string filename)
{
    // ONLY FOR USE IN THE STATE CONSTRUCTOR !!!

    cout << endl << "CState::load" << endl;

    ifstream ofile;
    string atomType;
    vec3 position, velocity;
    char Comment[256];

    ofile.open(filename.c_str());
    if (ofile.is_open())
    {
        cout << "Managed to open the file!" << endl;

        if (ofile.good())
        {
            ofile >> nAtoms;
            ofile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            ofile.getline(Comment, 256);
        }
        nMovingAtoms = nAtoms;
        cout << "nAtoms = " << nAtoms << endl;
        cout << "Comment: " << Comment << endl;

        atoms.reserve(nAtoms);

        while (ofile.good() && !ofile.eof())
        {
            ofile >> atomType;
            for (int dim = 0; dim < 3; dim++)
                ofile >> position(dim);
            position /= L0;
            for (int dim = 0; dim < 3; dim++)
                ofile >> velocity(dim);
            velocity *= (t0/L0);

            // ignore until end of line
            ofile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            // creating the atom (pointer) and push to vector
            CAtom* atom = new CAtom(position, velocity);
            nMatrixAtoms += atom->setAtomType(atomType);
            atoms.push_back(atom);
        }
        cout << "Exiting CState::load, returning 0" << endl << endl;
        nMovingAtoms = nAtoms - nMatrixAtoms;

        return 0;
    }
    else
    {
        cout << "! Unable to open file, returning 1 !" << endl;
        return 1;
    }

//    cout << "Exiting CState::load" << endl << endl;
}

void CState::move(const double dt, const bool statistics)
{
    vec3 position, velocity, force;
    vec3 newposition, newvelocity;
    vec3 vdt2;
//    mat vHalf(3, nAtoms);
    mat vHalf = zeros<mat>(3, nAtoms);
    ivec3 iBoundaryCrossings;
    vec3 fBoundaryCrossings;

    // first we find the new positions for all the atoms, and v(t + dt/2)
    for (int i = 0; i < nAtoms; i++)
    {
        if (atoms[i]->matrixAtom)
        {
            // vHalf.col(i) = zeros<mat>(3, 1); (done in initialization)
            continue;
        }

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

    // flow thing
    gravity();

    // then we find the new velocities for all the atoms, and move them forwards in time
    for (int i = 0; i < nAtoms; i++)
    {
        // checking that none of the (new) forces are larger than a set cutoff
        // this is to avoid infinite velocities etc...
//        atoms[i]->checkForceCutoff();

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
//    #pragma omp parallel for
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

void CState::gravity()
{
//    const double Fx = 0.1*E0/L0;
//    vec3 forcevec;
//    forcevec << Fx << 0.0 << 0.0;

    vec3 forcevec;
    forcevec << 0.1*E0/L0 << 0.0 << 0.0;

    for (int i = 0; i < nAtoms; i++)
        atoms[i]->addToNewForce(forcevec);
}

void CState::berendsen(const double Tbath, const double T, const double tt)
{
    double gamma = sqrt(1 + tt*(Tbath/T - 1));
    vec3 velocity;

    if(!isfinite(gamma)) {
//        cout << endl;
//        cout << "! You're trying to use the berendsen thermostat without having sampled ";
//        cout << "the temperature with the statistics sampler. This will give ";
//        cout << "infinite velocities...! Exiting now." << endl << endl;
//        exit(1);

        cout << endl;
        cout << "! The berendsen thermostat reports an infinite gamma-factor. ";
        cout << "This is usually caused by not sampling the temperature, but ";
        cout << "could also be something else. NOT using the thermostat this ";
        cout << "timestep. !" << endl;
        return;
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

int CState::getnMovingAtoms() const
{
    return nMovingAtoms;
}

int CState::getnMatrixAtoms() const
{
    return nMatrixAtoms;
}

int CState::getnBoxes() const
{
    return nBoxes;
}

vec3 CState::getSize() const
{
    return size;
}

/* Inlined in header file*/

//const CAtom &CState::getAtom(int i) const
//{
//    return *atoms[i];
//}

/* Unused stuff */

//vector<CAtom*> CState::getAtoms() const
//{
//    return atoms;
//}

//vector<CBox*> CState::getBoxes() const
//{
//    return boxes;
//}

//CBox CState::getBox(int i) const
//{
//    return *boxes[i];
//}
