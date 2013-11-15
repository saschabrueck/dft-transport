#include <math.h>   // sqrt()
#include "AtomStrain.H"

/************************************************************************************************/
AtomStrain::AtomStrain(WireGenerator *wire, Material *mat)
{
    // Lattice constants
    a0 = wire->nwire->a0;
    ideal_a0 = mat->ideal_a0;
    a0_scale = ideal_a0/a0;

    // Simulation size
    Lx = wire->lenx;
    Ly = wire->leny;
    Lz = wire->lenz;
    Lx_half = 0.5*Lx;
    Ly_half = 0.5*Ly;
    Lz_half = 0.5*Lz;

    /* XXX This is zincblende specific.  If arbitrary lattices are implemented, we should 
     * probably get this as a parameter from wire.
     */
    Nneighbors = 4;

    /* Set the number of atoms involved in the strain calculation.  For now, this includes the
     * atoms in the channel (Layer_Matrix), and the contacts (Boundary[0,1]).  This can be 
     * changed to include any arbitrary domain (e.g. a large strain domain around a smaller
     * electronic structure domain in a quantum dot).
     */
    Natoms = wire->No_Atom + wire->Boundary[0]->NA + wire->Boundary[1]->NA;

    /* Pointers to the atomic positions.  This is complicated, because in a nanowire simulation,
     * the atomic coordinates are stored in different structures.  The channel atoms are stored
     * in wire->Layer_Matrix and the contacts are stored in wire->Boundary[0,1].  One possible 
     * approach is to copy all the data from Layer_Matrix and Boundary[*] into a separate array
     * to do the strain calculation.  This effectively doubles the storage requirements for the
     * atomic coordinates.  The approach below uses the AtomPtr struct, defined in AtomStrain.H,
     * to point to the beginning of the coordinates in their current location (i.e. it might 
     * point to an address in Layer_Matrix or Boundary[0,1]).  The assumption is that the 
     * coordinates are stored in the same manner as Layer_Matrix or Boundary[0,1]; namely that
     *  *p     is the x coordinate,
     *  *(p+1) is the y coordinate,
     *  *(p+2) is the z coordinate,
     *  *(p+3) is the type.
     *  *(p+4+N) is the index of a nearest neighbor atom in the same structure, for 
     *           N = 0, ..., Nneighbors-1.
     * The nbrs[] field of the AtomPtr structure re-calculates the indices of the nearest
     * neighbor atoms, shifting them from their values in the original structure (Layer_Matrix,
     * Boundary[*]) to the appropriate value in the Atoms[] array.
     *
     * This approach requires less memory than copying all the coordinate and neighbor data into
     * a new array, at the cost of some complexity.  It does require more memory than merely 
     * pointing to the atomic coordinates in place.  This method can also be extended to arbitrary
     * systems (e.g. to consider a huge strain domain around a smaller electronic structure 
     * domain in a quanum dot, for example, allocate the surrounding strain domain and populate 
     * it.  Then make the Atoms[] array point to the appropriate coordinates in either the 
     * surrounding strain domain or the electronic structure domain; the strain calculation below
     * is unaffected, because it uses only the pointers in the Atoms[] array.)  Thus, only this
     * constructor requires modification to allow strain calculation for different types of systems.
     */
    Atoms = new AtomPtr [Natoms];

    int i, j, k, m;
    int base;

    for (i = j = 0; j < wire->Boundary[0]->NA; i+=4+Nneighbors, j++)
    {   // i loops over indices in Boundary[0]->UNN*, j over indices in Atoms
        Atoms[j].p = &wire->Boundary[0]->UNNm[i];
        Atoms[j].nbrs = new int [Nneighbors];

        for (k = 0; k < Nneighbors; k++)
        {
            m = (int)wire->Boundary[0]->UNN[i+4+k];
            if ((-1 == m) && (1 == wire->nwire->strain_bc))  // periodic in transport direction
            {
               if (i < 0.5*wire->Boundary[0]->NA)
               {   // index in Boundary[1]->UNNp
                   Atoms[j].nbrs[k] = wire->Boundary[0]->NA + wire->No_Atom + 
                                      wire->Boundary[0]->UNNm[i+4+k]-1;
               }
               else
               {   // index in Layer_Matrix
                   Atoms[j].nbrs[k] = wire->Boundary[0]->NA + wire->Layer_Matrix[i+4+k]-1;
               }
            }
            else if (m > i)  // avoid m == 0 (dangling bond) and m < i (double counting)
                Atoms[j].nbrs[k] = m-1;  // index in Boundary[0]->UNNm
            else
                Atoms[j].nbrs[k] = NO_NEIGHBOR;
        }
    }
    base = j;   // wire->Boundary[0]->NA
    for (i = 0; j < base+wire->No_Atom; i+= 4+Nneighbors, j++)
    {   // i loops over indices in Layer_Matrix, j continues over indices in Atoms
        Atoms[j].p = &wire->Layer_Matrix[i];
        Atoms[j].nbrs = new int [Nneighbors];

        for (k = 0; k < Nneighbors; k++)
        {
            m = (int)wire->Layer_Matrix[i+4+k];
            if ((-1 == m) && (1 == wire->nwire->strain_bc))  // periodic in transport direction
            {  
                if (i > 0.5*wire->No_Atom)
                {   // index in Boundary[1]->UNNp
                    int off = wire->Boundary[1]->NA - (wire->No_Atom - (j-base));

                    Atoms[j].nbrs[k] = base + wire->No_Atom + 
                                       wire->Boundary[1]->UNNp[off*(4+Nneighbors)+4+k]-1;
                }
                else
                {
                    /* Avoid indices (<j) in Boundary[0]->UNNm (double counting):
                     * if (i < 0.5*wire->No_Atom)  // index in Boundary[0]->UNNm
                     *     Atoms[j].nbrs[k] = wire->Boundary[0]->UNNm[i+4+k]-1;
                     */
                    Atoms[j].nbrs[k] = NO_NEIGHBOR;
                }
            }
            else if (m > i)  // avoid m == 0 (dangling bond) and m < i (double counting)
                Atoms[j].nbrs[k] = base+m-1;   // index in Layer_Matrix
            else
                Atoms[j].nbrs[k] = NO_NEIGHBOR;
        }
    }
    base = j;  // wire->Boundary[0]->NA + wire->No_Atom
    for (i = 0; j < base+wire->Boundary[1]->NA; i+= 4+Nneighbors, j++)
    {   // i loops over indices in Boundary[1]->UNN*, j continues over indices in Atoms
        Atoms[j].p = &wire->Boundary[1]->UNNp[i];
        Atoms[j].nbrs = new int [Nneighbors];

        for (k = 0; k < Nneighbors; k++)
        {
            m = (int)wire->Boundary[1]->UNN[i+4+k];
            /* Don't need to consider any periodic connections (-1 == m) here: connections from 
             * Layer_Matrix have already been counted, as have connections to Boundary[0]->UNNm.
             */
            if (m > i)  // avoid m == 0 (dangling bond) and m < i (double counting)
                Atoms[j].nbrs[k] = base+m-1;
            else
                Atoms[j].nbrs[k] = NO_NEIGHBOR;
        }
    }
}


/************************************************************************************************/
AtomStrain::~AtomStrain()
{

    delete[] Atoms;

}


/************************************************************************************************/
/* Return the unstrained (equilibrium, perfect crystal) distance between the passed atoms 
 * (atom indices), with the vector components stored in v.
 * 
 * The logic behind this method is pretty convoluted: NEMO 3D defines the unstrained 
 * coordinates in terms of an ideal lattice constant (ideal_a0 here, unstrnd_cubic_cell_length 
 * in the N3D input deck) and the strained coordinates in terms of a user-specified lattice 
 * constant (a0 here), plus strain displacement.  Doing so allows the strained lattice constants
 * to change.   The coordinates pointed to by Atoms[] are calculated using the user-specified a0,
 * so if a0 != ideal_a0, then the Atoms[] coordinates correspond to the "strained" coordinates
 * of NEMO 3D (minus strain displacements).  If a0 == ideal_a0, the Atoms[] coordinates are
 * identical to the "unstrained" coordinates of NEMO 3D.  This function returns the unstrained
 * distance in the same manner as NEMO 3D.
 */
double 
AtomStrain::get_unstrained_distance(double *v, int atom1, int atom2)
{
    v[0] = (*(Atoms[atom2].p)   - *(Atoms[atom1].p))   * a0_scale;
    v[1] = (*(Atoms[atom2].p+1) - *(Atoms[atom1].p+1)) * a0_scale;
    v[2] = (*(Atoms[atom2].p+2) - *(Atoms[atom1].p+2)) * a0_scale;

    /* Here we assume that if dim is greater than half the dimension, we need to shift the
     * the coordinates to account for periodic boundarys.  The check of the BCs occurs at a
     * higher level and results in the recognition (or not) of nearest neighbors.  If atoms on
     * the edges are _not_ periodic neighbors, they won't be recognized as neighbors, and this
     * function won't be called.
     */
#define PERIODIC_SHIFT(dim, Lhalf, Lfull)  \
    if (dim > Lhalf)  {dim -= Lfull;} else if (dim < -Lhalf)  {dim += Lfull;}
    PERIODIC_SHIFT(v[0], Lx_half, Lx);
    PERIODIC_SHIFT(v[1], Ly_half, Ly);
    PERIODIC_SHIFT(v[2], Lz_half, Lz);
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}


/************************************************************************************************/
VFFStrain::VFFStrain(WireGenerator *wire, Material *mat) : AtomStrain(wire, mat)
{
    if (wire->nwire->VFF_alpha > 0.0)  // specified in input deck
        alpha = wire->nwire->VFF_alpha;
    else
        alpha = mat->alpha;
    if (wire->nwire->VFF_beta > 0.0)   // specified in input deck
        beta = wire->nwire->VFF_beta;
    else
        beta = mat->beta;

    displen = 3*Natoms;
    disp = new double[displen];

    for (int i = 0; i < displen; i++)
       disp[i] = 0.0;
}


/************************************************************************************************/
VFFStrain::~VFFStrain()
{
    int i, j;

    // Add final displacements to atomic coordinates
    for (i = 0; i < Natoms; i++)
    {
        j = 3*i;
        *(Atoms[i].p)   += disp[j];
        *(Atoms[i].p+1) += disp[j+1];
        *(Atoms[i].p+2) += disp[j+2];
    }
    delete[] disp;

}

/************************************************************************************************/
/* Return the VFF alpha (bond length distortion) parameter for the atoms a1 and a2 (inidces in 
 * Atoms[]).  
 * XXX Currently OMEN supports only 1 material; if multiple materials are supported in the future,
 * this may need to be an array of alphas between all possible atom types.
 */
double VFFStrain::get_alpha(int a1, int a2)
{
   return alpha;
}

/************************************************************************************************/
/* Return the VFF beta (bond angle distortion) parameter for the atoms a1 and a2 (inidces in 
 * Atoms[]).  
 * XXX Currently OMEN supports only 1 material; if multiple materials are supported in the future,
 * this may need to be an array of betas between all possible atom types.
 */
double VFFStrain::get_beta(int a1, int a2)
{
   return beta;
}

/************************************************************************************************/
/* TODO Calculate the strain energy.  This is needed only for maccheckgrad(), and is not used
 * in this strain calculation.
 */
double VFFStrain::energy(double *x)
{
    return 1.0;
}

/************************************************************************************************/
/* Calculate the strain energy gradient.  x[] holds the current displacements, and the new 
 * gradient is stored in gradx[].
 */
void VFFStrain::gradient(double *x, double *gradx)
{
    int i, j, k, n;
    int ioff, joff, koff;
    double dxj, dyj, dzj;
    double dxk, dyk, dzk;
    double d1, d2;
    double b1, b2;
    double z, pref;
    double x0j[3];
    double x0k[3];

    for (i = 0; i < displen; i++)
        gradx[i] = 0.0;

    // Loop over all atoms
    for (i = 0; i < Natoms; i++)
    {
        ioff = 3*i;  // first index of ith atom in x[]
        // Loop over all neighbors
        for (j = 0; j < Nneighbors; j++)
        {
            n = Atoms[i].nbrs[j];
            if (NO_NEIGHBOR != n)
            {   // Bond length distortion
                joff = 3*n;   // first index of nth atom in x[]
                dxj = (x[joff]   - x[ioff])   + (*(Atoms[n].p)   - *(Atoms[i].p));
                dyj = (x[joff+1] - x[ioff+1]) + (*(Atoms[n].p+1) - *(Atoms[i].p+1));
                dzj = (x[joff+2] - x[ioff+2]) + (*(Atoms[n].p+2) - *(Atoms[i].p+2));
                d1 = get_unstrained_distance(x0j, i, n);
                pref = 2.0*get_alpha(i, n)*(dxj*dxj+dyj*dyj+dzj*dzj-d1*d1)/(d1*d1);
                z = pref*dxj;
                gradx[ioff]   -= z;
                gradx[joff]   += z;
                z = pref*dyj;
                gradx[ioff+1] -= z;
                gradx[joff+1] += z;
                z = pref*dzj;
                gradx[ioff+2] -= z;
                gradx[joff+2] += z;

                // Bond angle distortion
                b1 = get_beta(i, n);
                for (k = j+1; k < Nneighbors; k++)
                {
                    if (Atoms[i].nbrs[k] > n)
                    {
                        n = Atoms[i].nbrs[k];
                        koff = 3*n;
                        dxk = (x[koff]   - x[ioff])   + (*(Atoms[n].p)   - *(Atoms[i].p));
                        dyk = (x[koff+1] - x[ioff+1]) + (*(Atoms[n].p+1) - *(Atoms[i].p+1));
                        dzk = (x[koff+2] - x[ioff+2]) + (*(Atoms[n].p+2) - *(Atoms[i].p+2));
                        d2 = get_unstrained_distance(x0k, i, n);
                        b2 = get_beta(i, n);
                        pref = 2.0*sqrt(b1*b2)*((dxj*dxk+dyj*dyk+dzj*dzk) - 
                               (x0j[0]*x0k[0]+x0j[1]*x0k[1]+x0j[2]*x0k[2]))/(d1*d2);
                        gradx[ioff]   -= pref*(dxk+dxj);
                        gradx[joff]   += pref*dxk;
                        gradx[koff]   += pref*dxj;
                        gradx[ioff+1] -= pref*(dyk+dyj);
                        gradx[joff+1] += pref*dyk;
                        gradx[koff+1] += pref*dyj;
                        gradx[ioff+2] -= pref*(dzk+dzj);
                        gradx[joff+2] += pref*dzk;
                        gradx[koff+2] += pref*dzj;
                    }
                }   // k
            }
        }   // j
    }   // i
}


/************************************************************************************************/
MacoptStrain::MacoptStrain(WireGenerator *wire, Material *mat) : 
    VFFStrain(wire, mat), Macopt(3*Natoms,    0,     1e-16,    1e4,    1)
//          XXX Macopt defaults:   N       verbose, tolerance, itmax, rich
{
}

/************************************************************************************************/
MacoptStrain::~MacoptStrain()
{
}

/************************************************************************************************/
void MacoptStrain::relax_atoms(void)
{
    macoptII(disp, displen);
}

/************************************************************************************************/
double MacoptStrain::func(double *x)
{
    return energy(x);
}

/************************************************************************************************/
void MacoptStrain::dfunc(double *x, double *gradx)
{
    gradient(x, gradx);
}

/* TODO 
 * doxygen-ize
 * parallelize: 3D decomposition and communication
 * larger domain, for long range strain effects
 * fix certain atoms, depending on BCs
 */
