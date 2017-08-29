/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __ATOMSTRAIN
#define __ATOMSTRAIN

#include "WireGenerator.H"  // includes Material.H
#include "Macopt.H"

/************************************************************************************************/
struct AtomPtr
{
    double *p;   // pointer to beginning of atomic coordinates, in some other data structure
    int *nbrs;   // array of size Nneighbors, holding indices of nearest neighbors in Atoms array
};

#define NO_NEIGHBOR  0

/************************************************************************************************/
class AtomStrain  // base class
{
    public:
        AtomStrain(WireGenerator *, Material *);
        ~AtomStrain();
        virtual void relax_atoms(void)=0;

    protected:
        double get_unstrained_distance(double*, int, int);

        int Natoms;           // total number of atoms in strain calculation
        int Nneighbors;       // number of nearest neighbors
        AtomPtr *Atoms;       // array of size Natoms
        double a0, ideal_a0;  // actual and ideal (unstrained) lattice constants
        double a0_scale;      // ratio of actual and ideal lattice constants
        double Lx, Lx_half;   // simulation domain size, and half size, in x
        double Ly, Ly_half;   // simulation domain size, and half size, in y
        double Lz, Lz_half;   // simulation domain size, and half size, in z
};

/************************************************************************************************/
class VFFStrain : public AtomStrain  // Use the Keating Valence Force Field method for strain
{
    public:
        VFFStrain(WireGenerator *, Material *);
        ~VFFStrain();

    protected:
        double get_alpha(int, int);
        double get_beta(int, int);
        double energy(double *);
        void gradient(double *, double *);

        double alpha;  // Keating strain parameter: bond length distortion
        double beta;   // Keating strain parameter: bond angle  distortion
        double *disp;  // displacement of each atomic coordinate
        int displen;   // size of disp (3*Natoms)
};

/************************************************************************************************/
class MacoptStrain : public VFFStrain, public Macopt   // Use the Macopt CG solver for VFF
{
    public:
        MacoptStrain(WireGenerator *, Material *);
        ~MacoptStrain();
        void relax_atoms(void);
        double func(double *);              // evaluate function
        void dfunc(double *, double *);     // evaluate gradient
};

#endif  // __ATOMSTRAIN
