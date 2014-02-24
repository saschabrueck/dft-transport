#ifndef _ARPACK
#define _ARPACK

#include <mpi.h>
#include <iostream>
#include "CSR.H"
#include "Blas.H"
#include "Utilities.H"

class Arpack{
public:

    Arpack();
    ~Arpack();
    bool eigs(CPX*,CPX*,CPX*,int,int,int,double&);
    bool eigs(CPX*,CPX*,double*,int,int,int,double&);

private:
	
};
    
#endif
