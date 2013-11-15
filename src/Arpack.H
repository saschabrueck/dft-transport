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
    bool eigs(CPX*,CPX*,CPX*,int,int);

private:
	
};
    
#endif
