#ifndef __INJECTION
#define __INJECTION

#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Types.H"
#include "Utilities.H"
#include "Blas.H"
#include "ScaLapack.H"

template<class T>
class Injection{
public:
    
    Injection(){};
    Injection(int,double){};

    virtual ~Injection(){};

    virtual void initialize(int,int,int) = 0;
    virtual void calc_kphase(TCSR<T>*,TCSR<T>*,TCSR<CPX>*,int,int,CPX*,CPX*,double*,\
			     int*,int*,int*,int*,int,int,int,MPI_Comm,int*) = 0;
    virtual void calc_kphase_sym(TCSR<T>*,TCSR<T>*,TCSR<CPX>*,int,int,CPX*,CPX*,CPX*,\
				 CPX*,double*,double*,int*,int*,int*,int*,int*,int*,\
				 int*,int*,int,int,MPI_Comm,int*) = 0;
    virtual void calc_kphase(TCSR<T>*,T*,TCSR<CPX>*,int,int,CPX*,CPX*,double*,\
			     int*,int*,int*,int*,int,int,int,MPI_Comm,int*) = 0;
    virtual void calc_kphase_sym(TCSR<T>*,T*,TCSR<CPX>*,int,int,CPX*,CPX*,CPX*,\
				 CPX*,double*,double*,int*,int*,int*,int*,int*,int*,\
				 int*,int*,int,int,MPI_Comm,int*) = 0;
};

#endif
