#ifndef __FERMI
#define __FERMI

#include <mpi.h>
#include "Utilities.H"
#include "Types.H"

class Fermi
{
public:
    
    Fermi();
    ~Fermi();

    double find_fermi(double,double,double*,double,int,double*,int,double*,int,int,double,\
		      double,double);
    void density(double*,double,double*,double,int,double*,int,double*,int,int,double,double,\
		 int,double,double,double);
    void derivate(double*,double,double*,double,int,double*,int,double*,int,int,double,double,\
		  int,double,double,double);

private:
    double Eps0,e,kB;

    double Newton(double,double,double,double*,double,int,double*,int,double*,int,int,double,\
		  double,double);
    double get_edge(double*,double,int,int,int);
};

#endif
