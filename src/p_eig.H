#ifndef __P_EIG
#define __P_EIG

#include "ScaLapack.H"

#define LOGCERR std::cerr << "ERROR IN LINE " << __LINE__ << " OF FILE " << __FILE__ << std::endl

int p_grid_desc_init(int&,int,int,int&,int&,int*,int*);
int p_eig(double*,double*,double*,int,MPI_Comm);
int p_eig(CPX*,CPX*,double*,int,MPI_Comm);
int p_inv(CPX*,int,MPI_Comm);
int p_lin(CPX*,CPX*,CPX*,int,int,MPI_Comm);

#endif

