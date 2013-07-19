/*!\class RGF RGF.H
\brief This class solves the ballistic NEGF equation (E-H-Sigma)*GR=I using the recursive Green's Function (RGF) algorithm described in J. Appl. Phys. 81, 7845 (1997). There is no parallelization. This class is derived from GFSolver.H.
*/

#ifndef _RGF
#define _RGF

#include "CSR.H"
#include "CSC.H"
#include "Types.H"
#include "Utilities.H"
#include "GFSolver.H"

class RGF : public GFSolver{
	
public:

    RGF(TCSR<CPX>* mat,int,MPI_Comm);

    virtual ~RGF();
    virtual void solve_equation(CPX*,CPX*,CPX*,CPX*,CPX*,CPX*,int*,int*,int,int,\
				int,int,int,int,int*,int*,int);
				
private:

    ///matrix E-H
    TCSR<CPX>* matrix;
    int findx,NPCS,TB;

    void extract_diag(CPX*,int*,int*,int*,int,int*,int);
    void extract_not_diag(TCSR<CPX>*,TCSC<CPX,int>*,int*,int*,int,int*,int);
    void extract_not_diag(TCSC<CPX,int>*,TCSR<CPX>*,int*,int*,int,int*,int);
    void init_green(CPX*,CPX*,int*,int,int*,int*,int,int,int*);
    void calc_gR(CPX*,CPX*,CPX*,int*,int*,int,int*,int,int*,int);
    void calc_SigR(CPX*,CPX*,CPX*,CPX*,TCSR<CPX>*,TCSC<CPX,int>*,int*,int*,int,int*,\
		   int,char*);
    void calc_GR(CPX*,CPX*,CPX*,CPX*,CPX*,TCSR<CPX>*,TCSC<CPX,int>*,int*,int*,int,
		 int*,int,int,char*);
    void write_matrix(const char*,CPX*,int,int);
    void transpose(CPX*,CPX*,int,int);
    void get_Gamma(CPX*,CPX*,int);
							
};

#endif
