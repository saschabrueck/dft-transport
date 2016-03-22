#ifndef _BANDED
#define _BANDED

#include "CSR.H"
#include "Types.H"
#include "LinearSolver.H"

template <class T>
class Banded : public LinearSolver<T>{
	
public:

    Banded(TCSR<T>*,MPI_Comm);
    Banded(TCSR<T>*,MPI_Comm,int);

    virtual ~Banded();
    virtual void prepare();
    virtual void prepare(int*,int*,int,int,int*,int);
    virtual void prepare_corner(CPX*,CPX*,int*,int*,int*,int,int,int*,int);
    virtual void solve_equation(T* res, T* rhs, int no_rhs);    
				
private:

    int icontxt;
    int pp;
    int bw;
    int NB;
    int size_tot;
    int size_csr_loc;
    int LAF;
    int descA[7];
    int *ipiv;
    CPX* A;
    CPX* AF;

};

/************************************************************************************************/

template <>
Banded<CPX>::Banded(TCSR<CPX>* mat,MPI_Comm solver_comm)
{
    pp = 1;
    int mpi_size;
    MPI_Comm_size(solver_comm,&mpi_size);

    icontxt = MPI_Comm_c2f(solver_comm);
    size_tot = mat->size_tot;
    size_csr_loc = mat->size;
    NB = int(ceil(double(size_tot)/double(mpi_size)));
    char gridr[1] = {'R'};
    Cblacs_gridinit(&icontxt,gridr,1,mpi_size);
    bw = mat->get_bandwidth(solver_comm);
    A = new CPX[(2*bw*pp+2*bw+1)*size_csr_loc];
    LAF = (NB+pp*bw)*2*bw+6*(bw+pp*bw)*(bw+pp*2*bw);
    AF = new CPX[LAF];
    ipiv = NULL;
    if (pp) ipiv = new int[NB+2*bw];
    descA[0]=501;
    descA[1]=icontxt;
    descA[2]=size_tot;
    descA[3]=NB;
    descA[4]=0;
    descA[5]=2*pp*bw+2*bw+1;
    mat->sparse_to_narrow_band(A,bw,pp);
}

/************************************************************************************************/

template <>
Banded<CPX>::~Banded()
{
    Cblacs_gridexit(icontxt);
    delete[] A;
    delete[] AF;
    if (pp) delete[] ipiv;
    
}

/************************************************************************************************/

template <>
void Banded<CPX>::prepare()
{
    int iinfo;
    if (pp) {
        int lwork = 1;
        CPX workq;
        c_pzgbtrf(size_tot,bw,bw,A,1,descA,ipiv,AF,LAF,&workq,lwork,&iinfo);
    } else {
        int lwork = -1;
        CPX workq;
        CPX *work = &workq;

        c_pzdbtrf(size_tot,bw,bw,A,1,descA,AF,LAF,work,lwork,&iinfo);
 
        lwork = int(real(workq));
        work = new CPX[lwork];
 
        c_pzdbtrf(size_tot,bw,bw,A,1,descA,AF,LAF,work,lwork,&iinfo);
 
        delete[] work;
    }
}

/************************************************************************************************/

template <class T>
void Banded<T>::prepare(int *Bmin,int *Bmax,int NBlock,int Bsize,int *orb_per_at,int tb)
{   
    prepare();
}

/************************************************************************************************/

template <class T>
void Banded<T>::prepare_corner(CPX *ML_ref,CPX *MR_ref,int *NBC,int *Bmin,int *Bmax,int NBlock,\
			      int Bsize,int *orb_per_at,int tb)
{
}

/************************************************************************************************/

template <>
void Banded<CPX>::solve_equation(CPX *res, CPX *arg_rhs, int no_rhs)
{
    if (no_rhs) {
        int iinfo;
        int descRHS[7];
        descRHS[0]=502;
        descRHS[1]=icontxt;
        descRHS[2]=size_tot;
        descRHS[3]=NB;
        descRHS[4]=0;
        descRHS[5]=size_csr_loc;

        c_zcopy(size_csr_loc*no_rhs,arg_rhs,1,res,1);
 
        int lwork = -1;
        CPX workq;
        CPX *work = &workq;

        if (pp) {
            c_pzgbtrs('N',size_tot,bw,bw,no_rhs,A,1,descA,ipiv,res,1,descRHS,AF,LAF,work,lwork,&iinfo);
        } else {
            c_pzdbtrs('N',size_tot,bw,bw,no_rhs,A,1,descA,res,1,descRHS,AF,LAF,work,lwork,&iinfo);
        }
 
        lwork = (2-pp)*int(real(workq));
        work = new CPX[lwork];
 
        if (pp) {
            c_pzgbtrs('N',size_tot,bw,bw,no_rhs,A,1,descA,ipiv,res,1,descRHS,AF,LAF,work,lwork,&iinfo);
        } else {
            c_pzdbtrs('N',size_tot,bw,bw,no_rhs,A,1,descA,res,1,descRHS,AF,LAF,work,lwork,&iinfo);
        }
 
        delete[] work;
    }
}

/************************************************************************************************/

#endif
