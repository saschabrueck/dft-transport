/*************************************************************************************************
To get more memory, replace zi and di by zl and dl, int *Wi by UF_long *Wi
and TCSC<T,int> *matrix by TCSC<T,UF_long> *matrix
Problem: results may be unstable
*************************************************************************************************/

#ifndef _UMFPACK
#define _UMFPACK

#include "CSR.H"
#include "CSC.H"
#include "Types.H"
#include "LinearSolver.H"
#include "umfpack.h"
#include "amd.h"

template <class T>
class Umfpack : public LinearSolver<T>{
	
public:

    Umfpack(TCSR<T>* mat,MPI_Comm);

    virtual ~Umfpack();
    virtual void prepare();
    virtual void prepare(int*,int*,int,int,int*,int);
    virtual void prepare_corner(CPX*,CPX*,int*,int*,int*,int,int,int*,int);
    virtual void solve_equation(T* res, T* rhs, int no_rhs);    
				
private:

    TCSC<T,UF_long>* matrix;

    int n;
    UF_long *Wi;
    double *W;
    bool factorized;
    double Info [UMFPACK_INFO];
    double Control [UMFPACK_CONTROL];
    void *Symbolic, *Numeric ;		
};

/************************************************************************************************/

template <class T>
Umfpack<T>::Umfpack(TCSR<T>* mat,MPI_Comm solver_comm)
{

    Control[UMFPACK_PRL]                 = 1;       //Printing level
    Control[UMFPACK_DENSE_ROW]           = 0.2;
    Control[UMFPACK_DENSE_COL]           = 0.2;
    Control[UMFPACK_PIVOT_TOLERANCE]     = 0.1;
    Control[UMFPACK_BLOCK_SIZE]          = 32;
    Control[UMFPACK_STRATEGY]            = 0;       //0 for automatic choice
    Control[UMFPACK_ALLOC_INIT]          = 0.7;
    Control[UMFPACK_IRSTEP]              = 0;       //No iterative refinement
    Control[UMFPACK_2BY2_TOLERANCE]      = 0.01;
    Control[UMFPACK_FIXQ]                = 0;
    Control[UMFPACK_AMD_DENSE]           = 10;
    Control[UMFPACK_SYM_PIVOT_TOLERANCE] = 0.001;
    Control[UMFPACK_SCALE]               = 1;
    Control[UMFPACK_FRONT_ALLOC_INIT]    = 0.5;
    Control[UMFPACK_DROPTOL]             = 0;
    Control[UMFPACK_AGGRESSIVE]          = 1;

    matrix                               = new TCSC<T,UF_long>(mat,0);
    
    n                                    = matrix->size;

    Wi                                   = new UF_long[n];
    W                                    = new double[4*n];
    
    factorized                           = false;
}

/************************************************************************************************/

template <>
inline Umfpack<CPX>::~Umfpack()
{

    umfpack_zl_free_symbolic (&Symbolic);
    umfpack_zl_free_numeric(&Numeric);
    
    delete[] Wi;
    delete[] W;

    delete matrix;
}

/************************************************************************************************/

template <>
inline Umfpack<double>::~Umfpack()
{

    umfpack_dl_free_symbolic (&Symbolic);
    umfpack_dl_free_numeric(&Numeric);
    
    delete[] Wi;
    delete[] W;
    delete matrix;
}

/************************************************************************************************/

template <>
inline void Umfpack<CPX>::prepare()
{
    int status;
    
    //symbolic factorization
    status=umfpack_zl_symbolic(n,n,matrix->edge_j,matrix->index_i,\
                               (double*)matrix->nnz,NULL,&Symbolic,Control,Info);

    if (status < 0)
    {
	umfpack_zl_report_info (Control, Info) ;
	umfpack_zl_report_status (Control, status) ;
	printf("umfpack_zl_symbolic failed\n");
    }
    
    // numeric factorization
    status=umfpack_zl_numeric(matrix->edge_j,matrix->index_i,(double*)matrix->nnz,NULL,\
			      Symbolic,&Numeric,Control,Info);
    
    if (status < 0)
    {
	umfpack_zl_report_info (Control, Info);
	umfpack_zl_report_status (Control, status);
	printf("umfpack_zl_numeric failed\n");
    }

    factorized = true;
}

/************************************************************************************************/

template <>
inline void Umfpack<double>::prepare()
{
    int status;
    
    //symbolic factorization
    status=umfpack_dl_symbolic(n,n,matrix->edge_j,matrix->index_i,\
                               matrix->nnz,&Symbolic,Control,Info);

    if (status < 0)
    {
	umfpack_dl_report_info (Control, Info);
	umfpack_dl_report_status (Control, status);
	printf("umfpack_dl_symbolic failed\n");
    }

    // numeric factorization
    status=umfpack_dl_numeric(matrix->edge_j,matrix->index_i,matrix->nnz,Symbolic,\
                              &Numeric,Control,Info);

    if (status < 0)
    {
	umfpack_dl_report_info (Control, Info);
	umfpack_dl_report_status (Control, status);
	printf("umfpack_dl_numeric failed\n");
    }

    factorized = true;
}

/************************************************************************************************/

template <class T>
void Umfpack<T>::prepare(int *Bmin,int *Bmax,int NBlock,int Bsize,int *orb_per_at,int tb)
{
    prepare();
}

/************************************************************************************************/

template <class T>
void Umfpack<T>::prepare_corner(CPX *ML_ref,CPX *MR_ref,int *NBC,int *Bmin,int *Bmax,int NBlock,\
				int Bsize,int *orb_per_at,int tb)
{
}

/************************************************************************************************/

template <>
inline void Umfpack<CPX>::solve_equation(CPX* res, CPX* arg_rhs, int no_rhs)
{   
    // solve Ax=b
    int IRHS,status;

    if(factorized==false){
        prepare();
    }
    
    for(IRHS=0;IRHS<no_rhs;IRHS++){
        status=umfpack_zl_wsolve(UMFPACK_A,matrix->edge_j,matrix->index_i,(double*)matrix->nnz,\
                                 NULL,(double*)&res[IRHS*n],NULL,(double*)&arg_rhs[IRHS*n],NULL,\
                                 Numeric,Control,Info,Wi,W);

        if (status < 0)
        {
            umfpack_zl_report_info (Control, Info);
            umfpack_zl_report_status (Control, status);
            printf("umfpack_zl_wsolve failed\n");
        }
        
    }

}

/************************************************************************************************/

template <>
inline void Umfpack<double>::solve_equation(double* res, double* arg_rhs, int no_rhs)
{
    
    // solve Ax=b
    int IRHS,status;

    if(factorized==false){
        prepare();
    }
    
    for(IRHS=0;IRHS<no_rhs;IRHS++){
        status=umfpack_dl_wsolve(UMFPACK_A,matrix->edge_j,matrix->index_i,matrix->nnz,\
                                 &res[IRHS*n],&arg_rhs[IRHS*n],Numeric,Control,Info,Wi,W);

        if (status < 0)
        {
            umfpack_dl_report_info (Control, Info);
            umfpack_dl_report_status (Control, status);
            printf("umfpack_dl_wsolve failed\n");
        }
        
    }
}

/************************************************************************************************/

#endif
