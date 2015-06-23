#include "c_pexsi_interface.h"
#include "CSR.H"

int densitypbc(TCSR<double> *KohnSham,TCSR<double> *Overlap,TCSR<double> *Ps,CPX energy,CPX weight,double dos,MPI_Comm matrix_comm)
{
    int matrix_procs,matrix_rank;
    MPI_Comm_size(matrix_comm,&matrix_procs);
    MPI_Comm_rank(matrix_comm,&matrix_rank);

    if (KohnSham->findx!=1 || Overlap->findx!=1) return (LOGCERR, EXIT_FAILURE);

    double *HS_nnz_inp = new double[2*Overlap->n_nonzeros]();
    double *HS_nnz_out = new double[2*Overlap->n_nonzeros]();

    c_dcopy(Overlap->n_nonzeros,Overlap->nnz,1,HS_nnz_inp,2);
    c_zscal(Overlap->n_nonzeros,-energy,(CPX*)HS_nnz_inp,1);
    c_daxpy(Overlap->n_nonzeros,1.0,KohnSham->nnz,1,HS_nnz_inp,2);

    int n_nonzeros_global;
    MPI_Allreduce(&Overlap->n_nonzeros,&n_nonzeros_global,1,MPI_INT,MPI_SUM,matrix_comm);
    int info;

    PPEXSIOptions  options;
    PPEXSISetDefaultOptions(&options);
    options.npSymbFact = matrix_procs;
    options.ordering = 0;
    options.verbosity = 0;
 
    PPEXSIPlan   plan;
    plan = PPEXSIPlanInitialize(matrix_comm,1,matrix_procs,-1,&info);
    if (info) return (LOGCERR, EXIT_FAILURE);
    PPEXSILoadRealSymmetricHSMatrix(plan,options,Overlap->size_tot,n_nonzeros_global,Overlap->n_nonzeros,Overlap->size,Overlap->edge_i,Overlap->index_j,HS_nnz_inp,1,NULL,&info);
    if (info) return (LOGCERR, EXIT_FAILURE);
    PPEXSISymbolicFactorizeComplexSymmetricMatrix(plan,options,&info);
    if (info) return (LOGCERR, EXIT_FAILURE);
    PPEXSISelInvComplexSymmetricMatrix(plan,options,HS_nnz_inp,HS_nnz_out,&info);
    if (info) return (LOGCERR, EXIT_FAILURE);
    PPEXSIPlanFinalize(plan,&info);
    if (info) return (LOGCERR, EXIT_FAILURE);

    delete[] HS_nnz_inp;
    c_zscal(Overlap->n_nonzeros,-weight/M_PI*CPX(0.0,1.0),(CPX*)HS_nnz_out,1);
    dos=c_ddot(Overlap->n_nonzeros,Overlap->nnz,1,HS_nnz_out,2);
    c_daxpy(Overlap->n_nonzeros,1.0,HS_nnz_out,2,Ps->nnz,1);
    delete[] HS_nnz_out;
    return 0;
}
