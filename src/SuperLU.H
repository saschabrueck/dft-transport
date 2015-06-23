#ifndef _SUPERLU
#define _SUPERLU

#include "CSR.H"
#include "Types.H"
#include "LinearSolver.H"
#include "superlu.h"

template <class T>
class SuperLU : public LinearSolver<T>{
	
public:

    SuperLU(TCSR<T>*,MPI_Comm);

    virtual ~SuperLU();
    virtual void prepare();
    virtual void prepare(int*,int*,int,int,int*,int);
    virtual void prepare_corner(CPX*,CPX*,int*,int*,int,int,int*,int);
    virtual void solve_equation(T* res, T* rhs, int no_rhs);    
				
private:

    int mpi_size,mpi_rank,fortran_index;
    int m_loc,n;
    MPI_Comm slu_comm;
    SuperMatrix A;
    superlu_options_t options;
    gridinfo_t grid;
    SuperLUStat_t stat;
    ScalePermstruct_t ScalePermstruct;
    LUstruct_t LUstruct;
    SOLVEstruct_t SOLVEstruct;

    void create_local_matrix(TCSR<T>*,SuperMatrix*);
};

/************************************************************************************************/

template <class T>
SuperLU<T>::SuperLU(TCSR<T>* mat,MPI_Comm solver_comm)
{
    slu_comm                  = solver_comm;
    
    MPI_Comm_size(slu_comm,&mpi_size);
    MPI_Comm_rank(slu_comm,&mpi_rank);

    set_default_options_dist(&options);
    options.IterRefine        = NOREFINE;
    options.PrintStat         = NO;
    
    fortran_index             = mat->findx;
    n                         = mat->size_tot;
    
    superlu_gridinit(slu_comm,mpi_size,1,&grid);
    create_local_matrix(mat,&A);
}

/************************************************************************************************/

template <>
SuperLU<CPX>::~SuperLU()
{
    
    PStatFree(&stat);
    Destroy_CompRowLoc_Matrix_dist(&A);
    ScalePermstructFree(&ScalePermstruct);
    Destroy_LU(n,&grid, &LUstruct);
    LUstructFree(&LUstruct);
    zSolveFinalize(&options, &SOLVEstruct);

    //superlu_gridexit(&grid);
}

/************************************************************************************************/
/*
template <>
SuperLU<double>::~SuperLU()
{
    
    PStatFree(&stat);
    Destroy_CompRowLoc_Matrix_dist(&A);
    ScalePermstructFree(&ScalePermstruct);
    Destroy_LU(n,&grid, &LUstruct);
    LUstructFree(&LUstruct);
    dSolveFinalize(&options, &SOLVEstruct);

    //superlu_gridexit(&grid);
}
*/
/************************************************************************************************/

template <>
void SuperLU<CPX>::prepare()
{
    int info;

    ScalePermstructInit(A.nrow,A.ncol,&ScalePermstruct);
    LUstructInit(A.nrow,A.ncol,&LUstruct);

    PStatInit(&stat);
    
    pzgssvx(&options,&A,&ScalePermstruct,NULL,m_loc,0,&grid,&LUstruct,&SOLVEstruct,\
            NULL,&stat,&info);

    options.Fact = FACTORED;
}

/************************************************************************************************/
/*
template <>
void SuperLU<double>::prepare()
{
    int info;

    ScalePermstructInit(A.nrow,A.ncol,&ScalePermstruct);
    LUstructInit(A.nrow,A.ncol,&LUstruct);

    PStatInit(&stat);
    
    pdgssvx(&options,&A,&ScalePermstruct,NULL,m_loc,0,&grid,&LUstruct,&SOLVEstruct,\
            NULL,&stat,&info);

    options.Fact = FACTORED;
}
*/
/************************************************************************************************/

template <class T>
void SuperLU<T>::prepare(int *Bmin,int *Bmax,int NBlock,int Bsize,int *orb_per_at,int tb)
{
    prepare();
}

/************************************************************************************************/

template <class T>
void SuperLU<T>::prepare_corner(CPX *ML_ref,CPX *MR_ref,int *Bmin,int *Bmax,int NBlock,int Bsize,\
				int *orb_per_at,int tb)
{
}

/************************************************************************************************/

template <>
void SuperLU<CPX>::solve_equation(CPX *res, CPX *arg_rhs, int no_rhs)
{
    int IRHS,IM,info;
    doublecomplex *rhs_local = new doublecomplex[no_rhs*2*m_loc];
    double *berr             = doubleMalloc_dist(no_rhs);
    
    for(IRHS=0;IRHS<no_rhs;IRHS++){
        for(IM=0;IM<m_loc;IM++){
            rhs_local[IRHS*m_loc+IM].r = real(arg_rhs[IRHS*m_loc+IM]);
            rhs_local[IRHS*m_loc+IM].i = imag(arg_rhs[IRHS*m_loc+IM]);
        }
    }
    
    pzgssvx(&options,&A,&ScalePermstruct,rhs_local,m_loc,no_rhs,&grid,&LUstruct,&SOLVEstruct,\
            berr,&stat,&info);
    
    PStatPrint(&options,&stat,&grid);

    for(IRHS=0;IRHS<no_rhs;IRHS++){
        for(IM=0;IM<m_loc;IM++){
            res[IRHS*m_loc+IM] = CPX(rhs_local[IRHS*m_loc+IM].r,\
				     rhs_local[IRHS*m_loc+IM].i);
        }
    }

    delete[] rhs_local;
    delete[] berr;
}

/************************************************************************************************/
/*
template <>
void SuperLU<double>::solve_equation(double *res, double *arg_rhs, int no_rhs)
{
    MPI_Status status;
    int IRHS,IM,info;
    double *rhs_local        = new double[no_rhs*2*m_loc];
    double *berr             = doubleMalloc_dist(no_rhs);
    
    for(IRHS=0;IRHS<no_rhs;IRHS++){
        for(IM=0;IM<m_loc;IM++){
            rhs_local[IRHS*m_loc+IM] = arg_rhs[IRHS*n+mpi_rank*m_loc_fst+IM];
        }
    }
    
    pdgssvx(&options,&A,&ScalePermstruct,rhs_local,m_loc,no_rhs,&grid,&LUstruct,&SOLVEstruct,\
            berr,&stat,&info);

    PStatPrint(&options,&stat,&grid);
    
    for(IM=0;IM<no_element;IM++){
        for(IRHS=0;IRHS<no_rhs;IRHS++){
            res[IRHS*m_loc+IM] = rhs_local[IRHS*m_loc+IM];
        }
    }

    delete[] rhs_local;
    delete[] berr;

}
*/
/************************************************************************************************/

template <>
void SuperLU<CPX>::create_local_matrix(TCSR<CPX> *matrix,SuperMatrix *SMat)
{
    int_t nnz_loc,fst_row,IZ,IR;
    doublecomplex *nnz_val_loc;
    int_t *rowptr,*colind;
    
    m_loc       = matrix->size;

    nnz_loc     = matrix->n_nonzeros;
    fst_row     = matrix->first_row;
    rowptr      = (int_t*)intMalloc_dist(m_loc+1);
    colind      = (int_t*)intMalloc_dist(nnz_loc);
    nnz_val_loc = (doublecomplex*)doublecomplexMalloc_dist(nnz_loc);
    
    for(IR=0;IR<=m_loc;IR++){
        rowptr[IR]        = matrix->edge_i[IR]-matrix->findx;
    }

    for(IZ=0;IZ<nnz_loc;IZ++){
        colind[IZ]        = matrix->index_j[IZ]-matrix->findx;
        nnz_val_loc[IZ].r = real(matrix->nnz[IZ]);
        nnz_val_loc[IZ].i = imag(matrix->nnz[IZ]);
    }
    
    zCreate_CompRowLoc_Matrix_dist(SMat,n,n,nnz_loc,m_loc,fst_row,nnz_val_loc,colind,rowptr,\
                                   SLU_NR_loc,SLU_Z,SLU_GE);
}

/************************************************************************************************/
/*
template <>
void SuperLU<CPX>::create_local_matrix(TCSR<CPX> *matrix,SuperMatrix *SMat)
{
    int_t nnz_loc,fst_row,ind_loc,IZ,IR;
    doublecomplex *nnz_val_loc;
    int_t *rowptr,*colind;
    
    m_loc     = matrix->size/mpi_size; 
    m_loc_fst = m_loc;
    
    if ((m_loc * mpi_size) != matrix->size) {
        m_loc     = m_loc+1;
        m_loc_fst = m_loc;
        if (mpi_rank == (mpi_size - 1))
            m_loc = matrix->size - m_loc_fst * (mpi_size - 1);
    }

    nnz_loc     = matrix->edge_i[mpi_rank*m_loc_fst+m_loc]-matrix->edge_i[mpi_rank*m_loc_fst];
    fst_row     = mpi_rank* m_loc_fst;
    ind_loc     = matrix->edge_i[mpi_rank*m_loc_fst]-matrix->findx;
    rowptr      = (int_t*)intMalloc_dist(m_loc+1);
    colind      = (int_t*)intMalloc_dist(nnz_loc);
    nnz_val_loc = (doublecomplex*)doublecomplexMalloc_dist(nnz_loc);
    
    for(IR=fst_row;IR<=fst_row+m_loc;IR++){
        rowptr[IR-fst_row]        = matrix->edge_i[IR]-matrix->edge_i[fst_row];
    }

    for(IZ=ind_loc;IZ<ind_loc+nnz_loc;IZ++){
        colind[IZ-ind_loc]        = matrix->index_j[IZ]-matrix->findx;
        nnz_val_loc[IZ-ind_loc].r = real(matrix->nnz[IZ]);
        nnz_val_loc[IZ-ind_loc].i = imag(matrix->nnz[IZ]);
    }
    
    zCreate_CompRowLoc_Matrix_dist(SMat,n,n,nnz_loc,m_loc,fst_row,nnz_val_loc,colind,rowptr,\
                                   SLU_NR_loc,SLU_Z,SLU_GE);
}
*/
/************************************************************************************************/
/*
template <>
void SuperLU<double>::create_local_matrix(TCSR<double> *matrix,SuperMatrix *SMat)
{
    int_t nnz_loc,fst_row,ind_loc,IZ,IR;
    double *nnz_val_loc;
    int_t *rowptr,*colind;
    
    m_loc     = matrix->size/mpi_size; 
    m_loc_fst = m_loc;
    
    if ((m_loc * mpi_size) != matrix->size) {
        m_loc     = m_loc+1;
        m_loc_fst = m_loc;
        if (mpi_rank == (mpi_size - 1))
            m_loc = matrix->size - m_loc_fst * (mpi_size - 1);
    }

    nnz_loc     = matrix->edge_i[mpi_rank*m_loc_fst+m_loc]-matrix->edge_i[mpi_rank*m_loc_fst];
    fst_row     = mpi_rank* m_loc_fst;
    ind_loc     = matrix->edge_i[mpi_rank*m_loc_fst]-matrix->findx;
    rowptr      = (int_t*)intMalloc_dist(m_loc+1);
    colind      = (int_t*)intMalloc_dist(nnz_loc);
    nnz_val_loc = (double*)doubleMalloc_dist(nnz_loc);
    
    for(IR=fst_row;IR<=fst_row+m_loc;IR++){
        rowptr[IR-fst_row]        = matrix->edge_i[IR]-matrix->edge_i[fst_row];
    }

    for(IZ=ind_loc;IZ<ind_loc+nnz_loc;IZ++){
        colind[IZ-ind_loc]        = matrix->index_j[IZ]-matrix->findx;
        nnz_val_loc[IZ-ind_loc]   = matrix->nnz[IZ];
    }
    
    dCreate_CompRowLoc_Matrix_dist(SMat,n,n,nnz_loc,m_loc,fst_row,nnz_val_loc,colind,rowptr,\
                                   SLU_NR_loc,SLU_D,SLU_GE);
}
*/
/************************************************************************************************/

#endif
