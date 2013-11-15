#ifndef _MUMPS
#define _MUMPS

#include "CSR.H"
#include "Types.H"
#include "LinearSolver.H"
#include "zmumps_c.h"

template <class T>
class MUMPS : public LinearSolver<T>{
	
public:

    MUMPS(TCSR<T>*,MPI_Comm);
    MUMPS(TCSR<T>*,MPI_Comm,int);

    virtual ~MUMPS();
    virtual void prepare();
    virtual void prepare(int*,int*,int,int,int*,int);
    virtual void solve_equation(T* res, T* rhs, int no_rhs);    
				
private:

    int mpi_size,mpi_rank,fortran_index,m_loc,fst_row;
    MPI_Comm slu_comm;
    ZMUMPS_STRUC_C id;

    int do_distribute;

    void create_local_matrix(TCSR<T>*,ZMUMPS_STRUC_C*);
};


extern "C"{
    void MUMPS_CALL zmumps_c(ZMUMPS_STRUC_C * zmumps_par);
}

/************************************************************************************************/
/*
template <class T>
MUMPS<T>::MUMPS(TCSR<T>* mat,MPI_Comm solver_comm)
{
    slu_comm                  = solver_comm;

    MPI_Comm_size(slu_comm,&mpi_size);
    MPI_Comm_rank(slu_comm,&mpi_rank);

    fortran_index             = 1;

    id.par                    = 1;
    id.sym                    = 0;
    id.comm_fortran           = slu_comm;

    id.job                    = -1;
    
    zmumps_c(&id);

    create_local_matrix(mat,&id);

}
*/
/************************************************************************************************/

template <class T>
MUMPS<T>::MUMPS(TCSR<T>* mat,MPI_Comm solver_comm,int pord=0)
{

    slu_comm                  = solver_comm;

    MPI_Comm_size(slu_comm,&mpi_size);
    MPI_Comm_rank(slu_comm,&mpi_rank);

    fortran_index             = 1;

    do_distribute             = 1;

    id.par                    = 1;
    id.sym                    = 0;
    id.comm_fortran           = slu_comm;

    id.job                    = -1;
    
    zmumps_c(&id);

    create_local_matrix(mat,&id);

    if(pord%2){
        id.icntl[6] = 4; //PORD reordering algorithm
    }

    if (pord/2){
        do_distribute=0;
    }

}

/************************************************************************************************/

template <>
MUMPS<CPX>::~MUMPS()
{
    id.job = -2;
    zmumps_c(&id);

    delete[] id.irn_loc;
    delete[] id.jcn_loc;
    delete[] id.a_loc;
    delete[] id.rhs;

}

/************************************************************************************************/

template <>
void MUMPS<CPX>::prepare()
{
    id.icntl[0]               = -1;
    id.icntl[1]               = -1;
    id.icntl[2]               = -1;
    id.icntl[3]               = -1;
    //id.icntl[6]               = 4; //PORD reordering algorithm
    id.icntl[17]              = 3;
    id.icntl[20]              = 0;
    
    id.job                    = 4;

    zmumps_c(&id);
}

/************************************************************************************************/

template <class T>
void MUMPS<T>::prepare(int *Bmin,int *Bmax,int NBlock,int Bsize,int *orb_per_at,int tb)
{   
    prepare();
}

/************************************************************************************************/

template <>
void MUMPS<CPX>::solve_equation(CPX *res, CPX *arg_rhs, int no_rhs)
{
    
    MPI_Status status;
    int IRHS,IM,i_proc,no_element,Itot;

    id.nrhs     = no_rhs;
    id.lrhs     = id.n;
    id.rhs      = new F_DOUBLE[no_rhs*2*id.n];
    //id.rhs      = new ZMUMPS_COMPLEX[no_rhs*2*id.n];
    id.job      = 3;
    
    if (do_distribute) {
        if(mpi_rank){
            MPI_Send(&m_loc,1,MPI_INT,0,0,slu_comm);
            MPI_Send(arg_rhs,no_rhs*m_loc,MPI_DOUBLE_COMPLEX,0,1,slu_comm);
        }else{
            Itot = 0;
            for(IM=0;IM<m_loc;IM++){
                for(IRHS=0;IRHS<no_rhs;IRHS++){
                    id.rhs[IRHS*id.n+Itot].r=real(arg_rhs[IRHS*m_loc+IM]);
                    id.rhs[IRHS*id.n+Itot].i=imag(arg_rhs[IRHS*m_loc+IM]);
                }
                Itot++;
            }
            for(i_proc=1;i_proc<mpi_size;i_proc++){
                MPI_Recv(&no_element,1,MPI_INT,i_proc,0,slu_comm,&status);
                MPI_Recv(arg_rhs,no_rhs*no_element,MPI_DOUBLE_COMPLEX,i_proc,1,slu_comm,&status);
                for(IM=0;IM<no_element;IM++){
                    for(IRHS=0;IRHS<no_rhs;IRHS++){
                        id.rhs[IRHS*id.n+Itot].r=real(arg_rhs[IRHS*no_element+IM]);
                        id.rhs[IRHS*id.n+Itot].i=imag(arg_rhs[IRHS*no_element+IM]);
                    }
                    Itot++;
                }
            }
        }
    }else{
        if(!mpi_rank){
            for(int IELE=0;IELE<no_rhs*id.n;IELE++){
                id.rhs[IELE].r=real(arg_rhs[IELE]);
                id.rhs[IELE].i=imag(arg_rhs[IELE]);
            }
        }
    }
    
    zmumps_c(&id);   

    MPI_Bcast(id.rhs,no_rhs*id.n,MPI_DOUBLE_COMPLEX,0,slu_comm);
 
    if (do_distribute) {
        for(IM=fst_row;IM<fst_row+m_loc;IM++){
            for(IRHS=0;IRHS<no_rhs;IRHS++){
                res[IRHS*m_loc+IM-fst_row] = \
                    CPX(id.rhs[IRHS*id.n+IM].r,id.rhs[IRHS*id.n+IM].i);
            }
        }
    }else{
        for(int IELE=0;IELE<no_rhs*id.n;IELE++){
            res[IELE] = CPX(id.rhs[IELE].r,id.rhs[IELE].i);
        }
    }
   
}

/************************************************************************************************/

template <>
void MUMPS<CPX>::create_local_matrix(TCSR<CPX> *matrix,ZMUMPS_STRUC_C *ZStruc)
{
    int IR,IC,index;

    ZStruc->n       = matrix->size_tot;
    ZStruc->nz_loc  = matrix->n_nonzeros;
    ZStruc->irn_loc = new int[matrix->n_nonzeros];
    ZStruc->jcn_loc = new int[matrix->n_nonzeros];
    ZStruc->a_loc   = new F_DOUBLE[matrix->n_nonzeros];
    //ZStruc->a_loc   = new ZMUMPS_COMPLEX[matrix->n_nonzeros];

    fst_row         = matrix->first_row;
    m_loc           = matrix->size;
        
    index           = 0;
    
    for(IR=0;IR<matrix->size;IR++){
        for(IC=matrix->edge_i[IR]-matrix->findx;IC<matrix->edge_i[IR+1]-matrix->findx;IC++){
            ZStruc->irn_loc[index] = fst_row+IR+fortran_index;
            ZStruc->jcn_loc[index] = matrix->index_j[index]-matrix->findx+fortran_index;
            ZStruc->a_loc[index].r = real(matrix->nnz[index]);
            ZStruc->a_loc[index].i = imag(matrix->nnz[index]);
            index++;
        }
    }

}

/************************************************************************************************/

#endif
