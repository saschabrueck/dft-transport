#ifndef _AZTEC
#define _AZTEC

#include "CSR.H"
#include "Types.H"
#include "LinearSolver.H"
#include "aztec.h"

template <class T>
class Aztec : public LinearSolver<T>{
	
public:

    Aztec(TCSR<T>* mat,MPI_Comm);

    virtual ~Aztec();
    virtual void prepare();
    virtual void prepare(int*,int*,int,int,int*,int);
    virtual void solve_equation(T* res, T* rhs, int no_rhs);
				
private:

    int      count_elements(int, int*, TCSR<T>*);
    void     create_matrix_row(int, int, TCSR<T>*, T*, int*);

    int      n,nnz,size,rank,first_row;
    int      proc_config[AZ_PROC_SIZE];/* Processor information.                */
    int      options[AZ_OPTIONS_SIZE]; /* Array used to select solver options.  */
    double   params[AZ_PARAMS_SIZE];   /* User selected solver paramters.       */
    int      *data_org;                /* Array to specify data layout          */
    double   status[AZ_STATUS_SIZE];   /* Information returned from AZ_solve(). */
    int      *update;                  /* vector elements updated on this node. */
    int	     *external;                /* vector elements needed by this node.  */
    int      *update_index;            /* ordering of update[] and external[]   */
    int      *extern_index;            /* locally on this processor.            */
    int      *bindx;                   /* Sparse matrix to be solved is stored  */
    T        *val;                     /* in these MSR arrays.                  */
    int      N_update;                 /* # of unknowns updated on this node    */
    MPI_Comm internal_comm;
							
};

/************************************************************************************************/

template <class T>
Aztec<T>::Aztec(TCSR<T>* mat,MPI_Comm solver_comm)
{   

    internal_comm = solver_comm;
    n             = mat->size_tot;
    first_row     = mat->first_row;
    
    MPI_Comm_size(internal_comm,&size);
    MPI_Comm_rank(internal_comm,&rank);
    
    AZ_set_proc_config(proc_config,internal_comm);
    AZ_read_update(&N_update,&update,proc_config,n,1,AZ_linear);
    
    nnz      = count_elements(N_update,update,mat);
    bindx    = new int[nnz+1];
    val      = new T[nnz+1];

    bindx[0] = N_update+1;
    
    for(int i=0;i<N_update;i++){
	create_matrix_row(update[i],i,mat,val,bindx);
    }
   
    AZ_transform(proc_config, &external, bindx, val, update, &update_index,
		 &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
		 AZ_MSR_MATRIX);
    
    AZ_defaults(options, params);
    options[AZ_solver]          = AZ_bicgstab;
    options[AZ_scaling]         = AZ_sym_row_sum;
    options[AZ_precond]         = AZ_dom_decomp;
    options[AZ_subdomain_solve] = AZ_ilut;
    options[AZ_ilut_fill]       = 2.0;
    options[AZ_max_iter]        = 1000;
    options[AZ_conv]            = AZ_rhs;
    options[AZ_output]          = AZ_none;
}

/************************************************************************************************/

template <class T>
Aztec<T>::~Aztec()
{

    delete[] val;
    delete[] bindx;
    delete[] update;
    delete[] external;
    delete[] update_index;
    delete[] extern_index;
    delete[] data_org;

}

/************************************************************************************************/

template <class T>
void Aztec<T>::prepare()
{
}

/************************************************************************************************/

template <class T>
void Aztec<T>::prepare(int *Bmin,int *Bmax,int NBlock,int Bsize,int *n_of_el,int tb)
{
}

/************************************************************************************************/

template <class T>
void Aztec<T>::solve_equation(T* res, T* arg_rhs, int no_rhs)
{
    MPI_Status mpi_status;
    T *b = new T[N_update];
    T *x = new T[N_update+data_org[AZ_N_external]];
   
    for(int i=0;i<N_update;i++){
	b[update_index[i]] = arg_rhs[update[i]-first_row];
	x[update_index[i]] = 0.0;
    }

    AZ_solve(x, b, options, params, NULL, bindx, NULL, NULL, NULL, val, data_org,
	     status, proc_config);

    if(!rank){
	for(int i=0;i<N_update;i++){
	    res[update[i]] = x[update_index[i]];
	}
	int *update_loc, *update_index_loc;
	T *x_loc;
	for(int i_proc=1;i_proc<size;i_proc++){
	    MPI_Recv(&N_update,1,MPI_INT,i_proc,0,internal_comm,&mpi_status);
	    update_loc       = new int[N_update];
	    update_index_loc = new int[N_update];
	    x_loc            = new T[N_update];

	    MPI_Recv(update_loc,N_update,MPI_INT,i_proc,1,internal_comm,&mpi_status);
	    MPI_Recv(update_index_loc,N_update,MPI_INT,i_proc,2,internal_comm,&mpi_status);
	    MPI_Recv(x_loc,N_update,MPI_DOUBLE,i_proc,3,internal_comm,&mpi_status);
	    for(int i=0;i<N_update;i++){
		res[update_loc[i]] = x_loc[update_index_loc[i]];
	    }

	    delete[] update_loc;
	    delete[] update_index_loc;
	    delete[] x_loc;

	}
    }else{
	MPI_Send(&N_update,1,MPI_INT,0,0,internal_comm);
	MPI_Send(update,N_update,MPI_INT,0,1,internal_comm);
	MPI_Send(update_index,N_update,MPI_INT,0,2,internal_comm);
	MPI_Send(x,N_update,MPI_DOUBLE,0,3,internal_comm);
    }

    MPI_Bcast(res,n,MPI_DOUBLE,0,internal_comm);

    delete[] b;
    delete[] x;

}

/************************************************************************************************/

template <class T>
int Aztec<T>::count_elements(int N_upd, int* upd, TCSR<T>* matrix)
{
    int element = 0;

    for(int i=0;i<N_upd;i++){
	element = element + matrix->edge_i[upd[i]+1-matrix->first_row] - \
	  matrix->edge_i[upd[i]-matrix->first_row];
    }

    return element;
}

/************************************************************************************************/

template <class T>
void Aztec<T>::create_matrix_row(int row, int location, TCSR<T>* matrix, T* value, int* bidx)
{
    int index = bidx[location];
    int start = matrix->edge_i[row-matrix->first_row]-matrix->findx;
    int stop  = matrix->edge_i[row+1-matrix->first_row]-matrix->findx;

    value[location] = 0.0;
    
    for(int i=start;i<stop;i++){
	if((matrix->index_j[i]-matrix->findx)==row){
	    value[location] = matrix->nnz[i];
	}else{
	    bidx[index]    = matrix->index_j[i]-matrix->findx;
	    value[index]   = matrix->nnz[i];
	    index++;
	}
    }
    
    bidx[location+1] = index;

}

/************************************************************************************************/

#endif


