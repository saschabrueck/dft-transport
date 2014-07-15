#ifndef __CONNECTION
#define __CONNECTION

#include <mpi.h>

template <typename T>
void MPI_Send_Sparse(TCSR<T> *M,int dest,int message,MPI_Comm comm);
template <typename T>
void MPI_Recv_Sparse(TCSR<T> *M,int source,int message,MPI_Comm comm);
template <typename T>
void MPI_Send_Full(T *M,int count,int dest,int message,MPI_Comm comm);
template <typename T>
void MPI_Isend_Full(T *M,int count,int dest,int message,MPI_Comm comm);
template <typename T>
void MPI_Recv_Full(T *M,int source,int message,MPI_Comm comm);
void MPI_Send_Res(CPX *M,int N,int Ntot,int no_rhs,int dest,MPI_Comm comm);
void MPI_Recv_Res(CPX *M,int source,MPI_Comm comm);
void MPI_Bcast_Sparse(TCSR<CPX> *M,int source,MPI_Comm comm);

/************************************************************************************************/

template <>
void MPI_Send_Sparse(TCSR<double> *M,int dest,int message,MPI_Comm comm)
{
    MPI_Send(&M->size,1,MPI_INT,dest,5*message,comm);
    MPI_Send(&M->n_nonzeros,1,MPI_INT,dest,5*message+1,comm);
    MPI_Send(M->edge_i,M->size+1,MPI_INT,dest,5*message+2,comm);
    MPI_Send(M->index_j,M->n_nonzeros,MPI_INT,dest,5*message+3,comm);
    MPI_Send(M->nnz,M->n_nonzeros,MPI_DOUBLE,dest,5*message+4,comm);
        
}

/************************************************************************************************/

template <>
void MPI_Send_Sparse(TCSR<CPX> *M,int dest,int message,MPI_Comm comm)
{
    MPI_Send(&M->size,1,MPI_INT,dest,5*message,comm);
    MPI_Send(&M->n_nonzeros,1,MPI_INT,dest,5*message+1,comm);
    MPI_Send(M->edge_i,M->size+1,MPI_INT,dest,5*message+2,comm);
    MPI_Send(M->index_j,M->n_nonzeros,MPI_INT,dest,5*message+3,comm);
    MPI_Send(M->nnz,M->n_nonzeros,MPI_DOUBLE_COMPLEX,dest,5*message+4,comm);
        
}

/************************************************************************************************/

template <>
void MPI_Recv_Sparse(TCSR<double> *M,int source,int message,MPI_Comm comm)
{
    MPI_Status status;

    MPI_Recv(&M->size,1,MPI_INT,source,5*message,comm,&status);
    MPI_Recv(&M->n_nonzeros,1,MPI_INT,source,5*message+1,comm,&status);
    MPI_Recv(M->edge_i,M->size+1,MPI_INT,source,5*message+2,comm,&status);
    MPI_Recv(M->index_j,M->n_nonzeros,MPI_INT,source,5*message+3,comm,&status);
    MPI_Recv(M->nnz,M->n_nonzeros,MPI_DOUBLE,source,5*message+4,comm,&status);
}

/************************************************************************************************/

template <>
void MPI_Recv_Sparse(TCSR<CPX> *M,int source,int message,MPI_Comm comm)
{
    MPI_Status status;
    
    MPI_Recv(&M->size,1,MPI_INT,source,5*message,comm,&status);
    MPI_Recv(&M->n_nonzeros,1,MPI_INT,source,5*message+1,comm,&status);
    MPI_Recv(M->edge_i,M->size+1,MPI_INT,source,5*message+2,comm,&status);
    MPI_Recv(M->index_j,M->n_nonzeros,MPI_INT,source,5*message+3,comm,&status);
    MPI_Recv(M->nnz,M->n_nonzeros,MPI_DOUBLE_COMPLEX,source,5*message+4,comm,&status);
}

/************************************************************************************************/

template <>
void MPI_Send_Full(double *M,int count,int dest,int message,MPI_Comm comm)
{
    MPI_Send(&count,1,MPI_INT,dest,2*message,comm);
    MPI_Send(M,count,MPI_DOUBLE,dest,2*message+1,comm);
        
}

/************************************************************************************************/

template <>
void MPI_Send_Full(CPX *M,int count,int dest,int message,MPI_Comm comm)
{
    MPI_Send(&count,1,MPI_INT,dest,2*message,comm);
    MPI_Send(M,count,MPI_DOUBLE_COMPLEX,dest,2*message+1,comm);
        
}

/************************************************************************************************/

template <>
void MPI_Isend_Full(double *M,int count,int dest,int message,MPI_Comm comm)
{
    MPI_Request request;

    MPI_Isend(&count,1,MPI_INT,dest,2*message,comm,&request);
    MPI_Isend(M,count,MPI_DOUBLE,dest,2*message+1,comm,&request);    
}

/************************************************************************************************/

template <>
void MPI_Isend_Full(CPX *M,int count,int dest,int message,MPI_Comm comm)
{
    MPI_Request request;

    MPI_Isend(&count,1,MPI_INT,dest,2*message,comm,&request);
    MPI_Isend(M,count,MPI_DOUBLE_COMPLEX,dest,2*message+1,comm,&request);        
}

/************************************************************************************************/

template <>
void MPI_Recv_Full(double *M,int source,int message,MPI_Comm comm)
{
    int count;
    MPI_Status status;

    MPI_Recv(&count,1,MPI_INT,source,2*message,comm,&status);
    MPI_Recv(M,count,MPI_DOUBLE,source,2*message+1,comm,&status);
}

/************************************************************************************************/

template <>
void MPI_Recv_Full(CPX *M,int source,int message,MPI_Comm comm)
{
    int count;
    MPI_Status status;
    
    MPI_Recv(&count,1,MPI_INT,source,2*message,comm,&status);
    MPI_Recv(M,count,MPI_DOUBLE_COMPLEX,source,2*message+1,comm,&status);
}

/************************************************************************************************/

void MPI_Send_Res(CPX *M,int N,int Ntot,int no_rhs,int dest,MPI_Comm comm)
{

    int IR;
    
    MPI_Send(&N,1,MPI_INT,dest,0,comm);
    MPI_Send(&no_rhs,1,MPI_INT,dest,1,comm);

    for(IR=0;IR<no_rhs;IR++){
	MPI_Send(&M[IR*Ntot],N,MPI_DOUBLE_COMPLEX,dest,2+IR,comm);
    }

}

/************************************************************************************************/

void MPI_Recv_Res(CPX *M,int source,MPI_Comm comm)
{
    MPI_Status status;
    int N,no_rhs,IR;

    MPI_Recv(&N,1,MPI_INT,source,0,comm,&status);
    MPI_Recv(&no_rhs,1,MPI_INT,source,1,comm,&status);

    for(IR=0;IR<no_rhs;IR++){
	MPI_Recv(&M[IR*N],N,MPI_DOUBLE_COMPLEX,source,2+IR,comm,&status);
    }

}

/************************************************************************************************/

inline void MPI_Bcast_Sparse(TCSR<CPX> *M,int source,MPI_Comm comm)
{
    MPI_Bcast(&M->size,1,MPI_INT,source,comm);
    MPI_Bcast(&M->n_nonzeros,1,MPI_INT,source,comm);
    MPI_Bcast(M->index_i,M->size,MPI_INT,source,comm);
    MPI_Bcast(M->edge_i,M->size+1,MPI_INT,source,comm);
    MPI_Bcast(M->index_j,M->n_nonzeros,MPI_INT,source,comm);
    MPI_Bcast(M->diag_pos,M->size,MPI_INT,source,comm);
    MPI_Bcast(M->nnz,M->n_nonzeros,MPI_DOUBLE_COMPLEX,source,comm);
}

/************************************************************************************************/

#endif
