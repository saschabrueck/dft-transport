/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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
template <typename T>
void MPI_Recv_Full(T *M,int count,int source,int message,MPI_Comm comm);
template <typename T>
void MPI_Dist_Sparse(TCSR<T> *M,int source,MPI_Comm comm);
template <typename T>
void MPI_Bcast_SpInd(TCSR<T> **M,int **M_index,int with_index,int trans_mpi_rank,MPI_Comm trans_comm);
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

template <>
void MPI_Recv_Full(double *M,int count,int source,int message,MPI_Comm comm)
{
    int var;
    MPI_Status status;
    
    MPI_Recv(&var,1,MPI_INT,source,2*message,comm,&status);
    MPI_Recv(M,count,MPI_DOUBLE,source,2*message+1,comm,&status);
}

/************************************************************************************************/

template <>
void MPI_Recv_Full(CPX *M,int count,int source,int message,MPI_Comm comm)
{
    int var;
    MPI_Status status;
    
    MPI_Recv(&var,1,MPI_INT,source,2*message,comm,&status);
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

template <>
void MPI_Dist_Sparse(TCSR<double> *M,int source,MPI_Comm comm)
{
    if(M->n_nonzeros>1){
        MPI_Bcast(M->index_i,M->size,MPI_INT,source,comm);
	MPI_Bcast(M->edge_i,M->size+1,MPI_INT,source,comm);
	MPI_Bcast(M->index_j,M->n_nonzeros,MPI_INT,source,comm);
	MPI_Bcast(M->diag_pos,M->size,MPI_INT,source,comm);
	MPI_Bcast(M->nnz,M->n_nonzeros,MPI_DOUBLE,source,comm);
    }
}

/************************************************************************************************/

template <>
void MPI_Dist_Sparse(TCSR<CPX> *M,int source,MPI_Comm comm)
{
    if(M->n_nonzeros>1){
        MPI_Bcast(M->index_i,M->size,MPI_INT,source,comm);
	MPI_Bcast(M->edge_i,M->size+1,MPI_INT,source,comm);
	MPI_Bcast(M->index_j,M->n_nonzeros,MPI_INT,source,comm);
	MPI_Bcast(M->diag_pos,M->size,MPI_INT,source,comm);
	MPI_Bcast(M->nnz,M->n_nonzeros,MPI_DOUBLE_COMPLEX,source,comm);
    }
}

/************************************************************************************************/

template <typename T>
void MPI_Bcast_SpInd(TCSR<T> **M,int **M_index,int with_index,int trans_mpi_rank,\
		     MPI_Comm trans_comm)
{

    int m_data[5];

    if(!trans_mpi_rank){
        m_data[0] = M[0]->size;
	m_data[1] = M[0]->n_nonzeros;
	m_data[2] = M[0]->findx;
	m_data[3] = M[0]->size_tot;
	m_data[4] = M[0]->first_row;
    }    
    MPI_Bcast(m_data,5,MPI_INT,0,trans_comm);

    if(trans_mpi_rank){
        *M       = new TCSR<T>(max(m_data[0],1),m_data[1],m_data[2]);
	*M_index = new int[max(m_data[1]*with_index,1)];
    }

    MPI_Dist_Sparse(*M,0,trans_comm);

    M[0]->size_tot  = m_data[3];
    M[0]->first_row = m_data[4];

    if(with_index){
	if(M[0]->n_nonzeros>1){
	    MPI_Bcast(M_index[0],M[0]->n_nonzeros,MPI_INT,0,trans_comm);
	}else{
	    M[0]->edge_i[0] = 0;
	    M[0]->edge_i[1] = 0;
	}
    }
}

/************************************************************************************************/

#endif
