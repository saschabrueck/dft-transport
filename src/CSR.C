using namespace std;

#include <math.h>
#include "CSR.H"

CSR::CSR(int N, int n_nnz)
{
    r_nnz      = new double[n_nnz];
    i_nnz      = new double[n_nnz];
    index_i    = new int[N];
    index_j    = new int[n_nnz];
    edge_i     = new int[N+1];
    diag_pos   = new int[N];

    size       = N;
    type       = 0;
    n_nonzeros = n_nnz;
}

/************************************************************************************************/

CSR::~CSR()
{
    delete [] r_nnz;
    delete [] i_nnz;
    delete [] index_i;
    delete [] index_j;
    delete [] edge_i;
    delete [] diag_pos;
}

/************************************************************************************************/

void CSR::update_diag(double *r_diag, double *i_diag)
{
    int i;

    for(i=0;i<size;i++){
        r_nnz[diag_pos[i]] = r_nnz[diag_pos[i]]+r_diag[i];
        i_nnz[diag_pos[i]] = i_nnz[diag_pos[i]]+i_diag[i];
    }
    
}

/************************************************************************************************/

void CSR::r_update_diag(double *r_diag)
{
    int i;

    for(i=0;i<size;i++){
        r_nnz[diag_pos[i]] = r_nnz[diag_pos[i]]+r_diag[i];
    }
    
}

/************************************************************************************************/

void CSR::i_update_diag(double *i_diag)
{
    int i;

    for(i=0;i<size;i++){
        i_nnz[diag_pos[i]] = i_nnz[diag_pos[i]]+i_diag[i];
    }
    
}

/************************************************************************************************/

void CSR::get_row_edge()
{
    int i;

    edge_i[0] = 0;
  
    for(i=0;i<size;i++){
        edge_i[i+1] = edge_i[i]+index_i[i];
    }
    
}

/************************************************************************************************/

template <>
void TCSR<CPX>::copy_contain(TCSR<double> *mat,double factor)
{

    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    n_nonzeros = mat->n_nonzeros;
    findx      = mat->findx;
    first_row  = mat->first_row;

    c_dcopy(n_nonzeros,mat->nnz,1,(double*)nnz,2);
    c_zscal(n_nonzeros,CPX(factor,0.0),nnz,1);
    c_icopy(size,mat->index_i,1,index_i,1);
    c_icopy(n_nonzeros,mat->index_j,1,index_j,1);
    c_icopy(size+1,mat->edge_i,1,edge_i,1);
    c_icopy(size,mat->diag_pos,1,diag_pos,1);
}

/************************************************************************************************/

template <>
void TCSR<double>::copy_contain(TCSR<double> *mat,double factor)
{

    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    n_nonzeros = mat->n_nonzeros;
    findx      = mat->findx;
    first_row  = mat->first_row;

    c_dcopy(n_nonzeros,mat->nnz,1,nnz,1);
    c_dscal(n_nonzeros,factor,nnz,1);
    c_icopy(size,mat->index_i,1,index_i,1);
    c_icopy(n_nonzeros,mat->index_j,1,index_j,1);
    c_icopy(size+1,mat->edge_i,1,edge_i,1);
    c_icopy(size,mat->diag_pos,1,diag_pos,1);
}

/************************************************************************************************/

template <>
void TCSR<CPX>::copy_contain(TCSR<CPX> *mat,double factor)
{

    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    n_nonzeros = mat->n_nonzeros;
    findx      = mat->findx;
    first_row  = mat->first_row;

    c_zcopy(n_nonzeros,mat->nnz,1,nnz,1);
    c_zscal(n_nonzeros,CPX(factor,0.0),nnz,1);
    c_icopy(size,mat->index_i,1,index_i,1);
    c_icopy(n_nonzeros,mat->index_j,1,index_j,1);
    c_icopy(size+1,mat->edge_i,1,edge_i,1);
    c_icopy(size,mat->diag_pos,1,diag_pos,1);
}

/************************************************************************************************/

template <>
void TCSR<double>::copy_contain(TCSR<CPX> *mat,double factor)
{

    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    n_nonzeros = mat->n_nonzeros;
    findx      = mat->findx;
    first_row  = mat->first_row;

    c_dcopy(n_nonzeros,(double*)mat->nnz,2,nnz,1);
    c_dscal(n_nonzeros,factor,nnz,1);
    c_icopy(size,mat->index_i,1,index_i,1);
    c_icopy(n_nonzeros,mat->index_j,1,index_j,1);
    c_icopy(size+1,mat->edge_i,1,edge_i,1);
    c_icopy(size,mat->diag_pos,1,diag_pos,1);
}

/************************************************************************************************/

template <>
void TCSR<double>::sparse_to_cmp_full(CPX *B,int nrow,int ncol)
{
    int i,j,index;

    init_variable(B,nrow*ncol);

    for(i=0;i<size;i++){
        for(j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            index             = index_j[j]-findx;
            B[i+index*nrow]   = CPX(nnz[j],0.0);
        }
    }
}

/************************************************************************************************/

template <>
void TCSR<CPX>::sparse_to_cmp_full(CPX *B,int nrow,int ncol)
{
}

/************************************************************************************************/

template <>
void TCSR<double>::cmp_full_to_sparse(CPX *B,int nrow,int ncol,CPX factor)
{
    int i,j;

    size       = nrow;
    n_nonzeros = 0;

    for(i=0;i<nrow;i++){

        index_i[i] = 0;

        for(j=0;j<ncol;j++){

        if(abs(B[i+j*nrow])>tollim){

                index_j[n_nonzeros] = j+findx;
                nnz[n_nonzeros]     = real(factor*B[i+j*nrow]);

                index_i[i]++;
                n_nonzeros++;
            }
        }
    }

    get_row_edge();
    get_diag_pos();
}

/************************************************************************************************/

template <>
void TCSR<CPX>::cmp_full_to_sparse(CPX *B,int nrow,int ncol,CPX factor)
{
    int i,j;

    size       = nrow;
    n_nonzeros = 0;

    for(i=0;i<nrow;i++){

        index_i[i] = 0;

        for(j=0;j<ncol;j++){

        if(abs(B[i+j*nrow])>tollim){

                index_j[n_nonzeros] = j+findx;
                nnz[n_nonzeros]     = factor*B[i+j*nrow];

                index_i[i]++;
                n_nonzeros++;
            }
        }
    }

    get_row_edge();
    get_diag_pos();
}
