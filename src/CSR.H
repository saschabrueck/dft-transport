#ifndef __CSR
#define __CSR

#include <omp.h>
#include <fstream>
#include <algorithm>
#include <complex>
#include "Blas.H"
#include "array_tools.H"

#define LOGCERR std::cerr << "ERROR IN LINE " << __LINE__ << " OF FILE " << __FILE__ << std::endl

struct CSR_Exception{
    CSR_Exception(const int line,const char* file) {std::cerr<<"Error in line "<<line<<" of file "<<file<<std::endl;}
};

class CSR{
public:
    
    CSR(int,int);
    ~CSR();
    void update_diag(double*,double*);
    void r_update_diag(double*);
    void i_update_diag(double*);
    void get_row_edge();

    int pos,size,n_nonzeros,type;
    double *r_nnz,*i_nnz;
    int *index_i,*index_j,*edge_i,*diag_pos;
    
private:
};

template <class T>
class TCSR{
public:
    
    TCSR(int,int,int);
    TCSR(TCSR<T>*);
    TCSR(const char*);
    TCSR(char*,int,int);
// sab begin
//create new matrix by cutting out from old one
    TCSR(TCSR<T>*,int,int,int,int);
// THIS CONSTRUCTOR WILL CREATE A MATRIX BY WEIGHTED SUMMING OF TWO OTHERS
    TCSR(T,TCSR<double>*,T,TCSR<double>*);
    TCSR(CPX,TCSR<CPX>*,CPX,TCSR<CPX>*);
// this constructs a matrix by summing of a vector of matrices
    TCSR(int,T*,TCSR<T>**,int=-1);
// this creates a matrix which is the input matrix multiplied by the sum
// of the input vector elements corresponding to the row and col number
    TCSR(TCSR<T>*,T*);
// and this is the gather
    TCSR(TCSR<T>*,MPI_Comm);
    TCSR(TCSR<T>*,int,MPI_Comm);
    TCSR(TCSR<T>*,int*,int,MPI_Comm);
    TCSR(TCSR<T>*,MPI_Comm,int,MPI_Comm*,MPI_Comm*);
//distribute matrix according to pattern
    TCSR(TCSR<T>*,TCSR<T>*,int,MPI_Comm);
// sab end
    ~TCSR();
    void change_contain(TCSR<T>*);
    void copy_index(TCSR<double>*);
    void copy_index(TCSR<CPX>*);
    void copy_contain(TCSR<double>*,double);
    void copy_contain(TCSR<CPX>*,double);
    void cpu_distribute(int*,int*,int,int,int);
    void update_diag(T*,T);
    void update_loc_diag(T*,T);
    void update_diag_single(T);
    void update_diag_single(T,int,int);
    void update_diag_pack(T*,T,int);
    void update_rdiag_pack(double*,int,int,double,int,int*);
    void get_row_edge();
    void get_diag_pos();
    void scale_element(int,int*,T);
    void mat_vec_mult(T*,T*,int);
    void mat_vec_mult(T*,T*,int,int);
    void trans_mat_vec_mult(T*,T*,int,int);
    void trans_mat_vec_mult_omp(T*,T*,int,int);
    void trans_mat_vec_mult(T*,T*,int,int,int);
    void spec_mat_vec_mult(T*,T*,int,int*,int*,int,int);
    void spec_mat_vec_mult(T*,T*,int,int*,int*,int,int,int*);
    void sparse_to_full(T*,int,int);
    void sparse_to_cmp_full(CPX*,int,int);
    void full_to_sparse(T*,int,int);
    void sparse_to_full_trans(T *,int ,int );
    void cmp_full_to_sparse(CPX*,int,int,CPX);
    void write_CSR(const char*);
    void write(const char*);
    void init_variable(int*,int);
    void init_variable(double*,int);
    void init_variable(CPX*,int);
    void add_to_corner(T*,T,int*,int);
    void change_findx(int);
    void get_full_index_i(int*);
    int get_msize(int,int,int,int*);
// sab begin
    void write_CSR_bin(const char*);
    void mat_vec_mult_add(T*,T*,int,T);
    void shift_resize(int,int,int,int);
    void replicate_cell(int,int,int,int,MPI_Comm);
    void fill_pbc_block(TCSR<T>*,int,int,int,MPI_Comm);
    int additionalentries(TCSR<T>*);
    void add(TCSR<T>*,T);
    void add_real(TCSR<CPX>*,CPX);
    void Bcast(int,MPI_Comm);
    void scatter(TCSR<T>*,int,MPI_Comm);
    void reduce(int,MPI_Comm);
    void reducescatter(TCSR<T>*,MPI_Comm);
    double psipsidagger(CPX*,int);
    void psipsidagger(CPX*,int,CPX);
    void psipsidagger(double*,int,double);
    void psipsidagger_transpose(double*,int,double);
    void psipsidagger(CPX*,CPX*,CPX*,int,int,int,CPX);
    void psipsidagger(CPX*,int,CPX,int*,double*);
    void psipsidagger(CPX*,int,CPX,int*,double*,double*,double);
    void psipsidagger(CPX*,CPX*,int,CPX,int*,double*);
    void atomdensity(TCSR<CPX>*,CPX,int*,double*);
    double collect_density(double*,double*,int*);
    void moveawaypbc(T*,T*,int,int);
    void settozeropbc(int,int);
    void removepbc(int,int);
    void tr_full_to_sparse(T*,int,int);
    void full_to_sparse(T*,int,int,int,int);
    void sparse_transpose(TCSR<T> *mat);
    void extract_lower_triangle(TCSR<T>*);
    void extract_upper_triangle(TCSR<T>*);
    void sparse_to_full(T*,int,int,int,int);
    void set_to_id();
// sab end
    void set_nnz(int,T*);
    void set_index_i(int,int*);
    void set_index_j(int,int*);
    void set_edge_i(int,int*);
    
    
    void get_dense_by_range(T* matrix_out, int start_row, int start_col, int end_row, int end_col);
    void get_dense_by_range_trans(T* matrix_out, int start_row, int start_col, int end_row, int end_col);

    int pos,size,size_tot,n_nonzeros,type,findx,first_row;
    T *nnz;
    int *index_i,*index_j,*edge_i,*diag_pos;
    
private:
    inline MPI_Datatype give_me_MPI_datatype (double*)
    {
        return MPI_DOUBLE;
    }
    inline MPI_Datatype give_me_MPI_datatype (CPX*)
    {
        return MPI_DOUBLE_COMPLEX;
    }
};

template <class T>
TCSR<T>::TCSR(int N, int n_nnz, int fortran_index)
{
    nnz        = new T[max(n_nnz,1)];
    index_i    = new int[max(N,1)];
    index_j    = new int[max(n_nnz,1)];
    edge_i     = new int[N+1];
    diag_pos   = new int[max(N,1)];

    size       = N;
    size_tot   = N;
    type       = 0;
    n_nonzeros = n_nnz;
    findx      = fortran_index;
    first_row  = 0;                                                    
}

/************************************************************************************************/

template <class T>
inline TCSR<T>::TCSR(TCSR<T> *mat)
{

    nnz        = new T[mat->n_nonzeros];
    index_i    = new int[mat->size];
    index_j    = new int[mat->n_nonzeros];
    edge_i     = new int[mat->size+1];
    diag_pos   = new int[mat->size];

    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    n_nonzeros = mat->n_nonzeros;
    findx      = mat->findx;
    first_row  = mat->first_row;
    
    c_tcopy(n_nonzeros,mat->nnz,1,nnz,1);
    c_icopy(size,mat->index_i,1,index_i,1);
    c_icopy(n_nonzeros,mat->index_j,1,index_j,1);
    c_icopy(size+1,mat->edge_i,1,edge_i,1);
    c_icopy(size,mat->diag_pos,1,diag_pos,1);
}

/******************** written by SAB ********************/

template <class T>
inline TCSR<T>::TCSR(TCSR<T> *mat,int istart, int ilength, int jstart, int jlength)
{
    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    findx      = mat->findx;
    first_row  = mat->first_row;

    int iend       = istart+ilength;
    int jend       = jstart+jlength;
    int loc_istart = max(0,min(size,istart-first_row));
    int loc_iend   = max(0,min(size,iend-first_row));

    n_nonzeros = 0;

    for (int i=loc_istart;i<loc_iend;i++) {
        for(int e=mat->edge_i[i]-mat->findx;e<mat->edge_i[i+1]-mat->findx;e++){
            if(mat->index_j[e]-mat->findx>=jstart && mat->index_j[e]-mat->findx<jend){
                n_nonzeros++;
            }
        }
    }

    nnz        = new T[max(n_nonzeros,1)];
    index_i    = new int[max(size,1)];
    index_j    = new int[max(n_nonzeros,1)];
    edge_i     = new int[size+1];
    diag_pos   = new int[max(size,1)];

    n_nonzeros = 0;
    init_variable(index_i,size);

    for (int i=loc_istart;i<loc_iend;i++) {
        for(int e=mat->edge_i[i]-mat->findx;e<mat->edge_i[i+1]-mat->findx;e++){
            if(mat->index_j[e]-mat->findx>=jstart && mat->index_j[e]-mat->findx<jend){
                index_i[i]++;
                index_j[n_nonzeros]=mat->index_j[e];
                nnz[n_nonzeros]=mat->nnz[e];
                n_nonzeros++;
            }
        }
    }

    get_row_edge();
    get_diag_pos();

    if (n_nonzeros!=edge_i[size]-findx) throw CSR_Exception(__LINE__,__FILE__);
}

/******************** written by SAB ********************/

template <>
inline TCSR<CPX>::TCSR(CPX alphafactor,TCSR<CPX> *A,CPX betafactor,TCSR<CPX> *B)
{
    if ( (size       = A->size      )!=B->size      ) throw CSR_Exception(__LINE__,__FILE__);
    if ( (size_tot   = A->size_tot  )!=B->size_tot  ) throw CSR_Exception(__LINE__,__FILE__);
    if ( (type       = A->type      )!=B->type      ) throw CSR_Exception(__LINE__,__FILE__);
    if ( (findx      = A->findx     )!=B->findx     ) throw CSR_Exception(__LINE__,__FILE__);
    if ( (first_row  = A->first_row )!=B->first_row ) throw CSR_Exception(__LINE__,__FILE__);

    int addit  = B->additionalentries(A);

    if (addit<0) throw CSR_Exception(__LINE__,__FILE__);
  
    n_nonzeros = addit+A->n_nonzeros;

    index_i    = new int[size];
    edge_i     = new int[size+1];
    nnz        = new CPX[n_nonzeros];
    index_j    = new int[n_nonzeros];
    diag_pos   = new int[size];

    int i,j,e;
    n_nonzeros=0;
    CPX *linefull=new CPX[size_tot];

    for (i=0;i<size;i++) {
        init_variable(linefull,size_tot);
        index_i[i]=0;
        for (e=A->edge_i[i]-A->findx;e<A->edge_i[i+1]-A->findx;e++) {
            linefull[A->index_j[e]-A->findx]+=alphafactor*A->nnz[e];
        }
        for (e=B->edge_i[i]-B->findx;e<B->edge_i[i+1]-B->findx;e++) {
            linefull[B->index_j[e]-B->findx]+=betafactor*B->nnz[e];
        }
        for (j=0;j<size_tot;j++) {
            if (abs(linefull[j])) {
                nnz[n_nonzeros]=linefull[j];
                index_j[n_nonzeros]=j+findx;
                index_i[i]++;
                n_nonzeros++;
            }
        }
    }

    delete[] linefull;

    get_row_edge();
    get_diag_pos();
}

/******************** written by SAB ********************/

template <class T>
inline TCSR<T>::TCSR(T alphafactor,TCSR<double> *A,T betafactor,TCSR<double> *B)
{
    if ( (size       = A->size      )!=B->size      ) throw CSR_Exception(__LINE__,__FILE__);
    if ( (size_tot   = A->size_tot  )!=B->size_tot  ) throw CSR_Exception(__LINE__,__FILE__);
    if ( (type       = A->type      )!=B->type      ) throw CSR_Exception(__LINE__,__FILE__);
    if ( (findx      = A->findx     )!=B->findx     ) throw CSR_Exception(__LINE__,__FILE__);
    if ( (first_row  = A->first_row )!=B->first_row ) throw CSR_Exception(__LINE__,__FILE__);

    int addit  = B->additionalentries(A);

    if (addit<0) throw CSR_Exception(__LINE__,__FILE__);
  
    n_nonzeros = addit+A->n_nonzeros;

    index_i    = new int[size];
    edge_i     = new int[size+1];
    nnz        = new T[max(n_nonzeros,1)];
    index_j    = new int[max(n_nonzeros,1)];
    diag_pos   = new int[size];

    int i,j,e;
    n_nonzeros=0;
    T *linefull=new T[size_tot];

    for (i=0;i<size;i++) {
        init_variable(linefull,size_tot);
        index_i[i]=0;
        for (e=A->edge_i[i]-A->findx;e<A->edge_i[i+1]-A->findx;e++) {
            linefull[A->index_j[e]-A->findx]+=alphafactor*A->nnz[e];
        }
        for (e=B->edge_i[i]-B->findx;e<B->edge_i[i+1]-B->findx;e++) {
            linefull[B->index_j[e]-B->findx]+=betafactor*B->nnz[e];
        }
        for (j=0;j<size_tot;j++) {
            if (abs(linefull[j])) {
                nnz[n_nonzeros]=linefull[j];
                index_j[n_nonzeros]=j+findx;
                index_i[i]++;
                n_nonzeros++;
            }
        }
    }

    delete[] linefull;

    get_row_edge();
    get_diag_pos();
}

/******************** written by SAB ********************/

template <class T>
inline TCSR<T>::TCSR(int N,T *factors,TCSR<T> **A,int jlength)
{
    int i,j,e,n;
    T factor = 1.0;

    size       = A[0]->size; 
    size_tot   = A[0]->size_tot;
    type       = A[0]->type;
    findx      = A[0]->findx;
    first_row  = A[0]->first_row;
    n_nonzeros = A[0]->n_nonzeros;
    for (n=1;n<N;n++) {
        if (size      != A[n]->size      ) throw CSR_Exception(__LINE__,__FILE__);
        if (size_tot  != A[n]->size_tot  ) throw CSR_Exception(__LINE__,__FILE__);
        if (type      != A[n]->type      ) throw CSR_Exception(__LINE__,__FILE__);
        if (findx     != A[n]->findx     ) throw CSR_Exception(__LINE__,__FILE__);
        if (first_row != A[n]->first_row ) throw CSR_Exception(__LINE__,__FILE__);
        n_nonzeros+=A[n]->n_nonzeros;
    }

    index_i  = new int[size];
    edge_i   = new int[size+1];
    diag_pos = new int[size];

    T *nnz_tmp       = new T[max(n_nonzeros,1)];
    int *index_j_tmp = new int[max(n_nonzeros,1)];

    n_nonzeros=0;

    if (jlength<0) jlength=size_tot;
    T *linefull = new T[jlength];

    for (i=0;i<size;i++) {
        init_variable(linefull,jlength);
        index_i[i]=0;
        for (n=0;n<N;n++) {
            if (factors) factor=factors[n];
            for (e=A[n]->edge_i[i]-A[n]->findx;e<A[n]->edge_i[i+1]-A[n]->findx;e++) {
                linefull[A[n]->index_j[e]-A[n]->findx]+=factor*A[n]->nnz[e];
            }
        }
        for (j=0;j<jlength;j++) {
            if (abs(linefull[j])) {
                nnz_tmp[n_nonzeros]=linefull[j];
                index_j_tmp[n_nonzeros]=j+findx;
                index_i[i]++;
                n_nonzeros++;
            }
        }
    }

    delete[] linefull;

    nnz = new T[max(n_nonzeros,1)];
    c_tcopy(n_nonzeros,nnz_tmp,1,nnz,1);
    delete[] nnz_tmp;
    index_j = new int[max(n_nonzeros,1)];
    c_icopy(n_nonzeros,index_j_tmp,1,index_j,1);
    delete[] index_j_tmp;

    get_row_edge();
    get_diag_pos();
}

/******************** written by SAB ********************/

template <class T>
inline TCSR<T>::TCSR(TCSR<T> *mat,T *vec)
{

    nnz        = new T[max(mat->n_nonzeros,1)];
    index_i    = new int[mat->size];
    index_j    = new int[max(mat->n_nonzeros,1)];
    edge_i     = new int[mat->size+1];
    diag_pos   = new int[mat->size];

    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    n_nonzeros = mat->n_nonzeros;
    findx      = mat->findx;
    first_row  = mat->first_row;
    
    c_icopy(size,mat->index_i,1,index_i,1);
    c_icopy(n_nonzeros,mat->index_j,1,index_j,1);
    c_icopy(size+1,mat->edge_i,1,edge_i,1);
    c_icopy(size,mat->diag_pos,1,diag_pos,1);

    for(int i=0;i<size;i++) {
        for(int e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            nnz[e]=(vec[i+first_row]+vec[index_j[e]-findx])/2*mat->nnz[e];
        }
    }

}

/******************** written by SAB ********************/

template <class T>
inline TCSR<T>::TCSR(TCSR<T> *matlocal,MPI_Comm my_comm)
{
    int iam,nprocs;
    MPI_Comm_size(my_comm,&nprocs);
    MPI_Comm_rank(my_comm,&iam);

    int loc_n_nonzeros=matlocal->n_nonzeros;
    int* dist_loc_n_nonzeros=new int[nprocs];
    MPI_Allgather(&loc_n_nonzeros,1,MPI_INT,dist_loc_n_nonzeros,1,MPI_INT,my_comm);
    int* displc_n_nonzeros=new int[nprocs+1];
    displc_n_nonzeros[0]=0;
    for (int iii=1;iii<nprocs+1;iii++)
        displc_n_nonzeros[iii]=displc_n_nonzeros[iii-1]+dist_loc_n_nonzeros[iii-1];

    int loc_size=matlocal->size;
    int* dist_loc_size=new int[nprocs];
    MPI_Allgather(&loc_size,1,MPI_INT,dist_loc_size,1,MPI_INT,my_comm);
    int* displc_size=new int[nprocs+1];
    displc_size[0]=0;
    for (int iii=1;iii<nprocs+1;iii++)
        displc_size[iii]=displc_size[iii-1]+dist_loc_size[iii-1];
    if (displc_size[iam]!=matlocal->first_row) throw CSR_Exception(__LINE__,__FILE__);

    n_nonzeros = displc_n_nonzeros[nprocs];
    size_tot   = matlocal->size_tot;
    size       = size_tot;
    findx      = matlocal->findx;
    first_row  = 0;
    type       = 0;
    nnz        = new T[max(n_nonzeros,1)];
    index_i    = new int[max(size,1)];
    index_j    = new int[max(n_nonzeros,1)];
    edge_i     = new int[size+1];
    diag_pos   = new int[max(size,1)];

    MPI_Allgatherv(matlocal->nnz,matlocal->n_nonzeros,give_me_MPI_datatype(nnz),nnz,dist_loc_n_nonzeros,displc_n_nonzeros,give_me_MPI_datatype(nnz),my_comm);
    MPI_Allgatherv(matlocal->index_j,matlocal->n_nonzeros,MPI_INT,index_j,dist_loc_n_nonzeros,displc_n_nonzeros,MPI_INT,my_comm);
    MPI_Allgatherv(matlocal->index_i,matlocal->size,MPI_INT,index_i,dist_loc_size,displc_size,MPI_INT,my_comm);

    delete[] dist_loc_n_nonzeros;
    delete[] displc_n_nonzeros;
    delete[] dist_loc_size;
    delete[] displc_size;

    get_row_edge();
    get_diag_pos();
}

/******************** written by SAB ********************/

template <class T>
inline TCSR<T>::TCSR(TCSR<T> *matlocal,int rootint,MPI_Comm my_comm)
{
    int iam,nprocs;
    MPI_Comm_size(my_comm,&nprocs);
    MPI_Comm_rank(my_comm,&iam);

    int loc_n_nonzeros=matlocal->n_nonzeros;
    int* dist_loc_n_nonzeros=new int[nprocs];
    MPI_Allgather(&loc_n_nonzeros,1,MPI_INT,dist_loc_n_nonzeros,1,MPI_INT,my_comm);
    int* displc_n_nonzeros=new int[nprocs+1]();
    for (int iii=1;iii<nprocs+1;iii++)
        displc_n_nonzeros[iii]=displc_n_nonzeros[iii-1]+dist_loc_n_nonzeros[iii-1];

    int loc_size=matlocal->size;
    int* dist_loc_size=new int[nprocs];
    MPI_Allgather(&loc_size,1,MPI_INT,dist_loc_size,1,MPI_INT,my_comm);
    int* displc_size=new int[nprocs+1]();
    for (int iii=1;iii<nprocs+1;iii++)
        displc_size[iii]=displc_size[iii-1]+dist_loc_size[iii-1];
    if (displc_size[iam]!=matlocal->first_row) throw CSR_Exception(__LINE__,__FILE__);

    if (iam==rootint) {
        n_nonzeros = displc_n_nonzeros[nprocs];
        size_tot   = matlocal->size_tot;
    } else {
        n_nonzeros = 0;
        size_tot   = 0;
    }
    size       = size_tot;
    findx      = matlocal->findx;
    first_row  = 0;
    type       = 0;
    nnz        = new T[max(n_nonzeros,1)];
    index_i    = new int[max(size,1)];
    index_j    = new int[max(n_nonzeros,1)];
    edge_i     = new int[size+1];
    diag_pos   = new int[max(size,1)];

    MPI_Gatherv(matlocal->nnz,matlocal->n_nonzeros,give_me_MPI_datatype(nnz),nnz,dist_loc_n_nonzeros,displc_n_nonzeros,give_me_MPI_datatype(nnz),rootint,my_comm);
    MPI_Gatherv(matlocal->index_j,matlocal->n_nonzeros,MPI_INT,index_j,dist_loc_n_nonzeros,displc_n_nonzeros,MPI_INT,rootint,my_comm);
    MPI_Gatherv(matlocal->index_i,matlocal->size,MPI_INT,index_i,dist_loc_size,displc_size,MPI_INT,rootint,my_comm);

    delete[] dist_loc_n_nonzeros;
    delete[] displc_n_nonzeros;
    delete[] dist_loc_size;
    delete[] displc_size;

    get_row_edge();
    get_diag_pos();
}

/******************** written by SAB ********************/

template <class T>
inline TCSR<T>::TCSR(TCSR<T> *matlocal,int *rootvec,int rootnum,MPI_Comm my_comm)
{
    int iam,nprocs;
    MPI_Comm_size(my_comm,&nprocs);
    MPI_Comm_rank(my_comm,&iam);

    int loc_n_nonzeros=matlocal->n_nonzeros;
    int* dist_loc_n_nonzeros=new int[nprocs];
    MPI_Allgather(&loc_n_nonzeros,1,MPI_INT,dist_loc_n_nonzeros,1,MPI_INT,my_comm);
    int* displc_n_nonzeros=new int[nprocs+1]();
    for (int iii=1;iii<nprocs+1;iii++)
        displc_n_nonzeros[iii]=displc_n_nonzeros[iii-1]+dist_loc_n_nonzeros[iii-1];

    int loc_size=matlocal->size;
    int* dist_loc_size=new int[nprocs];
    MPI_Allgather(&loc_size,1,MPI_INT,dist_loc_size,1,MPI_INT,my_comm);
    int* displc_size=new int[nprocs+1]();
    for (int iii=1;iii<nprocs+1;iii++)
        displc_size[iii]=displc_size[iii-1]+dist_loc_size[iii-1];
    if (displc_size[iam]!=matlocal->first_row) throw CSR_Exception(__LINE__,__FILE__);

    n_nonzeros = 0;
    size_tot   = 0;
    for (int iroot=0;iroot<rootnum;iroot++){
        if (iam==rootvec[iroot]){
            n_nonzeros = displc_n_nonzeros[nprocs];
            size_tot   = matlocal->size_tot;
        }
    }
    size       = size_tot;
    findx      = matlocal->findx;
    first_row  = 0;
    type       = 0;
    nnz        = new T[max(n_nonzeros,1)];
    index_i    = new int[max(size,1)];
    index_j    = new int[max(n_nonzeros,1)];
    edge_i     = new int[size+1];
    diag_pos   = new int[max(size,1)];

    MPI_Gatherv(matlocal->nnz,matlocal->n_nonzeros,give_me_MPI_datatype(nnz),nnz,dist_loc_n_nonzeros,displc_n_nonzeros,give_me_MPI_datatype(nnz),rootvec[0],my_comm);
    MPI_Gatherv(matlocal->index_j,matlocal->n_nonzeros,MPI_INT,index_j,dist_loc_n_nonzeros,displc_n_nonzeros,MPI_INT,rootvec[0],my_comm);
    MPI_Gatherv(matlocal->index_i,matlocal->size,MPI_INT,index_i,dist_loc_size,displc_size,MPI_INT,rootvec[0],my_comm);

/*
//Bcast segfaults on pilatus
    MPI_Group my_group,sub_group;
    MPI_Comm_group(my_comm,&my_group);
    MPI_Group_incl(my_group,rootnum,rootvec,&sub_group);
    MPI_Comm sub_comm;
    MPI_Comm_create(my_comm,sub_group,&sub_comm);
    MPI_Group_free(&my_group);
    MPI_Group_free(&sub_group);
    MPI_Bcast(nnz,n_nonzeros,give_me_MPI_datatype(nnz),rootvec[0],sub_comm);
    MPI_Bcast(index_j,n_nonzeros,MPI_INT,rootvec[0],sub_comm);
    MPI_Bcast(index_i,size,MPI_INT,rootvec[0],sub_comm);
    MPI_Comm_free(&sub_comm);
*/

// /*
    struct tags {
        enum {
            nnz,
            index_j,
            index_i
        };
    };

    MPI_Status status;
    if (iam==rootvec[0]){
        for (int iroot=1;iroot<rootnum;iroot++){
            MPI_Send(nnz,n_nonzeros,give_me_MPI_datatype(nnz),rootvec[iroot],tags::nnz,my_comm);
            MPI_Send(index_j,n_nonzeros,MPI_INT,rootvec[iroot],tags::index_j,my_comm);
            MPI_Send(index_i,size,MPI_INT,rootvec[iroot],tags::index_i,my_comm);
        }
    }else{
        for (int iroot=1;iroot<rootnum;iroot++){
            if (iam==rootvec[iroot]){
                MPI_Recv(nnz,n_nonzeros,give_me_MPI_datatype(nnz),rootvec[0],tags::nnz,my_comm,&status);
                MPI_Recv(index_j,n_nonzeros,MPI_INT,rootvec[0],tags::index_j,my_comm,&status);
                MPI_Recv(index_i,size,MPI_INT,rootvec[0],tags::index_i,my_comm,&status);
            }
        }
    }
// */

    delete[] dist_loc_n_nonzeros;
    delete[] displc_n_nonzeros;
    delete[] dist_loc_size;
    delete[] displc_size;

    get_row_edge();
    get_diag_pos();
}

/******************** written by SAB ********************/

template <class T>
inline TCSR<T>::TCSR(TCSR<T> *matlocal,MPI_Comm my_comm_in,int size_my_comm_out,MPI_Comm *my_comm_out,MPI_Comm *equal_rank_comm)
{
    struct tags {
        enum {
            first_row,
            nnz,
            index_j,
            index_i
        };
    };

    int iam,nprocs;
    MPI_Comm_size(my_comm_in,&nprocs);
    MPI_Comm_rank(my_comm_in,&iam);
    if (nprocs % size_my_comm_out) throw CSR_Exception(__LINE__,__FILE__);
    int nmat = nprocs / size_my_comm_out;
    int color = iam / size_my_comm_out;
    MPI_Comm_split(my_comm_in,color,iam,my_comm_out);
    int rank_my_comm_out;
    MPI_Comm_rank(*my_comm_out,&rank_my_comm_out);
    MPI_Comm_split(my_comm_in,rank_my_comm_out,iam,equal_rank_comm);

    int loc_n_nonzeros=matlocal->n_nonzeros;
    int* dist_loc_n_nonzeros=new int[nprocs];
    MPI_Allgather(&loc_n_nonzeros,1,MPI_INT,dist_loc_n_nonzeros,1,MPI_INT,my_comm_in);
    int* displc_n_nonzeros=new int[nmat+1]();
    for (int irecv=0;irecv<nmat;irecv++) {
        int from_node = (irecv+iam*nmat)%nprocs;
        displc_n_nonzeros[irecv+1]=displc_n_nonzeros[irecv]+dist_loc_n_nonzeros[from_node];
    }

    int loc_size=matlocal->size;
    int* dist_loc_size=new int[nprocs];
    MPI_Allgather(&loc_size,1,MPI_INT,dist_loc_size,1,MPI_INT,my_comm_in);
    int* displc_size=new int[nmat+1]();
    for (int irecv=0;irecv<nmat;irecv++) {
        int from_node = (irecv+iam*nmat)%nprocs;
        displc_size[irecv+1]=displc_size[irecv]+dist_loc_size[from_node];
    }

    MPI_Request reqs;
    MPI_Status status;
    MPI_Irecv(&first_row,1,MPI_INT,(iam*nmat)%(nmat*size_my_comm_out),tags::first_row,my_comm_in,&reqs);
    if (!(iam%nmat))
        for (int isend=0;isend<nmat;isend++)
            MPI_Send(&matlocal->first_row,1,MPI_INT,iam/nmat+isend*size_my_comm_out,tags::first_row,my_comm_in);
    MPI_Wait(&reqs,&status);

    n_nonzeros = displc_n_nonzeros[nmat];
    size_tot   = matlocal->size_tot;
    size       = displc_size[nmat];
    findx      = matlocal->findx;
    type       = 0;
    nnz        = new T[max(n_nonzeros,1)];
    index_i    = new int[max(size,1)];
    index_j    = new int[max(n_nonzeros,1)];
    edge_i     = new int[size+1];
    diag_pos   = new int[max(size,1)];

    MPI_Status* nnz_status=new MPI_Status[nmat];
    MPI_Status* ind_j_status=new MPI_Status[nmat];
    MPI_Status* ind_i_status=new MPI_Status[nmat];
    MPI_Request* nnz_reqs=new MPI_Request[nmat];
    MPI_Request* ind_j_reqs=new MPI_Request[nmat];
    MPI_Request* ind_i_reqs=new MPI_Request[nmat];

    for (int irecv=0;irecv<nmat;irecv++) {
        int from_node = (irecv+iam*nmat)%nprocs;
        MPI_Irecv(&nnz[displc_n_nonzeros[irecv]],dist_loc_n_nonzeros[from_node],give_me_MPI_datatype(nnz),from_node,tags::nnz,my_comm_in,&nnz_reqs[irecv]);
        MPI_Irecv(&index_j[displc_n_nonzeros[irecv]],dist_loc_n_nonzeros[from_node],MPI_INT,from_node,tags::index_j,my_comm_in,&ind_j_reqs[irecv]);
        MPI_Irecv(&index_i[displc_size[irecv]],dist_loc_size[from_node],MPI_INT,from_node,tags::index_i,my_comm_in,&ind_i_reqs[irecv]);
    }

    for (int isend=0;isend<nmat;isend++) {
        int to_node = iam/nmat+isend*size_my_comm_out;
        MPI_Send(matlocal->nnz,matlocal->n_nonzeros,give_me_MPI_datatype(nnz),to_node,tags::nnz,my_comm_in);
        MPI_Send(matlocal->index_j,matlocal->n_nonzeros,MPI_INT,to_node,tags::index_j,my_comm_in);
        MPI_Send(matlocal->index_i,matlocal->size,MPI_INT,to_node,tags::index_i,my_comm_in);
    }

    MPI_Waitall(nmat,nnz_reqs,nnz_status);
    MPI_Waitall(nmat,ind_j_reqs,ind_j_status);
    MPI_Waitall(nmat,ind_i_reqs,ind_i_status);

    delete[] nnz_status;
    delete[] ind_j_status;
    delete[] ind_i_status;
    delete[] nnz_reqs;
    delete[] ind_j_reqs;
    delete[] ind_i_reqs;

    delete[] dist_loc_n_nonzeros;
    delete[] displc_n_nonzeros;
    delete[] dist_loc_size;
    delete[] displc_size;

    get_row_edge();
    get_diag_pos();
}

/******************** written by SAB ********************/

template <class T>
inline TCSR<T>::TCSR(TCSR<T> *mat,TCSR<T> *dep_mat,int rootrank,MPI_Comm matrix_comm)
{
    int iam,nprocs;
    MPI_Comm_size(matrix_comm,&nprocs);
    MPI_Comm_rank(matrix_comm,&iam);

    size_tot   = mat->size_tot;
    type       = mat->type;
    findx      = mat->findx;

    if (iam==rootrank) {
        if (size_tot != dep_mat->size_tot ) throw CSR_Exception(__LINE__,__FILE__);
        if (type     != dep_mat->type     ) throw CSR_Exception(__LINE__,__FILE__);
        if (findx    != dep_mat->findx    ) throw CSR_Exception(__LINE__,__FILE__);
    }

    size       = mat->size;
    first_row  = mat->first_row;
    
    int *first_rows = new int[nprocs+1];
    MPI_Allgather(&first_row,1,MPI_INT,first_rows,1,MPI_INT,matrix_comm);
    first_rows[nprocs]=size_tot;
    int *sizes = new int[nprocs];
    for (int irank=0;irank<nprocs;irank++) {
        sizes[irank]=first_rows[irank+1]-first_rows[irank];
    }
    if (size!=sizes[iam]) throw CSR_Exception(__LINE__,__FILE__);

    int *n_nonzeros_array = new int[nprocs];
    if (iam==rootrank) {
        for (int irank=0;irank<nprocs;irank++) {
            n_nonzeros_array[irank]=dep_mat->edge_i[first_rows[irank+1]]-dep_mat->edge_i[first_rows[irank]];
        }
    }
    MPI_Bcast(n_nonzeros_array,nprocs,MPI_INT,rootrank,matrix_comm);
    int *n_nonzeros_disp = new int[nprocs+1];
    n_nonzeros_disp[0]=0;
    for (int irank=0;irank<nprocs;irank++) {
        n_nonzeros_disp[irank+1]=n_nonzeros_disp[irank]+n_nonzeros_array[irank];
    }
    n_nonzeros=n_nonzeros_array[iam];

    nnz        = new T[max(n_nonzeros,1)];
    index_i    = new int[size];
    index_j    = new int[max(n_nonzeros,1)];
    edge_i     = new int[size+1];
    diag_pos   = new int[size];

    init_variable(index_i,size);
    if (!n_nonzeros) {
        init_variable(nnz,1);
        init_variable(index_j,1);
    }

    T *nnz_ptr = NULL;
    int *indj_ptr = NULL;
    int *indi_ptr = NULL;
    if (iam==rootrank) {
        nnz_ptr=dep_mat->nnz;
        indj_ptr=dep_mat->index_j;
        indi_ptr=dep_mat->index_i;
    }

    MPI_Scatterv(nnz_ptr,n_nonzeros_array,n_nonzeros_disp,give_me_MPI_datatype(nnz),nnz,n_nonzeros,give_me_MPI_datatype(nnz),rootrank,matrix_comm);
    MPI_Scatterv(indj_ptr,n_nonzeros_array,n_nonzeros_disp,MPI_INT,index_j,n_nonzeros,MPI_INT,rootrank,matrix_comm);
    MPI_Scatterv(indi_ptr,sizes,first_rows,MPI_INT,index_i,size,MPI_INT,rootrank,matrix_comm);

    delete[] first_rows;
    delete[] sizes;
    delete[] n_nonzeros_array;
    delete[] n_nonzeros_disp;

    get_row_edge();
    get_diag_pos();
}

/************************************************************************************************/

template <>
inline TCSR<CPX>::TCSR(const char *filename)
{

    int i, ind_i;
    FILE *F = fopen(filename, "r");
    double real_i, imag_i;

    fscanf(F,"%i",&size);
    fscanf(F,"%i",&n_nonzeros);
    fscanf(F,"%i",&findx);

    nnz        = new CPX[n_nonzeros];
    index_i    = new int[size];
    index_j    = new int[n_nonzeros];
    edge_i     = new int[size+1];
    diag_pos   = new int[size];

    size_tot   = size;
    type       = 0;
    first_row  = 0;

    for(i=0;i<size;i++) index_i[i] = 0;

    for(i=0;i<n_nonzeros;i++){
        fscanf(F,"%i",&ind_i);
        fscanf(F,"%i",index_j+i);
        fscanf(F,"%lg",&real_i);
        fscanf(F,"%lg",&imag_i);
        nnz[i] = CPX(real_i,imag_i);

        if(ind_i==index_j[i]){
            diag_pos[ind_i-findx] = i;
        }
        index_i[ind_i-findx]++;
    }
    get_row_edge();
    fclose(F);
}

/***********************************************************************************************/

template <>
inline TCSR<double>::TCSR(const char *filename)
{

    int i,ind_i;
    FILE *F = fopen(filename,"r");
   
    fscanf(F,"%i",&size);
    fscanf(F,"%i",&n_nonzeros);
    fscanf(F,"%i",&findx);
    
    nnz        = new double[n_nonzeros];
    index_i    = new int[size];                                           
    index_j    = new int[n_nonzeros];                                       
    edge_i     = new int[size+1];                                         
    diag_pos   = new int[size];                                           

    size_tot   = size;                                                    
    type       = 0;                                        
    first_row  = 0;                                                    

    for(i=0;i<size;i++) index_i[i] = 0;
   
    for(i=0;i<n_nonzeros;i++){
        fscanf(F,"%i",&ind_i);
	      fscanf(F,"%i",index_j+i);
	      fscanf(F,"%lg",nnz+i);
	      if(ind_i==index_j[i]){
	          diag_pos[ind_i-findx] = i;
	      }
	      index_i[ind_i-findx]++;
    }
    get_row_edge();
    fclose(F);

}

/***********************************************************************************************/

template <>
inline TCSR<CPX>::TCSR(char *filename,int mrank,int msize)
{
}

/***********************************************************************************************/

template <>
inline TCSR<double>::TCSR(char *filename,int msize,int mrank)
{

    int i,i0,ind_i,var_i,n_tot;
    double var_d;
    FILE *F = fopen(filename,"r");

    fscanf(F,"%i",&size_tot);
    fscanf(F,"%i",&n_tot);
    fscanf(F,"%i",&findx);

    cpu_distribute(&first_row,&size,size_tot,msize,mrank);

    n_nonzeros = 0;
    for(i=0;i<n_tot;i++){
        fscanf(F,"%i",&ind_i);
	fscanf(F,"%i",&var_i);
	fscanf(F,"%lg",&var_d);
	if(((ind_i-findx)>=first_row)&&((ind_i-findx)<first_row+size)){
	    n_nonzeros++;
	}
    }
    fclose(F);
    
    nnz        = new double[n_nonzeros];
    index_i    = new int[size];                                           
    index_j    = new int[n_nonzeros];                                       
    edge_i     = new int[size+1];                                         
    diag_pos   = new int[size];                                           
    
    type       = 0;                                                                                         

    for(i=0;i<size;i++) index_i[i] = 0;

    F = fopen(filename,"r");

    fscanf(F,"%i",&var_i);
    fscanf(F,"%i",&var_i);
    fscanf(F,"%i",&var_i);
   
    i0=0;
    for(i=0;i<n_tot;i++){
        fscanf(F,"%i",&ind_i);
	fscanf(F,"%i",&var_i);
	fscanf(F,"%lg",&var_d);
	if(((ind_i-findx)>=first_row)&&((ind_i-findx)<first_row+size)){
	    index_j[i0] = var_i;
	    nnz[i0]     = var_d;
	    if(ind_i==index_j[i0]){
	        diag_pos[ind_i-findx-first_row] = i0;
	    }
	    index_i[ind_i-findx-first_row]++;
	    i0++;
	}
    }
    get_row_edge();
    fclose(F);

}

/***********************************************************************************************/

template <class T>
TCSR<T>::~TCSR()
{

    delete[] nnz;
    delete[] index_i;
    delete[] index_j;
    delete[] edge_i;
    delete[] diag_pos;

}

/************************************************************************************************/

template <class T>
void TCSR<T>::change_contain(TCSR<T> *mat)
{
    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    n_nonzeros = mat->n_nonzeros;
    findx      = mat->findx;
    first_row  = mat->first_row;

    c_tcopy(n_nonzeros,mat->nnz,1,nnz,1);
    c_icopy(size,mat->index_i,1,index_i,1);
    c_icopy(n_nonzeros,mat->index_j,1,index_j,1);
    c_icopy(size+1,mat->edge_i,1,edge_i,1);
    c_icopy(size,mat->diag_pos,1,diag_pos,1);
}

/************************************************************************************************/

template <class T>
void TCSR<T>::copy_index(TCSR<double> *mat)
{
    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    n_nonzeros = mat->n_nonzeros;
    findx      = mat->findx;
    first_row  = mat->first_row;

    c_icopy(size,mat->index_i,1,index_i,1);
    c_icopy(n_nonzeros,mat->index_j,1,index_j,1);
    c_icopy(size+1,mat->edge_i,1,edge_i,1);
    c_icopy(size,mat->diag_pos,1,diag_pos,1);
}

/************************************************************************************************/

template <class T>
void TCSR<T>::copy_index(TCSR<CPX> *mat)
{
    size       = mat->size;
    size_tot   = mat->size_tot;
    type       = mat->type;
    n_nonzeros = mat->n_nonzeros;
    findx      = mat->findx;
    first_row  = mat->first_row;

    c_icopy(size,mat->index_i,1,index_i,1);
    c_icopy(n_nonzeros,mat->index_j,1,index_j,1);
    c_icopy(size+1,mat->edge_i,1,edge_i,1);
    c_icopy(size,mat->diag_pos,1,diag_pos,1);
}
/************************************************************************************************/
/*
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
*/
/************************************************************************************************/
/*
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
*/
/************************************************************************************************/
/*
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
*/
/************************************************************************************************/
/*
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
*/
/************************************************************************************************/

template <class T>
void TCSR<T>::cpu_distribute(int *Istart,int *Nlocal,int N,int msize,int mrank)
{
    int max_rank,NCPU;

    NCPU             = int(ceil((double)N/msize));
    max_rank         = N-(NCPU-1)*msize;
    *Nlocal          = NCPU-1;

    if(mrank<max_rank){*Nlocal = NCPU;}

    *Istart          = min(mrank,max_rank)*NCPU+max(mrank-max_rank,0)*(NCPU-1);
}

/************************************************************************************************/

template <class T>
void TCSR<T>::update_diag(T *diag,T factor)
{
    int i;

    for(i=0;i<size;i++){
        if(diag_pos[i]>=0){
	    nnz[diag_pos[i]] = nnz[diag_pos[i]]+factor*diag[first_row+i];
	}else{
	    printf("No diagonal at index %i (%i)\n",i,i+first_row);
	    abort();
	}
    }
    
}

/************************************************************************************************/

template <class T>
void TCSR<T>::update_loc_diag(T *diag,T factor)
{
    int i;

    for(i=0;i<size;i++){
        if(diag_pos[i]>=0){
	    nnz[diag_pos[i]] = nnz[diag_pos[i]]+factor*diag[i];
	}else{
	    printf("No diagonal at index %i\n",i);
	    abort();
	}
    }
    
}

/************************************************************************************************/

template <class T>
void TCSR<T>::update_diag_single(T diag_element)
{
    int i;

    for(i=0;i<size;i++){
        if(diag_pos[i]>=0){
	    nnz[diag_pos[i]] = nnz[diag_pos[i]]+diag_element;
	}else{
	    printf("No diagonal at index %i\n",i);
	    abort();
	}
    }
    
}

/************************************************************************************************/

template <class T>
void TCSR<T>::update_diag_single(T diag_element,int start,int stop)
{
    int i;

    for(i=max(first_row,start);i<min(first_row+size,stop);i++){
        if(diag_pos[i-first_row]>=0){
	    nnz[diag_pos[i-first_row]] = nnz[diag_pos[i-first_row]]+diag_element;
	}else{
	  printf("No diagonal at index %i (%i)\n",i-first_row,i);
	    abort();
	}
    }
    
}

/************************************************************************************************/

template <class T>
void TCSR<T>::update_diag_pack(T *diag,T factor,int N)
{
    int i,j;

    for(i=0;i<size/N;i++){
        for(j=0;j<N;j++){
	    if(diag_pos[i*N+j]>=0){
	        nnz[diag_pos[i*N+j]] = nnz[diag_pos[i*N+j]]+factor*diag[first_row/N+i];
	    }else{
	        printf("No diagonal at index %i (%i)\n",i*N+j,i*N+j+first_row);
		abort();
	    }
        }
    }
    
}

/************************************************************************************************/

template <>
inline void TCSR<double>::update_rdiag_pack(double *diag,int start,int stop,double factor,\
					    int N,int *n_of_elements)
{
    int i,j,ind;

    for(i=start;i<stop;i++){
        ind = get_msize(start,i-1,N,n_of_elements);
        for(j=0;j<(n_of_elements[i+1]-n_of_elements[i])*N/10;j++){
	    if(diag_pos[ind+j]>=0){
	        nnz[diag_pos[ind+j]] = nnz[diag_pos[ind+j]]+factor*diag[i];
	    }else{
	        printf("No diagonal at index %i (%i)\n",ind+j,i);
		abort();
	    }
        }
    }
    
}

/************************************************************************************************/

template <>
inline void TCSR<CPX>::update_rdiag_pack(double *diag,int start,int stop,double factor,int N,
					 int *n_of_elements)
{
    int i,j,ind;

    for(i=start;i<stop;i++){
        ind = get_msize(start,i-1,N,n_of_elements);
        for(j=0;j<(n_of_elements[i+1]-n_of_elements[i])*N/10;j++){
	    if(diag_pos[ind+j]>=0){
	        nnz[diag_pos[ind+j]] = nnz[diag_pos[ind+j]]+CPX(factor*diag[i],0.0);
	    }else{
	        printf("No diagonal at index %i (%i)\n",ind+j,i);
		abort();
	    }
        }
    }
    
}

/************************************************************************************************/

template <class T>
void TCSR<T>::get_row_edge()
{
    int i;

    edge_i[0] = findx;
  
    for(i=0;i<size;i++){
        edge_i[i+1] = edge_i[i]+index_i[i];
    }
    
}

/************************************************************************************************/

template <class T>
void TCSR<T>::get_diag_pos()
{
    int i,j;
    int found_diag;

    edge_i[0] = findx;
  
    for(i=0;i<size;i++){
        found_diag = 0;
        for(j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
	    if((index_j[j]-findx)==(i+first_row)){
	        diag_pos[i] = j;
		found_diag  = 1;
	    }
	}
	if(!found_diag){
	    diag_pos[i] = -1;
	}
    }
    
}

/************************************************************************************************/

template <class T>
void TCSR<T>::scale_element(int N, int* dat_ind, T scale_fact)
{
    int i;

    for(i=0;i<N;i++){
        nnz[dat_ind[i]] = scale_fact*nnz[dat_ind[i]];
    }
}


/************************************************************************************************/

template <>
inline void TCSR<double>::mat_vec_mult(double *in, double *out, int nrhs)
{
    int i,j,until;
    double *vec = new double[nrhs];
    
    for(i=0;i<size;i++){
        init_variable(vec,nrhs);
        until = edge_i[i+1]-edge_i[i];
        for(j=0;j<until;j++){
            c_daxpy(nrhs,nnz[edge_i[i]-findx+j],&in[index_j[edge_i[i]-findx+j]-findx],\
                    size_tot,vec,1);
        }
        c_dcopy(nrhs,vec,1,&out[i],size);
    }
    delete[] vec;

}

/************************************************************************************************/

template <>
inline void TCSR<CPX>::mat_vec_mult(CPX *in, CPX *out, int nrhs)
{
    int i,j,until;
    CPX *vec = new CPX[nrhs];
    
    for(i=0;i<size;i++){
        init_variable(vec,nrhs);
        until = edge_i[i+1]-edge_i[i];
        for(j=0;j<until;j++){
            c_zaxpy(nrhs,nnz[edge_i[i]-findx+j],&in[index_j[edge_i[i]-findx+j]-findx],\
                    size_tot,vec,1);
        }
        c_zcopy(nrhs,vec,1,&out[i],size);
    }
    delete[] vec;

}

/************************************************************************************************/

template <>
inline void TCSR<CPX>::mat_vec_mult(CPX *in, CPX *out, int nrhs, int NR)
{
    int i,j,until;
    CPX *vec = new CPX[nrhs];
    
    for(i=0;i<size;i++){
        init_variable(vec,nrhs);
        until = edge_i[i+1]-edge_i[i];
        for(j=0;j<until;j++){
            c_zaxpy(nrhs,nnz[edge_i[i]-findx+j],&in[index_j[edge_i[i]-findx+j]-findx],\
                    NR,vec,1);
        }
        c_zcopy(nrhs,vec,1,&out[i],size);
    }
    delete[] vec;

}

/************************************************************************************************/

/*This function calculate sparse_mat*in with the help of the transpose of in tr_in (more
  efficient when blas is used)*/
template <>
inline void TCSR<CPX>::trans_mat_vec_mult(CPX *tr_in, CPX *out, int nrhs, int NR)
{
    int i,j;
    CPX *vec = new CPX[nrhs];
    
    for(i=0;i<size;i++){
        init_variable(vec,nrhs);
        for(j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            c_zaxpy(nrhs,nnz[j],&tr_in[nrhs*(index_j[j]-findx)],1,vec,1);
        }
        c_zcopy(nrhs,vec,1,&out[i],size);
    }

    delete[] vec;

}

/************************************************************************************************/

/*This function calculate sparse_mat*in with the help of the transpose of in tr_in (more
  efficient when blas is used)*/
template <>
inline void TCSR<CPX>::trans_mat_vec_mult_omp(CPX *tr_in, CPX *out, int nrhs, int NR)
{
    int nthr=omp_get_max_threads();
    int i;//DECLARE J PRIVATE TO THREAD
    CPX *vec = new CPX[nrhs*nthr];

#pragma omp parallel for schedule(dynamic)
    for(i=0;i<size;i++){
        int ithr=omp_get_thread_num();
        init_variable(&vec[nrhs*ithr],nrhs);
        for(int j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            c_zaxpy(nrhs,nnz[j],&tr_in[nrhs*(index_j[j]-findx)],1,&vec[nrhs*ithr],1);
        }
        c_zcopy(nrhs,&vec[nrhs*ithr],1,&out[i],size);
    }

    delete[] vec;

}

/************************************************************************************************/

/*This function calculate sparse_mat*in with the help of the transpose of in tr_in (more
  efficient when blas is used)*/
template <>
inline void TCSR<CPX>::trans_mat_vec_mult(CPX *tr_in, CPX *out, int nrhs, int NR,int NLine)
{
    int i,j;
    CPX *vec = new CPX[nrhs];
    
    for(i=0;i<NLine;i++){
        init_variable(vec,nrhs);
        for(j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            c_zaxpy(nrhs,nnz[j],&tr_in[nrhs*(index_j[j]-findx)],1,vec,1);
        }
        c_zcopy(nrhs,vec,1,&out[i],size);
    }

    delete[] vec;

}

/************************************************************************************************/

template <>
inline void TCSR<double>::mat_vec_mult(double *in, double *out, int nrhs, int NR)
{
    int i,j;
    double *vec = new double[nrhs];
    
    for(i=0;i<size;i++){
        init_variable(vec,nrhs);
        for(j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            c_daxpy(nrhs,nnz[j],&in[index_j[j]-findx],NR,vec,1);
        }
        c_dcopy(nrhs,vec,1,&out[i],size);
    }

    delete[] vec;

}

/************************************************************************************************/

/*This function calculate sparse_mat*in with the help of the transpose of in tr_in (more
  efficient when blas is used)*/
template <>
inline void TCSR<double>::trans_mat_vec_mult(double *tr_in, double *out, int nrhs, int NR)
{
    int i,j;
    double *vec = new double[nrhs];
    
    for(i=0;i<size;i++){
        init_variable(vec,nrhs);
        for(j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            c_daxpy(nrhs,nnz[j],&tr_in[nrhs*(index_j[j]-findx)],1,vec,1);
        }
        c_dcopy(nrhs,vec,1,&out[i],size);
    }

    delete[] vec;

}

/************************************************************************************************/

/*This function calculate sparse_mat*in with the help of the transpose of in tr_in (more
  efficient when blas is used)*/
template <>
inline void TCSR<double>::trans_mat_vec_mult_omp(double *tr_in, double *out, int nrhs, int NR)
{
    int nthr=omp_get_max_threads();
    int i;//DECLARE J PRIVATE TO THREAD
    double *vec = new double[nrhs*nthr];

#pragma omp parallel for schedule(dynamic)
    for(i=0;i<size;i++){
        int ithr=omp_get_thread_num();
        init_variable(&vec[nrhs*ithr],nrhs);
        for(int j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            c_daxpy(nrhs,nnz[j],&tr_in[nrhs*(index_j[j]-findx)],1,&vec[nrhs*ithr],1);
        }
        c_dcopy(nrhs,&vec[nrhs*ithr],1,&out[i],size);
    }

    delete[] vec;

}

/************************************************************************************************/

/*This function calculate sparse_mat*in with the help of the transpose of in tr_in (more
  efficient when blas is used)*/
template <>
inline void TCSR<double>::trans_mat_vec_mult(double *tr_in, double *out, int nrhs, int NR, int NLine)
{
    int i,j;
    double *vec = new double[nrhs];
    
    for(i=0;i<NLine;i++){
        init_variable(vec,nrhs);
        for(j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            c_daxpy(nrhs,nnz[j],&tr_in[nrhs*(index_j[j]-findx)],1,vec,1);
        }
        c_dcopy(nrhs,vec,1,&out[i],size);
    }

    delete[] vec;

}

/************************************************************************************************/

template <class T>
inline void TCSR<T>::spec_mat_vec_mult(T *in, T *out, int in_fst_row, int *low, int *high, \
				       int NS, int fac)
{
    int i,j,start,stop,IS,until,in_index;

    init_variable(out,size);

    for(IS=0;IS<NS-1;IS++){
        start = max(low[IS]*fac-first_row,0);
        stop  = min((high[IS]+1)*fac-first_row,size);
        for(i=start;i<stop;i++){
            until = edge_i[i+1]-edge_i[i];
            for(j=0;j<until;j++){
                if((index_j[edge_i[i]-findx+j]-findx)>=low[IS+1]*fac){
		    in_index = index_j[edge_i[i]-findx+j]-findx-in_fst_row;
		    out[i]   = out[i]+nnz[edge_i[i]-findx+j]*in[in_index];
                }
            }
        }
    }
}

/************************************************************************************************/

template <class T>
inline void TCSR<T>::spec_mat_vec_mult(T *in, T *out, int in_fst_row, int *low, int *high, \
				       int NS, int fac, int *size_of_el)
{
    int i,j,start,stop,IS,until,in_index;

    init_variable(out,size);

    for(IS=0;IS<NS-1;IS++){
        start = max(get_msize(0,low[IS]-1,fac,size_of_el)-first_row,0);
	stop  = min(get_msize(0,high[IS],fac,size_of_el)-first_row,size);
	for(i=start;i<stop;i++){
	    until = edge_i[i+1]-edge_i[i];
	    for(j=0;j<until;j++){
	        if((index_j[edge_i[i]-findx+j]-findx)>=get_msize(0,low[IS+1]-1,fac,size_of_el)){
		    in_index = index_j[edge_i[i]-findx+j]-findx-in_fst_row;
		    out[i]   = out[i]+nnz[edge_i[i]-findx+j]*in[in_index];
		}
	    }
	}
    }
}

/************************************************************************************************/
/*
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
*/
/************************************************************************************************/
/*
template <>
void TCSR<CPX>::sparse_to_cmp_full(CPX *B,int nrow,int ncol)
{
}
*/
/************************************************************************************************/
/*
template <class T>
void TCSR<T>::sparse_to_full(T *B,int nrow,int ncol)
{
    int i,j,index;

    init_variable(B,nrow*ncol);
    
    for(i=0;i<size;i++){
        for(j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            index             = index_j[j]-findx;
            B[i+index*nrow]   = nnz[j];
        }
    }
}
*/
/*SAB********************************************************************************************/

template <class T>
void TCSR<T>::sparse_to_full(T *B,int nrow,int ncol)
{
    int i,j,e;

    init_variable(B,nrow*ncol);
    for(i=0;i<min(nrow,size);i++){
        for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
            if ((j=index_j[e]-findx)<ncol){
                B[i+j*nrow] = nnz[e];
            }
        }
    }
}

/*SAB********************************************************************************************/

template <class T>
void TCSR<T>::sparse_to_full(T *B,int nrow,int ncol,int ipos, int jpos)
{
    int i,j,e;

    init_variable(B,nrow*ncol);
    for(i=ipos;i<min(ipos+nrow,size);i++){
        for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
            j=index_j[e]-findx-jpos;
            if (j>=0 && j<ncol){
                B[i-ipos+j*nrow] = nnz[e];
            }
        }
    }
}

/************************************************************************************************/

template <>
inline void TCSR<double>::mat_vec_mult_add(double *in, double *out, int nrhs, double fac)
{
    int i,j,until;
    double *vec = new double[nrhs];
    
    for(i=0;i<size;i++){
        init_variable(vec,nrhs);
        until = edge_i[i+1]-edge_i[i];
        for(j=0;j<until;j++){
            c_daxpy(nrhs,nnz[edge_i[i]-findx+j],&in[index_j[edge_i[i]-findx+j]-findx],\
                    size_tot,vec,1);
        }
        c_daxpy(nrhs,fac,vec,1,&out[i],size);
    }
    delete[] vec;

}

/************************************************************************************************/

template <>
inline void TCSR<CPX>::mat_vec_mult_add(CPX *in, CPX *out, int nrhs, CPX fac)
{
    int i,j,until;
    CPX *vec = new CPX[nrhs];
    
    for(i=0;i<size;i++){
        init_variable(vec,nrhs);
        until = edge_i[i+1]-edge_i[i];
        for(j=0;j<until;j++){
            c_zaxpy(nrhs,nnz[edge_i[i]-findx+j],&in[index_j[edge_i[i]-findx+j]-findx],\
                    size_tot,vec,1);
        }
        c_zaxpy(nrhs,fac,vec,1,&out[i],size);
    }
    delete[] vec;

}

/************************************************************************************************/

template <class T>
void TCSR<T>::shift_resize(int istart, int ilength, int jstart, int jlength)
{
    if (first_row)           throw CSR_Exception(__LINE__,__FILE__);
    if (istart<0)            throw CSR_Exception(__LINE__,__FILE__);
    if (istart+ilength>size) throw CSR_Exception(__LINE__,__FILE__);

    for (int e=0;e<n_nonzeros;e++){
        index_j[e]-=jstart;
        if (index_j[e]-findx<0)        throw CSR_Exception(__LINE__,__FILE__);
        if (index_j[e]-findx>=jlength) throw CSR_Exception(__LINE__,__FILE__);
    }

    for (int i=0;i<ilength;i++){
        index_i[i]=index_i[i+istart];
    }

    size = ilength;
    size_tot = ilength;

    get_row_edge();
    get_diag_pos();
}

/************************************************************************************************/

template <class T>
void TCSR<T>::replicate_cell(int ndof,int bandw,int contact,int strinum,MPI_Comm my_comm)
{
    int i,e;
    T* linefull = new T[(2*bandw+1)*ndof];
    if (contact==1) {
        TCSR<T> *Pscut = new TCSR<T>(this,strinum*ndof,ndof,0,size_tot*ndof);
        TCSR<T> *Pcell = new TCSR<T>(Pscut,my_comm);
        delete Pscut;
        for (int i_elem=strinum*ndof;i_elem<(strinum+1)*ndof;i_elem++) {
            for (e=Pcell->edge_i[i_elem]-Pcell->findx;e<Pcell->edge_i[i_elem+1]-Pcell->findx;e++) {
                Pcell->index_j[e]=(Pcell->index_j[e]+(bandw-strinum)*ndof+size_tot)%size_tot;
            }
        }
        Pcell->shift_resize(strinum*ndof,ndof,0,(2*bandw+1)*ndof);
        for (i=0;i<size;i++) {
            int i_band=(i+first_row)/ndof;
            int i_elem=(i+first_row)%ndof;
            if (i_band<bandw) {
                init_variable(linefull,(2*bandw+1)*ndof);
                for (e=Pcell->edge_i[i_elem]-Pcell->findx;e<Pcell->edge_i[i_elem+1]-Pcell->findx;e++) {
                    linefull[Pcell->index_j[e]-Pcell->findx]=Pcell->nnz[e];
                }
                for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
                    nnz[e]=linefull[(index_j[e]-findx+(bandw-i_band)*ndof+size_tot)%size_tot];
                }
            }
        }
        delete Pcell;
    } else if (contact==2) {
        TCSR<T> *Pscut = new TCSR<T>(this,size_tot-(strinum+1)*ndof,ndof,0,size_tot*ndof);
        TCSR<T> *Pcell = new TCSR<T>(Pscut,my_comm);
        delete Pscut;
        for (int i_elem=size_tot-(strinum+1)*ndof;i_elem<size_tot-strinum*ndof;i_elem++) {
            for (e=Pcell->edge_i[i_elem]-Pcell->findx;e<Pcell->edge_i[i_elem+1]-Pcell->findx;e++) {
                Pcell->index_j[e]=(Pcell->index_j[e]-(bandw-strinum)*ndof+size_tot)%size_tot;
            }
        }
        Pcell->shift_resize(size_tot-(strinum+1)*ndof,ndof,size_tot-(2*bandw+1)*ndof,(2*bandw+1)*ndof);
        for (i=0;i<size;i++) {
            int i_band=(i+first_row)/ndof;
            int i_elem=(i+first_row)%ndof;
            if (i_band>size_tot/ndof-bandw-1) {
                init_variable(linefull,(2*bandw+1)*ndof);
                for (e=Pcell->edge_i[i_elem]-Pcell->findx;e<Pcell->edge_i[i_elem+1]-Pcell->findx;e++) {
                    linefull[Pcell->index_j[e]-Pcell->findx]=Pcell->nnz[e];
                }
                for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
                    nnz[e]=linefull[(index_j[e]-findx+(bandw-i_band)*ndof+size_tot)%size_tot];
                }
            }
        }
        delete Pcell;
    }
    delete[] linefull;
}

/************************************************************************************************/

template <class T>
void TCSR<T>::fill_pbc_block(TCSR<T> *A,int ndof,int bandw,int contact,MPI_Comm my_comm)
{
    if (size      != A->size      ) throw CSR_Exception(__LINE__,__FILE__);
    if (size_tot  != A->size_tot  ) throw CSR_Exception(__LINE__,__FILE__);
    if (type      != A->type      ) throw CSR_Exception(__LINE__,__FILE__);
    if (first_row != A->first_row ) throw CSR_Exception(__LINE__,__FILE__);

    int i,e;
    T* linefull = new T[size_tot];
    if (contact==1) {
        TCSR<T> **Pcell = new TCSR<T>*[bandw];
        for (int i_band=0;i_band<bandw;i_band++) {
            TCSR<T> *Pscut = new TCSR<T>(A,(i_band+1)*ndof,(bandw-i_band)*ndof,0,ndof);
            Pcell[i_band] = new TCSR<T>(Pscut,my_comm);
            delete Pscut;
            Pcell[i_band]->shift_resize((i_band+1)*ndof,bandw*ndof,-A->size_tot+(i_band+1)*ndof,size_tot);
        }
        TCSR<T> *Pblock = new TCSR<T>(bandw,NULL,Pcell,size_tot);
        for (int i_band=0;i_band<bandw;i_band++) {
            delete Pcell[i_band];
        }
        delete[] Pcell;
        for (i=0;i<Pblock->size;i++) {
            init_variable(linefull,size_tot);
            for (e=Pblock->edge_i[i]-Pblock->findx;e<Pblock->edge_i[i+1]-Pblock->findx;e++) {
                linefull[Pblock->index_j[e]-Pblock->findx]=Pblock->nnz[e];
            }
            for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
                nnz[e]=linefull[index_j[e]-findx];
            }
        }
        delete Pblock;
    } else if (contact==2) {
        TCSR<T> **Pcell = new TCSR<T>*[bandw];
        for (int i_band=0;i_band<bandw;i_band++) {
            TCSR<T> *Pscut = new TCSR<T>(A,size_tot-(bandw+1+i_band)*ndof,bandw*ndof,size_tot-ndof,ndof);
            Pcell[i_band] = new TCSR<T>(Pscut,my_comm);
            delete Pscut;
            Pcell[i_band]->shift_resize(size_tot-(bandw+1+i_band)*ndof,bandw*ndof,size_tot-(i_band+1)*ndof,size_tot);
        }
        TCSR<T> *Pblock = new TCSR<T>(bandw,NULL,Pcell,size_tot);
        for (int i_band=0;i_band<bandw;i_band++) {
            delete Pcell[i_band];
        }
        delete[] Pcell;
        for (i=0;i<Pblock->size;i++) {
            init_variable(linefull,size_tot);
            for (e=Pblock->edge_i[i]-Pblock->findx;e<Pblock->edge_i[i+1]-Pblock->findx;e++) {
                linefull[Pblock->index_j[e]-Pblock->findx]=Pblock->nnz[e];
            }
            for (e=edge_i[i+size_tot-bandw*ndof]-findx;e<edge_i[i+size_tot-bandw*ndof+1]-findx;e++) {
                nnz[e]=linefull[index_j[e]-findx];
            }
        }
        delete Pblock;
    }
    delete[] linefull;
}

/************************************************************************************************/

template <class T>
int TCSR<T>::additionalentries(TCSR<T> *mat)
{
    int additional=0;
    int i,edgy1,edgy2,foundit;
    for (i=0;i<size;i++) {
        for (edgy1=edge_i[i]-findx;edgy1<edge_i[i+1]-findx;edgy1++) {
            foundit=0;
            for (edgy2=mat->edge_i[i]-mat->findx;edgy2<mat->edge_i[i+1]-mat->findx;edgy2++) {
                if (index_j[edgy1]==mat->index_j[edgy2]) {
                    foundit++;
                }
            }
            if (!foundit) {
                additional++;
            } else if (foundit>1) {
                return -1;
            }
        }
    }
    return additional;
}

/************************************************************************************************/

template <class T>
void TCSR<T>::add(TCSR<T> *A,T factor)
{
    if (size      != A->size      ) throw CSR_Exception(__LINE__,__FILE__);
    if (size_tot  != A->size_tot  ) throw CSR_Exception(__LINE__,__FILE__);
    if (type      != A->type      ) throw CSR_Exception(__LINE__,__FILE__);
    if (first_row != A->first_row ) throw CSR_Exception(__LINE__,__FILE__);

    int i,e;
    T *linefull=new T[size_tot];

    for (i=0;i<size;i++) {
        init_variable(linefull,size_tot);
        for (e=A->edge_i[i]-A->findx;e<A->edge_i[i+1]-A->findx;e++) {
            linefull[A->index_j[e]-A->findx]=factor*A->nnz[e];
        }
        for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            nnz[e]+=linefull[index_j[e]-findx];
        }
    }
    delete[] linefull;
}

/************************************************************************************************/

template <>
inline void TCSR<double>::add_real(TCSR<CPX> *A,CPX factor)
{
    if (size      != A->size      ) throw CSR_Exception(__LINE__,__FILE__);
    if (size_tot  != A->size_tot  ) throw CSR_Exception(__LINE__,__FILE__);
    if (type      != A->type      ) throw CSR_Exception(__LINE__,__FILE__);
    if (first_row != A->first_row ) throw CSR_Exception(__LINE__,__FILE__);

    int i,e;
    double *linefull=new double[size_tot];

    for (i=0;i<size;i++) {
        init_variable(linefull,size_tot);
        for (e=A->edge_i[i]-A->findx;e<A->edge_i[i+1]-A->findx;e++) {
            linefull[A->index_j[e]-A->findx]=real(factor*A->nnz[e]);
        }
        for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            nnz[e]+=linefull[index_j[e]-findx];
        }
    }
    delete[] linefull;
}

/************************************************************************************************/

template <class T>
void TCSR<T>::Bcast(int rootrank,MPI_Comm my_comm)
{
    MPI_Bcast(&size,1,MPI_INT,rootrank,my_comm);
    MPI_Bcast(&size_tot,1,MPI_INT,rootrank,my_comm);
    MPI_Bcast(&type,1,MPI_INT,rootrank,my_comm);
    MPI_Bcast(&n_nonzeros,1,MPI_INT,rootrank,my_comm);
    MPI_Bcast(&findx,1,MPI_INT,rootrank,my_comm);
    MPI_Bcast(&first_row,1,MPI_INT,rootrank,my_comm);
    
    MPI_Bcast(nnz,n_nonzeros,give_me_MPI_datatype(nnz),rootrank,my_comm);
    MPI_Bcast(index_i,size,MPI_INT,rootrank,my_comm);
    MPI_Bcast(index_j,n_nonzeros,MPI_INT,rootrank,my_comm);
    MPI_Bcast(edge_i,size+1,MPI_INT,rootrank,my_comm);
    MPI_Bcast(diag_pos,size,MPI_INT,rootrank,my_comm);
}

/************************************************************************************************/

template <class T>
void TCSR<T>::scatter(TCSR<T>*matdist,int root,MPI_Comm my_comm)
{
    int iam,nprocs;
    MPI_Comm_size(my_comm,&nprocs);
    MPI_Comm_rank(my_comm,&iam);

    int *scounts=new int[nprocs];
    int sloc=matdist->n_nonzeros;
    MPI_Allgather(&sloc,1,MPI_INT,scounts,1,MPI_INT,my_comm);

    int *first_rows=new int[nprocs];
    int frloc=matdist->first_row;
    MPI_Allgather(&frloc,1,MPI_INT,first_rows,1,MPI_INT,my_comm);

    int rows_start;
    int rows_end;
    if (iam==root){
        rows_start=first_row;
        rows_end=first_row+size;
    }
    MPI_Bcast(&rows_start,1,MPI_INT,root,my_comm);
    MPI_Bcast(&rows_end,1,MPI_INT,root,my_comm);

    for (int IR=0;IR<nprocs;IR++){
        int row=first_rows[IR];
        if (row<rows_start || row>=rows_end){
            scounts[IR]=0;
        }
    }

    int *displc=new int[nprocs];
    displc[0]=0;
    for (int IP=1;IP<nprocs;IP++){
        displc[IP]=displc[IP-1]+scounts[IP-1];
    }

    MPI_Scatterv(nnz,scounts,displc,give_me_MPI_datatype(nnz),matdist->nnz,scounts[iam],give_me_MPI_datatype(nnz),root,my_comm);

    delete[] scounts;
    delete[] first_rows;
    delete[] displc;
}

/************************************************************************************************/

template <class T>
void TCSR<T>::reduce(int root,MPI_Comm my_comm)
{
    T *nnz_tmp=NULL;

    int iam,nprocs;
    MPI_Comm_size(my_comm,&nprocs);
    MPI_Comm_rank(my_comm,&iam);

    if (iam==root){
        nnz_tmp=new T[n_nonzeros];
    }

    if (n_nonzeros){
        MPI_Reduce(nnz,nnz_tmp,n_nonzeros,give_me_MPI_datatype(nnz),MPI_SUM,root,my_comm);
    }

    if (iam==root){
        c_tcopy(n_nonzeros,nnz_tmp,1,nnz,1);
        delete[] nnz_tmp;
    }
}

/************************************************************************************************/

template <class T>
void TCSR<T>::reducescatter(TCSR<T>*matdist,MPI_Comm my_comm)
{
    int iam,nprocs;
    MPI_Comm_size(my_comm,&nprocs);
    MPI_Comm_rank(my_comm,&iam);

    int *recvcounts=new int[nprocs];
    int recvloc=matdist->n_nonzeros;
    MPI_Allgather(&recvloc,1,MPI_INT,recvcounts,1,MPI_INT,my_comm);

    MPI_Reduce_scatter(nnz,matdist->nnz,recvcounts,give_me_MPI_datatype(nnz),MPI_SUM,my_comm);

    delete[] recvcounts;
}

/************************************************************************************************/

template <>
inline double TCSR<double>::psipsidagger(CPX* psi,int nkval)
{
    int i,e;
    double result=0;
    for(i=0;i<size;i++){
        for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
            result+=nnz[e]*real(c_zdotc(nkval,&psi[index_j[e]-findx],size_tot,&psi[first_row+i],size_tot));
        }
    }
    return result;
}

/************************************************************************************************/

template <>
inline void TCSR<double>::psipsidagger(CPX* psi,int nkval,CPX factor)
{
    int i,e;
    for(i=0;i<size;i++){
        for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
            nnz[e]+=real(factor*c_zdotc(nkval,&psi[index_j[e]-findx],size_tot,&psi[first_row+i],size_tot));
        }
    }
}

/************************************************************************************************/

template <>
inline void TCSR<double>::psipsidagger(double* psi,int nkval,double factor)
{
    int i,e;
    for(i=0;i<size;i++){
        for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
            nnz[e]+=factor*c_ddot(nkval,&psi[index_j[e]-findx],size_tot,&psi[first_row+i],size_tot);
        }
    }
}

/************************************************************************************************/

template <>
inline void TCSR<double>::psipsidagger_transpose(double* psi,int nkval,double factor)
{
    int i,e;
    for(i=0;i<size;i++){
        for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
            nnz[e]+=factor*c_ddot(nkval,&psi[(index_j[e]-findx)*size],1,&psi[i*size],1);
        }
    }
}

/************************************************************************************************/

template <>
inline void TCSR<double>::psipsidagger(CPX* psi,CPX* psil,CPX* psir,int nkval,int ndof,int bw,CPX factor)
{
    int i,e;
    for(i=0;i<size;i++) {
        for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            if ( (index_j[e]-findx)/ndof - (i+first_row)/ndof > bw ) {
                nnz[e]+=real(factor*c_zdotc(nkval,&psil[index_j[e]-findx],size_tot,&psi[first_row+i],size_tot));
            } else if ( (index_j[e]-findx)/ndof - (i+first_row)/ndof < -bw ) {
                nnz[e]+=real(factor*c_zdotc(nkval,&psir[index_j[e]-findx],size_tot,&psi[first_row+i],size_tot));
            } else {
                nnz[e]+=real(factor*c_zdotc(nkval,&psi[index_j[e]-findx],size_tot,&psi[first_row+i],size_tot));
            }
        }
    }
}

/************************************************************************************************/

template <>
inline void TCSR<double>::psipsidagger(CPX* psi,int nkval,CPX factor,int* atom_of_bf,double* result)
{
    int i,e;
    for(i=0;i<size;i++){
        for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
            result[atom_of_bf[first_row+i]]+=nnz[e]*real(factor*c_zdotc(nkval,&psi[index_j[e]-findx],size_tot,&psi[first_row+i],size_tot));
        }
    }
}

/************************************************************************************************/

template <>
inline void TCSR<double>::psipsidagger(CPX* psi,int nkval,CPX factor,int* atom_of_bf,double* result,double* Vatom,double energy)
{
    int i,e;
    for(i=0;i<size;i++){
        if (energy>Vatom[atom_of_bf[first_row+i]]) {
            for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
                result[atom_of_bf[first_row+i]]+=nnz[e]*real(factor*c_zdotc(nkval,&psi[index_j[e]-findx],size_tot,&psi[first_row+i],size_tot));
            }
        }
    }
}

/************************************************************************************************/

template <>
inline void TCSR<double>::psipsidagger(CPX* psi,CPX* psi2,int nkval,CPX factor,int* atom_of_bf,double* result)
{
    int i,e;
    for(i=0;i<size;i++){
        for(e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
            result[atom_of_bf[first_row+i]]+=nnz[e]*real(factor*c_zdotc(nkval,&psi2[index_j[e]-findx],size_tot,&psi[first_row+i],size_tot));
        }
    }
}

/************************************************************************************************/

template <>
inline void TCSR<double>::atomdensity(TCSR<CPX>* A,CPX factor,int* atom_of_bf,double* result)
{
    if (size      != A->size      ) throw CSR_Exception(__LINE__,__FILE__);
    if (size_tot  != A->size_tot  ) throw CSR_Exception(__LINE__,__FILE__);
    if (type      != A->type      ) throw CSR_Exception(__LINE__,__FILE__);
    if (first_row != A->first_row ) throw CSR_Exception(__LINE__,__FILE__);

    int i,e;
    double *linefull=new double[size_tot];

    for (i=0;i<size;i++) {
        init_variable(linefull,size_tot);
        for (e=A->edge_i[i]-A->findx;e<A->edge_i[i+1]-A->findx;e++) {
            linefull[A->index_j[e]-A->findx]=real(factor*A->nnz[e]);
        }
        for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            result[atom_of_bf[first_row+i]]+=nnz[e]*linefull[index_j[e]-findx];
        }
    }
    delete[] linefull;
}

/************************************************************************************************/

template <class T>
double TCSR<T>::collect_density(double* xyz_grid,double* xyz_atom,int* atom_of_bf)
{
    double basis_exp_C[]={3.35827350,0.76975780,0.20003810};
//    int fac[]={1, 1, 2, 6, 24, 120, 720, 5040};
    double basis_coeff_C[]={1.768061,6.480152,6.480152,6.480152,0.585701,1.027739,1.027739,1.027739,0.213179,0.190691,0.190691,0.190691};
    double br=18.8971616;

    double rho_grid=0.0;
    for (int i=first_row;i<first_row+size;i++) {
        double xi=br*(xyz_atom[3*atom_of_bf[i]+0]-xyz_grid[0]);
        double yi=br*(xyz_atom[3*atom_of_bf[i]+1]-xyz_grid[1]);
        double zi=br*(xyz_atom[3*atom_of_bf[i]+2]-xyz_grid[2]);
        int bf_wrt_at_i=i%12;
        double exp_i=basis_exp_C[bf_wrt_at_i/4];
        int l_of_bf_i=bf_wrt_at_i%4;
        int lxi=0;
        int lyi=0;
        int lzi=0;
        if (l_of_bf_i==1) lyi=1;
        if (l_of_bf_i==2) lzi=1;
        if (l_of_bf_i==3) lxi=1;
//        double norm_i=1.0/sqrt( (fac[2*lxi]*fac[2*lyi]*fac[2*lzi]*pow(M_PI,1.5)) / (pow(2,2*(lxi+lyi+lzi))*fac[lxi]*fac[lyi]*fac[lzi]*pow(exp_i,lxi+lyi+lzi+1.5)) );//DOESNT WORK
//        double norm_i=pow(2,lxi+lyi+lzi)*pow(2.0/M_PI,0.75)*pow(exp_i,(lxi+lyi+lzi+3.0/2.0)/2.0);//MAYBE WORKS
        double norm_i=basis_coeff_C[bf_wrt_at_i];
        double gto_i=norm_i*pow(xi,lxi)*pow(yi,lyi)*pow(zi,lzi)*exp(-exp_i*(xi*xi+yi*yi+zi*zi));
        if (gto_i>1.0E-3) for (int e=edge_i[i-first_row]-findx;e<edge_i[i-first_row+1]-findx;e++) {
            int j=index_j[e]-findx;

            int bf_wrt_at_j=j%12;
            double exp_j=basis_exp_C[bf_wrt_at_j/4];
            int l_of_bf_j=bf_wrt_at_j%4;
            int lxj=0;
            int lyj=0;
            int lzj=0;
            if (l_of_bf_j==1) lyj=1;
            if (l_of_bf_j==2) lzj=1;
            if (l_of_bf_j==3) lxj=1;
            double xj=br*(xyz_atom[3*atom_of_bf[j]+0]-xyz_grid[0]);
            double yj=br*(xyz_atom[3*atom_of_bf[j]+1]-xyz_grid[1]);
            double zj=br*(xyz_atom[3*atom_of_bf[j]+2]-xyz_grid[2]);
//            double norm_j=1.0/sqrt( (fac[2*lxj]*fac[2*lyj]*fac[2*lzj]*pow(M_PI,1.5)) / (pow(2,2*(lxj+lyj+lzj))*fac[lxj]*fac[lyj]*fac[lzj]*pow(exp_j,lxj+lyj+lzj+1.5)) );
//            double norm_j=pow(2,lxj+lyj+lzj)*pow(2.0/M_PI,0.75)*pow(exp_j,(lxj+lyj+lzj+3.0/2.0)/2.0);
            double norm_j=basis_coeff_C[bf_wrt_at_j];
            double gto_j=norm_j*pow(xj,lxj)*pow(yj,lyj)*pow(zj,lzj)*exp(-exp_j*(xj*xj+yj*yj+zj*zj));

//            if (lxi+lyi+lzi+lyj+lzj+lxj==0)
            rho_grid+=nnz[e]*gto_i*gto_j;
        }
    }
    return rho_grid;

}

/************************************************************************************************/

template <class T>
void TCSR<T>::moveawaypbc(T *left, T *right, int bw, int ndof)
{
    int i,e;
    int bwndof=bw*ndof;
    T zero=(T)0;
    init_variable(left,size_tot*size_tot);
    init_variable(right,size_tot*size_tot);

    for (i=0;i<bwndof;i++) {
        for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            if (index_j[e]-findx-i>=bwndof) {
                right[i+(index_j[e]-findx)*size_tot]=nnz[e];
                nnz[e]=zero;
            }
        }
    }

    for (i=size_tot-bwndof;i<size_tot;i++) {
        for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            if (index_j[e]-findx-i<=-bwndof) {
                left[i+(index_j[e]-findx)*size_tot]=nnz[e];
                nnz[e]=zero;
            }
        }
    }
}

/************************************************************************************************/

template <class T>
void TCSR<T>::settozeropbc(int bw, int ndof)
{
    int i,e;
    T zero=(T)0;

    for (i=0;i<max(0,min(size,bw*ndof-first_row));i++) {
        for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            if ( (index_j[e]-findx)/ndof - (i+first_row)/ndof > bw ) {
                nnz[e]=zero;
            }
        }
    }

    for (i=max(0,min(size,size_tot-bw*ndof-first_row));i<size;i++) {
        for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            if ( (index_j[e]-findx)/ndof - (i+first_row)/ndof < -bw ) {
                nnz[e]=zero;
            }
        }
    }
}

/************************************************************************************************/

template <class T>
void TCSR<T>::removepbc(int bw, int ndof)
{
    int i,e;

    for (i=0;i<max(0,min(size,bw*ndof-first_row));i++) {
        for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            if ( (index_j[e]-findx)/ndof - (i+first_row)/ndof > bw ) {
                index_j[e]=-1;
            }
        }
    }

    for (i=max(0,min(size,size_tot-bw*ndof-first_row));i<size;i++) {
        for (e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++) {
            if ( (index_j[e]-findx)/ndof - (i+first_row)/ndof < -bw ) {
                index_j[e]=-1;
            }
        }
    }

    n_nonzeros = 0;
    init_variable(index_i,size);

    for (int i=0;i<size;i++) {
        for(int e=edge_i[i]-findx;e<edge_i[i+1]-findx;e++){
            if(index_j[e]>=0){
                index_i[i]++;
                index_j[n_nonzeros]=index_j[e];
                nnz[n_nonzeros]=nnz[e];
                n_nonzeros++;
            }
        }
    }

    get_row_edge();
    get_diag_pos();

    if (n_nonzeros!=edge_i[size]-findx) throw CSR_Exception(__LINE__,__FILE__);
}

/************************************************************************************************/

template <class T>
void TCSR<T>::tr_full_to_sparse(T *B,int nrow,int ncol)
{
    int i,j;

    size       = nrow;
    n_nonzeros = 0;

    for(i=0;i<nrow;i++){

        index_i[i] = 0;

        for(j=0;j<ncol;j++){

            if(abs(B[i+j*nrow])>tollim){

                index_j[n_nonzeros] = j+findx;
                nnz[n_nonzeros]     = B[i*nrow+j];

                index_i[i]++;
                n_nonzeros++;
            }
        }
    }
    get_row_edge();
    get_diag_pos();
}

/************************************************************************************************/

template <class T>
void TCSR<T>::extract_lower_triangle(TCSR<T> *mat)
{
    findx=mat->findx;
    int i,e;
    n_nonzeros=0;
    for(i=0;i<size;i++){
        index_i[i]=0;
        for(e=mat->edge_i[i]-mat->findx;e<mat->edge_i[i+1]-mat->findx;e++){
            if(mat->index_j[e]-mat->findx<=i){
                index_i[i]++;
                index_j[n_nonzeros]=mat->index_j[e];
                nnz[n_nonzeros]=mat->nnz[e];
                n_nonzeros++;
            }
        }
    }
    get_row_edge();
    get_diag_pos();

    if (n_nonzeros!=edge_i[size]-findx) throw CSR_Exception(__LINE__,__FILE__);
}

/************************************************************************************************/

template <class T>
void TCSR<T>::extract_upper_triangle(TCSR<T> *mat)
{
    findx=mat->findx;
    int i,e;
    n_nonzeros=0;
    for(i=0;i<size;i++){
        index_i[i]=0;
        for(e=mat->edge_i[i]-mat->findx;e<mat->edge_i[i+1]-mat->findx;e++){
            if(mat->index_j[e]-mat->findx>=i){
                index_i[i]++;
                index_j[n_nonzeros]=mat->index_j[e];
                nnz[n_nonzeros]=mat->nnz[e];
                n_nonzeros++;
            }
        }
    }
    get_row_edge();
    get_diag_pos();

    if (n_nonzeros!=edge_i[size]-findx) throw CSR_Exception(__LINE__,__FILE__);
}

/************************************************************************************************/

template <class T>
void TCSR<T>::full_to_sparse(T *B,int nrow,int ncol,int ipos,int jpos)
{
    int i,j;

    n_nonzeros = 0;
    init_variable(index_i,size);

    int imax=min(ipos+nrow,size)-ipos;    
    for(i=0;i<imax;i++){

        for(j=0;j<ncol;j++){

	    if(abs(B[i+j*nrow])>tollim){

                index_j[n_nonzeros] = jpos+j+findx;
                nnz[n_nonzeros]     = B[i+j*nrow];

                index_i[i+ipos]++;
                n_nonzeros++;
            }
        }
    }
    get_row_edge();
    get_diag_pos();
}

/************************************************************************************************/

template <class T>
void TCSR<T>::sparse_transpose(TCSR<T> *mat)
{
    int i, j, e;
    int fin  = mat->findx;
    int fout = findx;

    n_nonzeros = mat->n_nonzeros;
    
    init_variable(index_i,size);
    for(i=0;i<size;i++){
        for(e=mat->edge_i[i]-fin;e<mat->edge_i[i+1]-fin;e++){
            index_i[mat->index_j[e]-fin]++;
        }
    }
    get_row_edge();

    int *ptr = new int[size];
    for(i=0;i<size;i++){
        ptr[i] = edge_i[i] - fout;
    }
    
    for(i=0;i<size;i++){
        for(e=mat->edge_i[i]-fin;e<mat->edge_i[i+1]-fin;e++){
            j = mat->index_j[e]-fin;
            index_j[ptr[j]] = i + fout;
            nnz[ptr[j]] = mat->nnz[e];
            ptr[j]++;
        }
    }

    delete[] ptr;

    get_diag_pos();
}

/************************************************************************************************/

template <class T>
void TCSR<T>::set_to_id()
{
    n_nonzeros=size;
    for (int i=0;i<size;i++) {
        nnz[i]=1.0;
        index_i[i]=1;
        index_j[i]=first_row+findx+i;
    }

    get_row_edge();
    get_diag_pos();
}

/*SAB*END****************************************************************************************/

template <class T>
void TCSR<T>::full_to_sparse(T *B,int nrow,int ncol)
{
    int i,j;

    size       = nrow;
    n_nonzeros = 0;
    
    for(i=0;i<nrow;i++){

        index_i[i] = 0;

        for(j=0;j<ncol;j++){

	    if(abs(B[i+j*nrow])>tollim){

                index_j[n_nonzeros] = j+findx;
                nnz[n_nonzeros]     = B[i+j*nrow];

                index_i[i]++;
                n_nonzeros++;
            }
        }
    }
    get_row_edge();
    get_diag_pos();
}

/************************************************************************************************/
/*
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
*/
/************************************************************************************************/
/*
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
*/
/************************************************************************************************/

template <>
inline void TCSR<double>::write_CSR_bin(const char* filename)
{
    int i,j,u=0;
    
    ofstream myfile;
    myfile.open (filename,ios::out | ios::binary);
    double dsize = (double) size;
    double dnnon = (double) n_nonzeros;
    double dfind = (double) findx;
    myfile.write((char*)&dsize,sizeof(double));
    myfile.write((char*)&dnnon,sizeof(double));
    myfile.write((char*)&dfind,sizeof(double));
    for(i=0;i<size;i++){                                              
        for(j=0;j<index_i[i];j++){                                    
            double ipos = (double) (i+findx);
            double jpos = (double) index_j[u];
            double rval = nnz[u];
            double ival = 0.0;
            myfile.write((char*)&ipos,sizeof(double));
            myfile.write((char*)&jpos,sizeof(double));
            myfile.write((char*)&rval,sizeof(double));
            myfile.write((char*)&ival,sizeof(double));
            u++;                                                      
        }                                                             
    }                                                                 
    myfile.close();                                                   

}

/************************************************************************************************/

template <>
inline void TCSR<CPX>::write_CSR_bin(const char* filename)
{
    int i,j,u=0;
    
    ofstream myfile;
    myfile.open (filename,ios::out | ios::binary);
    double dsize = (double) size;
    double dnnon = (double) n_nonzeros;
    double dfind = (double) findx;
    myfile.write((char*)&dsize,sizeof(double));
    myfile.write((char*)&dnnon,sizeof(double));
    myfile.write((char*)&dfind,sizeof(double));
    for(i=0;i<size;i++){                                              
        for(j=0;j<index_i[i];j++){                                    
            double ipos = (double) (i+findx);
            double jpos = (double) index_j[u];
            double rval = real(nnz[u]);
            double ival = imag(nnz[u]);
            myfile.write((char*)&ipos,sizeof(double));
            myfile.write((char*)&jpos,sizeof(double));
            myfile.write((char*)&rval,sizeof(double));
            myfile.write((char*)&ival,sizeof(double));
            u++;                                                      
        }                                                             
    }                                                                 
    myfile.close();                                                   

}

/************************************************************************************************/

template <>
inline void TCSR<double>::write_CSR(const char* filename)
{
    int i,j,u=0;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    myfile<<size<<"\n";
    myfile<<n_nonzeros<<"\n";
    myfile<<findx<<"\n";
    for(i=0;i<size;i++){                                              
        for(j=0;j<index_i[i];j++){                                    
            myfile<<i+findx<<" "<<index_j[u]<<" "<<nnz[u]<<"\n";  
            u++;                                                      
        }                                                             
    }                                                                 
    myfile.close();                                                   

}

/************************************************************************************************/

template <>
inline void TCSR<CPX>::write_CSR(const char* filename)
{
    
    int i,j,u=0;                                                        
                                                                      
    ofstream myfile;                                                  
    myfile.open (filename);                                           
    myfile.precision(8);                                              
    myfile<<size<<"\n";                                               
    myfile<<n_nonzeros<<"\n";                                         
    myfile<<findx<<"\n";                                              
    for(i=0;i<size;i++){                                              
        for(j=0;j<index_i[i];j++){                                    
	  myfile<<i+findx<<" "<<index_j[u]<<" "<<real(nnz[u])<<" "<<imag(nnz[u])<<"\n";      
            u++;                                                      
        }                                                             
    }                                                                 
    myfile.close();                                                   
    
}

/************************************************************************************************/

template <>
inline void TCSR<double>::write(const char* filename)
{
    int u=0,i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    for(i=0;i<size;i++){
        for(j=0;j<index_i[i];j++){
            myfile<<i+1<<" "<<index_j[u]+1-findx<<" "<<nnz[u]<<"\n";
            u++;
        }
    }
    myfile.close();                                                   
}

/************************************************************************************************/

template <>
inline void TCSR<CPX>::write(const char* filename)
{
    int u=0,i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    for(i=0;i<size;i++){
        for(j=0;j<index_i[i];j++){
            myfile<<i+1<<" "<<index_j[u]+1-findx<<" "<<real(nnz[u])<<" "\
                  <<imag(nnz[u])<<"\n";
            u++;
        }
    }
    myfile.close();
}

/************************************************************************************************/

template <class T>
void TCSR<T>::init_variable(int *var,int N)
{
    for(int i=0;i<N;i++){
        var[i] = 0;
    }
}

/************************************************************************************************/

template <class T>
void TCSR<T>::init_variable(double *var,int N)
{
    for(int i=0;i<N;i++){
        var[i] = 0.0;
    }
}

/************************************************************************************************/

template <class T>
void TCSR<T>::init_variable(CPX *var,int N)
{
    for(int i=0;i<N;i++){
        var[i] = CPX(0.0,0.0);
    }
}

/************************************************************************************************/

template <class T>
void TCSR<T>::add_to_corner(T *V,T sign,int *ind,int N)
{
    int i,j;
    
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            nnz[ind[i*N+j]] = nnz[ind[i*N+j]]+sign*V[i+j*N];
        }
    }
}

/************************************************************************************************/

template<class T>
void TCSR<T>::change_findx(int fortran_index)
{
    int fold = findx;
    int fnew = fortran_index;
    int i;
    
    for(i=0;i<=size;i++){edge_i[i]=edge_i[i]-fold+fnew;}
    for(i=0;i<n_nonzeros;i++){index_j[i]=index_j[i]-fold+fnew;}

    findx = fortran_index;
}

/************************************************************************************************/

template<class T>
void TCSR<T>::get_full_index_i(int *full_index_i)
{
    int i,j;

    for(i=0;i<size;i++){
        for(j=edge_i[i]-findx;j<edge_i[i+1]-findx;j++){
            full_index_i[j] = i+findx;
        }
    }
}

/************************************************************************************************/

template<class T>
int TCSR<T>::get_msize(int start,int stop,int tb,int* n_of_el)
{
  /*
    int IM,msize,factor;

    factor = 1;
    if(tb==20) factor = 2;

    msize = 0;
    for(IM=start;IM<=stop;IM++){
        msize = msize+factor*n_of_el[IM];
    }

    return msize;
  */
    return (tb/10)*(n_of_el[stop+1]-n_of_el[start]);
}
/************************************************************************************************/

template<class T>
void TCSR<T>::set_nnz(int n_nnz, T *nnz_arr)
{
    c_tcopy(n_nnz,nnz_arr,1,nnz,1);
}

template<class T>
void TCSR<T>::set_index_i(int N, int *index_i_arr)
{
    c_icopy(N,index_i_arr,1,index_i,1);
}

template<class T>
void TCSR<T>::set_index_j(int n_nnz, int *index_j_arr)
{
    c_icopy(n_nnz,index_j_arr,1,index_j,1);
}

template<class T>
void TCSR<T>::set_edge_i(int N, int *edge_i_arr)
{
    c_icopy(N+1,edge_i_arr,1,edge_i,1);
}

#endif
