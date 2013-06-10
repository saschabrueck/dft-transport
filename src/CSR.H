#ifndef __CSR
#define __CSR

#include <fstream.h>
#include "Blas.H"

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
    TCSR(char*);
    TCSR(char*,int,int);
    ~TCSR();
    void change_contain(TCSR<T>*);
    void copy_index(TCSR<T>*);
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
    void trans_mat_vec_mult(T*,T*,int,int,int);
    void spec_mat_vec_mult(T*,T*,int,int*,int*,int,int);
    void spec_mat_vec_mult(T*,T*,int,int*,int*,int,int,int*);
    void sparse_to_full(T*,int,int);
    void sparse_to_cmp_full(CPX*,int,int);
    void full_to_sparse(T*,int,int);
    void cmp_full_to_sparse(CPX*,int,int,CPX);
    void write_CSR(const char*);
    void write(const char*);
    void init_variable(double*,int);
    void init_variable(CPX*,int);
    void add_to_corner(T*,T,int*,int);
    void change_findx(int);
    void get_full_index_i(int*);
    int get_msize(int,int,int,int*);

    int pos,size,size_tot,n_nonzeros,type,findx,first_row;
    T *nnz;
    int *index_i,*index_j,*edge_i,*diag_pos;
    
private:
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

/************************************************************************************************/

template <>
inline TCSR<CPX>::TCSR(char *filename)
{
}

/***********************************************************************************************/

template <>
inline TCSR<double>::TCSR(char *filename)
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
void TCSR<T>::copy_index(TCSR<T> *mat)
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

template <class T>
void TCSR<T>::cpu_distribute(int *Istart,int *Nlocal,int N,int msize,int mrank)
{
    int max_rank,NCPU;

    NCPU             = ceil((double)N/msize);
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

/************************************************************************************************/

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

#endif
