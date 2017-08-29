/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __CSC
#define __CSC

#include "CSR.H"
#include "Types.H"

template <class T, class W>
class TCSC{
public:
    
    TCSC(int,int,int);
    TCSC(TCSR<T>*,int);
    ~TCSC();
    
    void write(const char*);
    void init_variable(T*,int);
    void get_column_edge();
    void vec_mat_mult(T*,T*,int);
    void vec_mat_mult_add(T*,T*,int);
    void convert_csr(TCSR<T>*);
    void convert_csr(TCSR<T>*,int,int);
    void imag_transpose(TCSR<double>*,CPX);
    void add_imag_transpose(TCSR<double>*,CPX);

    int size,n_nonzeros,findx;
    T *nnz;
    W *index_i,*index_j,*edge_j;
};


template <class T, class W>
TCSC<T,W>::TCSC(int N, int n_nnz, int fortran_index)
{

    nnz        = new T[n_nnz];
    index_i    = new W[n_nnz];
    index_j    = new W[N];
    edge_j     = new W[N+1];

    size       = N;
    n_nonzeros = n_nnz;
    findx      = fortran_index;
}

/************************************************************************************************/

template <class T,class W>
TCSC<T,W>::TCSC(TCSR<T> *CSRMat,int fortran_index)
{
    size       = CSRMat->size;
    n_nonzeros = CSRMat->n_nonzeros;
    findx      = fortran_index;
    
    nnz        = new T[n_nonzeros];
    index_i    = new W[n_nonzeros];
    index_j    = new W[size];
    edge_j     = new W[size+1];

    convert_csr(CSRMat);
}

/************************************************************************************************/

template <class T,class W>
void TCSC<T,W>::convert_csr(TCSR<T> *CSRMat)
{
    int i, j, val;
    W *rowptr;
    int col, rowptr_col;
    int fin  = CSRMat->findx;
    int fout = findx;
    
    rowptr   = new W[size];

    for(i=0;i<size;i++){rowptr[i] = 0;}

    for(i=0;i<size;i++){
        for(j=CSRMat->edge_i[i]-fin;j<CSRMat->edge_i[i+1]-fin;j++){
            rowptr[CSRMat->index_j[j]-fin]++;
        }
    }

    val = 0;
    for(i=0;i<size;i++){
        edge_j[i] = val + fout;
        val       = val + rowptr[i];
        rowptr[i] = edge_j[i] - fout;
    }
    edge_j[size] = val + fout;
    
    for(i=0;i<size;i++){
        for(j=CSRMat->edge_i[i]-fin;j<CSRMat->edge_i[i+1]-fin;j++){
            col = CSRMat->index_j[j]-fin;
            rowptr_col = rowptr[col];
            index_i[rowptr_col] = i + fout;
            nnz[rowptr_col] = CSRMat->nnz[j];
            rowptr[col]++;
        }
    }

    delete[] rowptr;

}
    
/************************************************************************************************/

template <class T,class W>
void TCSC<T,W>::convert_csr(TCSR<T> *CSRMat,int nrow,int ncol)
{
    int i, j, val;
    W *rowptr;
    int col, rowptr_col;
    int fin  = CSRMat->findx;
    int fout = findx;
    size     = ncol;
    
    rowptr   = new W[size];

    for(i=0;i<size;i++){rowptr[i] = 0;}

    for(i=0;i<nrow;i++){
        for(j=CSRMat->edge_i[i]-fin;j<CSRMat->edge_i[i+1]-fin;j++){
            rowptr[CSRMat->index_j[j]-fin]++;
        }
    }

    val = 0;
    for(i=0;i<size;i++){
        edge_j[i] = val + fout;
        val       = val + rowptr[i];
        rowptr[i] = edge_j[i] - fout;
    }
    edge_j[size] = val + fout;
    
    for(i=0;i<nrow;i++){
        for(j=CSRMat->edge_i[i]-fin;j<CSRMat->edge_i[i+1]-fin;j++){
            col = CSRMat->index_j[j]-fin;
            rowptr_col = rowptr[col];
            index_i[rowptr_col] = i + fout;
            nnz[rowptr_col] = CSRMat->nnz[j];
            rowptr[col]++;
        }
    }

    delete[] rowptr;

}
    
/************************************************************************************************/

template <class T,class W>
TCSC<T,W>::~TCSC()
{
    delete[] nnz;
    delete[] index_i;
    delete[] index_j;
    delete[] edge_j;

}

/************************************************************************************************/

template <class T,class W>
void TCSC<T,W>::init_variable(T *var,int N)
{
    for(int i=0;i<N;i++){
        var[i]=(T)0;
    }
}

/************************************************************************************************/

template <>
inline void TCSC<double,int>::write(const char* filename)
{
    int i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    for(j=0;j<size;j++){
        for(i=edge_j[j]-findx;i<edge_j[j+1]-findx;i++){
            myfile<<index_i[i]+1-findx<<" "<<j+1<<" "<<nnz[i]<<"\n";
        }
    }
    myfile.close();
    
}

/************************************************************************************************/

template <>
inline void TCSC<double,long int>::write(const char* filename)
{
    int i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    for(j=0;j<size;j++){
        for(i=edge_j[j]-findx;i<edge_j[j+1]-findx;i++){
            myfile<<index_i[i]+1-findx<<" "<<j+1<<" "<<nnz[i]<<"\n";
        }
    }
    myfile.close();
    
}

/************************************************************************************************/

template <>
inline void TCSC<CPX,int>::write(const char* filename)
{
    int i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    for(j=0;j<size;j++){
        for(i=edge_j[j]-findx;i<edge_j[j+1]-findx;i++){
            myfile<<index_i[i]+1-findx<<" "<<j+1<<" "<<real(nnz[i])<<" "\
                  <<imag(nnz[i])<<"\n";
        }
    }
    myfile.close();
    
}

/************************************************************************************************/

template <>
inline void TCSC<CPX,long int>::write(const char* filename)
{
    int i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    for(j=0;j<size;j++){
        for(i=edge_j[j]-findx;i<edge_j[j+1]-findx;i++){
            myfile<<index_i[i]+1-findx<<" "<<j+1<<" "<<real(nnz[i])<<" "\
                  <<imag(nnz[i])<<"\n";
        }
    }
    myfile.close();
    
}

/************************************************************************************************/

template <class T,class W>
void TCSC<T,W>::get_column_edge()
{
    int j;

    edge_j[0] = findx;

    for(j=0;j<size;j++){
        edge_j[j+1] = edge_j[j]+index_j[j];
    }
    
}

/************************************************************************************************/

template <>
inline void TCSC<CPX,int>::vec_mat_mult(CPX *in, CPX *out, int NR)
{
    int i,j;
    CPX *vec = new CPX[NR];

    
    for(j=0;j<size;j++){
        init_variable(vec,NR);
        for(i=edge_j[j]-findx;i<(edge_j[j+1]-findx);i++){
            c_zaxpy(NR,nnz[i],&in[(index_i[i]-findx)*NR],1,vec,1);
        }
        c_zcopy(NR,vec,1,&out[j*NR],1);
    }
    
    delete[] vec;

}

/************************************************************************************************/

template <>
inline void TCSC<double,int>::vec_mat_mult(double *in, double *out, int NR)
{
    int i,j;
    double *vec = new double[NR];
    
    for(j=0;j<size;j++){
        init_variable(vec,NR);
        for(i=edge_j[j]-findx;i<(edge_j[j+1]-findx);i++){
            c_daxpy(NR,nnz[i],&in[(index_i[i]-findx)*NR],1,vec,1);
        }
        c_dcopy(NR,vec,1,&out[j*NR],1);
    }
    
    delete[] vec;

}

/************************************************************************************************/

template <>
inline void TCSC<CPX,int>::vec_mat_mult_add(CPX *in,CPX *out,int NR)
{
    int i,j;
    CPX *vec = new CPX[NR];
    
    for(j=0;j<size;j++){
        init_variable(vec,NR);
        for(i=edge_j[j]-findx;i<(edge_j[j+1]-findx);i++){
            c_zaxpy(NR,nnz[i],&in[(index_i[i]-findx)*NR],1,vec,1);
        }
	c_zaxpy(NR,CPX(1.0,0.0),vec,1,&out[j*NR],1);
    }
    
    delete[] vec;

}

/************************************************************************************************/

template <>
inline void TCSC<double,int>::vec_mat_mult_add(double *in,double *out,int NR)
{
    int i,j;
    double *vec = new double[NR];
    
    for(j=0;j<size;j++){
        init_variable(vec,NR);
        for(i=edge_j[j]-findx;i<(edge_j[j+1]-findx);i++){
            c_daxpy(NR,nnz[i],&in[(index_i[i]-findx)*NR],1,vec,1);
        }
	c_daxpy(NR,1.0,vec,1,&out[j*NR],1);
    }
    
    delete[] vec;

}

/************************************************************************************************/

template <>
inline void TCSC<CPX,int>::imag_transpose(TCSR<double> *A,CPX factor)
{
    if (size != A->size) throw CSR_Exception(__LINE__,__FILE__);

    int i,e;
    double *linefull=new double[size];

    for (i=0;i<size;i++) {
        for (int ii=0;ii<size;ii++) linefull[ii]=0.0;
        for (e=edge_j[i]-findx;e<edge_j[i+1]-findx;e++) {
            linefull[index_i[e]-findx]=imag(factor*nnz[e]);
        }
        for (e=A->edge_i[i]-A->findx;e<A->edge_i[i+1]-A->findx;e++) {
            A->nnz[e]=linefull[A->index_j[e]-A->findx];
        }
    }
    delete[] linefull;
}

/************************************************************************************************/

template <>
inline void TCSC<CPX,int>::add_imag_transpose(TCSR<double> *A,CPX factor)
{
    if (size != A->size) throw CSR_Exception(__LINE__,__FILE__);

    int i,e;
    double *linefull=new double[size];

    for (i=0;i<size;i++) {
        for (int ii=0;ii<size;ii++) linefull[ii]=0.0;
        for (e=edge_j[i]-findx;e<edge_j[i+1]-findx;e++) {
            linefull[index_i[e]-findx]=imag(factor*nnz[e]);
        }
        for (e=A->edge_i[i]-A->findx;e<A->edge_i[i+1]-A->findx;e++) {
            A->nnz[e]+=linefull[A->index_j[e]-A->findx];
        }
    }
    delete[] linefull;
}

/************************************************************************************************/

#endif
