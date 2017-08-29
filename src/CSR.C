/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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

void CSR::write(const char* filename)
{
    int u=0,i,j;
    
    ofstream myfile;
    myfile.open (filename);
    myfile.precision(8);
    for(i=0;i<size;i++){
        for(j=0;j<index_i[i];j++){
	    myfile<<i+1<<" "<<index_j[u]+1<<" "<<r_nnz[u]<<" "<<i_nnz[u]<<"\n";
            u++;
        }
    }
    myfile.close();                                                   
}

/************************************************************************************************/
