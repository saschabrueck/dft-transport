/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FULLINVERT
#define _FULLINVERT

#include "CSR.H"

class FullInvert{
	
public:

    FullInvert(TCSR<CPX>*,TCSR<double>*,CPX,MPI_Comm);
    ~FullInvert();

};

/************************************************************************************************/

FullInvert::FullInvert(TCSR<CPX>* mat,TCSR<double>* P,CPX factor,MPI_Comm solver_comm)
{

    int size_tot;
    int size_csr_loc;
    int icontxt;
    int icontxt_csr;
    int nbl;
    int NB_csr;
    int descAloc[9];
    int *ipiv;
    CPX *Aloc;
    int mpi_size;
    MPI_Comm_size(solver_comm,&mpi_size);

    icontxt = MPI_Comm_c2f(solver_comm);
    icontxt_csr = MPI_Comm_c2f(solver_comm);
    size_tot = mat->size_tot;
    size_csr_loc = mat->size;
    NB_csr = int(ceil(double(size_tot)/double(mpi_size)));
    CPX* Acsr = new CPX[size_csr_loc*size_tot];
    mat->sparse_to_full(Acsr,size_csr_loc,size_tot);
    char gridr[1] = {'R'};
    Cblacs_gridinit(&icontxt_csr,gridr,mpi_size,1);
    int descAcsr[9];
    int iinfo;
    c_descinit(descAcsr,size_tot,size_tot,NB_csr,size_tot,0,0,icontxt_csr,size_csr_loc,&iinfo);
//    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    int nprowcol[2]={0,0};
    MPI_Dims_create(mpi_size,2,nprowcol);
    int nprow = nprowcol[0];
    int npcol = nprowcol[1];
    int myrow, mycol;
    Cblacs_gridinit(&icontxt,gridr,nprow,npcol);
    Cblacs_gridinfo(icontxt,&nprow,&npcol,&myrow,&mycol);
    int nbl_in               = 64;
    int row_per_processor    = size_tot/nprow;
    int block_per_rprocessor = int(ceil(double(row_per_processor)/nbl_in));
    int mb                   = row_per_processor/block_per_rprocessor;
    int col_per_processor    = size_tot/npcol;
    int block_per_cprocessor = int(ceil(double(col_per_processor)/nbl_in));
    int nb                   = col_per_processor/block_per_cprocessor;
    nbl                      = (mb+nb)/2;//SOME ROUTINES REQUIRE MB==NB

    int rloc,cloc;
    rloc = max(1,c_numroc(size_tot,nbl,myrow,0,nprow));
    cloc = c_numroc(size_tot,nbl,mycol,0,npcol);
    c_descinit(descAloc,size_tot,size_tot,nbl,nbl,0,0,icontxt,rloc,&iinfo);
    Aloc = new CPX[rloc*cloc];

    c_pzgemr2d(size_tot,size_tot,Acsr,1,1,descAcsr,Aloc,1,1,descAloc,icontxt);

    ipiv = new int[rloc+descAloc[4]];
    c_pzgetrf(size_tot,size_tot,Aloc,1,1,descAloc,ipiv,&iinfo);
//    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    int lwork  = -1;
    int liwork = -1;
    CPX workq;
    int iworkq;
    CPX *work  = &workq;
    int *iwork = &iworkq;

    c_pzgetri(size_tot,Aloc,1,1,descAloc,ipiv,work,lwork,iwork,liwork,&iinfo);
//    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    lwork  = int(real(workq));
    liwork = iworkq;
    work   = new CPX[lwork];
    iwork  = new int[liwork];

    c_pzgetri(size_tot,Aloc,1,1,descAloc,ipiv,work,lwork,iwork,liwork,&iinfo);
//    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    delete[] work;
    delete[] iwork;
    delete[] ipiv;

    c_pzgemr2d(size_tot,size_tot,Aloc,1,1,descAloc,Acsr,1,1,descAcsr,icontxt);
    delete[] Aloc;

    for (int i=0;i<P->size;i++) {
        for (int e=P->edge_i[i]-P->findx;e<P->edge_i[i+1]-P->findx;e++) {
            P->nnz[e]+=real(factor*Acsr[i+size_csr_loc*(P->index_j[e]-P->findx)]);
        }
    }
    delete[] Acsr;

    Cblacs_gridexit(icontxt);
    Cblacs_gridexit(icontxt_csr);
}

/************************************************************************************************/

FullInvert::~FullInvert()
{
}

/************************************************************************************************/

#endif
