/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FULL
#define _FULL

#include "CSR.H"
#include "Types.H"
#include "LinearSolver.H"

/*! \class Full
 *  
 *  \brief Interface for dense distributed solver
 *  
 */
template <class T>
class Full : public LinearSolver<T>{
	
public:

    Full(TCSR<T>*,MPI_Comm);

    virtual ~Full();
    virtual void prepare();
    virtual void prepare(int*,int*,int,int,int*,int);
    virtual void prepare_corner(CPX*,CPX*,int*,int*,int*,int,int,int*,int);
    virtual void solve_equation(T* res, T* rhs, int no_rhs);    
				
private:

    int size_tot;
    int size_csr_loc;
    int icontxt;
    int icontxt_csr;
    int nbl;
    int NB_csr;
    int descAloc[9];
    int *ipiv;
    CPX *Aloc;
};

/************************************************************************************************/

template <>
Full<CPX>::Full(TCSR<CPX>* mat,MPI_Comm solver_comm)
{

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

    delete[] Acsr;

    ipiv = new int[rloc+descAloc[4]];

}

/************************************************************************************************/

template <>
Full<CPX>::~Full()
{
    Cblacs_gridexit(icontxt);
    Cblacs_gridexit(icontxt_csr);
    delete[] Aloc;
    delete[] ipiv;
}

/************************************************************************************************/

template <>
void Full<CPX>::prepare()
{
    int iinfo;
    c_pzgetrf(size_tot,size_tot,Aloc,1,1,descAloc,ipiv,&iinfo);
}

/************************************************************************************************/

template <class T>
void Full<T>::prepare(int *Bmin,int *Bmax,int NBlock,int Bsize,int *orb_per_at,int tb)
{   
    prepare();
}

/************************************************************************************************/

template <class T>
void Full<T>::prepare_corner(CPX *ML_ref,CPX *MR_ref,int *NBC,int *Bmin,int *Bmax,int NBlock,\
			      int Bsize,int *orb_per_at,int tb)
{
}

/************************************************************************************************/

template <>
void Full<CPX>::solve_equation(CPX *res, CPX *arg_rhs, int no_rhs)
{
    if (no_rhs) {
        int iinfo;
        int descRHScsr[9];
        c_descinit(descRHScsr,size_tot,no_rhs,NB_csr,no_rhs,0,0,icontxt_csr,size_csr_loc,&iinfo);
        int nprow, npcol, myrow, mycol;
        Cblacs_gridinfo(icontxt,&nprow,&npcol,&myrow,&mycol);
        int rloc,cloc;
        rloc = max(1,c_numroc(size_tot,nbl,myrow,0,nprow));
        cloc = c_numroc(no_rhs,nbl,mycol,0,npcol);
        int descRHSloc[9];
        c_descinit(descRHSloc,size_tot,no_rhs,nbl,nbl,0,0,icontxt,rloc,&iinfo);
        CPX *RHSloc = new CPX[rloc*cloc];

        c_pzgemr2d(size_tot,no_rhs,arg_rhs,1,1,descRHScsr,RHSloc,1,1,descRHSloc,icontxt);

        c_pzgetrs('N',size_tot,no_rhs,Aloc,1,1,descAloc,ipiv,RHSloc,1,1,descRHSloc,&iinfo);

        c_pzgemr2d(size_tot,no_rhs,RHSloc,1,1,descRHSloc,res,1,1,descRHScsr,icontxt);

        delete[] RHSloc;
    }
}

/************************************************************************************************/

#endif
