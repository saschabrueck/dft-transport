#ifndef __P_EIG
#define __P_EIG

#include "ScaLapack.H"

inline int p_eig(double *KSfull,double *OVfull,double *eigval,int nvec,MPI_Comm p_eig_comm)
{
    int iam, nprocs;
    MPI_Comm_size(p_eig_comm,&nprocs);
    MPI_Comm_rank(p_eig_comm,&iam);
    int icontxt=MPI_Comm_c2f(p_eig_comm);

    int nprow, npcol, myrow, mycol;
    nprow=int(sqrt(double(nprocs)));
    npcol=nprocs/nprow;
    while (npcol*nprow!=nprocs) {
        nprow++;
        npcol=nprocs/nprow;
    } // end while
    char gridr[1] = {'R'};
    Cblacs_gridinit(&icontxt,gridr,nprow,npcol);
    Cblacs_gridinfo(icontxt,&nprow,&npcol,&myrow,&mycol);
    int row_per_processor    = nvec/nprow;
    int block_per_rprocessor = int(ceil(double(row_per_processor)/64));
    int mb                   = row_per_processor/block_per_rprocessor;
    int col_per_processor    = nvec/npcol;
    int block_per_cprocessor = int(ceil(double(col_per_processor)/64));
    int nb                   = col_per_processor/block_per_cprocessor;
    int nbl                  = (mb+nb)/2;

    int rloc      = max(1,c_numroc(nvec,nbl,myrow,0,nprow));
    int cloc      = c_numroc(nvec,nbl,mycol,0,npcol);

    int iinfo;
    int descKSfull[9],descOVfull[9];
    c_descinit(descKSfull,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
    c_descinit(descOVfull,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
    int descKS[9],descOV[9],descZ[9];
    c_descinit(descKS,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);
    c_descinit(descOV,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);
    c_descinit(descZ,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);

    double *KSloc = new double[rloc*cloc];
    double *OVloc = new double[rloc*cloc];
    double *Zloc  = new double[rloc*cloc];

    c_pdgeadd('N',nvec,nvec,1.0,KSfull,1,1,descKSfull,0.0,KSloc,1,1,descKS);
    c_pdgeadd('N',nvec,nvec,1.0,OVfull,1,1,descOVfull,0.0,OVloc,1,1,descOV);

    double eps_eigval_degen=1.0E-6;
    double abstol = 2.0*c_pdlamch(icontxt,'S');
    double workytest[3];
    int iworkytest[3];
    int mout, nzout;
    int *ifail      = new int[nvec];
    int *iclu       = new int[2*nprocs];
    double *dgap    = new double[nprocs];
    if (nprocs<1) return (LOGCERR, EXIT_FAILURE);  
    c_pdsygvx(1,'V','A','U',nvec,KSloc,1,1,descKS,OVloc,1,1,descOV,\
              0.0,0.0,1,1,abstol,&mout,&nzout,eigval,eps_eigval_degen,Zloc,1,1,descZ,\
              workytest,-1,iworkytest,-1,ifail,iclu,dgap,&iinfo);
    int lworky    = int(workytest[0])+nvec;
    int liworky   = iworkytest[0];
    double *worky = new double[lworky];
    int *iworky   = new int[liworky];
    c_pdsygvx(1,'V','A','U',nvec,KSloc,1,1,descKS,OVloc,1,1,descOV,\
              0.0,0.0,1,1,abstol,&mout,&nzout,eigval,eps_eigval_degen,Zloc,1,1,descZ,\
              worky,lworky,iworky,liworky,ifail,iclu,dgap,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    delete[] worky;
    delete[] iworky;
    delete[] ifail;
    delete[] iclu;
    delete[] dgap;
    delete[] OVloc;
    delete[] KSloc;

    c_pdgeadd('N',nvec,nvec,1.0,Zloc,1,1,descZ,0.0,KSfull,1,1,descKSfull);

    delete[] Zloc;

    Cblacs_gridexit(icontxt);

    return 0;
}

inline int p_eig(CPX *KSfull,CPX *OVfull,double *eigval,int nvec,MPI_Comm p_eig_comm)
{
    int iam, nprocs;
    MPI_Comm_size(p_eig_comm,&nprocs);
    MPI_Comm_rank(p_eig_comm,&iam);
    int icontxt=MPI_Comm_c2f(p_eig_comm);

    int nprow, npcol, myrow, mycol;
    nprow=int(sqrt(double(nprocs)));
    npcol=nprocs/nprow;
    while (npcol*nprow!=nprocs) {
        nprow++;
        npcol=nprocs/nprow;
    } // end while
    char gridr[1] = {'R'};
    Cblacs_gridinit(&icontxt,gridr,nprow,npcol);
    Cblacs_gridinfo(icontxt,&nprow,&npcol,&myrow,&mycol);
    int row_per_processor    = nvec/nprow;
    int block_per_rprocessor = int(ceil(double(row_per_processor)/64));
    int mb                   = row_per_processor/block_per_rprocessor;
    int col_per_processor    = nvec/npcol;
    int block_per_cprocessor = int(ceil(double(col_per_processor)/64));
    int nb                   = col_per_processor/block_per_cprocessor;
    int nbl                  = (mb+nb)/2;

    int rloc      = max(1,c_numroc(nvec,nbl,myrow,0,nprow));
    int cloc      = c_numroc(nvec,nbl,mycol,0,npcol);

    int iinfo;
    int descKSfull[9],descOVfull[9];
    c_descinit(descKSfull,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
    c_descinit(descOVfull,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
    int descKS[9],descOV[9],descZ[9];
    c_descinit(descKS,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);
    c_descinit(descOV,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);
    c_descinit(descZ,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);

    CPX *KSloc = new CPX[rloc*cloc];
    CPX *OVloc = new CPX[rloc*cloc];
    CPX *Zloc  = new CPX[rloc*cloc];

    c_pzgeadd('N',nvec,nvec,CPX(1.0,0.0),KSfull,1,1,descKSfull,CPX(0.0,0.0),KSloc,1,1,descKS);
    c_pzgeadd('N',nvec,nvec,CPX(1.0,0.0),OVfull,1,1,descOVfull,CPX(0.0,0.0),OVloc,1,1,descOV);

    double eps_eigval_degen=1.0E-6;
    double abstol = 2.0*c_pdlamch(icontxt,'S');
    CPX workytest[3];
    double rworkytest[3];
    int iworkytest[3];
    int mout, nzout;
    int *ifail      = new int[nvec];
    int *iclu       = new int[2*nprocs];
    double *dgap    = new double[nprocs];
    if (nprocs<1) return (LOGCERR, EXIT_FAILURE);  
    c_pzhegvx(1,'V','A','U',nvec,KSloc,1,1,descKS,OVloc,1,1,descOV,\
              0.0,0.0,1,1,abstol,&mout,&nzout,eigval,eps_eigval_degen,Zloc,1,1,descZ,\
              workytest,-1,rworkytest,-1,iworkytest,-1,ifail,iclu,dgap,&iinfo);
    int lworky      = int(real(workytest[0]));
    int lrworky     = int(rworkytest[0])+3*nvec;
    int liworky     = iworkytest[0];
    CPX *worky      = new CPX[lworky];
    double *rworky  = new double[lrworky];
    int *iworky     = new int[liworky];
    c_pzhegvx(1,'V','A','U',nvec,KSloc,1,1,descKS,OVloc,1,1,descOV,\
              0.0,0.0,1,1,abstol,&mout,&nzout,eigval,eps_eigval_degen,Zloc,1,1,descZ,\
              worky,lworky,rworky,lrworky,iworky,liworky,ifail,iclu,dgap,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    delete[] worky;
    delete[] rworky;
    delete[] iworky;
    delete[] ifail;
    delete[] iclu;
    delete[] dgap;
    delete[] OVloc;
    delete[] KSloc;

    c_pzgeadd('N',nvec,nvec,1.0,Zloc,1,1,descZ,0.0,KSfull,1,1,descKSfull);

    delete[] Zloc;

    Cblacs_gridexit(icontxt);

    return 0;
}

#endif

