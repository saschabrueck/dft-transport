#include "ParallelEig.H"

int p_grid_desc_init(int &icontxt,int nprocs,int nvec,int &rloc,int &cloc,int *descfull,int *descloc)
{
    int nprowcol[2]={0,0};
    MPI_Dims_create(nprocs,2,nprowcol);
    int nprow = nprowcol[0];
    int npcol = nprowcol[1];
    int myrow, mycol;
    char gridr[1] = {'R'};
    Cblacs_gridinit(&icontxt,gridr,nprow,npcol);
    Cblacs_gridinfo(icontxt,&nprow,&npcol,&myrow,&mycol);
    int nbl_in               = 64;
    int row_per_processor    = nvec/nprow;
    int block_per_rprocessor = int(ceil(double(row_per_processor)/nbl_in));
    int mb                   = row_per_processor/block_per_rprocessor;
    int col_per_processor    = nvec/npcol;
    int block_per_cprocessor = int(ceil(double(col_per_processor)/nbl_in));
    int nb                   = col_per_processor/block_per_cprocessor;
    int nbl                  = (mb+nb)/2;//SOME ROUTINES REQUIRE MB==NB

    rloc = max(1,c_numroc(nvec,nbl,myrow,0,nprow));
    cloc = c_numroc(nvec,nbl,mycol,0,npcol);

    int iinfo;
    c_descinit(descfull,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    c_descinit(descloc,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    return 0;
}

int p_eig(double *KSfull,double *OVfull,double *eigval,int nvec,MPI_Comm p_eig_comm)
{
    int nprocs;
    MPI_Comm_size(p_eig_comm,&nprocs);
    int icontxt=MPI_Comm_c2f(p_eig_comm);
    int rloc,cloc;
    int descfull[9];
    int descloc[9];
    if (p_grid_desc_init(icontxt,nprocs,nvec,rloc,cloc,descfull,descloc)) return (LOGCERR, EXIT_FAILURE);

    double *KSloc = new double[rloc*cloc];
    double *OVloc = new double[rloc*cloc];
    double *Zloc  = new double[rloc*cloc];

    c_pdgeadd('N',nvec,nvec,1.0,KSfull,1,1,descfull,0.0,KSloc,1,1,descloc);
    c_pdgeadd('N',nvec,nvec,1.0,OVfull,1,1,descfull,0.0,OVloc,1,1,descloc);

    int fac_degen=max(0,int(floor(double(nvec)/sqrt(double(nprocs))))-1);
    double eps_eigval_degen=1.0E-6;
    double abstol = 2.0*c_pdlamch(icontxt,'S');
    double workytest[3];
    int iworkytest[3];
    int mout, nzout, iinfo;
    int *ifail      = new int[nvec];
    int *iclu       = new int[2*nprocs];
    double *dgap    = new double[nprocs];
    if (nprocs<1) return (LOGCERR, EXIT_FAILURE);  
    c_pdsygvx(1,'V','A','U',nvec,KSloc,1,1,descloc,OVloc,1,1,descloc,\
              0.0,0.0,1,1,abstol,&mout,&nzout,eigval,eps_eigval_degen,Zloc,1,1,descloc,\
              workytest,-1,iworkytest,-1,ifail,iclu,dgap,&iinfo);
    int lworky    = int(workytest[0])+fac_degen*nvec;
    int liworky   = iworkytest[0];
    double *worky = new double[lworky];
    int *iworky   = new int[liworky];
    c_pdsygvx(1,'V','A','U',nvec,KSloc,1,1,descloc,OVloc,1,1,descloc,\
              0.0,0.0,1,1,abstol,&mout,&nzout,eigval,eps_eigval_degen,Zloc,1,1,descloc,\
              worky,lworky,iworky,liworky,ifail,iclu,dgap,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    delete[] worky;
    delete[] iworky;
    delete[] ifail;
    delete[] iclu;
    delete[] dgap;
    delete[] OVloc;
    delete[] KSloc;

    c_pdgeadd('N',nvec,nvec,1.0,Zloc,1,1,descloc,0.0,KSfull,1,1,descfull);

    delete[] Zloc;

    Cblacs_gridexit(icontxt);

    return 0;
}

int p_eig(CPX *KSfull,CPX *OVfull,double *eigval,int nvec,MPI_Comm p_eig_comm)
{
    int iam, nprocs;
    MPI_Comm_size(p_eig_comm,&nprocs);
    MPI_Comm_rank(p_eig_comm,&iam);
    int icontxt=MPI_Comm_c2f(p_eig_comm);
    int rloc,cloc;
    int descfull[9];
    int descloc[9];
    if (p_grid_desc_init(icontxt,nprocs,nvec,rloc,cloc,descfull,descloc)) return (LOGCERR, EXIT_FAILURE);

    CPX *KSloc = new CPX[rloc*cloc];
    CPX *OVloc = new CPX[rloc*cloc];
    CPX *Zloc  = new CPX[rloc*cloc];

    c_pzgeadd('N',nvec,nvec,CPX(1.0,0.0),KSfull,1,1,descfull,CPX(0.0,0.0),KSloc,1,1,descloc);
    c_pzgeadd('N',nvec,nvec,CPX(1.0,0.0),OVfull,1,1,descfull,CPX(0.0,0.0),OVloc,1,1,descloc);

    int fac_degen=max(0,int(floor(double(nvec)/sqrt(double(nprocs))))-1);
    double eps_eigval_degen=1.0E-6;
    double abstol = 2.0*c_pdlamch(icontxt,'S');
    CPX workytest[3];
    double rworkytest[3];
    int iworkytest[3];
    int mout, nzout, iinfo;
    int *ifail      = new int[nvec];
    int *iclu       = new int[2*nprocs];
    double *dgap    = new double[nprocs];
    if (nprocs<1) return (LOGCERR, EXIT_FAILURE);  
    c_pzhegvx(1,'V','A','U',nvec,KSloc,1,1,descloc,OVloc,1,1,descloc,\
              0.0,0.0,1,1,abstol,&mout,&nzout,eigval,eps_eigval_degen,Zloc,1,1,descloc,\
              workytest,-1,rworkytest,-1,iworkytest,-1,ifail,iclu,dgap,&iinfo);
    int lworky      = int(real(workytest[0]));
    int lrworky     = int(rworkytest[0])+fac_degen*nvec;
    int liworky     = iworkytest[0];
    CPX *worky      = new CPX[lworky];
    double *rworky  = new double[lrworky];
    int *iworky     = new int[liworky];
    c_pzhegvx(1,'V','A','U',nvec,KSloc,1,1,descloc,OVloc,1,1,descloc,\
              0.0,0.0,1,1,abstol,&mout,&nzout,eigval,eps_eigval_degen,Zloc,1,1,descloc,\
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

    c_pzgeadd('N',nvec,nvec,1.0,Zloc,1,1,descloc,0.0,KSfull,1,1,descfull);

    delete[] Zloc;

    Cblacs_gridexit(icontxt);

    return 0;
}

int p_inv(CPX *Afull,int nvec,MPI_Comm p_eig_comm)
{
    int nprocs;
    MPI_Comm_size(p_eig_comm,&nprocs);
    int icontxt=MPI_Comm_c2f(p_eig_comm);
    int rloc,cloc;
    int descfull[9];
    int descloc[9];
    if (p_grid_desc_init(icontxt,nprocs,nvec,rloc,cloc,descfull,descloc)) return (LOGCERR, EXIT_FAILURE);

    CPX *Aloc = new CPX[rloc*cloc];

    c_pzgeadd('N',nvec,nvec,1.0,Afull,1,1,descfull,0.0,Aloc,1,1,descloc);

    int iinfo;
    int *ipiv = new int[rloc+descloc[4]];

    c_pzgetrf(nvec,nvec,Aloc,1,1,descloc,ipiv,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    int lwork  = -1;
    int liwork = -1;
    CPX workq;
    int iworkq;
    CPX *work  = &workq;
    int *iwork = &iworkq;

    c_pzgetri(nvec,Aloc,1,1,descloc,ipiv,work,lwork,iwork,liwork,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    lwork  = int(real(workq));
    liwork = iworkq;
    work   = new CPX[lwork];
    iwork  = new int[liwork];

    c_pzgetri(nvec,Aloc,1,1,descloc,ipiv,work,lwork,iwork,liwork,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    delete[] work;
    delete[] iwork;
    delete[] ipiv;

    c_pzgeadd('N',nvec,nvec,1.0,Aloc,1,1,descloc,0.0,Afull,1,1,descfull);

    delete[] Aloc;

    Cblacs_gridexit(icontxt);

    return 0;
}

int p_inv(CPX *Afull,CPX *Bfull,int nvec,CPX z,MPI_Comm p_eig_comm)
{
    int nprocs;
    MPI_Comm_size(p_eig_comm,&nprocs);
    int icontxt=MPI_Comm_c2f(p_eig_comm);
    int rloc,cloc;
    int descfull[9];
    int descloc[9];
    if (p_grid_desc_init(icontxt,nprocs,nvec,rloc,cloc,descfull,descloc)) return (LOGCERR, EXIT_FAILURE);

    CPX *Aloc = new CPX[rloc*cloc];

    c_pzgeadd('N',nvec,nvec,1.0,Afull,1,1,descfull,0.0,Aloc,1,1,descloc);

    int iinfo;
    int *ipiv = new int[rloc+descloc[4]];

    c_pzgetrf(nvec,nvec,Aloc,1,1,descloc,ipiv,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    int lwork  = -1;
    int liwork = -1;
    CPX workq;
    int iworkq;
    CPX *work  = &workq;
    int *iwork = &iworkq;

    c_pzgetri(nvec,Aloc,1,1,descloc,ipiv,work,lwork,iwork,liwork,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    lwork  = int(real(workq));
    liwork = iworkq;
    work   = new CPX[lwork];
    iwork  = new int[liwork];

    c_pzgetri(nvec,Aloc,1,1,descloc,ipiv,work,lwork,iwork,liwork,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    delete[] work;
    delete[] iwork;
    delete[] ipiv;

    c_pzgeadd('N',nvec,nvec,   1.0,Aloc,1,1,descloc,0.0,Afull,1,1,descfull);
    c_pzgeadd('T',nvec,nvec,  -1.0,Aloc,1,1,descloc,1.0,Afull,1,1,descfull);
    c_pzgeadd('N',nvec,nvec,     z,Aloc,1,1,descloc,0.0,Bfull,1,1,descfull);
    c_pzgeadd('T',nvec,nvec,-1.0/z,Aloc,1,1,descloc,1.0,Bfull,1,1,descfull);

    delete[] Aloc;

    Cblacs_gridexit(icontxt);

    return 0;
}

int p_lin(CPX *Afull,CPX *RHSfull,CPX *SOLfull,int nvec,int nrhs,MPI_Comm p_eig_comm)
{
    int nprocs;
    MPI_Comm_size(p_eig_comm,&nprocs);
    int icontxt=MPI_Comm_c2f(p_eig_comm);
    int rloc,cloc;
    int descfull[9];
    int descloc[9];
    if (p_grid_desc_init(icontxt,nprocs,nvec,rloc,cloc,descfull,descloc)) return (LOGCERR, EXIT_FAILURE);

    CPX *Aloc   = new CPX[rloc*cloc];
    CPX *RHSloc = new CPX[rloc*cloc];

    c_pzgeadd('N',nvec,nvec,1.0,Afull  ,1,1,descfull,0.0,Aloc  ,1,1,descloc);
    c_pzgeadd('N',nvec,nrhs,1.0,RHSfull,1,1,descfull,0.0,RHSloc,1,1,descloc);

    int iinfo;
    int *ipiv = new int[rloc+descloc[4]];

    c_pzgetrf(nvec,nvec,Aloc,1,1,descloc,ipiv,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    c_pzgetrs('N',nvec,nrhs,Aloc,1,1,descloc,ipiv,RHSloc,1,1,descloc,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    delete[] Aloc;
    delete[] ipiv;

    c_pzgeadd('N',nvec,nrhs,1.0,RHSloc,1,1,descloc,0.0,SOLfull,1,1,descfull);

    delete[] RHSloc;

    Cblacs_gridexit(icontxt);

    return 0;
}

int p_lin(CPX *Afull,CPX *RHSfull,int nvec,int nrhs,CPX z,MPI_Comm p_eig_comm)
{
    int nprocs;
    MPI_Comm_size(p_eig_comm,&nprocs);
    int icontxt=MPI_Comm_c2f(p_eig_comm);
    int rloc,cloc;
    int descfull[9];
    int descloc[9];
    if (p_grid_desc_init(icontxt,nprocs,nvec,rloc,cloc,descfull,descloc)) return (LOGCERR, EXIT_FAILURE);

    CPX *Aloc   = new CPX[rloc*cloc];
    CPX *RHSloc = new CPX[rloc*cloc];
    CPX *RSSloc = new CPX[rloc*cloc];

    c_pzgeadd('N',nvec,nvec,1.0,Afull  ,1,1,descfull,0.0,Aloc  ,1,1,descloc);
    c_pzgeadd('N',nvec,nrhs,1.0,RHSfull,1,1,descfull,0.0,RHSloc,1,1,descloc);
    c_pzgeadd('N',nvec,nrhs,1.0,RHSfull,1,1,descfull,0.0,RSSloc,1,1,descloc);

    int iinfo;
    int *ipiv = new int[rloc+descloc[4]];

    c_pzgetrf(nvec,nvec,Aloc,1,1,descloc,ipiv,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    c_pzgetrs('N',nvec,nrhs,Aloc,1,1,descloc,ipiv,RHSloc,1,1,descloc,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    c_pzgetrs('T',nvec,nrhs,Aloc,1,1,descloc,ipiv,RSSloc,1,1,descloc,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);

    delete[] Aloc;
    delete[] ipiv;

    c_pzgeadd('N',nvec,nrhs,   1.0,RHSloc,1,1,descloc,0.0,  Afull,1,1,descfull);
    c_pzgeadd('N',nvec,nrhs,     z,RHSloc,1,1,descloc,0.0,RHSfull,1,1,descfull);

    c_pzgeadd('N',nvec,nrhs,  -1.0,RSSloc,1,1,descloc,1.0,  Afull,1,1,descfull);
    c_pzgeadd('N',nvec,nrhs,-1.0/z,RSSloc,1,1,descloc,1.0,RHSfull,1,1,descfull);

    delete[] RHSloc;
    delete[] RSSloc;

    Cblacs_gridexit(icontxt);

    return 0;
}
