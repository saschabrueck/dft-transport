#include "DiagScaLapack.H"
#include "p_eig.H"

int diagscalapack(TCSR<double> *Overlap,TCSR<double> *KohnSham,transport_parameters *parameters_transport)
{
    int mpi_size, mpi_rank;
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);

    int nocc=parameters_transport->n_occ;

    int icontxt = MPI_Comm_c2f(MPI_COMM_WORLD);
    int icontxt_csr = MPI_Comm_c2f(MPI_COMM_WORLD);
    int size_tot = KohnSham->size_tot;
    int size_csr_loc = KohnSham->size;
    int NB_csr = int(ceil(double(size_tot)/double(mpi_size)));
    char gridr[1] = {'R'};
    Cblacs_gridinit(&icontxt_csr,gridr,mpi_size,1);
    int descAcsr[9];
    int descloc[9];
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
    int nbl                  = (mb+nb)/2;//SOME ROUTINES REQUIRE MB==NB

    int rloc,cloc;
    rloc = max(1,c_numroc(size_tot,nbl,myrow,0,nprow));
    cloc = c_numroc(size_tot,nbl,mycol,0,npcol);
    c_descinit(descloc,size_tot,size_tot,nbl,nbl,0,0,icontxt,rloc,&iinfo);

    double *OVloc = new double[rloc*cloc];
    double *KSloc = new double[rloc*cloc];
    double* Acsr = new double[size_csr_loc*size_tot];
    Overlap->sparse_to_full(Acsr,size_csr_loc,size_tot);
    c_pdgemr2d(size_tot,size_tot,Acsr,1,1,descAcsr,OVloc,1,1,descloc,icontxt);
    KohnSham->sparse_to_full(Acsr,size_csr_loc,size_tot);
    c_pdgemr2d(size_tot,size_tot,Acsr,1,1,descAcsr,KSloc,1,1,descloc,icontxt);
    delete[] Acsr;

    double *eigval = new double[size_tot];
    double *Zloc  = new double[rloc*cloc];
    int nvec=size_tot;
    int nprocs=mpi_size;

    int fac_degen=max(0,int(floor(double(nvec)/sqrt(double(nprocs))))-1);
    double eps_eigval_degen=1.0E-6;
    double abstol = 2.0*c_pdlamch(icontxt,'S');
    double workytest[3];
    int iworkytest[3];
    int mout, nzout;
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

/*
    c_dscal(OverlapCollect->n_nonzeros,0.0,OverlapCollect->nnz,1);
if (!iam) sabtime=get_time(0.0);
//    if (!iam) OverlapCollect->psipsidagger(KSfull,nocc,1.0); THIS TAKES REALLY A LOT OF TIME
    if (!iam) full_transpose(nocc,nvec,KSfull,OVfull);
    if (!iam) OverlapCollect->psipsidagger_transpose(OVfull,nocc,1.0);
if (!iam) cout << "Time for Dens " << get_time(sabtime) << endl;
*/

    delete[] Zloc;

    if (!mpi_rank) for (int i_out=0;i_out<nocc;i_out++) cout << eigval[i_out] << "\n";

    delete[] eigval;

    return 0;
}

