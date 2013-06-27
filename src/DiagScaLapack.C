#include "ScaLapack.H"
#include "array_tools.H"

int diagscalapack(TCSR<double> *Overlap,TCSR<double> *KohnSham,TCSR<double> *P_Matrix,int nocc)
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);

    double sabtime;
    double *OVfull, *KSfull;
    int iinfo;

    int nvec=KohnSham->size_tot;
    if (nvec!=Overlap->size_tot) return (cerr<<__LINE__<<endl, EXIT_FAILURE);

    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
    if (!iam) {
        KSfull = new double[nvec*nvec];
        KohnShamCollect->sparse_to_full(KSfull,nvec,nvec);
    } // end if
    delete KohnShamCollect;

    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
    if (!iam) {
        OVfull = new double[nvec*nvec];
        OverlapCollect ->sparse_to_full(OVfull,nvec,nvec);
    } // end if

    double *eigval = new double [nvec];
   
    int nprow, npcol, myrow, mycol;
    int rloc, cloc;
    double *KSloc, *OVloc, *Zloc;
    nprow=int(sqrt(double(nprocs)));
    npcol=nprocs/nprow;
    while (npcol*nprow!=nprocs) {
        nprow++;
        npcol=nprocs/nprow;
    } // end while
    if (!iam) cout <<"Procs "<<nprocs<<" Rows "<<nprow<<" Cols "<<npcol<<endl;
//    int icontxt=int(MPI_COMM_WORLD);
    int icontxt=MPI_Comm_c2f(MPI_COMM_WORLD);
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

    rloc     = max(1,c_numroc(nvec,nbl,myrow,0,nprow));
    cloc     = c_numroc(nvec,nbl,mycol,0,npcol);
    KSloc    = new double[rloc*cloc];
    OVloc    = new double[rloc*cloc];
    Zloc     = new double[rloc*cloc];

    int descKSfull[9],descOVfull[9];
    c_descinit(descKSfull,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
    c_descinit(descOVfull,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
    int descKS[9],descOV[9],descZ[9];
    c_descinit(descKS,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);
    c_descinit(descOV,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);
    c_descinit(descZ,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);

    c_pdgeadd('N',nvec,nvec,1.0,KSfull,1,1,descKSfull,0.0,KSloc,1,1,descKS);
    c_pdgeadd('N',nvec,nvec,1.0,OVfull,1,1,descOVfull,0.0,OVloc,1,1,descOV);

    if (!iam) {
        delete[] KSfull;
        delete[] OVfull;
    }

    int *iworkytest = new int[2];
    double *workytest= new double[2];
    int mout, nzout;
    int *ifail      = new int[nvec];
    int *iclu       = new int[2*nprocs];
    double *dgap    = new double[nprocs];
    if (nprocs<1) return (cerr<<__LINE__<<endl, EXIT_FAILURE);  
    c_pdsygvx(1,'V','A','U',nvec,KSloc,1,1,descKS,OVloc,1,1,descOV,\
              0.0,0.0,1,1,0.0,&mout,&nzout,eigval,0.0,Zloc,1,1,descZ,\
              workytest,-1,iworkytest,-1,ifail,iclu,dgap,&iinfo);
    int lworky=int(workytest[0]);
    int liworky=iworkytest[0];
    double *worky  = new double[lworky];
    int *iworky = new int[liworky];
    delete[] workytest;
    delete[] iworkytest;
    if (!iam) sabtime=get_time(0.0);
    c_pdsygvx(1,'V','A','U',nvec,KSloc,1,1,descKS,OVloc,1,1,descOV,\
              0.0,0.0,1,1,0.0,&mout,&nzout,eigval,0.0,Zloc,1,1,descZ,\
              worky,lworky,iworky,liworky,ifail,iclu,dgap,&iinfo);
    if (iinfo) return (cerr<<__LINE__<<endl, EXIT_FAILURE);
    if (!iam) cout << "Time after Diag " << get_time(sabtime) << endl;
    delete[] worky;
    delete[] iworky;
    delete[] KSloc;
    delete[] ifail;
    delete[] iclu;
    delete[] dgap;

    delete[] OVloc;

    double *Zfull;
    int descZfull[9];
    c_descinit(descZfull,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
    if (!iam) Zfull = new double[nvec*nvec];
    c_pdgeadd('N',nvec,nvec,1.0,Zloc,1,1,descZ,0.0,Zfull,1,1,descZfull);
    double *Zfulltranspose;
    if (!iam) Zfulltranspose = new double [nvec*nvec];
    if (!iam) full_transpose(nvec,nvec,Zfull,Zfulltranspose);

// THIS IS OF COURSE UNNECESSARILY COMPLICATED BUT IT IS TO TRY MPI ROUTINES
    if (!iam) sabtime=get_time(0.0);
    TCSR<double> *Density;
    if (!iam)
        Density = new TCSR<double>(OverlapCollect,Zfulltranspose,1.0/nprocs,nocc);
    else
        Density = new TCSR<double>(OverlapCollect->size,OverlapCollect->n_nonzeros,OverlapCollect->findx);
    delete OverlapCollect;
    if (!iam) cout << "Time after Dens " << get_time(sabtime) << endl;
    if (!iam) sabtime=get_time(0.0);
    MPI_Bcast(Density->nnz,Density->n_nonzeros,MPI_DOUBLE,0,MPI_COMM_WORLD);
    Density->reducescatter(P_Matrix,MPI_COMM_WORLD);

    delete[] eigval;

    Cblacs_gridexit(icontxt);

    return 0;
}

