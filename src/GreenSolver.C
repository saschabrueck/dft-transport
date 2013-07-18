#include "GreenSolver.H"

// Globals
PARAM *parameter;
WireStructure *nanowire;
ENERGY *En;
VOLTAGE *voltage;

int greensolver(TCSR<double> *Overlap,TCSR<double> *KohnSham,TCSR<double> *P_Matrix)
{
//    MPI_Comm sample_comm;
//    int size,rank,IT,IC;
//    int cond1,cond2,cond3,cond4;
//    int sample_id;

//    MPI_Init(&argc,&argv);
//    MPI_Comm_size(MPI_COMM_WORLD,&size);
//    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

//    #ifdef MPI_TIMING
//    // This is for the use of mpi timing
//    mpi_omen_id = rank;
//    #endif // MPI_TIMING

//    if(!rank){
//      time_t rawtime;
//      struct tm *timeinfo;

//      time(&rawtime);
//      timeinfo = localtime(&rawtime);
//      printf ("The current date/time is: %s\n",asctime(timeinfo));
//    }

    ofstream myfile;
//    ifstream myinpfile;
    int iam, nprocs;
    // THIS ROUTINE INCLUDES MPI INIT SIZE AND RANK

    Cblacs_pinfo(&iam,&nprocs);
 //  MPI::Intercomm Comm;
 //  Comm = MPI::COMM_WORLD;
    //iam=Comm.Get_rank();
//nprocs= Comm.Get_size();

// SKIP DIAG CODE
//if (1) {


/*
   Vector1D<int> row_block_size, col_block_size, local_rows, row_dist, col_dist, nblkrows_local_all;
   int sendbuf;
   int num_nodes, rank;
   int* recvbuf;

   row_block_size.assign(S.row_blk_size, S.row_blk_size + S.nblkrows_total);
   col_block_size.assign(S.col_blk_size, S.col_blk_size + S.nblkcols_total);
   local_rows.assign(S.local_rows, S.local_rows + S.nblkrows_local);
   row_dist.assign(S.row_dist, S.row_dist + S.nblkrows_total);
   col_dist.assign(S.col_dist, S.col_dist + S.nblkcols_total);

   num_nodes = Comm.Get_size();
   rank = Comm.Get_rank();

   recvbuf = new int[num_nodes];
   std::fill_n(recvbuf, num_nodes, 0);
   sendbuf = S.nblkrows_local;
   Comm.Allgather(&sendbuf, 1, MPI::INT, recvbuf, 1, MPI::INT);
   nblkrows_local_all.assign(recvbuf, recvbuf+num_nodes);

   TCSR<double> *Overlap, *KohnSham;
   Overlap = new TCSR<double>(S.fullmatrix_nrows, S.n_nze, 0);
   KohnSham = new TCSR<double>(KS.fullmatrix_nrows, KS.n_nze, 0);
   cDBCSR_to_CSRmatrix(S, Overlap);
   cDBCSR_to_CSRmatrix(KS, KohnSham);
   Overlap->change_findx(1);
   KohnSham->change_findx(1);

cout<<"Rank "<<rank<<" nblkrows_local "<<S.nblkrows_local<<"   ";
for (int iii=0;iii<S.nblkrows_local;iii++) cout<<S.local_rows[iii]*S.row_blk_size[iii]<<"   ";
cout<<endl;
*/
//   Overlap->write_CSR("csr.overlap");
//   KohnSham->write_CSR("csr.kohnsham");

/*    char *sabfileks = new char[255];
    char *sabfileov = new char[255];
    double sabtime;
    TCSR<double> *KohnSham, *Overlap;
    double *OVfull, *KSfull;
    int iinfo,nvec,nocc,ii;

    sabfileks="matrices.KOHN-SHAM";
    sabfileov="matrices.OVERLAP";
    sabtime=get_time(0.0);

    if (!iam) {
        KohnSham = new TCSR<double>(sabfileks);//,size,rank);
        Overlap  = new TCSR<double>(sabfileov);//,size,rank);
        ifstream occfile("nocc");
        occfile >> nocc;
        occfile.close();
        cout << "Read Time " << get_time(sabtime) << endl;
        nvec=KohnSham->size_tot;
        if (nvec!=Overlap->size_tot) LOGCERR, exit(EXIT_FAILURE);
        KSfull = new double[nvec*nvec];
        OVfull = new double[nvec*nvec];
        KohnSham->sparse_to_full(KSfull,nvec,nvec);
        Overlap ->sparse_to_full(OVfull,nvec,nvec);
        delete KohnSham;
        delete Overlap;
    } // end if
*/
/*
    int iinfo,nvec,nocc=5,ii;
    double sabtime;
    double *OVfull, *KSfull;
        nvec=KohnSham->size_tot;
        if (nvec!=Overlap->size_tot) LOGCERR, exit(EXIT_FAILURE);
        KSfull = new double[nvec*nvec];
        OVfull = new double[nvec*nvec];
        KohnSham->sparse_to_full(KSfull,nvec,nvec);
        Overlap ->sparse_to_full(OVfull,nvec,nvec);

//    MPI_Bcast(&nvec,1,MPI_INT,0,MPI_COMM_WORLD);
//    MPI_Bcast(&nocc,1,MPI_INT,0,MPI_COMM_WORLD);

     
//    MPI_Status mpistatus;
//    if (!iam) {
//        cout << KSfull[0] << " sended" << endl;
//        for (int iproc=1;iproc<nprocs;iproc++) MPI_Send(KSfull,nvec*nvec,MPI_INT,iproc,11,MPI_COMM_WORLD);
//    } // end if
//    else {
//        MPI_Recv(KSfull,nvec*nvec,MPI_INT,0,11,MPI_COMM_WORLD,&mpistatus);
//        cout << KSfull[0] << " from " << iam << endl;
//    } // end if
////    MPI_Barrier(MPI_COMM_WORLD);

//  int nmat=nvec*(nvec+1)/2;
//  double *KStri = new double[nmat];
//  double *OVtri = new double[nmat];

//    int ik=0;
//    for(ii=0;ii<nvec;ii++){
//        for(ij=0;ij<=ii;ij++){
//            OVtri[ik]=OVfull[ii*nvec+ij];
//            KStri[ik]=KSfull[ii*nvec+ij];
//            ik++;
//        }//end for ij
//    }//end for ii

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
    int icontxt=int(MPI_COMM_WORLD);
    Cblacs_gridinit(&icontxt,"R",nprow,npcol);
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

    int *iworkytest = new int[2];
    double *workytest= new double[2];
    int mout, nzout;
    int *ifail      = new int[nvec];
    int *iclu       = new int[2*nprocs];
    double *dgap    = new double[nprocs];
    if (nprocs<1) LOGCERR, exit(EXIT_FAILURE);  
    if (!iam) cout<<"Start diagonalization"<<endl;
    c_pdsygvx(1,'V','A','U',nvec,KSloc,1,1,descKS,OVloc,1,1,descOV,\
              0.0,0.0,1,1,0.0,&mout,&nzout,eigval,0.0,Zloc,1,1,descZ,\
              workytest,-1,iworkytest,-1,ifail,iclu,dgap,&iinfo);
    int lworky=int(workytest[0]);
    int liworky=iworkytest[0];
    double *worky  = new double[lworky];
    int *iworky = new int[liworky];
    delete[] workytest;
    delete[] iworkytest;
    c_pdsygvx(1,'V','A','U',nvec,KSloc,1,1,descKS,OVloc,1,1,descOV,\
              0.0,0.0,1,1,0.0,&mout,&nzout,eigval,0.0,Zloc,1,1,descZ,\
              worky,lworky,iworky,liworky,ifail,iclu,dgap,&iinfo);
    if (iinfo) {
  LOGCERR;
cout<<"info "<<iinfo<<endl;  
 exit(EXIT_FAILURE);}
    if (!iam) cout << "Time after diag " << get_time(sabtime) << endl;
    delete[] worky;
    delete[] iworky;
    delete[] KSloc;
    delete[] ifail;
    delete[] iclu;
    delete[] dgap;
    double *OVdec;
    int descOVdec[9];
    c_descinit(descOVdec,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
    if (!iam) OVdec = new double[nvec*nvec];
    c_pdgeadd('N',nvec,nvec,1.0,Zloc,1,1,descOV,0.0,OVdec,1,1,descOVdec);
    delete[] OVloc;
    // Make use of the decomposed Overlap matrix later
//hey, here i save z which are the coefficients on ovdec!!!
    if (!iam) {
        myfile.open ("eigval");
        myfile.precision(12);
        for(ii=0;ii<nvec;ii++) myfile<<eigval[ii]<<endl;
        myfile.close();
    } // end if

//    double *Densloc = new double[rloc*cloc];
//    int descD[9];
//    c_descinit(descD,nvec,nvec,nbl,nbl,0,0,icontxt,rloc,&iinfo);
//    c_pdgemm('N','T',nvec,nvec,nocc,2.0,Zloc,1,1,descZ,Zloc,1,1,descZ,0.0,Densloc,1,1,descD);
//    delete[] Zloc;

//    double *Densfull;
//    int descDensfull[9];
//    c_descinit(descDensfull,nvec,nvec,nvec,nvec,0,0,icontxt,nvec,&iinfo);
//    if (!iam) Densfull = new double[nvec*nvec];
//    c_pdgeadd('N',nvec,nvec,1.0,Densloc,1,1,descD,0.0,Densfull,1,1,descDensfull);
//    if (!iam) cout << "Time after mult " << get_time(sabtime) << endl;
//    delete[] Densloc;

//    if (!iam) {
//        double trace = c_ddot(nvec*nvec,Densfull,1,OVfull,1);
//        cout <<"Trace "<< trace <<endl;

//        TCSR<double> *Density;
//        Density = new TCSR<double>(nvec,nvec*nvec,1);
//        Density->full_to_sparse(Densfull,nvec,nvec);
//        delete[] Densfull;
//        Density->write_CSR("matrices.DENSITY");
//        delete Density;
//    } // end if

    if (!iam) {
        delete[] KSfull;
        delete[] OVfull;
//        delete[] OVdec;
    }

    delete[] eigval;

    Cblacs_gridexit(icontxt);
        TCSR<double> *Ps = new TCSR<double>(Overlap,OVdec,1.0,nocc);
        Ps->change_findx(0);
   CSRmatrix_to_cDBCSR(Ps, *P, row_block_size, col_block_size, row_dist, col_dist, local_rows, nblkrows_local_all);
cout<<"whatever"<<endl;

    //Cblacs_exit(1); doest exist because mpi finalize is sufficient
    //AND I NEED TO FINALIZE BECAUSE ELSE ONE PROCESS TERMINATES THE
    //OTHERS WHEN FINISHED

} //END IF 0 WHICH MEANS THAT IT SHOULD SKIP THE CODE UNTIL NOW
*/

    double d_one=1.0;
    double d_zer=0.0;
    CPX z_one=CPX(d_one,d_zer);
    CPX z_zer=CPX(d_zer,d_zer);
    double sabtime;
    int iinfo=0;
// get parameters
    int bandwidth;
    int ndof,ndofsq,ndofsqbandwidth;
    double evoltfactor;
    double colzerothr;
    bandwidth = 2;
    ndof = 624;
    evoltfactor = 27.2107;
    colzerothr = 1.0e-12;
    double energy;
    double energystart;
//    double energyend;
    double eps_limit;
    double eps_decay;
    energystart = 1.2;
//    energyend = 1.2;
    eps_limit = 1.0e-6;
    eps_decay = 1.0e-6;
    int inj_sign=-1;

// SHOULB BE GONE
 //   string inpfilename;
  //  ifstream paramfile(argv[1]);
   // paramfile >> inpfilename >> bandwidth >> ndof >> evoltfactor >> colzerothr >> eps_limit >> eps_decay >> energystart >> energyend;
   // paramfile.close();


//    double energystep;
    if (bandwidth<1) LOGCERR, exit(EXIT_FAILURE);
//    if (nprocs==1)
        energy=energystart;
//    else {
//        energystep=(energyend-energystart)/(nprocs-1);
//        energy=energystart+iam*energystep;
//    }
    double energyimag=d_zer;
    int complexenergypoint=1;
    ndofsq=ndof*ndof;
    ndofsqbandwidth=ndofsq*(2*bandwidth+1);
    int ntriblock=bandwidth*ndof;
    int triblocksize=ntriblock*ntriblock;
/*
   Vector1D<int> row_block_size, col_block_size, local_rows, row_dist, col_dist, nblkrows_local_all;
   int sendbuf;
   MPI::Intercomm Comm;
   int num_nodes, rank;
   int* recvbuf;

   row_block_size.assign(S.row_blk_size, S.row_blk_size + S.nblkrows_total);
   col_block_size.assign(S.col_blk_size, S.col_blk_size + S.nblkcols_total);
   local_rows.assign(S.local_rows, S.local_rows + S.nblkrows_local);
   row_dist.assign(S.row_dist, S.row_dist + S.nblkrows_total);
   col_dist.assign(S.col_dist, S.col_dist + S.nblkcols_total);

   Comm = MPI::COMM_WORLD;
   num_nodes = Comm.Get_size();
   rank = Comm.Get_rank();

   recvbuf = new int[num_nodes];
   std::fill_n(recvbuf, num_nodes, 0);
   sendbuf = S.nblkrows_local;
   Comm.Allgather(&sendbuf, 1, MPI::INT, recvbuf, 1, MPI::INT);
   nblkrows_local_all.assign(recvbuf, recvbuf+num_nodes);

   TCSR<double> *Overlap, *KohnSham;
   Overlap = new TCSR<double>(S.fullmatrix_nrows, S.n_nze, 1);
   KohnSham = new TCSR<double>(KS.fullmatrix_nrows, KS.n_nze, 1);

   cDBCSR_to_CSRmatrix(S, Overlap);
   cDBCSR_to_CSRmatrix(KS, KohnSham);
*/
//   Overlap->write_CSR("/home/seyedb/interface/csr.overlap");
//   KohnSham->write_CSR("/home/seyedb/interface/csr.kohnsham");

//   CSRmatrix_to_cDBCSR(Overlap, *P, row_block_size, col_block_size, row_dist, col_dist, local_rows, nblkrows_local_all);

//    TCSR<double> *KohnSham, *Overlap;




    TCSR<CPX> *SumHamC;

//    if (KohnSham->size_tot!=Overlap->size_tot) LOGCERR, exit(EXIT_FAILURE);  delete this line
// WARNING BECAUSE OF PBC IN HAMILTONIAN IT CANNOT DETERMINE BANDWIDTH
//    bandwidth=max(KohnSham->get_bandwidth(ndof),Overlap->get_bandwidth(ndof));
        //TCSR<CPX> *SumHamC = new TCSR<CPX>(SumHam->size_tot,SumHam->n_nonzeros,SumHam->findx);
        //SumHamC->copy_contain(SumHam,d_one);
    SumHamC = new TCSR<CPX>(CPX(evoltfactor,d_zer),KohnSham,-CPX(energy,energyimag),Overlap);
    int ncells=SumHamC->size_tot/ndof;										//COPY TO THIS PLACE
    if (ndof*ncells!=SumHamC->size_tot) LOGCERR, exit(EXIT_FAILURE); 				//COPY TO THIS PLACE 
//    delete KohnSham; //delete Overlap later

    CPX* KScpx;
    CPX* KScollect;
    CPX* KScpx_loc;
    int nrowsofkscpx=max(0,min(SumHamC->size,ndof-SumHamC->first_row));
    int KScpx_loc_size=nrowsofkscpx*ndof*(bandwidth+1);
    if (!iam) KScollect=new CPX[ndofsqbandwidth];
    if (nrowsofkscpx>0) {
        KScpx_loc=new CPX[KScpx_loc_size];
    }
    SumHamC->contactunitcell(KScpx_loc,nrowsofkscpx,ndof,bandwidth,1);
    int* KScollectcount=new int[nprocs];
    MPI::COMM_WORLD.Allgather(&KScpx_loc_size,1,MPI::INT,KScollectcount,1,MPI::INT);
    int* displcKScoll = new int[nprocs]();
    for (int iii=1;iii<nprocs;iii++) 
        displcKScoll[iii]=displcKScoll[iii-1]+KScollectcount[iii-1];
    MPI_Gatherv(KScpx_loc,KScpx_loc_size,MPI_DOUBLE_COMPLEX,KScollect,KScollectcount,displcKScoll,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
    if (nrowsofkscpx>0) {
        delete[] KScpx_loc;
    }
    if (!iam) {
        KScpx=new CPX[ndofsqbandwidth];
        full_transpose(ndof,ndof*(bandwidth+1),KScollect,&KScpx[bandwidth*ndof*ndof]);
        delete[] KScollect;
        for (int ibw=1;ibw<=bandwidth;ibw++)
            full_conjugate_transpose(ndof,ndof,&KScpx[(bandwidth+ibw)*ndof*ndof],&KScpx[(bandwidth-ibw)*ndof*ndof]);
    }

    double* KSfull;
    CPX* H0cpx;
    CPX* H1cpx;
    CPX* mats;
    CPX* matb;
    if (!iam) {
        KSfull=new double[ndofsqbandwidth];
        c_dcopy(ndofsqbandwidth,(double*)KScpx,2,KSfull,1);
// assemble tridiagonalblocks
        H0cpx=new CPX[triblocksize];
        H1cpx=new CPX[triblocksize]();
        for (int ibandwidth=0;ibandwidth<bandwidth;ibandwidth++)
            for (int jbandwidth=0;jbandwidth<bandwidth;jbandwidth++)
                for (int jdof=0;jdof<ndof;jdof++)
                    c_zcopy(ndof,&KScpx[(bandwidth-ibandwidth+jbandwidth)*ndofsq+jdof*ndof],1,&H0cpx[(bandwidth*(jdof+jbandwidth*ndof)+ibandwidth)*ndof],1);
        for (int ibandwidth=0;ibandwidth<bandwidth;ibandwidth++)
            for (int jbandwidth=0;jbandwidth<=ibandwidth;jbandwidth++)
                for (int jdof=0;jdof<ndof;jdof++)
                    c_zcopy(ndof,&KScpx[(2*bandwidth-ibandwidth+jbandwidth)*ndofsq+jdof*ndof],1,&H1cpx[(bandwidth*(jdof+jbandwidth*ndof)+ibandwidth)*ndof],1);

/****** HERE STARTS THE NEW METHOD ******/

        mats=new CPX[ndofsq];
        c_zcopy(ndofsq,KScpx,1,mats,1);
        for (int ibandwidth=1;ibandwidth<2*bandwidth+1;ibandwidth++)
            c_zaxpy(ndofsq,d_one,&KScpx[ibandwidth*ndofsq],1,mats,1);//easier is possible because mats is symmetric
        matb=new CPX[2*bandwidth*ndofsq];//lowest block of right hand side matrix
// i think this is right for the general case
        c_zcopy(ndofsq,KScpx,1,matb,1);
        for (int iband=1;iband<=2*(bandwidth-1);iband++) {
            c_zcopy(ndofsq,&matb[(iband-1)*ndofsq],1,&matb[iband*ndofsq],1);
            c_zaxpy(ndofsq,d_one,&KScpx[iband*ndofsq],1,&matb[iband*ndofsq],1);
        }
        c_zcopy(ndofsq,&KScpx[2*bandwidth*ndofsq],1,&matb[(2*bandwidth-1)*ndofsq],1);
        c_zscal(ndofsq,-z_one,&matb[(2*bandwidth-1)*ndofsq],1);
    }
// now do the inversion
    int *pivarrayn;
    double workyytest;
    if (!iam) {
        pivarrayn=new int[ndof];
        if (complexenergypoint) {
// the matrix mats is not hermitian for complex energy point i think
            sabtime=get_time(d_zer);
            c_zgetrf(ndof,ndof,mats,ndof,pivarrayn,&iinfo);
            if (iinfo) LOGCERR, exit(EXIT_FAILURE); 
            c_zgetrs('N',ndof,ndof*2*bandwidth,mats,ndof,pivarrayn,matb,ndof,&iinfo);
            if (iinfo) LOGCERR, exit(EXIT_FAILURE); 
            cout << "TIME FOR CPX INVERSION " << get_time(sabtime) << endl;
        } else {
            sabtime=get_time(d_zer);
            double *matsreal=new double[ndofsq];
            double *matbreal=new double[2*bandwidth*ndofsq];
            c_dcopy(ndofsq,(double*)mats,2,matsreal,1);
            c_dcopy(2*bandwidth*ndofsq,(double*)matb,2,matbreal,1);
            c_dsysv('U',ndof,ndof*2*bandwidth,matsreal,ndof,pivarrayn,matbreal,ndof,&workyytest,-1,&iinfo);
            int lworkyn=int(workyytest);
            double *workyn=new double[lworkyn];
            c_dsysv('U',ndof,ndof*2*bandwidth,matsreal,ndof,pivarrayn,matbreal,ndof,workyn,lworkyn,&iinfo);
            if (iinfo) LOGCERR, exit(EXIT_FAILURE); 
            delete[] workyn;
            delete[] matsreal;
            for (int imat=0;imat<2*bandwidth*ndofsq;imat++) matb[imat]=CPX(matbreal[imat],d_zer);
            delete[] matbreal;
            cout << "TIME FOR REAL INVERSION " << get_time(sabtime) << endl;
        }
        delete[] pivarrayn;
        delete[] mats;
    }
    CPX* sigmal;
    CPX* sigmar;
    CPX* matctri;
    CPX* presigmal;
    CPX* presigmar;
    CPX* H1cpxhc;
    int nprotra;
    CPX* Vtra;
    CPX* Vref;
    CPX* lambdavec;
    int* protravec;
    int* prorefvec;
    double* veltra;
    double* velref;
    CPX* injl;
    CPX* injr;
    CPX* injl1;
    CPX* injr1;
    if (!iam) {
// replicate matrix and add identity
        CPX *matn=new CPX[4*bandwidth*bandwidth*ndofsq];
        for (int icolumns=0;icolumns<2*bandwidth*ndof;icolumns++)
            for (int ibandwidth=0;ibandwidth<2*bandwidth;ibandwidth++)
                c_zcopy(ndof,&matb[icolumns*ndof],1,&matn[(icolumns*2*bandwidth+ibandwidth)*ndof],1);
        for (int ibandwidth=0;ibandwidth<2*bandwidth-1;ibandwidth++)
            for (int idiag=0;idiag<(2*bandwidth-1-ibandwidth)*ndof;idiag++)
                matn[(idiag+ibandwidth*ndof)*2*bandwidth*ndof+idiag]-=z_one;
        delete[] matb;
// identify and remove zero columns
        sabtime=get_time(d_zer);
        int *indnzcolvecn=new int[2*ntriblock];
        int *indrzcolvecn=new int[2*ntriblock];
        int nindnzcoln=0;
        int nindzerocoln=0;
        int nindrzcoln=0;
        for (int itriblock=0;itriblock<2*ntriblock;itriblock++)
            if (abs(matn[(itriblock+1)*2*ntriblock-1])>colzerothr)
                indnzcolvecn[nindnzcoln++]=itriblock;
            else {
                double sumcoln=c_dzasum(2*ntriblock,&matn[itriblock*2*ntriblock],1);
                CPX matndiag=matn[itriblock*2*ntriblock+itriblock];
                int asumnotzero=(sumcoln>colzerothr*2*ntriblock);
                int asumnotone=(sumcoln>1.0+colzerothr*2*ntriblock || sumcoln<1.0-colzerothr*2*ntriblock);
                int diagnotmone=(real(matndiag)>-d_one+colzerothr*2*ntriblock || real(matndiag)<-d_one-colzerothr*2*ntriblock);
                if ( asumnotzero && (asumnotone || diagnotmone) )
                    indnzcolvecn[nindnzcoln++]=itriblock;
                else {
                    nindzerocoln++;
                    if (itriblock>=bandwidth*ndof && itriblock<(bandwidth+1)*ndof)
                        indrzcolvecn[nindrzcoln++]=itriblock;
                }
            }
        if (nindzerocoln+nindnzcoln!=2*ntriblock) LOGCERR, exit(EXIT_FAILURE); 
        CPX *matmn=new CPX[nindnzcoln*nindnzcoln];
        CPX *matmr=new CPX[nindnzcoln*nindrzcoln];
        for (int jindnzcol=0;jindnzcol<nindnzcoln;jindnzcol++)
            for (int iindnzcol=0;iindnzcol<nindnzcoln;iindnzcol++)
                matmn[jindnzcol*nindnzcoln+iindnzcol]=matn[indnzcolvecn[jindnzcol]*2*ntriblock+indnzcolvecn[iindnzcol]];
        for (int jindnzcol=0;jindnzcol<nindnzcoln;jindnzcol++)
            for (int iindnzcol=0;iindnzcol<nindrzcoln;iindnzcol++)
                matmr[jindnzcol*nindrzcoln+iindnzcol]=matn[indnzcolvecn[jindnzcol]*2*ntriblock+indrzcolvecn[iindnzcol]];
        delete[] matn;
        cout << "TIME FOR REMOVING ZERO COLUMNS " << get_time(sabtime) << endl;
        int blockbwpos=-1;
        int blockbwpoe=-1;
        for (int ipos=0;ipos<nindnzcoln;ipos++) {
            if (indnzcolvecn[ipos]>=bandwidth*ndof && blockbwpos==-1) blockbwpos=ipos;
            if (indnzcolvecn[ipos]<=((bandwidth+1)*ndof-1)) blockbwpoe=ipos;
        }
        if (blockbwpos<0 || blockbwpoe<0) LOGCERR, exit(EXIT_FAILURE); 
// get eigenvalues and eigenvectors
        double *eigvalreal=new double[nindnzcoln];
        double *eigvalimag=new double[nindnzcoln];
        CPX *eigveccfull=new CPX[nindnzcoln*nindnzcoln];
        if (complexenergypoint) {
            CPX cdummy;
            CPX workyctest;
            CPX *eigvalcpx=new CPX[nindnzcoln];
            double *workdouble=new double[2*nindnzcoln];
            c_zgeev('N','V',nindnzcoln,matmn,nindnzcoln,eigvalcpx,&cdummy,1,eigveccfull,nindnzcoln,&workyctest,-1,workdouble,&iinfo);
            int lworkyc=int(real(workyctest));
            CPX *workyc=new CPX[lworkyc];
            c_zgeev('N','V',nindnzcoln,matmn,nindnzcoln,eigvalcpx,&cdummy,1,eigveccfull,nindnzcoln,workyc,lworkyc,workdouble,&iinfo);
            if (iinfo) LOGCERR, exit(EXIT_FAILURE);
            cout << "TIME FOR CPX DIAGONALIZATION " << get_time(sabtime) << endl;
            delete[] workdouble;
            delete[] workyc;
            for (int ieigval=0;ieigval<nindnzcoln;ieigval++) {
                eigvalreal[ieigval]=real(eigvalcpx[ieigval]);
                eigvalimag[ieigval]=imag(eigvalcpx[ieigval]);
            }
            delete[] eigvalcpx; 
        } else {
            double *eigvec=new double[nindnzcoln*nindnzcoln];
            double dddummy;
            double *matmnreal=new double[nindnzcoln*nindnzcoln];
            sabtime=get_time(d_zer);
            c_dcopy(nindnzcoln*nindnzcoln,(double*)matmn,2,matmnreal,1);
            c_dgeev('N','V',nindnzcoln,matmnreal,nindnzcoln,eigvalreal,eigvalimag,&dddummy,1,eigvec,nindnzcoln,&workyytest,-1,&iinfo);
            int lworky=int(workyytest);
            double *worky=new double[lworky];
            c_dgeev('N','V',nindnzcoln,matmnreal,nindnzcoln,eigvalreal,eigvalimag,&dddummy,1,eigvec,nindnzcoln,worky,lworky,&iinfo);
            if (iinfo) LOGCERR, exit(EXIT_FAILURE); 
            cout << "TIME FOR REAL DIAGONALIZATION " << get_time(sabtime) << endl;
            delete[] worky;
            delete[] matmnreal;
            int iwhile=0;
            while (iwhile<nindnzcoln) {
                if (!eigvalimag[iwhile]) {
                    for (int ivecelements=0;ivecelements<nindnzcoln;ivecelements++) {
                        eigveccfull[iwhile*nindnzcoln+ivecelements]=CPX(eigvec[iwhile*nindnzcoln+ivecelements],d_zer);
                    }
                    iwhile++;
                }
                else if (eigvalreal[iwhile]==eigvalreal[iwhile+1]) {
                    for (int ivecelements=0;ivecelements<nindnzcoln;ivecelements++) {
                        eigveccfull[ iwhile   *nindnzcoln+ivecelements]=CPX(eigvec[iwhile*nindnzcoln+ivecelements], eigvec[(iwhile+1)*nindnzcoln+ivecelements]);
                        eigveccfull[(iwhile+1)*nindnzcoln+ivecelements]=CPX(eigvec[iwhile*nindnzcoln+ivecelements],-eigvec[(iwhile+1)*nindnzcoln+ivecelements]);
                     }
                    iwhile+=2;
                }
                else LOGCERR, exit(EXIT_FAILURE);
            } // END WHILE
            delete[] eigvec;
        }
        delete[] matmn;
        CPX *eigvecc=new CPX[ndof*nindnzcoln];
        if (blockbwpoe-blockbwpos==ndof-1) {
            c_zlacpy('A',ndof,nindnzcoln,&eigveccfull[blockbwpos],nindnzcoln,eigvecc,ndof);
        } else {
// RECONSTRUCT
//add this rightshift NO JUST PUT ALL OF IT IN A SPECIAL ROUTINE
        sabtime=get_time(d_zer);
        CPX *eigveccrec=new CPX[nindnzcoln*nindrzcoln];
        c_zgemm('N','N',nindrzcoln,nindnzcoln,nindnzcoln,z_one,matmr,nindrzcoln,eigveccfull,nindnzcoln,z_zer,eigveccrec,nindrzcoln);
        for (int iind=0;iind<nindnzcoln;iind++)
            c_zscal(nindrzcoln,z_one/CPX(eigvalreal[iind],eigvalimag[iind]),&eigveccrec[nindrzcoln*iind],1);
        for (int iind=blockbwpos;iind<=blockbwpoe;iind++)
            c_zcopy(nindnzcoln,&eigveccfull[iind],nindnzcoln,&eigvecc[indnzcolvecn[iind]-bandwidth*ndof],ndof);
        for (int iind=0;iind<nindrzcoln;iind++)
            c_zcopy(nindnzcoln,&eigveccrec[iind],nindrzcoln,&eigvecc[indrzcolvecn[iind]-bandwidth*ndof],ndof);
        delete[] eigveccrec;
        cout << "TIME FOR RECONSTRUCTION " << get_time(sabtime) << endl;
    }
    delete[] eigveccfull;
    delete[] matmr;
// DETERMINE TYPE OF EIGENVALUE/VECTOR
    lambdavec=new CPX[nindnzcoln];
    int *dectravec=new int[nindnzcoln];
    int *decrefvec=new int[nindnzcoln];
    protravec=new int[nindnzcoln];
    prorefvec=new int[nindnzcoln];
    veltra=new double[nindnzcoln];
    velref=new double[nindnzcoln];
    int ndectra=0;
    int ndecref=0;
    nprotra=0;
    int nproref=0;
    CPX *matcdof=new CPX[ndofsq];
    CPX *vecout=new CPX[ndof];
    for (int iindnzcoln=0;iindnzcoln<nindnzcoln;iindnzcoln++) {
        double eigr=eigvalreal[iindnzcoln];
        double eigi=eigvalimag[iindnzcoln];
        if ( (abs(eigr)>eps_limit || abs(eigi)>eps_limit) && (abs(eigr+1)>eps_limit || abs(eigi)>eps_limit) ) {
           CPX lambda=z_one/(z_one/CPX(eigr,eigi)+z_one);
           if ((abs(lambda)>d_one+eps_decay && inj_sign>0) || (abs(lambda)<d_one-eps_decay && inj_sign<0))
               dectravec[ndectra++]=iindnzcoln;
           else if ((abs(lambda)<d_one-eps_decay && inj_sign>0) || (abs(lambda)>d_one+eps_decay && inj_sign<0))
               decrefvec[ndecref++]=iindnzcoln;
           else {
               c_zcopy(ndofsq,&KScpx[(bandwidth+1)*ndofsq],1,matcdof,1);
               for (int ibandw=2;ibandw<=bandwidth;ibandw++)
                   c_zaxpy(ndofsq,CPX(ibandw,d_zer)*pow(lambda,-(ibandw-1)),&KScpx[(bandwidth+ibandw)*ndofsq],1,matcdof,1);
               c_zgemv('N',ndof,ndof,z_one,matcdof,ndof,&eigvecc[iindnzcoln*ndof],1,z_zer,vecout,1);
               double velnum=-2.0*imag(z_one/lambda*c_zdotc(ndof,&eigvecc[iindnzcoln*ndof],1,vecout,1));
// the velocity is this numerator divided by "bandwidth*C'*(s_0+sum_i=1^bandwidth(lambda^i s_i+lambda^-i s_-i))*C"
// but as i think this number is always positive i didnt implement it yet
               if (velnum*inj_sign>0) {
                   veltra[nprotra]=abs(velnum);
                   protravec[nprotra++]=iindnzcoln;
                } else if (velnum*inj_sign<0) {
                   velref[nproref]=abs(velnum);
                   prorefvec[nproref++]=iindnzcoln;
                } 
                else LOGCERR, exit(EXIT_FAILURE);
           } // END IF decaying or propagating
           lambdavec[iindnzcoln]=lambda;
        } // END IF k not infinite
    } // END FOR
    delete[] eigvalreal;
    delete[] eigvalimag;
    delete[] matcdof;
    delete[] vecout;
    if (nprotra!=nproref) LOGCERR, exit(EXIT_FAILURE); 
    if (ndectra!=ndecref) ndectra=min(ndectra,ndecref);//return (LOGCERR, EXIT_FAILURE);
// FOR ABOVE CASE I NEED TO IMPLEMENT SORTING OF LAMBDA? OR JUST THROW AWAY SOME VECS OF THE LIST SETTING NDEC TO MIN OF PRO AND REF?
// ONCE I REMOVE ONE OF THE DIRECTIONS I WILL ALSO REMOVE THE LINE ABOVE
    int neigbas=nprotra+ndectra;
    Vtra=new CPX[bandwidth*ndof*neigbas];
    Vref=new CPX[bandwidth*ndof*neigbas];
    for (int ipro=0;ipro<nprotra;ipro++) {
        c_zcopy(ndof,&eigvecc[ndof*protravec[ipro]],1,&Vtra[ndof*bandwidth*ipro+ndof*(bandwidth-1)],1);
        for (int ibandw=2;ibandw<=bandwidth;ibandw++) {
            c_zcopy(ndof,&Vtra[ndof*bandwidth*ipro+ndof*(bandwidth+1-ibandw)],1,&Vtra[ndof*bandwidth*ipro+ndof*(bandwidth-ibandw)],1);
            c_zscal(ndof,lambdavec[protravec[ipro]],&Vtra[ndof*bandwidth*ipro+ndof*(bandwidth-ibandw)],1);
        }
    }
    for (int idec=0;idec<ndectra;idec++) {
        c_zcopy(ndof,&eigvecc[ndof*dectravec[idec]],1,&Vtra[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth-1)],1);
        for (int ibandw=2;ibandw<=bandwidth;ibandw++) {
            c_zcopy(ndof,&Vtra[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth+1-ibandw)],1,&Vtra[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth-ibandw)],1);
            c_zscal(ndof,lambdavec[dectravec[idec]],&Vtra[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth-ibandw)],1);
        }
    }
    for (int ipro=0;ipro<nprotra;ipro++) {
        c_zcopy(ndof,&eigvecc[ndof*prorefvec[ipro]],1,&Vref[ndof*bandwidth*ipro+ndof*(bandwidth-1)],1);
        for (int ibandw=2;ibandw<=bandwidth;ibandw++) {
            c_zcopy(ndof,&Vref[ndof*bandwidth*ipro+ndof*(bandwidth+1-ibandw)],1,&Vref[ndof*bandwidth*ipro+ndof*(bandwidth-ibandw)],1);
            c_zscal(ndof,lambdavec[prorefvec[ipro]],&Vref[ndof*bandwidth*ipro+ndof*(bandwidth-ibandw)],1);
        }
    }
    for (int idec=0;idec<ndectra;idec++) {
        c_zcopy(ndof,&eigvecc[ndof*decrefvec[idec]],1,&Vref[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth-1)],1);
        for (int ibandw=2;ibandw<=bandwidth;ibandw++) {
            c_zcopy(ndof,&Vref[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth+1-ibandw)],1,&Vref[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth-ibandw)],1);
            c_zscal(ndof,lambdavec[decrefvec[idec]],&Vref[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth-ibandw)],1);
        }
    }
    CPX normy;
    for (int ivec=0;ivec<neigbas;ivec++) {
        normy=CPX(d_one/c_dznrm2(bandwidth*ndof,&Vtra[ivec*bandwidth*ndof],1),d_zer);
        c_zscal(bandwidth*ndof,normy,&Vtra[ivec*bandwidth*ndof],1);
        veltra[ivec]*=real(normy)*real(normy);
        velref[ivec]*=real(normy)*real(normy);
        normy=CPX(d_one/c_dznrm2(bandwidth*ndof,&Vref[ivec*bandwidth*ndof],1),d_zer);
        c_zscal(bandwidth*ndof,normy,&Vref[ivec*bandwidth*ndof],1);
    }
    delete[] eigvecc;
// inverse g snake WARNING THE pow OF lambda IS DONE HERE
    CPX *matcpx=new CPX[ntriblock*neigbas];
    CPX *invgls=new CPX[neigbas*neigbas];
    CPX *invgrs=new CPX[neigbas*neigbas];
    sabtime=get_time(d_zer);
    //left
    c_zgemm('C','N',ntriblock,neigbas,ntriblock,z_one,H1cpx,ntriblock,Vtra,ntriblock,z_zer,matcpx,ntriblock);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vtra,ntriblock,matcpx,ntriblock,z_zer,invgls,neigbas);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[protravec[ieigbas]],+bandwidth),&invgls[ieigbas*neigbas],1);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[dectravec[ieigbas-nprotra]],+bandwidth),&invgls[ieigbas*neigbas],1);
    c_zgemm('N','N',ntriblock,neigbas,ntriblock,z_one,H0cpx,ntriblock,Vtra,ntriblock,z_zer,matcpx,ntriblock);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vtra,ntriblock,matcpx,ntriblock,z_one,invgls,neigbas);
    //right
    c_zgemm('N','N',ntriblock,neigbas,ntriblock,z_one,H1cpx,ntriblock,Vref,ntriblock,z_zer,matcpx,ntriblock);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_zer,invgrs,neigbas);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[prorefvec[ieigbas]],-bandwidth),&invgrs[ieigbas*neigbas],1);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[decrefvec[ieigbas-nprotra]],-bandwidth),&invgrs[ieigbas*neigbas],1);
    c_zgemm('N','N',ntriblock,neigbas,ntriblock,z_one,H0cpx,ntriblock,Vref,ntriblock,z_zer,matcpx,ntriblock);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_one,invgrs,neigbas);
    cout << "TIME FOR SOME MATRIX MATRIX MULTIPLICATIONS " << get_time(sabtime) << endl;
    if (0) {
    CPX *KSeig=new CPX[neigbas*neigbas*(2*bandwidth+1)];
    CPX *invgtmp=new CPX[neigbas*neigbas];
    sabtime=get_time(d_zer);
// KSeig is KSfull in V-base
    for (int iband=bandwidth;iband<(2*bandwidth+1);iband++) {
        c_zgemm('N','N',ndof,neigbas,ndof,z_one,&KScpx[ndofsq*iband],ndof,Vtra,ntriblock,z_zer,matcpx,ntriblock);
        c_zgemm('C','N',neigbas,neigbas,ndof,z_one,Vtra,ntriblock,matcpx,ntriblock,z_zer,&KSeig[neigbas*neigbas*iband],neigbas);
    }
    for (int ibandwidth=1;ibandwidth<=bandwidth;ibandwidth++)
        for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
            for (int jeigbas=0;jeigbas<neigbas;jeigbas++)
                KSeig[(bandwidth-ibandwidth)*neigbas*neigbas+ieigbas*neigbas+jeigbas]=conj(KSeig[(bandwidth+ibandwidth)*neigbas*neigbas+jeigbas*neigbas+ieigbas]);
// invgls=h0
    c_zcopy(neigbas*neigbas,&KSeig[neigbas*neigbas*bandwidth],1,invgls,1);
// invgls+=h1*diag(lambda**-1)+c.c.
    c_zcopy(neigbas*neigbas,&KSeig[neigbas*neigbas*(bandwidth+1)],1,invgtmp,1);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[protravec[ieigbas]],-1),&invgtmp[ieigbas*neigbas],1);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[dectravec[ieigbas-nprotra]],-1),&invgtmp[ieigbas*neigbas],1);
    c_zaxpy(neigbas*neigbas,z_one,invgtmp,1,invgls,1);
    for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
        for (int jeigbas=0;jeigbas<neigbas;jeigbas++)
            invgls[jeigbas*neigbas+ieigbas]+=conj(invgtmp[ieigbas*neigbas+jeigbas]);
// invgls+=diag(lambda**-1)*h0*diag(lambda**-1)
    c_zcopy(neigbas*neigbas,&KSeig[neigbas*neigbas*bandwidth],1,invgtmp,1);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[protravec[ieigbas]],-1),&invgtmp[ieigbas*neigbas],1);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[dectravec[ieigbas-nprotra]],-1),&invgtmp[ieigbas*neigbas],1);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,conj(pow(lambdavec[protravec[ieigbas]],-1)),&invgtmp[ieigbas],neigbas);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,conj(pow(lambdavec[dectravec[ieigbas-nprotra]],-1)),&invgtmp[ieigbas],neigbas);
    c_zaxpy(neigbas*neigbas,z_one,invgtmp,1,invgls,1);
// invgls+=h2'*diag(lambda**2)
    c_zcopy(neigbas*neigbas,&KSeig[neigbas*neigbas*(bandwidth-2)],1,invgtmp,1);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[protravec[ieigbas]],+bandwidth),&invgtmp[ieigbas*neigbas],1);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[dectravec[ieigbas-nprotra]],+bandwidth),&invgtmp[ieigbas*neigbas],1);
    c_zaxpy(neigbas*neigbas,z_one,invgtmp,1,invgls,1);
// invgls+=diag(lambda**-1)*h2'*diag(lambda**2)*diag(lambda**-1)
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[protravec[ieigbas]],-1),&invgtmp[ieigbas*neigbas],1);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[dectravec[ieigbas-nprotra]],-1),&invgtmp[ieigbas*neigbas],1);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,conj(pow(lambdavec[protravec[ieigbas]],-1)),&invgtmp[ieigbas],neigbas);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,conj(pow(lambdavec[dectravec[ieigbas-nprotra]],-1)),&invgtmp[ieigbas],neigbas);
    c_zaxpy(neigbas*neigbas,z_one,invgtmp,1,invgls,1);
// invgls+=h1'*diag(lambda**-1)
    c_zcopy(neigbas*neigbas,&KSeig[neigbas*neigbas*(bandwidth-1)],1,invgtmp,1);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[protravec[ieigbas]],+bandwidth-1),&invgtmp[ieigbas*neigbas],1);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[dectravec[ieigbas-nprotra]],+bandwidth-1),&invgtmp[ieigbas*neigbas],1);
    c_zaxpy(neigbas*neigbas,z_one,invgtmp,1,invgls,1);
// now for grs
    if (complexenergypoint||1)
        for (int iband=bandwidth;iband<(2*bandwidth+1);iband++) {
            c_zgemm('N','N',ndof,neigbas,ndof,z_one,&KScpx[ndofsq*iband],ndof,Vref,ntriblock,z_zer,matcpx,ntriblock);
            c_zgemm('C','N',neigbas,neigbas,ndof,z_one,Vref,ntriblock,matcpx,ntriblock,z_zer,&KSeig[neigbas*neigbas*iband],neigbas);
        }
    else {
// this does not really save time, at least not for the CNT where I tested it
        double *Vreal=new double[bandwidth*ndof*neigbas];
        double *Vimag=new double[bandwidth*ndof*neigbas];
        double *matdb=new double[bandwidth*ndof*neigbas];
        double *KStmp=new double[neigbas*neigbas];
        for (int ii=0;ii<bandwidth*ndof*neigbas;ii++) {
            Vreal[ii]=real(Vref[ii]);
            Vimag[ii]=imag(Vref[ii]);
        }
        for (int iband=bandwidth;iband<(2*bandwidth+1);iband++) {
            c_dgemm('N','N',ndof,neigbas,ndof,d_one,&KSfull[ndofsq*iband],ndof,Vreal,ntriblock,d_zer,matdb,ntriblock);
            c_dgemm('T','N',neigbas,neigbas,ndof,d_one,Vreal,ntriblock,matdb,ntriblock,d_zer,KStmp,neigbas);
            for (int imat=0;imat<neigbas*neigbas;imat++)
                KSeig[neigbas*neigbas*iband+imat]=CPX(KStmp[imat],d_zer);
            c_dgemm('T','N',neigbas,neigbas,ndof,d_one,Vimag,ntriblock,matdb,ntriblock,d_zer,KStmp,neigbas);
            for (int imat=0;imat<neigbas*neigbas;imat++)
                KSeig[neigbas*neigbas*iband+imat]-=CPX(d_zer,KStmp[imat]);
            c_dgemm('N','N',ndof,neigbas,ndof,d_one,&KSfull[ndofsq*iband],ndof,Vimag,ntriblock,d_zer,matdb,ntriblock);
            c_dgemm('T','N',neigbas,neigbas,ndof,d_one,Vreal,ntriblock,matdb,ntriblock,d_zer,KStmp,neigbas);
            for (int imat=0;imat<neigbas*neigbas;imat++)
                KSeig[neigbas*neigbas*iband+imat]+=CPX(d_zer,KStmp[imat]);
            c_dgemm('T','N',neigbas,neigbas,ndof,d_one,Vimag,ntriblock,matdb,ntriblock,d_zer,KStmp,neigbas);
            for (int imat=0;imat<neigbas*neigbas;imat++)
                KSeig[neigbas*neigbas*iband+imat]+=CPX(KStmp[imat],d_zer);
        }
        delete[] Vreal;
        delete[] Vimag;
        delete[] matdb;
        delete[] KStmp;
    }
// scale h_i with diag(lambda**-i)
    for (int ibandwidth=1;ibandwidth<=bandwidth;ibandwidth++) {
        for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
            c_zscal(neigbas,pow(lambdavec[prorefvec[ieigbas]],-ibandwidth),&KSeig[(bandwidth+ibandwidth)*neigbas*neigbas+ieigbas*neigbas],1);
        for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
            c_zscal(neigbas,pow(lambdavec[decrefvec[ieigbas-nprotra]],-ibandwidth),&KSeig[(bandwidth+ibandwidth)*neigbas*neigbas+ieigbas*neigbas],1);
    }
// do the h.c. again CAN BE SIMPLIFIED WHEN DISCARDING gls 
// CAN EVEN BE LEFT OUT BECAUSE WE NEVER USE THE h.c. MATRICES BUT ALWAYS DO h.c. WHEN WE NEED IT
    for (int ibandwidth=1;ibandwidth<=bandwidth;ibandwidth++)
        for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
            for (int jeigbas=0;jeigbas<neigbas;jeigbas++)
                KSeig[(bandwidth-ibandwidth)*neigbas*neigbas+ieigbas*neigbas+jeigbas]=conj(KSeig[(bandwidth+ibandwidth)*neigbas*neigbas+jeigbas*neigbas+ieigbas]);
// first for h0
    c_zcopy(neigbas*neigbas,&KSeig[neigbas*neigbas*bandwidth],1,invgrs,1);
    if (bandwidth>1) {                                                    // ODER MACHE ICH HIER EINE AEUSSERE LOOP UEBER IBANDW?
        c_zcopy(neigbas*neigbas,&KSeig[neigbas*neigbas*bandwidth],1,invgtmp,1);
        for (int ibandwidth=2;ibandwidth<=bandwidth;ibandwidth++) {
            for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
                c_zscal(neigbas,pow(lambdavec[prorefvec[ieigbas]],-1),&invgtmp[ieigbas*neigbas],1);
            for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
                c_zscal(neigbas,pow(lambdavec[decrefvec[ieigbas-nprotra]],-1),&invgtmp[ieigbas*neigbas],1);
            for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
                c_zscal(neigbas,conj(pow(lambdavec[prorefvec[ieigbas]],-1)),&invgtmp[ieigbas],neigbas);
            for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
                c_zscal(neigbas,conj(pow(lambdavec[decrefvec[ieigbas-nprotra]],-1)),&invgtmp[ieigbas],neigbas);
            c_zaxpy(neigbas*neigbas,z_one,invgtmp,1,invgrs,1);
        }
    }
// then h1,2,...
    for (int jbandwidth=1;jbandwidth<=bandwidth;jbandwidth++) {
        c_zcopy(neigbas*neigbas,&KSeig[(bandwidth+jbandwidth)*neigbas*neigbas],1,invgtmp,1);
        c_zaxpy(neigbas*neigbas,z_one,invgtmp,1,invgrs,1);
// add h.c. 
        if (jbandwidth+1<=bandwidth) {
            for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
                for (int jeigbas=0;jeigbas<neigbas;jeigbas++)
                    invgrs[jeigbas*neigbas+ieigbas]+=conj(invgtmp[ieigbas*neigbas+jeigbas]);
        }
        for (int ibandwidth=2;ibandwidth<=bandwidth;ibandwidth++) {
            for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
                c_zscal(neigbas,pow(lambdavec[prorefvec[ieigbas]],-1),&invgtmp[ieigbas*neigbas],1);
            for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
                c_zscal(neigbas,pow(lambdavec[decrefvec[ieigbas-nprotra]],-1),&invgtmp[ieigbas*neigbas],1);
            for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
                c_zscal(neigbas,conj(pow(lambdavec[prorefvec[ieigbas]],-1)),&invgtmp[ieigbas],neigbas);
            for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
                c_zscal(neigbas,conj(pow(lambdavec[decrefvec[ieigbas-nprotra]],-1)),&invgtmp[ieigbas],neigbas);
            c_zaxpy(neigbas*neigbas,z_one,invgtmp,1,invgrs,1);
            if (jbandwidth+ibandwidth<=bandwidth) {
                for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
                    for (int jeigbas=0;jeigbas<neigbas;jeigbas++)
                        invgrs[jeigbas*neigbas+ieigbas]+=conj(invgtmp[ieigbas*neigbas+jeigbas]);
            } //end if
        } // end for ibandwidth
    } //end for jbandwidth
    cout << "TIME FOR NEW MATRIX MATRIX STUFF " << get_time(sabtime) << endl;
    delete[] KSeig;
    delete[] invgtmp;
    } // END IF 0 (SKIPPED)
// inversion
    int *pivarrayg=new int[neigbas];
    sigmal=new CPX[triblocksize];
    full_conjugate_transpose(neigbas,ntriblock,Vtra,matcpx);
    sabtime=get_time(d_zer);
    c_zgetrf(neigbas,neigbas,invgls,neigbas,pivarrayg,&iinfo);
    if (iinfo) LOGCERR, exit(EXIT_FAILURE);  
    c_zgetrs('N',neigbas,ntriblock,invgls,neigbas,pivarrayg,matcpx,neigbas,&iinfo);
    if (iinfo) LOGCERR, exit(EXIT_FAILURE); 
    cout << "TIME FOR gls INVERSION " << get_time(sabtime) << endl;
    c_zgemm('N','N',ntriblock,ntriblock,neigbas,z_one,Vtra,ntriblock,matcpx,neigbas,z_zer,sigmal,ntriblock);
    sigmar=new CPX[triblocksize];
    full_conjugate_transpose(neigbas,ntriblock,Vref,matcpx);
    sabtime=get_time(d_zer);
    c_zgetrf(neigbas,neigbas,invgrs,neigbas,pivarrayg,&iinfo);
    if (iinfo) LOGCERR, exit(EXIT_FAILURE); 
    c_zgetrs('N',neigbas,ntriblock,invgrs,neigbas,pivarrayg,matcpx,neigbas,&iinfo);
    if (iinfo) LOGCERR, exit(EXIT_FAILURE); 
    cout << "TIME FOR grs INVERSION " << get_time(sabtime) << endl;
    c_zgemm('N','N',ntriblock,ntriblock,neigbas,z_one,Vref,ntriblock,matcpx,neigbas,z_zer,sigmar,ntriblock);
    delete[] pivarrayg;
    delete[] invgls;
    delete[] invgrs;
    delete[] matcpx;
    presigmal= new CPX[triblocksize];
    presigmar= new CPX[triblocksize];
    matctri= new CPX[triblocksize];
    H1cpxhc = new CPX[triblocksize];
    full_conjugate_transpose(ntriblock,ntriblock,H1cpx,H1cpxhc);
    sabtime=get_time(d_zer);
    if (complexenergypoint||1) {
        TCSR<CPX> *tauhc = new TCSR<CPX>(ntriblock,triblocksize,SumHamC->findx);
        tauhc->full_to_sparse(H1cpxhc,ntriblock,ntriblock);
        full_transpose(ntriblock,ntriblock,sigmal,matctri);
        tauhc->trans_mat_vec_mult(matctri,presigmal,ntriblock,1);
        full_conjugate_transpose(ntriblock,ntriblock,presigmal,matctri);
        full_transpose(ntriblock,ntriblock,matctri,sigmal);
        // i think i could replace the two lines above with the two below, but need to test it with complex tau
        //c_zcopy(triblocksize,presigmal,1,sigmal,1);
        //c_dscal(triblocksize,-d_one,(double*)sigmal+1,2);
        tauhc->trans_mat_vec_mult(sigmal,matctri,ntriblock,1);
        full_conjugate_transpose(ntriblock,ntriblock,matctri,sigmal);
        TCSR<CPX> *tau = new TCSR<CPX>(ntriblock,triblocksize,SumHamC->findx);
        tau->full_to_sparse(H1cpx,ntriblock,ntriblock);
        full_transpose(ntriblock,ntriblock,sigmar,matctri);
        tau->trans_mat_vec_mult(matctri,presigmar,ntriblock,1);
        full_conjugate_transpose(ntriblock,ntriblock,presigmar,matctri);
        full_transpose(ntriblock,ntriblock,matctri,sigmar);
        tau->trans_mat_vec_mult(sigmar,matctri,ntriblock,1);
        full_conjugate_transpose(ntriblock,ntriblock,matctri,sigmar);
    } else {
        double *sigmad= new double[triblocksize];
        double *sigmai= new double[triblocksize];
        double *sigmat= new double[triblocksize];
        double *H1    = new double[triblocksize];
        c_dcopy(ndofsqbandwidth,(double*)H1cpx,2,H1,1);
//
        TCSR<double> *taureal = new TCSR<double>(ntriblock,triblocksize,SumHamC->findx);
        taureal->full_to_sparse(H1,ntriblock,ntriblock);
    sabtime=get_time(d_zer);
        taureal->mat_vec_mult(sigmad,sigmat,ntriblock);
    cout << "SPARSE MATRIX MATRIX MULTIPLICATION FOR REAL SIGMA " << get_time(sabtime) << endl;
    sabtime=get_time(d_zer);
        taureal->trans_mat_vec_mult(sigmad,sigmat,ntriblock,1);
    cout << "SPARSE MATRIX MATRIX MULTIPLICATION FOR REAL SIGMA TRANSP " << get_time(sabtime) << endl;
//
        for (int imat=0;imat<triblocksize;imat++) {
            sigmad[imat]=real(sigmal[imat]);
            sigmai[imat]=imag(sigmal[imat]);
        }
        c_dgemm('T','N',ntriblock,ntriblock,ntriblock,d_one,H1,ntriblock,sigmad,ntriblock,d_zer,sigmat,ntriblock);
        c_dgemm('N','N',ntriblock,ntriblock,ntriblock,d_one,sigmat,ntriblock,H1,ntriblock,d_zer,sigmad,ntriblock);
        for (int imat=0;imat<triblocksize;imat++) presigmal[imat]=CPX(sigmat[imat],d_zer);
        c_dgemm('T','N',ntriblock,ntriblock,ntriblock,d_one,H1,ntriblock,sigmai,ntriblock,d_zer,sigmat,ntriblock);
        for (int imat=0;imat<triblocksize;imat++) presigmal[imat]+=CPX(d_zer,sigmat[imat]);
        c_dgemm('N','N',ntriblock,ntriblock,ntriblock,d_one,sigmat,ntriblock,H1,ntriblock,d_zer,sigmai,ntriblock);
        for (int imat=0;imat<triblocksize;imat++) sigmal[imat]=CPX(sigmad[imat],sigmai[imat]);
        for (int imat=0;imat<triblocksize;imat++) {
            sigmad[imat]=real(sigmar[imat]);
            sigmai[imat]=imag(sigmar[imat]);
        }
        c_dgemm('N','N',ntriblock,ntriblock,ntriblock,d_one,H1,ntriblock,sigmad,ntriblock,d_zer,sigmat,ntriblock);
        c_dgemm('N','T',ntriblock,ntriblock,ntriblock,d_one,sigmat,ntriblock,H1,ntriblock,d_zer,sigmad,ntriblock);
        for (int imat=0;imat<triblocksize;imat++) presigmar[imat]=CPX(sigmat[imat],d_zer);
        c_dgemm('N','N',ntriblock,ntriblock,ntriblock,d_one,H1,ntriblock,sigmai,ntriblock,d_zer,sigmat,ntriblock);
        for (int imat=0;imat<triblocksize;imat++) presigmar[imat]+=CPX(d_zer,sigmat[imat]);
        c_dgemm('N','T',ntriblock,ntriblock,ntriblock,d_one,sigmat,ntriblock,H1,ntriblock,d_zer,sigmai,ntriblock);
        for (int imat=0;imat<triblocksize;imat++) sigmar[imat]=CPX(sigmad[imat],sigmai[imat]);
        delete[] sigmad;
        delete[] sigmai;
        delete[] sigmat;
        delete[] H1;
    }
    cout << "MATRIX MATRIX MULTIPLICATIONS FOR SIGMA " << get_time(sabtime) << endl;
    } // END IF MASTER
    int method=2;
    if (method==1) {
      if (!iam) {
        c_zaxpy(triblocksize,-z_one,sigmal,1,H0cpx,1);
        c_zaxpy(triblocksize,-z_one,sigmar,1,H0cpx,1);
        int *pivarrayd= new int[ntriblock];
        sabtime=get_time(d_zer);
        c_zgetrf(ntriblock,ntriblock,H0cpx,ntriblock,pivarrayd,&iinfo);
        if (iinfo) LOGCERR, exit(EXIT_FAILURE);
        CPX workytestc;
        c_zgetri(ntriblock,H0cpx,ntriblock,pivarrayd,&workytestc,-1,&iinfo);
        int lworkyd=int(real(workytestc));
        CPX *workyd=new CPX[lworkyd];
        c_zgetri(ntriblock,H0cpx,ntriblock,pivarrayd,workyd,lworkyd,&iinfo);
        if (iinfo) LOGCERR, exit(EXIT_FAILURE); 
        cout << "TIME FOR G INVERSION " << get_time(sabtime) << endl;
        delete[] pivarrayd;
        delete[] workyd;
//transm=-trace((sigmar-sigmar')*G*(sigmal-sigmal')*G'); i think its better not to use imag() because of symmetry error
        sabtime=get_time(d_zer);
        for (int jtri=0;jtri<ntriblock;jtri++)
            for (int itri=0;itri<ntriblock;itri++)
                matctri[itri*ntriblock+jtri]=sigmal[itri*ntriblock+jtri]-conj(sigmal[jtri*ntriblock+itri]);
        c_zgemm('N','C',ntriblock,ntriblock,ntriblock,z_one,matctri,ntriblock,H0cpx,ntriblock,z_zer,sigmal,ntriblock);
        for (int jtri=0;jtri<ntriblock;jtri++)
            for (int itri=0;itri<ntriblock;itri++)
                matctri[itri*ntriblock+jtri]=sigmar[itri*ntriblock+jtri]-conj(sigmar[jtri*ntriblock+itri]);
        c_zgemm('N','N',ntriblock,ntriblock,ntriblock,z_one,matctri,ntriblock,H0cpx,ntriblock,z_zer,sigmar,ntriblock);
        CPX trace=0;
        for (int jtri=0;jtri<ntriblock;jtri++)
            for (int itri=0;itri<ntriblock;itri++)
                trace+=-sigmar[itri*ntriblock+jtri]*sigmal[jtri*ntriblock+itri];
        cout << "MATRIX MATRIX MULTIPLICATIONS AND TRACE FOR TRANSMISSION " << get_time(sabtime) << endl;
        cout << "Energy " << energy << " Transmission " << real(trace) << endl;
     } //END IF MASTER
    } else if (method==2) {
        TCSR<CPX> *LeftCorner;
        TCSR<CPX> *RightCorner;
//        int nleftcorner;
//        int nrightcorner;
     if (!iam) {
        sabtime=get_time(d_zer);
        CPX *Vtracp=new CPX[ntriblock*nprotra];
        CPX *Vrefcp=new CPX[ntriblock*nprotra];
        for (int ivec=0;ivec<ntriblock*nprotra;ivec++) {
            Vtracp[ivec]=conj(Vtra[ivec]);
            Vrefcp[ivec]=conj(Vref[ivec]);
        }
        injl=new CPX[ntriblock*nprotra];
        injr=new CPX[ntriblock*nprotra];
        c_zgemm('C','N',ntriblock,nprotra,ntriblock,z_one,H1cpx,ntriblock,Vtracp,ntriblock,z_zer,injl,ntriblock);
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,z_one,H1cpx,ntriblock,Vrefcp,ntriblock,z_zer,injr,ntriblock);
        injl1=new CPX[ntriblock*nprotra];
        injr1=new CPX[ntriblock*nprotra];
        c_zcopy(ntriblock*nprotra,injl,1,injl1,1);
        c_zcopy(ntriblock*nprotra,injr,1,injr1,1);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,pow(lambdavec[protravec[ipro]],-bandwidth),&injl1[ipro*ntriblock],1);
            c_zscal(ntriblock,pow(lambdavec[prorefvec[ipro]],+bandwidth),&injr1[ipro*ntriblock],1);
        }
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,z_one,H0cpx,ntriblock,Vtracp,ntriblock,z_one,injl1,ntriblock);
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,z_one,H0cpx,ntriblock,Vrefcp,ntriblock,z_one,injr1,ntriblock);
// add I1 to I0 to form I
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,-z_one,presigmal,ntriblock,injl1,ntriblock,z_one,injl,ntriblock);
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,-z_one,presigmar,ntriblock,injr1,ntriblock,z_one,injr,ntriblock);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,CPX(d_one/sqrt(2*M_PI*veltra[ipro]),d_zer),&injl[ipro*ntriblock],1);
            c_zscal(ntriblock,CPX(d_one/sqrt(2*M_PI*velref[ipro]),d_zer),&injr[ipro*ntriblock],1);
        }
        cout << "MATRIX MATRIX MULTIPLICATIONS FOR INJECTION " << get_time(sabtime) << endl;
        full_transpose(nprotra,ntriblock,injl,injl1);
        full_transpose(nprotra,ntriblock,injr,injr1);
//        delete[] injl1;
//        delete[] injr1;
        delete[] presigmal;
        delete[] presigmar;
        } //end if master
// now add sigma and solve
// i think i did it this weird way because i thought that the full matrix needs to be square
// and as the zero elements will be deleted anyway, and the size only taken up to the edge_i
// of ntriblock and not two ntriblock, there is no problem, just a bit of extra memory
//BUT I SHOULD CHANGE THIS AND IGNORE IT ANYWAY BUT DISTRIBUTE THE MATRICES IN FULL FORMAT
//THEY ARE ALMOST FULL ANYWAY
        CPX *Hlcpx;
        CPX *Hrcpx;
        CPX *Hlcpxtr;
        CPX *Hrcpxtr;
        if (!iam) {
        Hlcpx = new CPX[4*triblocksize]();
        Hrcpx = new CPX[4*triblocksize]();
        c_zaxpy(triblocksize,-z_one,sigmal,1,H0cpx,1);
        c_zlacpy('A',ntriblock,ntriblock,H0cpx,ntriblock,Hlcpx,2*ntriblock);
        c_zlacpy('A',ntriblock,ntriblock,H1cpx,ntriblock,&Hlcpx[2*triblocksize],2*ntriblock);
//        LeftCorner = new TCSR<CPX>(2*ntriblock,2*triblocksize,SumHamC->findx);
//        LeftCorner->cmp_full_to_sparse(Hlcpx,2*ntriblock,2*ntriblock,z_one);
//        nleftcorner=LeftCorner->edge_i[ntriblock]-LeftCorner->findx;
        c_zaxpy(triblocksize,+z_one,sigmal,1,H0cpx,1);
        c_zaxpy(triblocksize,-z_one,sigmar,1,H0cpx,1);
        c_zlacpy('A',ntriblock,ntriblock,H0cpx,ntriblock,&Hrcpx[2*triblocksize],2*ntriblock);
        c_zlacpy('A',ntriblock,ntriblock,H1cpxhc,ntriblock,Hrcpx,2*ntriblock);
//        RightCorner = new TCSR<CPX>(2*ntriblock,2*triblocksize,SumHamC->findx);
//        RightCorner->cmp_full_to_sparse(Hrcpx,2*ntriblock,2*ntriblock,z_one);
//        nrightcorner=RightCorner->edge_i[ntriblock]-RightCorner->findx;
        delete[] H0cpx;
        delete[] H1cpx;
        delete[] H1cpxhc;
            delete[] sigmal;
            delete[] sigmar;
         Hlcpxtr = new CPX[4*triblocksize]();
         Hrcpxtr = new CPX[4*triblocksize]();
         full_transpose(2*ntriblock,2*ntriblock,Hlcpx,Hlcpxtr);
         full_transpose(2*ntriblock,2*ntriblock,Hrcpx,Hrcpxtr);
        delete[] Hlcpx;
        delete[] Hrcpx;
        }// end if master
// copy corners to complete sparse matrix
    int nrowsofleftcorner=max(0,min(SumHamC->size,ntriblock-SumHamC->first_row));
    int nrowsofrightcorner=max(0,min(SumHamC->size,ntriblock-SumHamC->size_tot+SumHamC->first_row+SumHamC->size));
    int nrowsofleftcornerelem=nrowsofleftcorner*2*ntriblock;
    int nrowsofrightcornerelem=nrowsofrightcorner*2*ntriblock;
    int* nrowsleftcornerarray=new int[nprocs];
    int* nrowsrightcornerarray=new int[nprocs];
    MPI_Allgather(&nrowsofleftcornerelem,1,MPI_INT,nrowsleftcornerarray,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&nrowsofrightcornerelem,1,MPI_INT,nrowsrightcornerarray,1,MPI_INT,MPI_COMM_WORLD);
    CPX* Hlcpx_local;
    CPX* Hrcpx_local;
    if (nrowsofleftcorner>0) {
        Hlcpx_local = new CPX[nrowsofleftcorner*2*ntriblock];
    }
    if (nrowsofrightcorner>0) {
        Hrcpx_local = new CPX[nrowsofrightcorner*2*ntriblock];
    }
    int* displcvecl = new int[nprocs]();
    int* displcvecr = new int[nprocs]();
    for (int iii=1;iii<nprocs;iii++) {
        displcvecl[iii]=displcvecl[iii-1]+nrowsleftcornerarray[iii-1];
        displcvecr[iii]=displcvecr[iii-1]+nrowsrightcornerarray[iii-1];
    }
    MPI_Scatterv(Hlcpxtr,nrowsleftcornerarray,displcvecl,MPI_DOUBLE_COMPLEX,Hlcpx_local,nrowsofleftcorner*2*ntriblock,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
    MPI_Scatterv(Hrcpxtr,nrowsrightcornerarray,displcvecr,MPI_DOUBLE_COMPLEX,Hrcpx_local,nrowsofrightcorner*2*ntriblock,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
    if (!iam) {
        delete[] Hlcpxtr;
        delete[] Hrcpxtr;
    }
        LeftCorner = new TCSR<CPX>(nrowsofleftcorner,nrowsofleftcorner*2*ntriblock,SumHamC->findx);
        LeftCorner->tr_full_to_sparse(Hlcpx_local,nrowsofleftcorner,2*ntriblock);
//        LeftCorner->first_row=displcvecl[iam];
        LeftCorner->first_row=SumHamC->first_row;
        LeftCorner->size_tot=SumHamC->size_tot;
        RightCorner = new TCSR<CPX>(nrowsofrightcorner,nrowsofrightcorner*2*ntriblock,SumHamC->findx);
        RightCorner->tr_full_to_sparse(Hrcpx_local,nrowsofrightcorner,2*ntriblock);
//        RightCorner->first_row=displcvecr[iam];
        RightCorner->first_row=SumHamC->first_row+SumHamC->size-nrowsofrightcorner;
        RightCorner->size_tot=SumHamC->size_tot;
        for (int iii=0;iii<RightCorner->n_nonzeros;iii++) RightCorner->index_j[iii]+=SumHamC->size_tot-2*ntriblock;
    if (nrowsofleftcorner>0) {
        delete[] Hlcpx_local;
    }
    if (nrowsofrightcorner>0) {
        delete[] Hrcpx_local;
    }
//        nrightcorner=RightCorner->edge_i[nrowsofrightcorner]-RightCorner->findx;i dont need this because nofnonzeros is set properly
/*    int* numofnnzleft;
    if (!iam) {
        int iposition=0;
        numofnnzleft=new int[nprocs];
        for (int iproc=0;iproc<nprocs;iproc++) {
            numofnnzleft[iproc]=LeftCorner->edge_i[iposition+1]-LeftCorner->edge_i[iposition];
            iposition+=nrowsleftcornerarray[iproc];
        }
    }
*/
        TCSR<CPX> *HamSig;
    if (nrowsofleftcorner==0 && nrowsofrightcorner==0) {
        HamSig=new TCSR<CPX>(SumHamC);
    } else if (nrowsofleftcorner==SumHamC->size && nrowsofrightcorner==0) {
        HamSig=new TCSR<CPX>(LeftCorner);
    } else if (nrowsofleftcorner==0 && nrowsofrightcorner==SumHamC->size) {
        HamSig=new TCSR<CPX>(RightCorner);
    } else { //or do i not use an else?
        int nhaminterior=SumHamC->edge_i[SumHamC->size-nrowsofrightcorner]-SumHamC->edge_i[nrowsofleftcorner];
        HamSig = new TCSR<CPX>(SumHamC->size,nhaminterior+LeftCorner->n_nonzeros+RightCorner->n_nonzeros,SumHamC->findx);
//warning for assemble i need index_i of all matrices DO I HAVE THAT ????????????????????????
        HamSig->assemble(LeftCorner,SumHamC,RightCorner,nrowsofleftcorner,nrowsofrightcorner);
    }
        delete LeftCorner;
        delete RightCorner;
        delete SumHamC;
// sparse solver
        CPX *Sol=new CPX[HamSig->size*2*nprotra];
        CPX *Inj=new CPX[HamSig->size*2*nprotra]();
        CPX *Injl=new CPX[HamSig->size*nprotra]();
        CPX *Injr=new CPX[HamSig->size*nprotra]();
// transposing makes scattering easier
//CHECK THIS IF IT IS RIGHT AND WORKS AFTER ADDING TRANSPOSING !!!
// especially the place at the target
    for (int iii=0;iii<nprocs;iii++) {
        nrowsleftcornerarray[iii]*=nprotra/(2*ntriblock);
        nrowsrightcornerarray[iii]*=nprotra/(2*ntriblock);
        displcvecl[iii]*=nprotra/(2*ntriblock);
        displcvecr[iii]*=nprotra/(2*ntriblock);
    }
    MPI_Scatterv(injl1,nrowsleftcornerarray,displcvecl,MPI_DOUBLE_COMPLEX,Injl,nrowsofleftcorner*nprotra,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
    MPI_Scatterv(injr1,nrowsrightcornerarray,displcvecr,MPI_DOUBLE_COMPLEX,&Injr[(HamSig->size-nrowsofrightcorner)*nprotra],nrowsofrightcorner*nprotra,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
    if (!iam) {
        delete[] injl;
        delete[] injr;
        delete[] injl1;
        delete[] injr1;
    }
    full_transpose(HamSig->size,nprotra,Injl,Inj);
    full_transpose(HamSig->size,nprotra,Injr,&Inj[HamSig->size*nprotra]);
    delete[] Injl;
    delete[] Injr;
    
//        c_zlacpy('A',ntriblock,nprotra,injl,ntriblock,Inj,HamSig->size_tot);
//        c_zlacpy('A',ntriblock,nprotra,injr,ntriblock,&Inj[HamSig->size_tot*(nprotra+1)-ntriblock],HamSig->size_tot);

        LinearSolver<CPX>* solver;
cout <<"INIT SOLVER"<<endl;
        solver = new Umfpack<CPX>(HamSig,MPI_COMM_WORLD);
        sabtime=get_time(d_zer);
cout <<"PREPARE SOLVER"<<endl;
        solver->prepare();
cout <<"SOLVE SOLVER"<<endl;
        solver->solve_equation(Sol,Inj,2*nprotra);
        cout << "TIME FOR WAVEFUNCTION SPARSE SOLVER WITH "<< ncells <<" UNIT CELLS " << get_time(sabtime) << endl;
        delete[] Inj;
        delete solver;
/*        delete HamSig;
// i think we need an allgather here
        CPX *Pmat = new CPX[ndof*ncells*ndof*ncells];
        CPX alphastep;
        if (nprocs==1)
            alphastep=z_one;
        else if (!iam || iam==nprocs-1)
            alphastep=CPX(energystep/2,d_zer);
        else 
            alphastep=CPX(energystep,d_zer);
        //sabtime=get_time(d_zer);
        //c_zgemm('N','C',ndof*ncells,ndof*ncells,2*nprotra,alphastep/(2*M_PI),Sol,ndof*ncells,Sol,ndof*ncells,z_zer,Pmat,ndof*ncells);
        //cout << "TIME FOR CONSTRUCTION OF FULL DENSITY MATRIX " << get_time(sabtime) << endl;
// transmission
        CPX *vecoutdof=new CPX[ndof];
        double transml=d_zer;
        for (int iband=1;iband<=bandwidth;iband++) {
            for (int ipro=0;ipro<nprotra;ipro++) {
                c_zgemv('N',ndof,ndof,z_one,&KScpx[ndofsq*(bandwidth+iband)],ndof,&Sol[ndof*(iband+ncells*ipro)],1,z_zer,vecoutdof,1);
                transml+=iband*4*M_PI*imag(c_zdotc(ndof,&Sol[ndof*ncells*ipro],1,vecoutdof,1));
            }
        }
        delete[] Sol;
        delete[] vecoutdof;
        cout << "Energy " << energy << " Transmission " << transml << endl;
*/
    } 
    else LOGCERR, exit(EXIT_FAILURE); 

/*
// SAB sets argc to 0 to avoid run
    argc=0;
    for(IT=1;IT<argc;IT++){
    
        init_parameters();

	double t0,t1;
	t0 = get_time(0.0);

        yyin = fopen(argv[IT],"r");
        yyrestart(yyin);
        yyparse();
        fclose(yyin);

	t0 = get_time(t0);
	MPI_Allreduce(&t0,&t1,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	if(!rank) printf("Time to parse the command file %lg (s)\n",t1);

        check_mpi(size,rank,parameter->CPU_per_vg_point,parameter->CPU_per_vd_point,\
		  parameter->CPU_per_kz_point,parameter->Nky,parameter->Nkz,\
		  parameter->CPU_ppoint);
    
        TransportBase *TR;

	cond1 = parameter->eta_res>0;
	cond2 = parameter->tb==20;
	cond3 = nanowire->NDim<=2;
	cond4 = nanowire->incoherent_injection;
        
        if(cond1||cond2||cond3||cond4){

            if(parameter->eta_res<1.0e-12) parameter->eta_res = 1.0e-12;

	    switch(parameter->transport_type){
	    case 0:
	        TR = new ElTransport<CPX>();
		break;
	    case 1:
	        TR = new PhTransport<CPX>();
		break;
	    }

        }else{

	    switch(parameter->transport_type){
	    case 0:
	        TR = new ElTransport<double>();
		break;
	    case 1:
	        TR = new PhTransport<double>();
		break;
	    }

        }

	if(parameter->CPU_per_sample==-1) parameter->CPU_per_sample = size;

	sample_id = (rank-rank%parameter->CPU_per_sample)/parameter->CPU_per_sample;

	MPI_Comm_split(MPI_COMM_WORLD,sample_id,rank,&sample_comm);
        
        TR->initialize(nanowire,parameter->mat_name,parameter->strain_model,parameter->mat_binary_x,\
		       parameter->lattype,parameter->tb,parameter->last_first,parameter->n_of_modes,\
		       parameter->Nk,parameter->bs_solver,parameter->Nky,parameter->Nkz,En,\
		       parameter->eta,parameter->eta_res,parameter->Temp,voltage,\
		       parameter->CPU_per_temp_point,parameter->CPU_per_vg_point,\
		       parameter->CPU_per_vd_point,parameter->CPU_per_kz_point,parameter->CPU_ppoint,\
		       parameter->CPU_per_wire,parameter->CPU_per_bc,parameter->NPROW,parameter->NPCOL,\
		       parameter->NPCS,parameter->spec_decomp,parameter->poisson_solver,\
		       parameter->poisson_criterion,parameter->poisson_iteration,\
		       parameter->max_proc_poisson,parameter->poisson_inner_criterion,\
		       parameter->poisson_inner_iteration,parameter->plot_all_k,sample_id,sample_comm);

        for(IC=0;IC<parameter->no_comm;IC++){
            TR->execute_task(parameter->command[IC],parameter->directory);
        }

	MPI_Comm_free(&sample_comm);

        delete TR;
        delete_parameters();

    }

    MPI_Finalize();
*/
    return 0;
}

/************************************************************************************************/
