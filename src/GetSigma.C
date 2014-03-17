#include <mpi.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include <valarray>
using namespace std;
#include "ScaLapack.H"
#include "InjectionFeast.H"
#include "Arpack.H"
#include "GetSigma.H"

BoundarySelfEnergy::BoundarySelfEnergy()
{
    sigmal = NULL;
    sigmar = NULL;
    spsigmaldist = NULL;
    spsigmardist = NULL;
    injl = NULL;
    injr = NULL;
    n_propagating=-1;

    master_rank=-1;;

    do_delete_H=0;
    do_delete_sigma=0;
    do_delete_spsigdist=0;
    do_delete_inj=0;
    do_delete_spainjdist=0;

    H0 = NULL;
    H1 = NULL;
    H1t = NULL;

    complexenergypoint=0;
    compute_inj=0;

    bandwidth=2;
    colzerothr=1.0E-12;
    eps_limit=1.0E-8;
    eps_decay=1.0E-4;
}

BoundarySelfEnergy::~BoundarySelfEnergy()
{
    if (do_delete_sigma) {
        delete[] sigmal;
        delete[] sigmar;
    }

    if (do_delete_inj) {
        delete[] injl;
        delete[] injr;
    }

    if (do_delete_spsigdist) {
        delete spsigmaldist;
        delete spsigmardist;
    }

    if (do_delete_spainjdist) {
        delete spainjldist;
        delete spainjrdist;
    }

    if (do_delete_H) {
        delete H0;
        delete H1;
        delete H1t;
    }
}

int BoundarySelfEnergy::Set_master(MPI_Comm matrix_comm,MPI_Comm boundary_comm)
{
    int matrix_rank, matrix_procs;
    MPI_Comm_size(matrix_comm,&matrix_procs);
    MPI_Comm_rank(matrix_comm,&matrix_rank);
    int boundary_rank, boundary_procs;
    MPI_Comm_size(boundary_comm,&boundary_procs);
    MPI_Comm_rank(boundary_comm,&boundary_rank);
    MPI_Group matrix_group;
    MPI_Group boundary_group;
    MPI_Comm_group(matrix_comm,&matrix_group);
    MPI_Comm_group(boundary_comm,&boundary_group);
    int *matrix_ranks_array = new int[matrix_procs];
    for (int i=0;i<matrix_procs;i++) matrix_ranks_array[i]=i;
    int *boundary_ranks_array = new int[matrix_procs];
    MPI_Group_translate_ranks(matrix_group,matrix_procs,matrix_ranks_array,boundary_group,boundary_ranks_array);
    master_rank=-1;
    for (int irank=0;irank<matrix_procs;irank++) {
        if (boundary_ranks_array[irank]==0 && master_rank==-1) {
            master_rank=matrix_ranks_array[irank];
        } else if (boundary_ranks_array[irank]==0 && master_rank!=-1) {
            return (LOGCERR, EXIT_FAILURE);
        }
    }
    delete[] matrix_ranks_array;
    delete[] boundary_ranks_array;
    return 0;
}

int BoundarySelfEnergy::Cutout(TCSR<CPX> *SumHamC,int contact,CPX energy,transport_methods::transport_method method,c_transport_type parameter_sab,MPI_Comm matrix_comm)
{
    energyp=energy;
    if (imag(energy)) complexenergypoint=1;
    if (method==transport_methods::WF) compute_inj=1;
    if (complexenergypoint && compute_inj) return (LOGCERR, EXIT_FAILURE);

    n_cells=parameter_sab.n_cells;
    ndof=SumHamC->size_tot/n_cells;
    int ntriblock=bandwidth*ndof;

    bandwidth=parameter_sab.bandwidth;
    colzerothr=parameter_sab.colzero_threshold;
    eps_limit=parameter_sab.eps_limit;
    eps_decay=parameter_sab.eps_decay;

    int iam,nprocs;
    MPI_Comm_size(matrix_comm,&nprocs);
    MPI_Comm_rank(matrix_comm,&iam);
    int *gathered_master_ranks = new int[nprocs];
    MPI_Allgather(&master_rank,1,MPI_INT,gathered_master_ranks,1,MPI_INT,matrix_comm);
    master_rank=-1;
    for (int irank=0;irank<nprocs;irank++) {
        if (gathered_master_ranks[irank]>=0) {
            master_rank=gathered_master_ranks[irank];
        }
    }
    delete[] gathered_master_ranks;

    if (contact==1) {
        TCSR<CPX> *H0cut = new TCSR<CPX>(SumHamC,0,ntriblock,0,ntriblock);
        TCSR<CPX> *H1cut = new TCSR<CPX>(SumHamC,0,ntriblock,ntriblock,ntriblock);
        H0 = new TCSR<CPX>(H0cut,master_rank,matrix_comm);
        H1 = new TCSR<CPX>(H1cut,master_rank,matrix_comm);
        delete H0cut;
        delete H1cut;
        if (iam==master_rank) {
            H0->shift_resize(0,ntriblock,0,ntriblock);
            H1->shift_resize(0,ntriblock,ntriblock,ntriblock);
            H1t = new TCSR<CPX>(H1->size,H1->n_nonzeros,H1->findx);
            H1t->sparse_transpose(H1);
        }
    } else if (contact==2) {
        TCSR<CPX> *H0cut = new TCSR<CPX>(SumHamC,SumHamC->size_tot-ntriblock,ntriblock,SumHamC->size_tot-ntriblock,ntriblock);
        TCSR<CPX> *H1tcut = new TCSR<CPX>(SumHamC,SumHamC->size_tot-ntriblock,ntriblock,SumHamC->size_tot-ntriblock-ntriblock,ntriblock);
        H0 = new TCSR<CPX>(H0cut,master_rank,matrix_comm);
        H1t = new TCSR<CPX>(H1tcut,master_rank,matrix_comm);
        delete H0cut;
        delete H1tcut;
        if (iam==master_rank) {
            H0->shift_resize(SumHamC->size_tot-ntriblock,ntriblock,SumHamC->size_tot-ntriblock,ntriblock);
            H1t->shift_resize(SumHamC->size_tot-ntriblock,ntriblock,SumHamC->size_tot-ntriblock-ntriblock,ntriblock);
            H1 = new TCSR<CPX>(H1t->size,H1t->n_nonzeros,H1t->findx);
            H1->sparse_transpose(H1t);
        }
    }
    do_delete_H=1;
    return 0;
}

void BoundarySelfEnergy::Distribute(TCSR<CPX> *SumHamC,MPI_Comm matrix_comm)
{
    int iam;
    MPI_Comm_rank(matrix_comm,&iam);
    int ntriblock=bandwidth*ndof;
    TCSR<CPX> *spsigmal = NULL;
    if (iam==master_rank) {
        spsigmal = new TCSR<CPX>(SumHamC->size_tot,ntriblock*ntriblock,SumHamC->findx);
        spsigmal->full_to_sparse(sigmal,ntriblock,ntriblock,0,0);
        delete[] sigmal;
        sigmal = NULL;
    }
    spsigmaldist = new TCSR<CPX>(SumHamC,spsigmal,master_rank,matrix_comm);
    if (iam==master_rank) delete spsigmal;
    TCSR<CPX> *spsigmar = NULL;
    if (iam==master_rank) {
        spsigmar = new TCSR<CPX>(SumHamC->size_tot,ntriblock*ntriblock,SumHamC->findx);
        spsigmar->full_to_sparse(sigmar,ntriblock,ntriblock,SumHamC->size_tot-ntriblock,SumHamC->size_tot-ntriblock);
        delete[] sigmar;
        sigmar = NULL;
    }
    spsigmardist = new TCSR<CPX>(SumHamC,spsigmar,master_rank,matrix_comm);
    if (iam==master_rank) delete spsigmar;
    do_delete_sigma=0;
    do_delete_spsigdist=1;

    if (compute_inj) {
        MPI_Bcast(&n_propagating,1,MPI_INT,master_rank,matrix_comm);
        TCSR<CPX> *spainjl = NULL;
        if (iam==master_rank) {
            spainjl = new TCSR<CPX>(SumHamC->size_tot,ntriblock*n_propagating,SumHamC->findx);
            spainjl->full_to_sparse(injl,ntriblock,n_propagating,0,0);
            delete[] injl;
            injl = NULL;
        }
        spainjldist = new TCSR<CPX>(SumHamC,spainjl,master_rank,matrix_comm);
        if (iam==master_rank) delete spainjl;
        TCSR<CPX> *spainjr = NULL;
        if (iam==master_rank) {
            spainjr = new TCSR<CPX>(SumHamC->size_tot,ntriblock*n_propagating,SumHamC->findx);
            spainjr->full_to_sparse(injr,ntriblock,n_propagating,SumHamC->size_tot-ntriblock,0);
            delete[] injr;
            injr = NULL;
        }
        spainjrdist = new TCSR<CPX>(SumHamC,spainjr,master_rank,matrix_comm);
        if (iam==master_rank) delete spainjr;
        do_delete_inj=0;
        do_delete_spainjdist=1;
        }
}

int BoundarySelfEnergy::Pascal(int r,int n)
{
    if( n == 0 )
        return 1;
    if( r == 0 || r == n )
        return 1;
    return Pascal( r - 1, n - 1 ) + Pascal( r, n - 1 );
}

int BoundarySelfEnergy::GetSigma(MPI_Comm boundary_comm)
{
    double d_one=1.0;
    double d_zer=0.0;
    CPX z_one=CPX(d_one,d_zer);
    CPX z_zer=CPX(d_zer,d_zer);
    double sabtime;
    int iinfo=0;
// inj_sign
    int inj_sign=-1;
// set parameters
    int ndofsq=ndof*ndof;
    int nblocksband=2*bandwidth+1;
    int ntriblock=bandwidth*ndof;
    int triblocksize=ntriblock*ntriblock;
    int boundary_rank;
    MPI_Comm_rank(boundary_comm,&boundary_rank);
// global for this class
    CPX *eigvalcpx;
    CPX *eigvecc;
    int neigval;
// copy block
    CPX* KScpx;
    if (boundary_rank==0) {
        KScpx=new CPX[ndofsq*nblocksband];
        H1t->sparse_to_full(KScpx,ndof,bandwidth*ndof);
        H0->sparse_to_full(&KScpx[bandwidth*ndofsq],ndof,bandwidth*ndof);
        H1->sparse_to_full(&KScpx[2*bandwidth*ndofsq],ndof,ndof);
    }
// after i tested if it makes a difference how i assemble the tridiagonalblocks i will remove everything there
// and build the full kscpx, which i need for my old method, by sparse_to_full out of the h0 and h1, and i build
// the sparse amat and bmat for feast out of h0 and h1 like i do below
/*
    CPX* H0cpx;
    CPX* H1cpx;
    CPX* H1cpxt;
    if (boundary_rank==0) {
        H0cpx=new CPX[triblocksize];
        H0->sparse_to_full(H0cpx,ntriblock,ntriblock);
        H1cpx=new CPX[triblocksize];
        H1->sparse_to_full(H1cpx,ntriblock,ntriblock);
        H1cpxt=new CPX[triblocksize];
        H1t->sparse_to_full(H1cpxt,ntriblock,ntriblock);
    }
*/
// assemble tridiagonalblocks, here just to see if there is an influence if i take the first unit cell and replicate or take the first two unit cells
/*
    if (boundary_rank==0) {
        for (int ibandwidth=0;ibandwidth<bandwidth;ibandwidth++)
            for (int jbandwidth=0;jbandwidth<bandwidth;jbandwidth++)
                for (int jdof=0;jdof<ndof;jdof++)
                    c_zcopy(ndof,&KScpx[(bandwidth-ibandwidth+jbandwidth)*ndofsq+jdof*ndof],1,&H0cpx[(bandwidth*(jdof+jbandwidth*ndof)+ibandwidth)*ndof],1);
        for (int iz=0;iz<triblocksize;iz++) H1cpx[iz]=z_zer;
        for (int ibandwidth=0;ibandwidth<bandwidth;ibandwidth++)
            for (int jbandwidth=0;jbandwidth<=ibandwidth;jbandwidth++)
                for (int jdof=0;jdof<ndof;jdof++)
                    c_zcopy(ndof,&KScpx[(2*bandwidth-ibandwidth+jbandwidth)*ndofsq+jdof*ndof],1,&H1cpx[(bandwidth*(jdof+jbandwidth*ndof)+ibandwidth)*ndof],1);
        full_transpose(ntriblock,ntriblock,H1cpx,H1cpxt);
//        delete H0;
//        delete H1;
//        delete H1t;
        H0->full_to_sparse(H0cpx,ntriblock,ntriblock);
        H1->full_to_sparse(H1cpx,ntriblock,ntriblock);
        H1t->full_to_sparse(H1cpxt,ntriblock,ntriblock);
//        delete[] H0cpx;
//        delete[] H1cpx;
//        delete[] H1cpxt;
    }
*/
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
int arpack=0;
// FEAST
if (0) {
    InjectionFeast<CPX> *k_inj = new InjectionFeast<CPX>();
    int neigfeast=200;
    k_inj->initialize(2*ntriblock,neigfeast);
    eigvalcpx=new CPX[neigfeast];
    CPX* eigvec=new CPX[2*ntriblock*neigfeast];
/*
    CPX *Amat=new CPX[2*ntriblock*2*ntriblock]();
    CPX *Bmat=new CPX[2*ntriblock*2*ntriblock]();
    for (int idiag=0;idiag<(2*bandwidth-1)*ndof;idiag++) {
        Amat[2*ntriblock*ndof+idiag+(2*ntriblock)*idiag]=z_one;
        Bmat[idiag+(2*ntriblock)*idiag]=z_one;
    }
    c_zlacpy('A',ndof,2*ntriblock,KScpx,ndof,&Amat[(2*bandwidth-1)*ndof],2*ntriblock);
    c_zlacpy('A',ndof,ndof,&KScpx[2*ntriblock*ndof],ndof,&Bmat[(2*ntriblock+1)*(2*bandwidth-1)*ndof],2*ntriblock);
    c_zscal(2*ntriblock*ndof,CPX(-1.0,0.0),&Bmat[2*ntriblock*(2*bandwidth-1)*ndof],1);
    TCSR<CPX> *spA = new TCSR<CPX>(2*ntriblock,2*triblocksize+ntriblock,1);
    spA->full_to_sparse(Amat,2*ntriblock,2*ntriblock);
    delete[] Amat;
    TCSR<CPX> *spB = new TCSR<CPX>(2*ntriblock,triblocksize+ntriblock,1);
    spB->full_to_sparse(Bmat,2*ntriblock,2*ntriblock);
    delete[] Bmat;
*/
    TCSR<CPX> *spA = NULL;
    TCSR<CPX> *spB = NULL;
    if (!boundary_rank) {
        TCSR<CPX> *H1u = new TCSR<CPX>(H1,0,ndof,0,ndof);
        c_zscal(H1u->n_nonzeros,-z_one,H1u->nnz,1);
        for (int i_ele=0;i_ele<H1u->n_nonzeros;i_ele++) {
            H1u->index_j[i_ele]+=2*ntriblock-ndof;
        }
        spB = new TCSR<CPX>(2*ntriblock,H1u->n_nonzeros+2*ntriblock-ndof,H1u->findx);
        for (int i_ele=0;i_ele<2*ntriblock-ndof;i_ele++) {
            spB->nnz[i_ele]=z_one;
            spB->index_j[i_ele]=i_ele+spB->findx;
            spB->index_i[i_ele]=1;
        }
        c_zcopy(H1u->n_nonzeros,H1u->nnz,1,&spB->nnz[2*ntriblock-ndof],1);
        c_icopy(H1u->n_nonzeros,H1u->index_j,1,&spB->index_j[2*ntriblock-ndof],1);
        c_icopy(ndof,H1u->index_i,1,&spB->index_i[2*ntriblock-ndof],1);
        delete H1u;
        spB->get_row_edge();
        spB->get_diag_pos();
        TCSR<CPX> *H1tu = new TCSR<CPX>(H1t,0,ndof,0,ntriblock);
        H1tu->size_tot=2*ntriblock;
        TCSR<CPX> *H0u = new TCSR<CPX>(H0,0,ndof,0,ntriblock);
        H0u->size_tot=2*ntriblock;
        for (int i_ele=0;i_ele<H0u->n_nonzeros;i_ele++) {
            H0u->index_j[i_ele]+=ntriblock;
        }
        TCSR<CPX> *Hu = new TCSR<CPX>(z_one,H1tu,z_one,H0u);
        delete H1tu;
        delete H0u;
        spA = new TCSR<CPX>(2*ntriblock,Hu->n_nonzeros+2*ntriblock-ndof,Hu->findx);
        for (int i_ele=0;i_ele<2*ntriblock-ndof;i_ele++) {
            spA->nnz[i_ele]=z_one;
            spA->index_j[i_ele]=i_ele+ndof+spA->findx;
            spA->index_i[i_ele]=1;
        }
        c_zcopy(Hu->n_nonzeros,Hu->nnz,1,&spA->nnz[2*ntriblock-ndof],1);
        c_icopy(Hu->n_nonzeros,Hu->index_j,1,&spA->index_j[2*ntriblock-ndof],1);
        c_icopy(ndof,Hu->index_i,1,&spA->index_i[2*ntriblock-ndof],1);
        delete Hu;
        spA->get_row_edge();
        spA->get_diag_pos();
    }
    sabtime=get_time(d_zer);
    k_inj->calc_kphase(spA,spB,2*ntriblock,eigvalcpx,eigvec,&neigval,boundary_comm,&iinfo);
    cout << "TIME FOR FEAST " << get_time(sabtime) << endl;
    delete spA;
    delete spB;
    delete k_inj;
    eigvecc=new CPX[ndof*neigval];//IST NEIGVAL AUCH RICHTIG? SIND DIE AUCH IN EIGVEC RICHTIG ANGEORDNET
    c_zlacpy('A',ndof,neigval,eigvec,2*ntriblock,eigvecc,ndof);
    delete[] eigvec;
//ARPACK
} else if (arpack) {
    sabtime=get_time(d_zer);
    double mmtime=d_zer;
    double artime=d_zer;
    double tmptime;
    Arpack* AR = new Arpack();
/*
if (worldrank==10) {
stringstream mysstreamy;
mysstreamy << "ksmat" << worldrank;
ofstream myfily(mysstreamy.str().c_str());
for (int imat=0;imat<ndofsq*5;imat++) myfily << real(KScpx[imat]) << endl;
myfily.close();
}
*/
// /*
    CPX shiftvec[] = {z_one,-z_one};
    int narpoints=2;
    int neigarpack=78;
    int neigarpack2=52;
    int neigarpack3=0;
// */
/*
    CPX shiftvec[] = {z_one,-z_one};
    int narpoints=2;
    int neigarpack=110;
    int neigarpack2=22;
    int neigarpack3=0;
*/
/*
    CPX shiftvec[] = {z_one,-z_one,CPX(d_zer,d_one),-CPX(d_zer,d_one)};
    int narpoints=4;
    int neigarpack=198;
    int neigarpack2=0;
    int neigarpack3=0;
*/
    eigvalcpx = new CPX[neigarpack*narpoints+neigarpack2+neigarpack3];
    CPX* eigvec = new CPX[2*ntriblock*(neigarpack*narpoints+neigarpack2+neigarpack3)];
    CPX* tmpmat = new CPX[ndofsq*nblocksband];
    double* inpmat = new double[ndofsq*nblocksband];
    double *trmat = new double[2*ntriblock*ndof];
    int *pivarraya=new int[ndof];
    for (int ishifts=0;ishifts<narpoints;ishifts++) {
        for (int imat=0;imat<ndofsq*nblocksband;imat++) tmpmat[imat]=z_zer;
        for (int iresult=0;iresult<nblocksband;iresult++) {
            for (int iksmat=0;iksmat<=iresult;iksmat++) {
                double pascalfactor=double(Pascal(nblocksband-1-iresult,nblocksband-1-iresult+iksmat));
                c_zaxpy(ndofsq,pascalfactor*pow(shiftvec[ishifts],iksmat),&KScpx[(nblocksband-1-iresult+iksmat)*ndofsq],1,&tmpmat[iresult*ndofsq],1);
            }
        }
        if (imag(shiftvec[ishifts])) {
            c_zgetrf(ndof,ndof,&tmpmat[2*bandwidth*ndofsq],ndof,pivarraya,&iinfo);
            if (iinfo) return (LOGCERR, EXIT_FAILURE);
            c_zgetrs('N',ndof,ndof*2*bandwidth,&tmpmat[2*bandwidth*ndofsq],ndof,pivarraya,tmpmat,ndof,&iinfo);
            if (iinfo) return (LOGCERR, EXIT_FAILURE);
            CPX *trmatc = new CPX[2*ntriblock*ndof];
            full_transpose(ntriblock*2,ndof,tmpmat,trmatc);
            //c_dcopy(2*ntriblock*ndof,(double*)trmatc,2,inpmat,1);
            tmptime=get_time(0.0);
            AR->eigs(&eigvec[ishifts*2*ntriblock*neigarpack],&eigvalcpx[ishifts*neigarpack],trmatc,ndof,2*ntriblock,neigarpack,mmtime);
            artime+=get_time(tmptime);
            delete[] trmatc;
        } else {
            c_dcopy(ndofsq*nblocksband,(double*)tmpmat,2,inpmat,1);
            c_dgetrf(ndof,ndof,&inpmat[2*bandwidth*ndofsq],ndof,pivarraya,&iinfo);
            if (iinfo) return (LOGCERR, EXIT_FAILURE);
            c_dgetrs('N',ndof,ndof*2*bandwidth,&inpmat[2*bandwidth*ndofsq],ndof,pivarraya,inpmat,ndof,&iinfo);
            if (iinfo) return (LOGCERR, EXIT_FAILURE);
            full_transpose(ntriblock*2,ndof,inpmat,trmat);
/*
stringstream mysstream;
mysstream << "arpackmat" << worldrank << "_" << ishifts << "new";
//ofstream myfile;
//myfile.open(mysstream.str().c_str());
ifstream my_file(mysstream.str().c_str());
int yesmyfile=1;
if (my_file) yesmyfile=0;
my_file.close();
if (yesmyfile) {
ofstream myfile(mysstream.str().c_str());
for (int imat=0;imat<ndofsq*nblocksband;imat++) myfile << real(tmpmat[imat]) << " " << imag(tmpmat[imat]) << endl;
myfile.close();
}
*/
            tmptime=get_time(0.0);
            AR->eigs(&eigvec[ishifts*2*ntriblock*neigarpack],&eigvalcpx[ishifts*neigarpack],trmat,ndof,2*ntriblock,neigarpack,mmtime);
            artime+=get_time(tmptime);
        }
        for (int ieigs=0;ieigs<neigarpack;ieigs++) {
            eigvalcpx[ishifts*neigarpack+ieigs]=z_one/eigvalcpx[ishifts*neigarpack+ieigs]+shiftvec[ishifts];
//}else{
//            eigvalcpx[ishifts*neigarpack+ieigs]=2.0*z_one/eigvalcpx[ishifts*neigarpack+ieigs];
//            eigvalcpx[ishifts*neigarpack+ieigs]+=sqrt(eigvalcpx[ishifts*neigarpack+ieigs]*eigvalcpx[ishifts*neigarpack+ieigs]-abs(shiftvec[ishifts])*abs(shiftvec[ishifts]));
//}
        }
    }
    delete[] tmpmat;
    delete[] inpmat;
    delete[] trmat;
    delete[] pivarraya;
//start arpack for square old method
/*
    ndof*=2;
    ndofsq*=4;
    bandwidth/=2;
    nblocksband=3;
    int *pivarraya2=new int[ndof];
    CPX* KScpxbig = new CPX[ndofsq*nblocksband];
    H0->sparse_to_full(&KScpxbig[ndofsq],ntriblock,ntriblock);
    H1->sparse_to_full(&KScpxbig[2*ndofsq],ntriblock,ntriblock);
    H1t->sparse_to_full(KScpxbig,ntriblock,ntriblock);
    CPX *tmpmat2 = new CPX[ndofsq*nblocksband]();
    for (int iresult=0;iresult<nblocksband;iresult++) {
        for (int iksmat=0;iksmat<=iresult;iksmat++) {
            double pascalfactor=double(Pascal(nblocksband-1-iresult,nblocksband-1-iresult+iksmat));
            c_zaxpy(ndofsq,pascalfactor*pow(-z_one,iksmat),&KScpxbig[(nblocksband-1-iresult+iksmat)*ndofsq],1,&tmpmat2[iresult*ndofsq],1);
        }
    }
    double *inpmat2 = new double[ndofsq*nblocksband];
    c_dcopy(ndofsq*nblocksband,(double*)tmpmat2,2,inpmat2,1);
    c_dgetrf(ndof,ndof,&inpmat2[2*bandwidth*ndofsq],ndof,pivarraya2,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    c_dgetrs('N',ndof,ndof*2*bandwidth,&inpmat2[2*bandwidth*ndofsq],ndof,pivarraya2,inpmat2,ndof,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    double *trmat2 = new double[ndofsq*nblocksband];
    full_transpose(ndof*2,ndof,inpmat2,trmat2);
//    CPX *trmatc2 = new CPX[2*ndof*ndof];
//    c_dcopy(2*ndof*ndof,trmat2,1,(double*)trmatc2,2);
    tmptime=get_time(0.0);
    AR->eigs(&eigvec[2*ntriblock*narpoints*neigarpack],&eigvalcpx[narpoints*neigarpack],trmat2,ndof,2*ntriblock,neigarpack2,mmtime);
    artime+=get_time(tmptime);
    delete[] pivarraya2;
    delete[] KScpxbig;
    bandwidth*=2;
    nblocksband=5;
    ndof/=2;
    ndofsq/=4;
//    delete[] trmatc2;
    delete[] tmpmat2;
    delete[] inpmat2;
    delete[] trmat2;
    CPX *poly = new CPX[nblocksband];
    CPX* vecoutdof = new CPX[ndof];
    CPX *rootmat = new CPX[4*4];
    double *workydouble=new double[2*4];
    int lworkycc=2*4;
    CPX *workycc=new CPX[lworkycc];
    for (int ieigs=0;ieigs<neigarpack2;ieigs++) {
        eigvalcpx[narpoints*neigarpack+ieigs]=sqrt(z_one/eigvalcpx[narpoints*neigarpack+ieigs]-z_one);
        for (int ip=0;ip<nblocksband;ip++) {
            c_zgemv('N',ndof,ndof,z_one,&KScpx[ndofsq*ip],ndof,&eigvec[(narpoints*neigarpack+ieigs)*2*ntriblock],1,z_zer,vecoutdof,1);
            poly[ip]=c_zdotc(ndof,&eigvec[(narpoints*neigarpack+ieigs)*2*ntriblock],1,vecoutdof,1);
        }
        for (int ip=0;ip<4*4;ip++) {
            rootmat[ip]=CPX(0.0,0.0);
        }
        for (int ip=1;ip<5;ip++) {
            rootmat[(ip-1)*4]=poly[4-ip]/poly[4];
        }
        for (int ip=1;ip<4;ip++) {
            rootmat[(ip-1)*4+ip]=CPX(1.0,0.0);
        }
        c_zgeev('N','N',4,rootmat,4,poly,NULL,1,NULL,1,workycc,lworkycc,workydouble,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        for (int ip=0;ip<4;ip++) {
            double eigcmpr = real(pow(eigvalcpx[narpoints*neigarpack+ieigs],2));
            double eigcmpi = imag(pow(eigvalcpx[narpoints*neigarpack+ieigs],2));
            double eigpolr = real(pow(poly[ip],2));
            double eigpoli = imag(pow(poly[ip],2));
            double tolcond = 1.0E-3;
            int conditionr = (eigpolr<eigcmpr+tolcond && eigpolr>eigcmpr-tolcond);
            int conditioni = (eigpoli<eigcmpi+tolcond && eigpoli>eigcmpi-tolcond);
            if (conditionr && conditioni) {
                if (real(poly[ip]/eigvalcpx[narpoints*neigarpack+ieigs])<0) {
//                  cout<<"SIGN CHANGE "<<poly[ip]<<" FOR "<<eigvalcpx[narpoints*neigarpack+ieigs]<<" IN "<<ieigs<<" ON "<<worldrank<<endl;
                    eigvalcpx[narpoints*neigarpack+ieigs]*=-1.0;
                }
                break;
            }
        }
    }
    delete[] workycc;
    delete[] workydouble;
    delete[] poly;
    delete[] vecoutdof;
    delete[] rootmat;
if (0) {
stringstream mysstream;
mysstream << "arpackmat" << worldrank;
//ofstream myfile;
//myfile.open(mysstream.str().c_str());
ifstream my_file(mysstream.str().c_str());
int yesmyfile=1;
if (my_file) yesmyfile=0;
my_file.close();
if (yesmyfile&&worldrank==26) {
ofstream myfile(mysstream.str().c_str());
for (int imat=0;imat<624*2*624*2*2;imat++) myfile << real(tmpmat2[imat]) << endl;
myfile.close();
}
}
//end arpack for square
*/
//start arpack for square new method
    if (neigarpack2) {
        double *matB = new double[1*ndofsq];
        double *matM = new double[4*ndofsq];
        double *matC = new double[8*ndofsq]();
        c_dcopy(ndofsq,(double*)&KScpx[2*ndofsq],2,matB,1);
        c_daxpy(ndofsq,-1.0,(double*)&KScpx[0*ndofsq],2,matB,1);
        c_daxpy(ndofsq,-1.0,(double*)&KScpx[4*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[0],2*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[2*ndofsq+ndof],2*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[3*ndofsq],2,matB,1);
        c_daxpy(ndofsq,-1.0,(double*)&KScpx[1*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[2*ndofsq],2*ndof);
        c_dscal(ndofsq,-1.0,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[ndof],2*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[0*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[0],2*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[2*ndofsq+ndof],2*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[1*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[2*ndofsq],2*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[3*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[4*ndofsq+ndof],2*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[4*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[4*ndofsq],2*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[6*ndofsq+ndof],2*ndof);
        delete[] matB;
        int *pivarrayM = new int[2*ndof];
        c_dgetrf(2*ndof,2*ndof,matM,2*ndof,pivarrayM,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_dgetrs('N',2*ndof,4*ndof,matM,2*ndof,pivarrayM,matC,2*ndof,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] pivarrayM;
        delete[] matM;
        tmptime=get_time(0.0);
        AR->eigs(&eigvec[2*ntriblock*narpoints*neigarpack],&eigvalcpx[narpoints*neigarpack],matC,2*ndof,2*ntriblock,neigarpack2,mmtime);
        artime+=get_time(tmptime);
        delete[] matC;
        for (int ieigs=0;ieigs<neigarpack2;ieigs++) {
            eigvalcpx[narpoints*neigarpack+ieigs]=sqrt(z_one/eigvalcpx[narpoints*neigarpack+ieigs]-z_one);
            CPX eigvalfromvec = eigvec[2*ntriblock*(narpoints*neigarpack+ieigs)+ndof]/eigvec[2*ntriblock*(narpoints*neigarpack+ieigs)];
            if (real(eigvalfromvec/eigvalcpx[narpoints*neigarpack+ieigs])<0) {
                eigvalcpx[narpoints*neigarpack+ieigs]*=-1.0;
            }
        }
    }
//end arpack for square new method
//start arpack for pow3
    if (neigarpack3) {
        double *matB = new double[ndofsq]();
        double *matM = new double[16*ndofsq]();
        double *matC = new double[16*ndofsq]();
        for (int idiag=0;idiag<ndof;idiag++) {
            matB[idiag+ndof*idiag]=1.0;
        }
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[0*4*ndofsq+3*ndof],4*ndof);
        for (int idiag=0;idiag<ndof;idiag++) {
            matB[idiag+ndof*idiag]=-1.0;
        }
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[3*4*ndofsq+3*ndof],4*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[0*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[0*4*ndofsq+0*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[1*4*ndofsq+1*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[2*4*ndofsq+2*ndof],4*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[1*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[0*4*ndofsq+2*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[1*4*ndofsq+0*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[2*4*ndofsq+1*ndof],4*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[2*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[0*4*ndofsq+1*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[1*4*ndofsq+2*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matC[2*4*ndofsq+0*ndof],4*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[3*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[0*4*ndofsq+0*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[1*4*ndofsq+1*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[2*4*ndofsq+2*ndof],4*ndof);
        c_dcopy(ndofsq,(double*)&KScpx[4*ndofsq],2,matB,1);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[1*4*ndofsq+0*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[2*4*ndofsq+1*ndof],4*ndof);
        c_dlacpy('A',ndof,ndof,matB,ndof,&matM[3*4*ndofsq+2*ndof],4*ndof);
        delete[] matB;
        c_daxpy(16*ndofsq,-1.0,matC,1,matM,1);
        int *pivarrayM = new int[4*ndof];
        c_dgetrf(4*ndof,4*ndof,matM,4*ndof,pivarrayM,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_dgetrs('N',4*ndof,4*ndof,matM,4*ndof,pivarrayM,matC,4*ndof,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] pivarrayM;
        delete[] matM;
/*
if (worldrank==46){
ofstream myfile("powthree46");
for (int ii=0;ii<16*ndofsq;ii++) myfile << matC[ii] << endl;
myfile.close();
}
MPI_Barrier(MPI_COMM_WORLD);
exit(0);
*/
        tmptime=get_time(0.0);
        AR->eigs(&eigvec[2*ntriblock*(narpoints*neigarpack+neigarpack2)],&eigvalcpx[narpoints*neigarpack+neigarpack2],matC,2*ntriblock,2*ntriblock,neigarpack2,mmtime);
        artime+=get_time(tmptime);
        delete[] matC;
        for (int ieigs=0;ieigs<neigarpack3;ieigs++) {
            eigvalcpx[narpoints*neigarpack+neigarpack2+ieigs]=pow(z_one/eigvalcpx[narpoints*neigarpack+neigarpack2+ieigs]-z_one,1.0/3.0);
            CPX eigvalfromvec = eigvec[2*ntriblock*(narpoints*neigarpack+neigarpack2+ieigs)]/eigvec[2*ntriblock*(narpoints*neigarpack+neigarpack2+ieigs)+ndof];
            if (real(eigvalfromvec/eigvalcpx[narpoints*neigarpack+neigarpack2+ieigs])<0) {
                if (imag(eigvalfromvec/eigvalcpx[narpoints*neigarpack+neigarpack2+ieigs])>0) {
                    eigvalcpx[narpoints*neigarpack+neigarpack2+ieigs]*=0.5*(-1.0+sqrt(3));
                } else {
                    eigvalcpx[narpoints*neigarpack+neigarpack2+ieigs]*=0.5*(-1.0-sqrt(3));
                }
            }
        }
    }
//end arpack for pow3
/*
stringstream mysstream;
mysstream << "eigvals" << worldrank;
ofstream myfile(mysstream.str().c_str());
for (int ieigout=0;ieigout<narpoints*neigarpack;ieigout++) myfile << real(eigvalcpx[ieigout]) << " " << imag(eigvalcpx[ieigout])  << endl;
myfile.close();
*/
    for (int ieigs=0;ieigs<narpoints*neigarpack+neigarpack2+neigarpack3;ieigs++) {
        eigvalcpx[ieigs]=z_one/(eigvalcpx[ieigs]-z_one);
    }
    neigval=narpoints*neigarpack+neigarpack2+neigarpack3;
    delete AR;
    eigvecc = new CPX[ndof*neigval];
    c_zlacpy('A',ndof,neigval,eigvec,2*ntriblock,eigvecc,ndof);
    delete[] eigvec;
    cout << worldrank << " NEEDED TIME FOR DIAGONALIZATION " << get_time(sabtime) << " INCLUDING ARPACK " << artime << " INCLUDING MATVECMULT "<< mmtime << endl;
// ALL EIGENVALUES
} else {
if (!boundary_rank) {
    CPX *mats=new CPX[ndofsq];
    c_zcopy(ndofsq,KScpx,1,mats,1);
    for (int ibandwidth=1;ibandwidth<nblocksband;ibandwidth++)
        c_zaxpy(ndofsq,d_one,&KScpx[ibandwidth*ndofsq],1,mats,1);//easier is possible because mats is symmetric
    CPX *matb=new CPX[2*bandwidth*ndofsq];//lowest block of right hand side matrix
// i think this is right for the general case
    c_zcopy(ndofsq,KScpx,1,matb,1);
    for (int iband=1;iband<=2*(bandwidth-1);iband++) {
        c_zcopy(ndofsq,&matb[(iband-1)*ndofsq],1,&matb[iband*ndofsq],1);
        c_zaxpy(ndofsq,d_one,&KScpx[iband*ndofsq],1,&matb[iband*ndofsq],1);
    }
    c_zcopy(ndofsq,&KScpx[2*bandwidth*ndofsq],1,&matb[(2*bandwidth-1)*ndofsq],1);
    c_zscal(ndofsq,-z_one,&matb[(2*bandwidth-1)*ndofsq],1);
// now do the inversion
    int *pivarrayn=new int[ndof];
    double workyytest;
    if (complexenergypoint) {
// the matrix mats is not hermitian for complex energy point i think
        sabtime=get_time(d_zer);
        c_zgetrf(ndof,ndof,mats,ndof,pivarrayn,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_zgetrs('N',ndof,ndof*2*bandwidth,mats,ndof,pivarrayn,matb,ndof,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
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
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] workyn;
        delete[] matsreal;
        for (int imat=0;imat<2*bandwidth*ndofsq;imat++) matb[imat]=CPX(matbreal[imat],d_zer);
        delete[] matbreal;
        cout << "TIME FOR REAL INVERSION " << get_time(sabtime) << endl;
    }
    delete[] pivarrayn;
    delete[] mats;
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
    if (nindzerocoln+nindnzcoln!=2*ntriblock) return (LOGCERR, EXIT_FAILURE);
    CPX *matmn=new CPX[nindnzcoln*nindnzcoln];
    CPX *matmr=new CPX[nindnzcoln*nindrzcoln];
    for (int jindnzcol=0;jindnzcol<nindnzcoln;jindnzcol++)
        for (int iindnzcol=0;iindnzcol<nindnzcoln;iindnzcol++)
            matmn[jindnzcol*nindnzcoln+iindnzcol]=matn[indnzcolvecn[jindnzcol]*2*ntriblock+indnzcolvecn[iindnzcol]];
    for (int jindnzcol=0;jindnzcol<nindnzcoln;jindnzcol++)
        for (int iindnzcol=0;iindnzcol<nindrzcoln;iindnzcol++)
            matmr[jindnzcol*nindrzcoln+iindnzcol]=matn[indnzcolvecn[jindnzcol]*2*ntriblock+indrzcolvecn[iindnzcol]];
    delete[] matn;
    int blockbwpos=-1;
    int blockbwpoe=-1;
    for (int ipos=0;ipos<nindnzcoln;ipos++) {
        if (indnzcolvecn[ipos]>=bandwidth*ndof && blockbwpos==-1) blockbwpos=ipos;
        if (indnzcolvecn[ipos]<=((bandwidth+1)*ndof-1)) blockbwpoe=ipos;
    }
    if (blockbwpos<0 || blockbwpoe<0) return (LOGCERR, EXIT_FAILURE);
// get eigenvalues and eigenvectors
    eigvalcpx=new CPX[nindnzcoln];
    CPX *eigveccfull=new CPX[nindnzcoln*nindnzcoln];
    if (complexenergypoint) {
        CPX cdummy;
        CPX workyctest;
        double *workdouble=new double[2*nindnzcoln];
        c_zgeev('N','V',nindnzcoln,matmn,nindnzcoln,eigvalcpx,&cdummy,1,eigveccfull,nindnzcoln,&workyctest,-1,workdouble,&iinfo);
        int lworkyc=int(real(workyctest));
        CPX *workyc=new CPX[lworkyc];
        c_zgeev('N','V',nindnzcoln,matmn,nindnzcoln,eigvalcpx,&cdummy,1,eigveccfull,nindnzcoln,workyc,lworkyc,workdouble,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        cout << "TIME FOR CPX DIAGONALIZATION " << get_time(sabtime) << endl;
        delete[] workdouble;
        delete[] workyc;
//        sabtime=get_time(d_zer);
//        Arpack* AR = new Arpack();
//        AR->eigs(eigveccfull,eigvalcpx,matmn,nindnzcoln,300);
//        cout << "TIME FOR ARPACK DIAGONALIZATION " << get_time(sabtime) << endl;
//        delete AR;
    } else {
        double *eigvalreal=new double[nindnzcoln];
        double *eigvalimag=new double[nindnzcoln];
        double *eigvec=new double[nindnzcoln*nindnzcoln];
        double dddummy;
        double *matmnreal=new double[nindnzcoln*nindnzcoln];
        sabtime=get_time(d_zer);
        c_dcopy(nindnzcoln*nindnzcoln,(double*)matmn,2,matmnreal,1);
        c_dgeev('N','V',nindnzcoln,matmnreal,nindnzcoln,eigvalreal,eigvalimag,&dddummy,1,eigvec,nindnzcoln,&workyytest,-1,&iinfo);
        int lworky=int(workyytest);
        double *worky=new double[lworky];
        c_dgeev('N','V',nindnzcoln,matmnreal,nindnzcoln,eigvalreal,eigvalimag,&dddummy,1,eigvec,nindnzcoln,worky,lworky,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
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
            else return (LOGCERR, EXIT_FAILURE);
        } // END WHILE
        delete[] eigvec;
        for (int ieigval=0;ieigval<nindnzcoln;ieigval++)
            eigvalcpx[ieigval]=CPX(eigvalreal[ieigval],eigvalimag[ieigval]);
        delete[] eigvalreal;
        delete[] eigvalimag;
    }
    delete[] matmn;
    eigvecc=new CPX[ndof*nindnzcoln];
    if (blockbwpoe-blockbwpos==ndof-1) {
        c_zlacpy('A',ndof,nindnzcoln,&eigveccfull[blockbwpos],nindnzcoln,eigvecc,ndof);
    } else {
// RECONSTRUCT
        sabtime=get_time(d_zer);
        CPX *eigveccrec=new CPX[nindnzcoln*nindrzcoln];
        c_zgemm('N','N',nindrzcoln,nindnzcoln,nindnzcoln,z_one,matmr,nindrzcoln,eigveccfull,nindnzcoln,z_zer,eigveccrec,nindrzcoln);
        for (int iind=0;iind<nindnzcoln;iind++)
            c_zscal(nindrzcoln,z_one/eigvalcpx[iind],&eigveccrec[nindrzcoln*iind],1);
        for (int iind=blockbwpos;iind<=blockbwpoe;iind++)
            c_zcopy(nindnzcoln,&eigveccfull[iind],nindnzcoln,&eigvecc[indnzcolvecn[iind]-bandwidth*ndof],ndof);
        for (int iind=0;iind<nindrzcoln;iind++)
            c_zcopy(nindnzcoln,&eigveccrec[iind],nindrzcoln,&eigvecc[indrzcolvecn[iind]-bandwidth*ndof],ndof);
        delete[] eigveccrec;
        cout << "TIME FOR RECONSTRUCTION " << get_time(sabtime) << endl;
    }
    delete[] indnzcolvecn;
    delete[] indrzcolvecn;
    delete[] eigveccfull;
    delete[] matmr;
    neigval=nindnzcoln;
} // END ONLY ON BOUNDARY MASTER
} // END METHOD FOR ALL EIGENVALUES
if (!boundary_rank) {
// DETERMINE TYPE OF EIGENVALUE/VECTOR
    CPX *lambdavec=new CPX[neigval];
    int *dectravec=new int[neigval];
    int *decrefvec=new int[neigval];
    int *protravec=new int[neigval];
    int *prorefvec=new int[neigval];
    double *veltra=new double[neigval];
    double *velref=new double[neigval];
    int ndectra=0;
    int ndecref=0;
    int nprotra=0;
    int nproref=0;
    CPX *matcdof=new CPX[ndofsq];
    CPX *vecout=new CPX[ndof];
    for (int iindnzcoln=0;iindnzcoln<neigval;iindnzcoln++) {
        double eigr=real(eigvalcpx[iindnzcoln]);
        double eigi=imag(eigvalcpx[iindnzcoln]);
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
                } else return (LOGCERR, EXIT_FAILURE);
           } // END IF decaying or propagating
           lambdavec[iindnzcoln]=lambda;
        } // END IF k not infinite
    } // END FOR
    delete[] eigvalcpx; 
    delete[] matcdof;
    delete[] vecout;
if (nprotra) { cout<<worldrank<<" AT "<<real(energyp)<<" NPROTRA "<<nprotra<<" NPROREF "<<nproref<<" NDECTRA "<<ndectra<<" NDECREF "<<ndecref<<endl; } else cout<<worldrank<<" FAILED ! ! ! ! !"<<endl;
for (int iout=0;iout<nprotra;iout++) cout<<worldrank<<" VALUE "<<lambdavec[protravec[iout]]<<endl;
#ifdef _ALLEIGVALS_WRITEOUT
ofstream myfile;
stringstream mysstream;
mysstream << "AllEigvals" << worldrank;
myfile.open(mysstream.str().c_str());
myfile.precision(15);
for (int iele=0;iele<neigval;iele++)
    myfile << real(lambdavec[iele]) << " " << imag(lambdavec[iele]) << endl;
myfile.close();
#endif
//if (!nprotra) {
//arpack=0;
//}
//MPI_Barrier(MPI_COMM_WORLD);
    if (nprotra!=nproref) return (LOGCERR, EXIT_FAILURE);
    if (ndectra!=ndecref) ndectra=min(ndectra,ndecref);//return (LOGCERR, EXIT_FAILURE);
// FOR ABOVE CASE I NEED TO IMPLEMENT SORTING OF LAMBDA? OR JUST THROW AWAY SOME VECS OF THE LIST SETTING NDEC TO MIN OF PRO AND REF?
// ONCE I REMOVE ONE OF THE DIRECTIONS I WILL ALSO REMOVE THE LINE ABOVE
    int neigbas=nprotra+ndectra;
    CPX *Vtra=new CPX[bandwidth*ndof*neigbas];
    CPX *Vref=new CPX[bandwidth*ndof*neigbas];
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
    CPX *VT=new CPX[neigbas*ntriblock];
    CPX *matcpx=new CPX[ntriblock*neigbas];
    CPX *invgls=new CPX[neigbas*neigbas];
    CPX *invgrs=new CPX[neigbas*neigbas];
    sabtime=get_time(d_zer);
    //left
//    c_zgemm('C','N',ntriblock,neigbas,ntriblock,z_one,H1cpx,ntriblock,Vtra,ntriblock,z_zer,matcpx,ntriblock);
// above its still a complex conjugate and no transpose for H1cpx, does that matter? it should make a difference...
    full_transpose(neigbas,ntriblock,Vtra,VT);
    H1t->trans_mat_vec_mult(VT,matcpx,neigbas,1);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vtra,ntriblock,matcpx,ntriblock,z_zer,invgls,neigbas);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[protravec[ieigbas]],+bandwidth),&invgls[ieigbas*neigbas],1);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[dectravec[ieigbas-nprotra]],+bandwidth),&invgls[ieigbas*neigbas],1);
//    c_zgemm('N','N',ntriblock,neigbas,ntriblock,z_one,H0cpx,ntriblock,Vtra,ntriblock,z_zer,matcpx,ntriblock);
    H0->trans_mat_vec_mult(VT,matcpx,neigbas,1);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vtra,ntriblock,matcpx,ntriblock,z_one,invgls,neigbas);
    //right
//    c_zgemm('N','N',ntriblock,neigbas,ntriblock,z_one,H1cpx,ntriblock,Vref,ntriblock,z_zer,matcpx,ntriblock);
    full_transpose(neigbas,ntriblock,Vref,VT);
    H1->trans_mat_vec_mult(VT,matcpx,neigbas,1);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_zer,invgrs,neigbas);
    for (int ieigbas=0;ieigbas<nprotra;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[prorefvec[ieigbas]],-bandwidth),&invgrs[ieigbas*neigbas],1);
    for (int ieigbas=nprotra;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdavec[decrefvec[ieigbas-nprotra]],-bandwidth),&invgrs[ieigbas*neigbas],1);
//    c_zgemm('N','N',ntriblock,neigbas,ntriblock,z_one,H0cpx,ntriblock,Vref,ntriblock,z_zer,matcpx,ntriblock);
    H0->trans_mat_vec_mult(VT,matcpx,neigbas,1);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_one,invgrs,neigbas);
    delete[] VT;
    cout << "TIME FOR SOME MATRIX MATRIX MULTIPLICATIONS " << get_time(sabtime) << endl;
    if (0) {
    CPX *KSeig=new CPX[neigbas*neigbas*nblocksband];
    CPX *invgtmp=new CPX[neigbas*neigbas];
    sabtime=get_time(d_zer);
// KSeig is KScpx in V-base
    for (int iband=bandwidth;iband<nblocksband;iband++) {
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
        for (int iband=bandwidth;iband<nblocksband;iband++) {
            c_zgemm('N','N',ndof,neigbas,ndof,z_one,&KScpx[ndofsq*iband],ndof,Vref,ntriblock,z_zer,matcpx,ntriblock);
            c_zgemm('C','N',neigbas,neigbas,ndof,z_one,Vref,ntriblock,matcpx,ntriblock,z_zer,&KSeig[neigbas*neigbas*iband],neigbas);
        }
    else {
// this does not really save time, at least not for the CNT where I tested it
        double *Vreal=new double[bandwidth*ndof*neigbas];
        double *Vimag=new double[bandwidth*ndof*neigbas];
        double *matdb=new double[bandwidth*ndof*neigbas];
        double *KStmp=new double[neigbas*neigbas];
        double* KSfull=new double[ndofsq*nblocksband];
        c_dcopy(ndofsq*nblocksband,(double*)KScpx,2,KSfull,1);
        for (int ii=0;ii<bandwidth*ndof*neigbas;ii++) {
            Vreal[ii]=real(Vref[ii]);
            Vimag[ii]=imag(Vref[ii]);
        }
        for (int iband=bandwidth;iband<nblocksband;iband++) {
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
        delete[] KSfull;
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
    delete[] KScpx;
// inversion
    int *pivarrayg=new int[neigbas];
    sigmal=new CPX[triblocksize];
    full_transpose(neigbas,ntriblock,Vtra,matcpx);
    c_dscal(neigbas*ntriblock,-d_one,((double*)matcpx)+1,2);
    double *worknorm=new double[neigbas];
    double anorm=c_zlange('1',neigbas,neigbas,invgls,neigbas,worknorm);
    delete[] worknorm;
    sabtime=get_time(d_zer);
    c_zgetrf(neigbas,neigbas,invgls,neigbas,pivarrayg,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    double rcond;
    CPX *cworkcond=new CPX[2*neigbas];
    double *dworkcond=new double[2*neigbas];
    c_zgecon('1',neigbas,invgls,neigbas,anorm,&rcond,cworkcond,dworkcond,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    delete[] cworkcond;
    delete[] dworkcond;
    cout << "Condition number is " << rcond << endl;
    if (rcond<numeric_limits<double>::epsilon()) cout << "Warning: Condition number below numerical precision" << endl;
    c_zgetrs('N',neigbas,ntriblock,invgls,neigbas,pivarrayg,matcpx,neigbas,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    cout << "TIME FOR gls INVERSION " << get_time(sabtime) << endl;
    c_zgemm('N','N',ntriblock,ntriblock,neigbas,z_one,Vtra,ntriblock,matcpx,neigbas,z_zer,sigmal,ntriblock);
    sigmar=new CPX[triblocksize];
    full_transpose(neigbas,ntriblock,Vref,matcpx);
    c_dscal(neigbas*ntriblock,-d_one,((double*)matcpx)+1,2);
    sabtime=get_time(d_zer);
    c_zgetrf(neigbas,neigbas,invgrs,neigbas,pivarrayg,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    c_zgetrs('N',neigbas,ntriblock,invgrs,neigbas,pivarrayg,matcpx,neigbas,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    cout << "TIME FOR grs INVERSION " << get_time(sabtime) << endl;
    c_zgemm('N','N',ntriblock,ntriblock,neigbas,z_one,Vref,ntriblock,matcpx,neigbas,z_zer,sigmar,ntriblock);
    delete[] pivarrayg;
    delete[] invgls;
    delete[] invgrs;
    delete[] matcpx;
    CPX *presigmal= new CPX[triblocksize];
    CPX *presigmar= new CPX[triblocksize];
    CPX *matctri = new CPX[triblocksize];
    sabtime=get_time(d_zer);
// i wouldnt need all the transpose if g00R, which is contained in sigma, was symmetric
    full_transpose(ntriblock,ntriblock,sigmal,matctri);
    H1t->trans_mat_vec_mult(matctri,presigmal,ntriblock,1);
    H1t->trans_mat_vec_mult(presigmal,matctri,ntriblock,1);
    full_transpose(ntriblock,ntriblock,matctri,sigmal);
    full_transpose(ntriblock,ntriblock,sigmar,matctri);
    H1->trans_mat_vec_mult(matctri,presigmar,ntriblock,1);
    H1->trans_mat_vec_mult(presigmar,matctri,ntriblock,1);
    full_transpose(ntriblock,ntriblock,matctri,sigmar);
    cout << "MATRIX MATRIX MULTIPLICATIONS FOR SIGMA " << get_time(sabtime) << endl;

    sabtime=get_time(d_zer);
    int inversion_with_sparse_mult=1;
    if (inversion_with_sparse_mult) {
        H0->sparse_to_full(matctri,ntriblock,ntriblock);
        c_zaxpy(triblocksize,-z_one,sigmal,1,matctri,1);
        int *pivarrays=new int[ntriblock];
        c_zgetrf(ntriblock,ntriblock,matctri,ntriblock,pivarrays,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        CPX nworks;
        c_zgetri(ntriblock,matctri,ntriblock,pivarrays,&nworks,-1,&iinfo);
        int lworks=int(real(nworks));
        CPX *works=new CPX[lworks];
        c_zgetri(ntriblock,matctri,ntriblock,pivarrays,works,lworks,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] works;
        full_transpose(ntriblock,ntriblock,matctri,sigmal);
        H1t->trans_mat_vec_mult(sigmal,presigmal,ntriblock,1);
        H1t->trans_mat_vec_mult(presigmal,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigmal);
        H0->sparse_to_full(matctri,ntriblock,ntriblock);
        c_zaxpy(triblocksize,-z_one,sigmar,1,matctri,1);
        c_zgetrf(ntriblock,ntriblock,matctri,ntriblock,pivarrays,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_zgetri(ntriblock,matctri,ntriblock,pivarrays,&nworks,-1,&iinfo);
        lworks=int(real(nworks));
        works=new CPX[lworks];
        c_zgetri(ntriblock,matctri,ntriblock,pivarrays,works,lworks,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] works;
        delete[] pivarrays;
        full_transpose(ntriblock,ntriblock,matctri,sigmar);
        H1->trans_mat_vec_mult(sigmar,presigmar,ntriblock,1);
        H1->trans_mat_vec_mult(presigmar,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigmar);
    } else {
        H0->sparse_to_full(matctri,ntriblock,ntriblock);
        c_zaxpy(triblocksize,-z_one,sigmal,1,matctri,1);
        int *pivarrays=new int[ntriblock];
        c_zgetrf(ntriblock,ntriblock,matctri,ntriblock,pivarrays,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        H1->sparse_to_full(sigmal,ntriblock,ntriblock);
        c_zgetrs('T',ntriblock,ntriblock,matctri,ntriblock,pivarrays,sigmal,ntriblock,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        full_transpose(ntriblock,ntriblock,sigmal,presigmal);
        H1t->trans_mat_vec_mult(presigmal,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigmal);
        H0->sparse_to_full(matctri,ntriblock,ntriblock);
        c_zaxpy(triblocksize,-z_one,sigmar,1,matctri,1);
        c_zgetrf(ntriblock,ntriblock,matctri,ntriblock,pivarrays,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        H1t->sparse_to_full(sigmar,ntriblock,ntriblock);
        c_zgetrs('T',ntriblock,ntriblock,matctri,ntriblock,pivarrays,sigmar,ntriblock,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] pivarrays;
        full_transpose(ntriblock,ntriblock,sigmar,presigmar);
        H1->trans_mat_vec_mult(presigmar,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigmar);
    }
    cout << "TIME FOR SYMMETRIZATION " << get_time(sabtime) << endl;

    delete[] matctri;

    if (compute_inj) {
        sabtime=get_time(d_zer);
        n_propagating=nprotra;
        c_dscal(nprotra*ntriblock,-d_one,((double*)Vtra)+1,2);
        c_dscal(nprotra*ntriblock,-d_one,((double*)Vref)+1,2);
        injl=new CPX[ntriblock*nprotra];
        injr=new CPX[ntriblock*nprotra];
        H1t->mat_vec_mult(Vtra,injl,nprotra);
        H1->mat_vec_mult(Vref,injr,nprotra);
        CPX *injl1=new CPX[ntriblock*nprotra];
        CPX *injr1=new CPX[ntriblock*nprotra];
        c_zcopy(ntriblock*nprotra,injl,1,injl1,1);
        c_zcopy(ntriblock*nprotra,injr,1,injr1,1);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,pow(lambdavec[protravec[ipro]],-bandwidth),&injl1[ipro*ntriblock],1);
            c_zscal(ntriblock,pow(lambdavec[prorefvec[ipro]],+bandwidth),&injr1[ipro*ntriblock],1);
        }
        H0->mat_vec_mult_add(Vtra,injl1,nprotra,z_one);
        H0->mat_vec_mult_add(Vref,injr1,nprotra,z_one);
// add I1 to I0 to form I
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,-z_one,presigmal,ntriblock,injl1,ntriblock,z_one,injl,ntriblock);
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,-z_one,presigmar,ntriblock,injr1,ntriblock,z_one,injr,ntriblock);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,CPX(d_one/sqrt(2*M_PI*veltra[ipro]),d_zer),&injl[ipro*ntriblock],1);
            c_zscal(ntriblock,CPX(d_one/sqrt(2*M_PI*velref[ipro]),d_zer),&injr[ipro*ntriblock],1);
        }
        cout << "MATRIX MATRIX MULTIPLICATIONS FOR INJECTION " << get_time(sabtime) << endl;
        delete[] injl1;
        delete[] injr1;
        do_delete_inj=1;
    }
    delete[] presigmal;
    delete[] presigmar;
    delete[] Vtra;
    delete[] Vref;
// lil arrays that do not take much memory
    delete[] lambdavec;
    delete[] dectravec;
    delete[] decrefvec;
    delete[] protravec;
    delete[] prorefvec;
    delete[] veltra;
    delete[] velref;
}//end if !boundary_rank
    do_delete_sigma=1;
    return 0;
}
