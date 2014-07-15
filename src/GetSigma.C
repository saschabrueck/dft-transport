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
#include "InjectionIEV.H"
#include "GetSigma.H"

BoundarySelfEnergy::BoundarySelfEnergy()
{
    sigma = NULL;
    spsigmadist = NULL;
    inj = NULL;
    lambdapro = NULL;
    n_propagating=-1;

    master_rank=-1;

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
        delete[] sigma;
    }

    if (do_delete_inj) {
        delete[] inj;
        delete[] lambdapro;
    }

    if (do_delete_spsigdist) {
        delete spsigmadist;
    }

    if (do_delete_spainjdist) {
        delete spainjdist;
        delete[] lambdapro;
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

int BoundarySelfEnergy::Cutout(TCSR<CPX> *SumHamC,int pcontact,CPX energy,transport_methods::transport_method method,c_transport_type parameter_sab,MPI_Comm matrix_comm)
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

    contact=pcontact;
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
    TCSR<CPX> *spsigma = NULL;
    if (iam==master_rank) {
        if (contact==1) {
            spsigma = new TCSR<CPX>(SumHamC->size_tot,ntriblock*ntriblock,SumHamC->findx);
            spsigma->full_to_sparse(sigma,ntriblock,ntriblock,0,0);
        } else if (contact==2) {
            spsigma = new TCSR<CPX>(SumHamC->size_tot,ntriblock*ntriblock,SumHamC->findx);
            spsigma->full_to_sparse(sigma,ntriblock,ntriblock,SumHamC->size_tot-ntriblock,SumHamC->size_tot-ntriblock);
        }
        delete[] sigma;
        sigma = NULL;
    }
    spsigmadist = new TCSR<CPX>(SumHamC,spsigma,master_rank,matrix_comm);
    if (iam==master_rank) delete spsigma;
    do_delete_sigma=0;
    do_delete_spsigdist=1;

    if (compute_inj) {
        MPI_Bcast(&n_propagating,1,MPI_INT,master_rank,matrix_comm);
        if (iam!=master_rank) {
            lambdapro=new CPX[n_propagating];
        }
        MPI_Bcast(&lambdapro,n_propagating,MPI_DOUBLE_COMPLEX,master_rank,matrix_comm);
        TCSR<CPX> *spainj = NULL;
        if (iam==master_rank) {
            if (contact==1) {
                spainj = new TCSR<CPX>(SumHamC->size_tot,ntriblock*n_propagating,SumHamC->findx);
                spainj->full_to_sparse(inj,ntriblock,n_propagating,0,0);
            } else if (contact==2) {
                spainj = new TCSR<CPX>(SumHamC->size_tot,ntriblock*n_propagating,SumHamC->findx);
                spainj->full_to_sparse(inj,ntriblock,n_propagating,SumHamC->size_tot-ntriblock,0);
            }
            delete[] inj;
            inj = NULL;
        }
        spainjdist = new TCSR<CPX>(SumHamC,spainj,master_rank,matrix_comm);
        if (iam==master_rank) delete spainj;
        do_delete_inj=0;
        do_delete_spainjdist=1;
        }
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
    int inj_sign;
    if (contact==1) {
        inj_sign=+1;
    } else if (contact==2) {
        inj_sign=-1;
    }
// set parameters
    int ndofsq=ndof*ndof;
    int nblocksband=2*bandwidth+1;
    int ntriblock=bandwidth*ndof;
    int triblocksize=ntriblock*ntriblock;
    int boundary_rank;
    MPI_Comm_rank(boundary_comm,&boundary_rank);
// global for this class
    CPX *Vtra=new CPX[2*ntriblock*2*ntriblock];
    CPX *Vref=new CPX[2*ntriblock*2*ntriblock];
    CPX *lambdatra=new CPX[2*ntriblock];
    CPX *lambdaref=new CPX[2*ntriblock];
    double *veltra=new double[2*ntriblock];
    double *velref=new double[2*ntriblock];
    int ndectra=0;
    int ndecref=0;
    int nprotra=0;
    int nproref=0;
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
int feast=0;
// FEAST
if (feast) {
    InjectionFeast<CPX> *k_inj = new InjectionFeast<CPX>(4,10);
    int neigfeast=45;
    k_inj->initialize(2*ntriblock,2*ntriblock,neigfeast);
    CPX* eigvalcpx=new CPX[neigfeast];
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
    int Ntr,Nref;
    int *ind_Ntr  = new int[2*ntriblock];
    int *ind_Nref = new int[2*ntriblock];
    double *vel_f = new double[2*ntriblock];
    k_inj->calc_kphase(spA,spB,H1,2*ntriblock,2*ntriblock,eigvalcpx,eigvec,vel_f,&Ntr,ind_Ntr,&Nref,ind_Nref,inj_sign,1,1,boundary_comm,&iinfo);
    cout << "TIME FOR FEAST " << get_time(sabtime) << endl;
    delete spA;
    delete spB;
    delete k_inj;
/*
    for(IV=0;IV<N;IV++){
        c_zcopy(inc,&Vin[ind[IV]*inc],1,&Vout[IV*inc],1);
    }
    for(IV=0;IV<N;IV++){
        c_zcopy(inc,&Vin[ind[IV]*inc],1,&Vout[IV*inc],1);
    }
*/
//oh no we dont know the number of propagating or do we? we have to fill out npro/tra/ref/dec and V and lambda and vel? and stuff
////the we can later also compare if it is worse if we use tra like in omen
    int neigval=Ntr+Nref;
    CPX *eigvecc=new CPX[ndof*neigval];//IST NEIGVAL AUCH RICHTIG? SIND DIE AUCH IN EIGVEC RICHTIG ANGEORDNET //ID O NOT NEED THSI HERE BECAUSE I PUT ON VEL DIRECTLY
    c_zlacpy('A',ndof,neigval,eigvec,2*ntriblock,eigvecc,ndof);
    delete[] ind_Ntr;
    delete[] ind_Nref;
    delete[] vel_f;
    delete[] eigvec;
} else {
    InjectionIEV inj_iev;
    inj_iev.Compute(KScpx,Vtra,Vref,lambdatra,lambdaref,ndectra,ndecref,nprotra,nproref,veltra,velref,ndof,bandwidth,inj_sign,complexenergypoint,colzerothr,eps_limit,eps_decay,energyp,boundary_comm);
}
    do_delete_sigma=1;
TCSR<CPX> *Hch;
if (contact==1){
Hch=H1;
H1=H1t;
H1t=Hch;
}
if (!boundary_rank) {
    delete[] KScpx;
// inverse g snake WARNING THE pow OF lambda IS DONE HERE
    if (nprotra!=nproref) return (LOGCERR, EXIT_FAILURE);
    if (ndectra!=ndecref) ndectra=min(ndectra,ndecref);
    int neigbas=nprotra+ndectra;
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
    for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdatra[ieigbas],+bandwidth),&invgls[ieigbas*neigbas],1);
//    c_zgemm('N','N',ntriblock,neigbas,ntriblock,z_one,H0cpx,ntriblock,Vtra,ntriblock,z_zer,matcpx,ntriblock);
    H0->trans_mat_vec_mult(VT,matcpx,neigbas,1);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vtra,ntriblock,matcpx,ntriblock,z_one,invgls,neigbas);
    //right
//    c_zgemm('N','N',ntriblock,neigbas,ntriblock,z_one,H1cpx,ntriblock,Vref,ntriblock,z_zer,matcpx,ntriblock);
    full_transpose(neigbas,ntriblock,Vref,VT);
    H1->trans_mat_vec_mult(VT,matcpx,neigbas,1);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_zer,invgrs,neigbas);
    for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
        c_zscal(neigbas,pow(lambdaref[ieigbas],-bandwidth),&invgrs[ieigbas*neigbas],1);
//    c_zgemm('N','N',ntriblock,neigbas,ntriblock,z_one,H0cpx,ntriblock,Vref,ntriblock,z_zer,matcpx,ntriblock);
    H0->trans_mat_vec_mult(VT,matcpx,neigbas,1);
    c_zgemm('C','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_one,invgrs,neigbas);
    delete[] VT;
    cout << "TIME FOR SOME MATRIX MATRIX MULTIPLICATIONS " << get_time(sabtime) << endl;
// inversion
    int *pivarrayg=new int[neigbas];
    CPX *sigmal=new CPX[triblocksize];
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
    CPX *sigmar=new CPX[triblocksize];
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
        CPX *injl=new CPX[ntriblock*nprotra];
        CPX *injr=new CPX[ntriblock*nprotra];
        CPX *lambdaprol=new CPX[nprotra];
        CPX *lambdapror=new CPX[nprotra];
        for (int ipro=0;ipro<nprotra;ipro++) {
            lambdaprol[ipro]=pow(lambdatra[ipro],-inj_sign);
            lambdapror[ipro]=pow(lambdaref[ipro],-inj_sign);
        }
        H1t->mat_vec_mult(Vtra,injl,nprotra);
        H1->mat_vec_mult(Vref,injr,nprotra);
        CPX *injl1=new CPX[ntriblock*nprotra];
        CPX *injr1=new CPX[ntriblock*nprotra];
        c_zcopy(ntriblock*nprotra,injl,1,injl1,1);
        c_zcopy(ntriblock*nprotra,injr,1,injr1,1);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,pow(lambdatra[ipro],+inj_sign*bandwidth),&injl1[ipro*ntriblock],1);
            c_zscal(ntriblock,pow(lambdaref[ipro],-inj_sign*bandwidth),&injr1[ipro*ntriblock],1);
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
        lambdapro=lambdapror;
        delete[] lambdaprol;
        inj=injr;
        delete[] injl;
        do_delete_inj=1;
    }
    sigma=sigmar;
    delete[] sigmal;
    delete[] presigmal;
    delete[] presigmar;
    delete[] Vtra;
    delete[] Vref;
// lil arrays that do not take much memory
    delete[] veltra;
    delete[] velref;
}//end if !boundary_rank
    do_delete_sigma=1;
    return 0;
}
