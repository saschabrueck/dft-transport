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
    inj = NULL;
    spsigmadist = NULL;
    spainjdist = NULL;
    lambdapro = NULL;
    n_propagating=-1;

    master_rank=-1;

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
    delete spsigmadist;

    if (compute_inj) {
        delete spainjdist;
        delete[] lambdapro;
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
    if (iam==master_rank) {
        delete spsigma;
        delete H0;
        delete H1;
        delete H1t;
    }

    if (compute_inj) {
        MPI_Bcast(&n_propagating,1,MPI_INT,master_rank,matrix_comm);
        if (iam!=master_rank) {
            lambdapro = new CPX[n_propagating];
        }
        MPI_Bcast(lambdapro,n_propagating,MPI_DOUBLE_COMPLEX,master_rank,matrix_comm);
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
    CPX *Vtra;
    CPX *Vref;
    CPX *lambdatra;
    CPX *lambdaref;
    double *veltra;
    double *velref;
    int ndectra=0;
    int ndecref=0;
    int nprotra=0;
    int nproref=0;
    if (boundary_rank==0) {
        Vtra=new CPX[ntriblock*ntriblock];
        Vref=new CPX[ntriblock*ntriblock];
        lambdatra=new CPX[ntriblock];
        lambdaref=new CPX[ntriblock];
        veltra=new double[ntriblock];
        velref=new double[ntriblock];
    }
// FEAST
//    if (!complexenergypoint) {
    if (0) {
        InjectionFeast<CPX> *k_inj = new InjectionFeast<CPX>(complexenergypoint,2*bandwidth,10);
        int neigfeast=45;
        k_inj->initialize(2*ntriblock,2*ntriblock,neigfeast);
        TCSR<CPX> *spA = NULL;
        TCSR<CPX> *spB = NULL;
        int Ntr;
        int Nref;
        int *ind_Ntr;
        int *ind_Nref;
        CPX* eigvalcpx;
        CPX* eigvec;
        if (!boundary_rank) {
            TCSR<CPX> *H1u = new TCSR<CPX>(H1t,0,ndof,0,ndof);
            c_zscal(H1u->n_nonzeros,-z_one,H1u->nnz,1);
            spB = new TCSR<CPX>(2*ntriblock,H1u->n_nonzeros+2*ntriblock-ndof,H1u->findx);
            c_zcopy(H1u->n_nonzeros,H1u->nnz,1,spB->nnz,1);
            c_icopy(H1u->n_nonzeros,H1u->index_j,1,spB->index_j,1);
            c_icopy(ndof,H1u->index_i,1,spB->index_i,1);
            for (int i_ele=0;i_ele<2*ntriblock-ndof;i_ele++) {
                spB->nnz[H1u->n_nonzeros+i_ele]=z_one;
                spB->index_j[H1u->n_nonzeros+i_ele]=i_ele+ndof+spB->findx;
                spB->index_i[ndof+i_ele]=1;
            }
            delete H1u;
            spB->get_row_edge();
            spB->get_diag_pos();
            TCSR<CPX> *H1tu = new TCSR<CPX>(H1t,0,ndof,ndof,ntriblock-ndof);
            for (int e=0;e<H1tu->n_nonzeros;e++) {
                H1tu->index_j[e]-=ndof;
            }
            H1tu->size_tot=2*ntriblock;
            TCSR<CPX> *H0u = new TCSR<CPX>(H0,0,ndof,0,ntriblock);
            for (int e=0;e<H0u->n_nonzeros;e++) {
                H0u->index_j[e]+=ntriblock-ndof;
            }
            H0u->size_tot=2*ntriblock;
            TCSR<CPX> *Hsu = new TCSR<CPX>(z_one,H1tu,z_one,H0u);
            delete H1tu;
            delete H0u;
            TCSR<CPX> *H1su = new TCSR<CPX>(H1,0,ndof,0,ndof);
            for (int e=0;e<H1su->n_nonzeros;e++) {
                H1su->index_j[e]+=2*ntriblock-ndof;
            }
            H1su->size_tot=2*ntriblock;
            TCSR<CPX> *Hu = new TCSR<CPX>(z_one,Hsu,z_one,H1su);
            delete Hsu;
            delete H1su;
            spA = new TCSR<CPX>(2*ntriblock,Hu->n_nonzeros+2*ntriblock-ndof,Hu->findx);
            c_zcopy(Hu->n_nonzeros,Hu->nnz,1,spA->nnz,1);
            c_icopy(Hu->n_nonzeros,Hu->index_j,1,spA->index_j,1);
            c_icopy(ndof,Hu->index_i,1,spA->index_i,1);
            for (int i_ele=0;i_ele<2*ntriblock-ndof;i_ele++) {
                spA->nnz[Hu->n_nonzeros+i_ele]=z_one;
                spA->index_j[Hu->n_nonzeros+i_ele]=i_ele+spA->findx;
                spA->index_i[ndof+i_ele]=1;
            }
            delete Hu;
            spA->get_row_edge();
            spA->get_diag_pos();
            ind_Ntr   = new int[ntriblock];
            ind_Nref  = new int[ntriblock];
            eigvec    = new CPX[2*ntriblock*ntriblock];
            eigvalcpx = new CPX[ntriblock];
        }
        sabtime=get_time(d_zer);
        k_inj->calc_kphase(spA,spB,H1,2*ntriblock,2*ntriblock,eigvalcpx,eigvec,velref,&Ntr,ind_Ntr,&Nref,ind_Nref,inj_sign,1,1,boundary_comm,&iinfo);
        cout << "TIME FOR FEAST " << get_time(sabtime) << endl;
        delete k_inj;
        if (!boundary_rank) {
            delete spA;
            delete spB;
            nprotra=Ntr;
            nproref=Ntr;
            ndectra=Nref-Ntr;
            ndecref=Nref-Ntr;
            for (int IV=0;IV<Nref;IV++) {
                c_zcopy(ntriblock,&eigvec[ind_Nref[IV]*2*ntriblock],1,&Vtra[IV*ntriblock],1);
            }
            for(int IV=0;IV<Nref;IV++){
                CPX lambdatmp;
                c_zcopy(1,&eigvalcpx[ind_Nref[IV]],1,&lambdatmp,1);
                lambdatra[IV]=exp(CPX(0.0,-1.0)*lambdatmp);
            }
            if (inj_sign==-1) {
                c_zcopy(Nref*ntriblock,Vtra,1,Vref,1);
                c_zcopy(Nref,lambdatra,1,lambdaref,1);
            } else {
                c_zcopy(Ntr*ntriblock,&Vtra[(Nref-Ntr)*ntriblock],1,Vref,1);
                c_zcopy((Nref-Ntr)*ntriblock,Vtra,1,&Vref[Ntr*ntriblock],1);
                c_zcopy(Ntr,&lambdatra[Nref-Ntr],1,lambdaref,1);
                c_zcopy(Nref-Ntr,lambdatra,1,&lambdaref[Ntr],1);
            }
            for (int IV=0;IV<Ntr;IV++) {
                c_zcopy(ntriblock,&eigvec[ind_Ntr[IV]*2*ntriblock],1,&Vtra[IV*ntriblock],1);
            }
            for(int IV=0;IV<Ntr;IV++){
                CPX lambdatmp;
                c_zcopy(1,&eigvalcpx[ind_Ntr[IV]],1,&lambdatmp,1);
                lambdatra[IV]=exp(CPX(0.0,-1.0)*lambdatmp);
            }
            c_dcopy(Ntr,velref,1,veltra,1);
            delete[] ind_Ntr;
            delete[] ind_Nref;
            delete[] eigvec;
            delete[] eigvalcpx;
        }
    } else {
        CPX* KScpx;
        if (!boundary_rank) {
            KScpx=new CPX[ndofsq*nblocksband];
            H1t->sparse_to_full(KScpx,ndof,bandwidth*ndof);
            H0->sparse_to_full(&KScpx[bandwidth*ndofsq],ndof,bandwidth*ndof);
            H1->sparse_to_full(&KScpx[2*bandwidth*ndofsq],ndof,ndof);
        }
        InjectionIEV inj_iev;
        inj_iev.Compute(KScpx,Vtra,Vref,lambdatra,lambdaref,ndectra,ndecref,nprotra,nproref,veltra,velref,ndof,bandwidth,inj_sign,complexenergypoint,colzerothr,eps_limit,eps_decay,energyp,boundary_comm);
        if (!boundary_rank) {
            delete[] KScpx;
        }
    }
TCSR<CPX> *Hch;
if (contact==1){
Hch=H1;
H1=H1t;
H1t=Hch;
}
if (!boundary_rank) {
    if (nprotra!=nproref) return (LOGCERR, EXIT_FAILURE);
    if (ndectra!=ndecref) ndectra=min(ndectra,ndecref);
    int neigbas=nprotra+ndectra;
    CPX *VT=new CPX[neigbas*ntriblock];
    CPX *matcpx=new CPX[ntriblock*neigbas];
    CPX *invgrs=new CPX[neigbas*neigbas];
    sabtime=get_time(d_zer);
    int dotheh0ph1 = 1;
    if (dotheh0ph1) {
        full_transpose(neigbas,ntriblock,Vref,VT);
        H1->trans_mat_vec_mult(VT,matcpx,neigbas,1);
        c_zgemm('T','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_zer,invgrs,neigbas);
        for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
            c_zscal(neigbas,pow(lambdaref[ieigbas],+inj_sign*bandwidth),&invgrs[ieigbas*neigbas],1);
        H0->trans_mat_vec_mult(VT,matcpx,neigbas,1);
        c_zgemm('T','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_one,invgrs,neigbas);
    } else {
        full_transpose(neigbas,ntriblock,Vref,VT);
        H1t->trans_mat_vec_mult(VT,matcpx,neigbas,1);
        c_zgemm('T','N',neigbas,neigbas,ntriblock,-z_one,Vref,ntriblock,matcpx,ntriblock,z_zer,invgrs,neigbas);
        for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
            c_zscal(neigbas,pow(lambdaref[ieigbas],-inj_sign*bandwidth),&invgrs[ieigbas*neigbas],1);
    }
    delete[] VT;
    delete[] matcpx;
    cout << "TIME FOR MATRIX MATRIX MULTIPLICATIONS FOR INVERSE OF G" << get_time(sabtime) << endl;
    sigma         = new CPX[triblocksize];
    CPX *presigma = new CPX[triblocksize];
    CPX *matctri  = new CPX[triblocksize];
    sabtime=get_time(d_zer);
    double *worknorm=new double[neigbas];
    double anorm=c_zlange('1',neigbas,neigbas,invgrs,neigbas,worknorm);
    delete[] worknorm;
    int *pivarrayg=new int[neigbas];
    c_zgetrf(neigbas,neigbas,invgrs,neigbas,pivarrayg,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    double rcond;
    CPX *cworkcond=new CPX[2*neigbas];
    double *dworkcond=new double[2*neigbas];
    c_zgecon('1',neigbas,invgrs,neigbas,anorm,&rcond,cworkcond,dworkcond,&iinfo);
    if (iinfo) return (LOGCERR, EXIT_FAILURE);
    delete[] cworkcond;
    delete[] dworkcond;
    if (rcond<numeric_limits<double>::epsilon()) return (LOGCERR, EXIT_FAILURE);
    CPX* RCOR = new CPX[ntriblock*neigbas];
    int inversion_small_multmult = 1;
    if (inversion_small_multmult) {
        full_transpose(neigbas,ntriblock,Vref,RCOR);
        c_zgetrs('N',neigbas,ntriblock,invgrs,neigbas,pivarrayg,RCOR,neigbas,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_zgemm('N','N',ntriblock,ntriblock,neigbas,z_one,Vref,ntriblock,RCOR,neigbas,z_zer,sigma,ntriblock);
// i wouldnt need all the transpose if g00R, which is contained in sigma, was symmetric
        full_transpose(ntriblock,ntriblock,sigma,matctri);
        H1->trans_mat_vec_mult(matctri,presigma,ntriblock,1);
        H1->trans_mat_vec_mult(presigma,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigma);
    } else {
        CPX* LCOR = new CPX[ntriblock*neigbas];
        H1->mat_vec_mult(Vref,LCOR,neigbas);
        full_transpose(neigbas,ntriblock,LCOR,RCOR);
        c_zgetrs('N',neigbas,ntriblock,invgrs,neigbas,pivarrayg,RCOR,neigbas,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_zgemm('N','N',ntriblock,ntriblock,neigbas,z_one,LCOR,ntriblock,RCOR,neigbas,z_zer,sigma,ntriblock);
        full_transpose(neigbas,ntriblock,Vref,RCOR);
        c_zgetrs('N',neigbas,ntriblock,invgrs,neigbas,pivarrayg,RCOR,neigbas,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_zgemm('N','N',ntriblock,ntriblock,neigbas,z_one,LCOR,ntriblock,RCOR,neigbas,z_zer,presigma,ntriblock);
        delete[] LCOR;
    }
    delete[] pivarrayg;
    delete[] RCOR;
    delete[] invgrs;
    cout << "TIME FOR INVERSION AND MATRIX MATRIX MULTIPLICATIONS FOR SIGMA " << get_time(sabtime) << endl;
/*
    sabtime=get_time(d_zer);
    int inversion_with_sparse_mult=1;
    if (inversion_with_sparse_mult) {
        H0->sparse_to_full(matctri,ntriblock,ntriblock);
        c_zaxpy(triblocksize,-z_one,sigma,1,matctri,1);
        c_zgetrf(ntriblock,ntriblock,matctri,ntriblock,pivarrays,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_zgetri(ntriblock,matctri,ntriblock,pivarrays,&nworks,-1,&iinfo);
        lworks=int(real(nworks));
        works=new CPX[lworks];
        c_zgetri(ntriblock,matctri,ntriblock,pivarrays,works,lworks,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] works;
        delete[] pivarrays;
        full_transpose(ntriblock,ntriblock,matctri,sigma);
        H1->trans_mat_vec_mult(sigma,presigma,ntriblock,1);
        H1->trans_mat_vec_mult(presigma,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigma);
    } else {
        H0->sparse_to_full(matctri,ntriblock,ntriblock);
        c_zaxpy(triblocksize,-z_one,sigma,1,matctri,1);
        c_zgetrf(ntriblock,ntriblock,matctri,ntriblock,pivarrays,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        H1t->sparse_to_full(sigma,ntriblock,ntriblock);
        c_zgetrs('T',ntriblock,ntriblock,matctri,ntriblock,pivarrays,sigma,ntriblock,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] pivarrays;
        full_transpose(ntriblock,ntriblock,sigma,presigma);
        H1->trans_mat_vec_mult(presigma,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigma);
    }
    cout << "TIME FOR SYMMETRIZATION " << get_time(sabtime) << endl;
*/
    delete[] matctri;

    if (compute_inj) {
        sabtime=get_time(d_zer);
        n_propagating=nprotra;
        c_dscal(nprotra*ntriblock,-d_one,((double*)Vtra)+1,2);
        c_dscal(nprotra*ntriblock,-d_one,((double*)Vref)+1,2);
        inj = new CPX[ntriblock*nprotra];
        lambdapro = new CPX[nprotra];
        for (int ipro=0;ipro<nprotra;ipro++) {
            lambdapro[ipro]=pow(lambdaref[ipro],-inj_sign);
//it seems to be exactly the same no matter what sign i use here
        }
        H1->mat_vec_mult(Vref,inj,nprotra);
        CPX *injr1=new CPX[ntriblock*nprotra];
        c_zcopy(ntriblock*nprotra,inj,1,injr1,1);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,pow(lambdaref[ipro],-inj_sign*bandwidth),&injr1[ipro*ntriblock],1);
        }
        H0->mat_vec_mult_add(Vref,injr1,nprotra,z_one);
// add I1 to I0 to form I
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,-z_one,presigma,ntriblock,injr1,ntriblock,z_one,inj,ntriblock);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,CPX(d_one/sqrt(2*M_PI*velref[ipro]),d_zer),&inj[ipro*ntriblock],1);
        }
        cout << "MATRIX MATRIX MULTIPLICATIONS FOR INJECTION " << get_time(sabtime) << endl;
        delete[] injr1;
    }
    delete[] presigma;
    delete[] Vtra;
    delete[] Vref;
    delete[] lambdatra;
    delete[] lambdaref;
    delete[] veltra;
    delete[] velref;
}//end if !boundary_rank
    return 0;
}
