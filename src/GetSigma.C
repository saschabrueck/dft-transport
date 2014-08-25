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

    spA = NULL;
    spB = NULL;

    compute_inj=0;
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

int BoundarySelfEnergy::Cutout(TCSR<CPX> *SumHamC,int pcontact,CPX penergy,transport_methods::transport_method method,c_transport_type parameter_sab,MPI_Comm matrix_comm)
{
    energy=penergy;
    if (method==transport_methods::WF) compute_inj=1;
    if (imag(energy) && compute_inj) return (LOGCERR, EXIT_FAILURE);

    contact=pcontact;
    if (contact==1) {
        inj_sign=+1;
    } else if (contact==2) {
        inj_sign=-1;
    }

    n_cells=parameter_sab.n_cells;
    ndof=SumHamC->size_tot/n_cells;
    bandwidth=parameter_sab.bandwidth;
    int ndofsq=ndof*ndof;
    int ntriblock=bandwidth*ndof;

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

    int i0,j0,i1,j1;
    if (contact==1) {
        i0=0;
        j0=i0;
        i1=i0;
        j1=j0+ntriblock;
    } else if (contact==2) {
        i0=SumHamC->size_tot-ntriblock;
        j0=i0;
        i1=i0;
        j1=j0-ntriblock;
    }
    TCSR<CPX> *H0cut = new TCSR<CPX>(SumHamC,i0,ntriblock,j0,ntriblock);
    TCSR<CPX> *H1cut = new TCSR<CPX>(SumHamC,i1,ntriblock,j1,ntriblock);
    H0 = new TCSR<CPX>(H0cut,master_rank,matrix_comm);
    H1 = new TCSR<CPX>(H1cut,master_rank,matrix_comm);
    delete H0cut;
    delete H1cut;
    if (iam==master_rank) {
        H0->shift_resize(i0,ntriblock,j0,ntriblock);
        H1->shift_resize(i1,ntriblock,j1,ntriblock);
        H1t = new TCSR<CPX>(H1->size,H1->n_nonzeros,H1->findx);
        H1t->sparse_transpose(H1);
        CPX *H = new CPX[ndofsq*(2*bandwidth+1)];
        if (contact==1) {
            H0->sparse_to_full(&H[bandwidth*ndofsq],ndof,ntriblock,0,0);
            H1->sparse_to_full(&H[2*bandwidth*ndofsq],ndof,ndof,0,0);
            for (int ibw=1;ibw<=bandwidth;ibw++)
                full_transpose(ndof,ndof,&H[(bandwidth+ibw)*ndofsq],&H[(bandwidth-ibw)*ndofsq]);
        } else if (contact==2) {
            H0->sparse_to_full(&H[ndofsq],ndof,ntriblock,ntriblock-ndof,0);
            H1->sparse_to_full(&H[0],ndof,ndof,ntriblock-ndof,ntriblock-ndof);
            for (int ibw=1;ibw<=bandwidth;ibw++)
                full_transpose(ndof,ndof,&H[(bandwidth-ibw)*ndofsq],&H[(bandwidth+ibw)*ndofsq]);
        }
        TCSR<CPX> *spAu = new TCSR<CPX>(2*ntriblock,2*bandwidth*ndofsq,SumHamC->findx);
        spAu->full_to_sparse(&H[ndofsq],ndof,2*bandwidth*ndof,0,0);
        TCSR<CPX> *spBu = new TCSR<CPX>(2*ntriblock,ndofsq,SumHamC->findx);
        spBu->full_to_sparse(H,ndof,ndof,0,0);
        delete[] H;
        spA = new TCSR<CPX>(2*ntriblock,spAu->n_nonzeros+2*ntriblock-ndof,spAu->findx);
        c_zcopy(spAu->n_nonzeros,spAu->nnz,1,spA->nnz,1);
        c_icopy(spAu->n_nonzeros,spAu->index_j,1,spA->index_j,1);
        c_icopy(ndof,spAu->index_i,1,spA->index_i,1);
        for (int i_ele=0;i_ele<2*ntriblock-ndof;i_ele++) {
            spA->nnz[spAu->n_nonzeros+i_ele]=CPX(1.0,0.0);
            spA->index_j[spAu->n_nonzeros+i_ele]=i_ele+spA->findx;
            spA->index_i[ndof+i_ele]=1;
        }
        delete spAu;
        spA->get_row_edge();
        spA->get_diag_pos();
        c_zscal(spBu->n_nonzeros,-CPX(1.0,0.0),spBu->nnz,1);
        spB = new TCSR<CPX>(2*ntriblock,spBu->n_nonzeros+2*ntriblock-ndof,spBu->findx);
        c_zcopy(spBu->n_nonzeros,spBu->nnz,1,spB->nnz,1);
        c_icopy(spBu->n_nonzeros,spBu->index_j,1,spB->index_j,1);
        c_icopy(ndof,spBu->index_i,1,spB->index_i,1);
        for (int i_ele=0;i_ele<2*ntriblock-ndof;i_ele++) {
            spB->nnz[spBu->n_nonzeros+i_ele]=CPX(1.0,0.0);
            spB->index_j[spBu->n_nonzeros+i_ele]=i_ele+ndof+spB->findx;
            spB->index_i[ndof+i_ele]=1;
        }
        delete spBu;
        spB->get_row_edge();
        spB->get_diag_pos();
    } else {
        delete H0;
        delete H1;
        H0 = NULL;
        H1 = NULL;
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
// set parameters
    int complexenergypoint=0;
    if (imag(energy)) complexenergypoint=1;
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
    int do_feast=0;
    if (do_feast) {
        InjectionFeast<CPX> *k_inj = new InjectionFeast<CPX>(complexenergypoint,2*bandwidth,10);
        int neigfeast=45;
        k_inj->initialize(2*ntriblock,2*ntriblock,neigfeast);
        int Ntr;
        int Nref;
        int *ind_Ntr;
        int *ind_Nref;
        CPX* eigvalcpx;
        CPX* eigvec;
        if (!boundary_rank) {
            ind_Ntr   = new int[ntriblock];
            ind_Nref  = new int[ntriblock];
            eigvec    = new CPX[2*ntriblock*ntriblock];
            eigvalcpx = new CPX[ntriblock];
        }
        sabtime=get_time(d_zer);
        TCSR<CPX> *H1v; //DO I NEED THIS AT ALL OR IS IT THE SAME AND I CAN JUST CHANGE DERIVATIVE SO THAT THERE THE INJ_SIGN IS INCLUDED
        if (inj_sign==1) {
            H1v=H1;
        } else if (inj_sign==-1) {
            H1v=H1t;
        }
        k_inj->calc_kphase(spA,spB,H1v,2*ntriblock,2*ntriblock,eigvalcpx,eigvec,velref,&Ntr,ind_Ntr,&Nref,ind_Nref,inj_sign,1,1,boundary_comm,&iinfo);
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
            KScpx = new CPX[ndofsq*nblocksband];
            spA->sparse_to_full(&KScpx[ndofsq],ndof,2*bandwidth*ndof);
            spB->sparse_to_full(KScpx,ndof,ndof);
            c_zscal(ndofsq,-z_one,KScpx,1);
            delete spA;
            delete spB;
        }
        InjectionIEV inj_iev;
        inj_iev.Compute(KScpx,Vtra,Vref,lambdatra,lambdaref,ndectra,ndecref,nprotra,nproref,veltra,velref,ndof,bandwidth,inj_sign,complexenergypoint,colzerothr,eps_limit,eps_decay,boundary_comm);
        if (!boundary_rank) {
            delete[] KScpx;
        }
    }
/*
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
if (!boundary_rank) cout<<worldrank<<" AT "<<real(energy)<<" NPROTRA "<<nprotra<<" NPROREF "<<nproref<<" NDECTRA "<<ndectra<<" NDECREF "<<ndecref<<" SIGN "<<inj_sign<<endl;
if (!boundary_rank) {
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
stringstream mysstream;
mysstream << "AllEigvals" << worldrank << "_" << inj_sign;
ofstream myfile(mysstream.str().c_str());
myfile.precision(15);
for (int iele=0;iele<nproref+ndecref;iele++)
    myfile << real(lambdaref[iele]) << " " << imag(lambdaref[iele]) << endl;
myfile.close();
}
*/
if (!boundary_rank) {
    if (nprotra!=nproref) return (LOGCERR, EXIT_FAILURE);
    if (ndectra!=ndecref) ndecref=min(ndectra,ndecref);//NEEDED?
    int neigbas=nproref+ndecref;
    CPX *VT=new CPX[neigbas*ntriblock];
    CPX *matcpx=new CPX[ntriblock*neigbas];
    CPX *invgrs=new CPX[neigbas*neigbas];
    sabtime=get_time(d_zer);
    int dotheh0ph1 = 0;
    if (dotheh0ph1) {
        full_transpose(neigbas,ntriblock,Vref,VT);
        H1t->trans_mat_vec_mult(VT,matcpx,neigbas,1);
        c_zgemm('T','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_zer,invgrs,neigbas);
        for (int ieigbas=0;ieigbas<neigbas;ieigbas++)
            c_zscal(neigbas,pow(lambdaref[ieigbas],+inj_sign*bandwidth),&invgrs[ieigbas*neigbas],1);
        H0->trans_mat_vec_mult(VT,matcpx,neigbas,1);
        c_zgemm('T','N',neigbas,neigbas,ntriblock,z_one,Vref,ntriblock,matcpx,ntriblock,z_one,invgrs,neigbas);
    } else {
        full_transpose(neigbas,ntriblock,Vref,VT);
        H1->trans_mat_vec_mult(VT,matcpx,neigbas,1);
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
        H1t->trans_mat_vec_mult(matctri,presigma,ntriblock,1);
        H1t->trans_mat_vec_mult(presigma,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigma);
    } else {
        CPX* LCOR = new CPX[ntriblock*neigbas];
        H1t->mat_vec_mult(Vref,LCOR,neigbas);
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

    sabtime=get_time(d_zer);
    int inversion_with_sparse_mult_symm=0;
    int linear_system_with_dense_mult_symm=0;
    if (inversion_with_sparse_mult_symm) {
        H0->sparse_to_full(matctri,ntriblock,ntriblock);
        c_zaxpy(triblocksize,-z_one,sigma,1,matctri,1);
        int *pivarrays=new int[ntriblock];
        c_zgetrf(ntriblock,ntriblock,matctri,ntriblock,pivarrays,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        CPX nworks;
        c_zgetri(ntriblock,matctri,ntriblock,pivarrays,&nworks,-1,&iinfo);
        int lworks=int(real(nworks));
        CPX* works=new CPX[lworks];
        c_zgetri(ntriblock,matctri,ntriblock,pivarrays,works,lworks,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] works;
        delete[] pivarrays;
        full_transpose(ntriblock,ntriblock,matctri,sigma);
        H1t->trans_mat_vec_mult(sigma,presigma,ntriblock,1);
        H1t->trans_mat_vec_mult(presigma,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigma);
    } else if (linear_system_with_dense_mult_symm) {
        H0->sparse_to_full(matctri,ntriblock,ntriblock);
        c_zaxpy(triblocksize,-z_one,sigma,1,matctri,1);
        int *pivarrays=new int[ntriblock];
        c_zgetrf(ntriblock,ntriblock,matctri,ntriblock,pivarrays,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        H1->sparse_to_full(sigma,ntriblock,ntriblock);
        c_zgetrs('T',ntriblock,ntriblock,matctri,ntriblock,pivarrays,sigma,ntriblock,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] pivarrays;
        full_transpose(ntriblock,ntriblock,sigma,presigma);
        H1t->trans_mat_vec_mult(presigma,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigma);
    }
    cout << "TIME FOR SYMMETRIZATION " << get_time(sabtime) << endl;

    delete[] matctri;

    if (compute_inj) {
        sabtime=get_time(d_zer);
        n_propagating=nprotra;
        c_dscal(nprotra*ntriblock,-d_one,((double*)Vref)+1,2);
        inj = new CPX[ntriblock*nprotra];
        lambdapro = new CPX[nprotra];
        for (int ipro=0;ipro<nprotra;ipro++) {
            lambdapro[ipro]=pow(lambdaref[ipro],-inj_sign);
        }
        int dotheh0ph1inj=0;
        if (dotheh0ph1inj) {
            H1t->mat_vec_mult(Vref,inj,nprotra);
//            H1t->mat_vec_mult(Vtra,inj,nprotra);
            CPX *injr1=new CPX[ntriblock*nprotra];
            c_zcopy(ntriblock*nprotra,inj,1,injr1,1);
            for (int ipro=0;ipro<nprotra;ipro++) {
                c_zscal(ntriblock,pow(lambdaref[ipro],-inj_sign*bandwidth),&injr1[ipro*ntriblock],1);
//                c_zscal(ntriblock,pow(lambdatra[ipro],inj_sign*bandwidth),&injr1[ipro*ntriblock],1);
            }
            H0->mat_vec_mult_add(Vref,injr1,nprotra,z_one);
//            H0->mat_vec_mult_add(Vtra,injr1,nprotra,z_one);
            c_zgemm('N','N',ntriblock,nprotra,ntriblock,-z_one,presigma,ntriblock,injr1,ntriblock,z_one,inj,ntriblock);
            delete[] injr1;
        } else {
            H1t->mat_vec_mult(Vref,inj,nprotra);
            for (int ipro=0;ipro<nprotra;ipro++) {
                c_zscal(ntriblock,pow(lambdaref[ipro],inj_sign*bandwidth),&Vref[ipro*ntriblock],1);
            }
            c_zgemm('N','N',ntriblock,nprotra,ntriblock,z_one,sigma,ntriblock,Vref,ntriblock,z_one,inj,ntriblock);
        }
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,CPX(d_one/sqrt(2*M_PI*velref[ipro]),d_zer),&inj[ipro*ntriblock],1);
//            c_zscal(ntriblock,CPX(d_one/sqrt(2*M_PI*veltra[ipro]),d_zer),&inj[ipro*ntriblock],1);
        }
        cout << "MATRIX MATRIX MULTIPLICATIONS FOR INJECTION " << get_time(sabtime) << endl;
    }
    delete[] presigma;
    delete[] Vtra;
    delete[] Vref;
    delete[] lambdatra;
    delete[] lambdaref;
    delete[] veltra;
    delete[] velref;
    delete H0;
    delete H1;
    delete H1t;
}//end if !boundary_rank
    return 0;
}
