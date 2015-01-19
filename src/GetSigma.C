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

int BoundarySelfEnergy::Cutout(TCSR<CPX> *SumHamC,contact_type pcontact,CPX penergy,transport_methods::transport_method method,MPI_Comm matrix_comm)
{
    energy=penergy;
    if (method==transport_methods::WF) compute_inj=1;
    if (imag(energy) && compute_inj) return (LOGCERR, EXIT_FAILURE);

    contact=pcontact;
    int start=contact.start;
    int inj_sign=contact.inj_sign;
    int ndof=contact.ndof;
    int bandwidth=contact.bandwidth;
    int ntriblock=bandwidth*ndof;

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

    TCSR<CPX> *H0cut = new TCSR<CPX>(SumHamC,start,ntriblock,start,ntriblock);
    TCSR<CPX> *H1cut = new TCSR<CPX>(SumHamC,start,ntriblock,start+inj_sign*ntriblock,ntriblock);
    H0 = new TCSR<CPX>(H0cut,master_rank,matrix_comm);
    H1 = new TCSR<CPX>(H1cut,master_rank,matrix_comm);
    delete H0cut;
    delete H1cut;
    if (iam==master_rank) {
        H0->shift_resize(start,ntriblock,start,ntriblock);
        H1->shift_resize(start,ntriblock,start+inj_sign*ntriblock,ntriblock);
        H1t = new TCSR<CPX>(H1->size,H1->n_nonzeros,H1->findx);
        H1t->sparse_transpose(H1);
        TCSR<CPX> **H = new TCSR<CPX>*[2*bandwidth+1];
        if (contact.inj_sign==+1) {
            for (int ibw=bandwidth;ibw<2*bandwidth;ibw++) {
                H[ibw] = new TCSR<CPX>(H0,0,ndof,(ibw-bandwidth)*ndof,ndof);
                H[ibw]->shift_resize(0,ndof,(ibw-bandwidth)*ndof,ndof);
            }
            H[2*bandwidth] = new TCSR<CPX>(H1,0,ndof,0,ndof);
            H[2*bandwidth]->shift_resize(0,ndof,0,ndof);
            for (int ibw=1;ibw<=bandwidth;ibw++) {
                H[bandwidth-ibw] = new TCSR<CPX>(H[bandwidth+ibw]);
                H[bandwidth-ibw]->sparse_transpose(H[bandwidth+ibw]);
            }
        } else if (contact.inj_sign==-1) {
            H[0] =  new TCSR<CPX>(H1,ntriblock-ndof,ndof,ntriblock-ndof,ndof);
            H[0]->shift_resize(ntriblock-ndof,ndof,ntriblock-ndof,ndof);
            for (int ibw=1;ibw<bandwidth+1;ibw++) {
                H[ibw] = new TCSR<CPX>(H0,ntriblock-ndof,ndof,(ibw-1)*ndof,ndof);
                H[ibw]->shift_resize(ntriblock-ndof,ndof,(ibw-1)*ndof,ndof);
            }
            for (int ibw=1;ibw<=bandwidth;ibw++) {
                H[bandwidth+ibw] = new TCSR<CPX>(H[bandwidth-ibw]);
                H[bandwidth+ibw]->sparse_transpose(H[bandwidth-ibw]);
            }
        }
        H[0]->shift_resize(0,ndof,0,2*ntriblock);
        TCSR<CPX> *spBu = new TCSR<CPX>(H[0]);
        for (int ibw=1;ibw<2*bandwidth+1;ibw++) {
            H[ibw]->shift_resize(0,ndof,-(ibw-1)*ndof,2*ntriblock);
        }
        TCSR<CPX> *spAu = new TCSR<CPX>(2*bandwidth,NULL,&H[1],2*ntriblock);
        for (int ibw=0;ibw<2*bandwidth+1;ibw++) {
            delete H[ibw];
        }
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
    int start=contact.start;
    int ndof=contact.ndof;
    int bandwidth=contact.bandwidth;
    int ntriblock=bandwidth*ndof;
    TCSR<CPX> *spsigma = NULL;
    if (iam==master_rank) {
        spsigma = new TCSR<CPX>(SumHamC->size_tot,ntriblock*ntriblock,SumHamC->findx);
        spsigma->full_to_sparse(sigma,ntriblock,ntriblock,start,start);
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
            spainj = new TCSR<CPX>(SumHamC->size_tot,ntriblock*n_propagating,SumHamC->findx);
            spainj->full_to_sparse(inj,ntriblock,n_propagating,start,0);
            delete[] inj;
            inj = NULL;
        }
        spainjdist = new TCSR<CPX>(SumHamC,spainj,master_rank,matrix_comm);
        if (iam==master_rank) delete spainj;
        }
}

int BoundarySelfEnergy::GetSigma(MPI_Comm boundary_comm,int evecpos,transport_parameters *parameter_sab)
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
    int inj_sign=contact.inj_sign;
    int ndof=contact.ndof;
    int bandwidth=contact.bandwidth;
    int ndofsq=ndof*ndof;
    int nblocksband=2*bandwidth+1;
    int ntriblock=bandwidth*ndof;
    int triblocksize=ntriblock*ntriblock;
    double colzerothr=parameter_sab->colzero_threshold;
    double eps_limit=parameter_sab->eps_limit;
    double eps_decay=parameter_sab->eps_decay;
    double eps_eigval_degen=1.0E-6;
    int boundary_rank;
    MPI_Comm_rank(boundary_comm,&boundary_rank);
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
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
    sabtime=get_time(d_zer);
// FEAST
    int do_feast=0;
    if (do_feast) {
        InjectionFeast<CPX> *k_inj = new InjectionFeast<CPX>(complexenergypoint,2*bandwidth,10,eps_decay);
        int neigfeast=50;
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
        TCSR<CPX> *H1v; //DO I NEED THIS AT ALL OR IS IT THE SAME AND I CAN JUST CHANGE DERIVATIVE SO THAT THERE THE INJ_SIGN IS INCLUDED
        if (inj_sign==1) {
            H1v=H1;
        } else if (inj_sign==-1) {
            H1v=H1t;
        }
        k_inj->calc_kphase(spA,spB,H1v,2*ntriblock,2*ntriblock,eigvalcpx,eigvec,velref,&Ntr,ind_Ntr,&Nref,ind_Nref,inj_sign,1,1,boundary_comm,&iinfo);
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
if (!worldrank) cout << "TIME FOR EIGENVALUE SOLVER " << get_time(sabtime) << endl;
 /*
if (!boundary_rank) cout<<worldrank<<" AT "<<real(energy)<<" NPROTRA "<<nprotra<<" NPROREF "<<nproref<<" NDECTRA "<<ndectra<<" NDECREF "<<ndecref<<" SIGN "<<inj_sign<<endl;
if (!boundary_rank) {
stringstream mysstream;
mysstream << "AllEigvals" << evecpos << "_" << inj_sign;
ofstream myfile(mysstream.str().c_str());
myfile.precision(15);
for (int iele=0;iele<nproref+ndecref;iele++) {
    CPX k_eigval=CPX(0.0,1.0)*log(lambdaref[iele]);
    myfile << real(k_eigval) << " " << imag(k_eigval) << endl;
}
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
if (!worldrank) cout << "TIME FOR MATRIX MATRIX MULTIPLICATIONS FOR INVERSE OF G " << get_time(sabtime) << endl;
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
    sigma = new CPX[triblocksize];
    CPX* RCOR = new CPX[ntriblock*neigbas];
    int inversion_small_multmult = 0;
    if (inversion_small_multmult) {
        full_transpose(neigbas,ntriblock,Vref,RCOR);
        c_zgetrs('N',neigbas,ntriblock,invgrs,neigbas,pivarrayg,RCOR,neigbas,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_zgemm('N','N',ntriblock,ntriblock,neigbas,z_one,Vref,ntriblock,RCOR,neigbas,z_zer,sigma,ntriblock);
// i wouldnt need all the transpose if g00R, which is contained in sigma, was symmetric
        CPX *matctri  = new CPX[triblocksize];
        full_transpose(ntriblock,ntriblock,sigma,matctri);
        H1t->trans_mat_vec_mult(matctri,sigma,ntriblock,1);
        H1t->trans_mat_vec_mult(sigma,matctri,ntriblock,1);
        full_transpose(ntriblock,ntriblock,matctri,sigma);
        delete[] matctri;
    } else {
        CPX* LCOR = new CPX[ntriblock*neigbas];
        H1t->mat_vec_mult(Vref,LCOR,neigbas);
        full_transpose(neigbas,ntriblock,LCOR,RCOR);
        c_zgetrs('N',neigbas,ntriblock,invgrs,neigbas,pivarrayg,RCOR,neigbas,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        c_zgemm('N','N',ntriblock,ntriblock,neigbas,z_one,LCOR,ntriblock,RCOR,neigbas,z_zer,sigma,ntriblock);
        delete[] LCOR;
    }
    delete[] pivarrayg;
    delete[] RCOR;
    delete[] invgrs;
if (!worldrank) cout << "TIME FOR INVERSION AND MATRIX MATRIX MULTIPLICATIONS FOR SIGMA " << get_time(sabtime) << endl;

// /*
    sabtime=get_time(d_zer);
    int inversion_with_sparse_mult_symm=0;
    int linear_system_with_dense_mult_symm=0;
    if (inversion_with_sparse_mult_symm) {
        CPX *matctri  = new CPX[triblocksize];
        CPX *presigma  = new CPX[triblocksize];
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
        int do_extra_symm=0;
        if (do_extra_symm) {
            c_zaxpy(ntriblock*ntriblock,z_one,matctri,1,sigma,1);
            c_zscal(ntriblock*ntriblock,z_one*0.5,sigma,1);
        }
        delete[] matctri;
        delete[] presigma;
    } else if (linear_system_with_dense_mult_symm) {
        CPX *matctri  = new CPX[triblocksize];
        CPX *presigma  = new CPX[triblocksize];
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
        int do_extra_symm=0;
        if (do_extra_symm) {
            c_zaxpy(ntriblock*ntriblock,z_one,matctri,1,sigma,1);
            c_zscal(ntriblock*ntriblock,z_one*0.5,sigma,1);
        }
        delete[] matctri;
        delete[] presigma;
    }
if (!worldrank) cout << "TIME FOR SYMMETRIZATION " << get_time(sabtime) << endl;
// */

    if (compute_inj) {
        sabtime=get_time(d_zer);
        n_propagating=nprotra;
        c_dscal(nprotra*ntriblock,-d_one,((double*)Vref)+1,2);
        for (int i=0;i<nprotra;i++) {
            int degeneracy=1;
            for (int j=0;j<i;j++) {
                if (abs(lambdaref[i]-lambdaref[j])<eps_eigval_degen && abs(velref[i]-velref[j])<eps_eigval_degen) {
                    degeneracy++;
                    CPX prod=c_zdotc(ntriblock,&Vref[i*ntriblock],1,&Vref[j*ntriblock],1);
                    c_zaxpy(ntriblock,-prod,&Vref[i*ntriblock],1,&Vref[j*ntriblock],1);
                    CPX norm=CPX(d_one/c_dznrm2(ntriblock,&Vref[j*ntriblock],1),d_zer);
                    c_zscal(ntriblock,norm,&Vref[j*ntriblock],1);
                }
            }
if (degeneracy==2) cout << "DEGENERATE ON " << evecpos << endl;
            if (degeneracy>2) return (LOGCERR, EXIT_FAILURE);
        }
        inj = new CPX[ntriblock*nprotra];
        lambdapro = new CPX[nprotra];
        for (int ipro=0;ipro<nprotra;ipro++) {
            lambdapro[ipro]=pow(lambdaref[ipro],+inj_sign);
        }
        H1t->mat_vec_mult(Vref,inj,nprotra);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,pow(lambdaref[ipro],inj_sign*bandwidth),&Vref[ipro*ntriblock],1);
        }
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,z_one,sigma,ntriblock,Vref,ntriblock,z_one,inj,ntriblock);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,CPX(d_one/sqrt(2*M_PI*velref[ipro]),d_zer),&inj[ipro*ntriblock],1);
        }
if (!worldrank) cout << "MATRIX MATRIX MULTIPLICATIONS FOR INJECTION " << get_time(sabtime) << endl;
    }
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
