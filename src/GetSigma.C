#include <mpi.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <valarray>
#include "ScaLapack.H"
#include "InjectionBeyn.H"
#include "InjectionIEV.H"
#include "GetSigma.H"

BoundarySelfEnergy::BoundarySelfEnergy()
{
    sigma = NULL;
    inj = NULL;
    gamma = NULL;
    spsigmadist = NULL;
    spainjdist = NULL;
    lambdapro = NULL;
    n_propagating=-1;
    rcond=1.0;
    n_dec=-1;
    eigval_degeneracy=-1;

    master_rank=-1;

    H0 = NULL;
    H1 = NULL;
    H1t = NULL;

    H = NULL;

    compute_inj=0;
    compute_gamma=0;
}

BoundarySelfEnergy::~BoundarySelfEnergy()
{
    if (sigma || inj || gamma || spsigmadist || spainjdist || lambdapro || H || H0 || H1 || H1t ) throw std::exception();
}

void BoundarySelfEnergy::create_M_matrix(CPX *M,int NR,TCSR<CPX> **A,int bw,CPX z)
{
    init_var(M,NR*NR);
 
    for(int IP=-bw;IP<bw+1;IP++){
        TCSR<CPX>* AA = A[IP+bw];
        int ffindx = AA->findx;
        CPX zz = pow(z,IP);
        for(int IR=0;IR<NR;IR++){
            for(int IJ=AA->edge_i[IR]-ffindx;IJ<AA->edge_i[IR+1]-ffindx;IJ++){
                M[IR+(AA->index_j[IJ]-ffindx)*NR] += zz*AA->nnz[IJ];
            }
        }
    }
}

void BoundarySelfEnergy::Deallocate_Sigma()
{
    delete spsigmadist;
    spsigmadist = NULL;
}

void BoundarySelfEnergy::Deallocate_Gamma()
{
    if (compute_gamma) {
        delete[] gamma;
        gamma = NULL;
    }
}

void BoundarySelfEnergy::Deallocate_Injection()
{
    if (compute_inj) {
        delete spainjdist;
        spainjdist = NULL;
        delete[] lambdapro;
        lambdapro = NULL;
    }
}

void BoundarySelfEnergy::Finalize()
{
    delete[] sigma;
    sigma = NULL;
    if (compute_inj) {
        delete[] inj;
        inj = NULL;
        delete[] lambdapro;
        lambdapro = NULL;
    }
    if (compute_gamma) {
        delete[] gamma;
        gamma = NULL;
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
/*
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
*/
    if (!boundary_rank) master_rank = matrix_rank;
    return 0;
}

int BoundarySelfEnergy::Cutout(TCSR<CPX> *SumHamC,contact_type pcontact,CPX penergy,transport_methods::transport_method_type transport_method,MPI_Comm matrix_comm)
{
    energy=penergy;
    if (transport_method==transport_methods::WF) compute_inj=1;
    if (transport_method==transport_methods::NEGF) compute_gamma=1;

    contact=pcontact;
    int inj_sign=contact.inj_sign;
    int ndof=contact.ndof;
    int bandwidth=contact.bandwidth;
    int ntriblock=bandwidth*ndof;
    int start=contact.start;
    if (inj_sign==-1) {
        start+=-(bandwidth-1)*ndof;
    }

    int iam,nprocs;
    MPI_Comm_size(matrix_comm,&nprocs);
    MPI_Comm_rank(matrix_comm,&iam);
/*
    int *gathered_master_ranks = new int[nprocs];
    MPI_Allgather(&master_rank,1,MPI_INT,gathered_master_ranks,1,MPI_INT,matrix_comm);
    master_rank=-1;
    for (int irank=0;irank<nprocs;irank++) {
        if (gathered_master_ranks[irank]>=0) {
            master_rank=gathered_master_ranks[irank];
        }
    }
    delete[] gathered_master_ranks;
*/
    int master_rank_local = master_rank;
    MPI_Allreduce(&master_rank_local,&master_rank,1,MPI_INT,MPI_MAX,matrix_comm);

    TCSR<CPX> *H0cut = new TCSR<CPX>(SumHamC,start,ntriblock,start,ntriblock);
    TCSR<CPX> *H1cut = new TCSR<CPX>(SumHamC,start,ntriblock,start+inj_sign*ntriblock,ntriblock);
    H0 = new TCSR<CPX>(H0cut,master_rank,matrix_comm);
    H1 = new TCSR<CPX>(H1cut,master_rank,matrix_comm);
    delete H0cut;
    delete H1cut;
    if (iam==master_rank) {
        H1t = new TCSR<CPX>(H1->size,H1->n_nonzeros,H1->findx);
        H1t->sparse_transpose(H1);
        H = new TCSR<CPX>*[2*bandwidth+1];
        if (contact.inj_sign==+1) {
            for (int ibw=bandwidth;ibw<2*bandwidth;ibw++) {
                H[ibw] = new TCSR<CPX>(H0,0,ndof,(ibw-bandwidth)*ndof,ndof);
            }
            H[2*bandwidth] = new TCSR<CPX>(H1,0,ndof,0,ndof);
            for (int ibw=1;ibw<=bandwidth;ibw++) {
                H[bandwidth-ibw] = new TCSR<CPX>(H[bandwidth+ibw]);
                H[bandwidth-ibw]->sparse_transpose(H[bandwidth+ibw]);
            }
        } else if (contact.inj_sign==-1) {
            H[0] =  new TCSR<CPX>(H1,ntriblock-ndof,ndof,ntriblock-ndof,ndof);
            for (int ibw=1;ibw<bandwidth+1;ibw++) {
                H[ibw] = new TCSR<CPX>(H0,ntriblock-ndof,ndof,(ibw-1)*ndof,ndof);
            }
            for (int ibw=1;ibw<=bandwidth;ibw++) {
                H[bandwidth+ibw] = new TCSR<CPX>(H[bandwidth-ibw]);
                H[bandwidth+ibw]->sparse_transpose(H[bandwidth-ibw]);
            }
        }
        int tau_from_ndof_blocks=0;
        if (tau_from_ndof_blocks) {
            CPX* Hfull = new CPX[3*ntriblock*ntriblock]();
            for (int ibw=0;ibw<bandwidth;ibw++) {
                for (int jbw=0;jbw<2*bandwidth+1;jbw++) {
                    for(int i=0;i<ndof;i++){
                        for(int e=H[jbw]->edge_i[i]-H[jbw]->findx;e<H[jbw]->edge_i[i+1]-H[jbw]->findx;e++){
                            int j=H[jbw]->index_j[e]-H[jbw]->findx;
                            Hfull[ndof*ibw+(ibw+jbw)*ndof*ntriblock+i+j*ntriblock]=H[jbw]->nnz[e];
                        }
                    }
                }
            }
            H1t->full_to_sparse(&Hfull[(1-inj_sign)*ntriblock*ntriblock],ntriblock,ntriblock);
            H0-> full_to_sparse(&Hfull[             ntriblock*ntriblock],ntriblock,ntriblock);
            H1-> full_to_sparse(&Hfull[(1+inj_sign)*ntriblock*ntriblock],ntriblock,ntriblock);
            delete[] Hfull;
        }
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
    int ndof=contact.ndof;
    int bandwidth=contact.bandwidth;
    int ntriblock=bandwidth*ndof;
    int start=contact.start;
    if (contact.inj_sign==-1) {
        start+=-(bandwidth-1)*ndof;
    }

    MPI_Bcast(&n_propagating,1,MPI_INT,master_rank,matrix_comm);
    MPI_Bcast(&rcond,1,MPI_DOUBLE,master_rank,matrix_comm);
    MPI_Bcast(&n_dec,1,MPI_INT,master_rank,matrix_comm);
    MPI_Bcast(&eigval_degeneracy,1,MPI_INT,master_rank,matrix_comm);

    if (compute_gamma) {
        gamma = new CPX[ntriblock*ntriblock];
        if (iam==master_rank) {
            full_transpose(ntriblock,ntriblock,sigma,gamma);
            c_dscal(ntriblock*ntriblock,-1.0,((double*)gamma)+1,2);
            c_zaxpy(ntriblock*ntriblock,CPX(-1.0,0.0),sigma,1,gamma,1);
            c_zscal(ntriblock*ntriblock,CPX(0.0,-1.0),gamma,1);
        }
        MPI_Bcast(gamma,ntriblock*ntriblock,MPI_DOUBLE_COMPLEX,master_rank,matrix_comm);
    }

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

int BoundarySelfEnergy::GetSigma(MPI_Comm boundary_comm,transport_parameters transport_params)
{
    if (imag(energy) && transport_params.n_points_inv) {
        if (GetSigmaInv(boundary_comm,transport_params)) return (LOGCERR, EXIT_FAILURE);
    } else {
        if (GetSigmaEig(boundary_comm,transport_params)) return (LOGCERR, EXIT_FAILURE);
    }
    return 0;
}

int BoundarySelfEnergy::GetSigmaInv(MPI_Comm boundary_comm,transport_parameters transport_params)
{
    int complexenergypoint=0;
    if (imag(energy)) complexenergypoint=1;
    int inj_sign=contact.inj_sign;
    int ndof=contact.ndof;
    int bandwidth=contact.bandwidth;
    int ndofsq=ndof*ndof;
    int nblocksband=2*bandwidth+1;
    int ntriblock=bandwidth*ndof;
    int triblocksize=ntriblock*ntriblock;
    int NK=transport_params.n_points_inv;
    double imag_shift=0.0;
    int boundary_rank,boundary_size;
    MPI_Comm_rank(boundary_comm,&boundary_rank);
    MPI_Comm_size(boundary_comm,&boundary_size);
    CPX **B=NULL;
    CPX *BB=NULL;
    CPX *M=NULL;
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
double sabtime=get_time(0.0);
    if (!boundary_rank) {
        B=new CPX*[4*bandwidth-1];
        for (int IB=0;IB<4*bandwidth-1;IB++) {
            B[IB]=new CPX[ndofsq]();
        }
        M=new CPX[ndofsq];
    }
    for (int IK=0;IK<NK;IK++) {
        CPX z = exp(CPX(0.0,2.0*(IK+0.5)*M_PI/double(NK)))+imag_shift*CPX(0.0,1.0);
        if (!boundary_rank) {
            create_M_matrix(M,ndof,H,bandwidth,z);
        }
        if (p_inv(M,ndof,boundary_comm)) return (LOGCERR, EXIT_FAILURE);
        if (!boundary_rank) {
            for (int IB=0;IB<4*bandwidth-1;IB++) {
                CPX expfac=pow(z,IB-(2*bandwidth-1))/double(NK);
                c_zaxpy(ndofsq,expfac,M,1,B[IB],1);
            }
        }
    }
    if (!boundary_rank) {
        delete[] M;
    }
if (!worldrank) cout << "TIME FOR SIGMA SOLVER INTEGRATION " << get_time(sabtime) << endl;
    if (!boundary_rank) {
        sigma = new CPX[triblocksize];
        BB = new CPX[triblocksize];
        for (int IB=0;IB<bandwidth;IB++) {
            for (int JB=0;JB<bandwidth;JB++) {
                c_zlacpy('A',ndof,ndof,B[IB-JB+2*bandwidth-1-inj_sign*bandwidth],ndof,&sigma[ndof*IB+ntriblock*ndof*JB],ntriblock);
            }
        }
        full_transpose(ntriblock,ntriblock,sigma,BB);
        for (int IB=0;IB<bandwidth;IB++) {
            for (int JB=0;JB<bandwidth;JB++) {
                c_zlacpy('A',ndof,ndof,B[IB-JB+2*bandwidth-1],ndof,&sigma[ndof*IB+ntriblock*ndof*JB],ntriblock);
            }
        }
    }
    if (p_lin(sigma,BB,BB,ntriblock,ntriblock,boundary_comm)) return (LOGCERR, EXIT_FAILURE);
    if (!boundary_rank) {
        H1t->trans_mat_vec_mult(BB,sigma,ntriblock,1);
        c_zscal(triblocksize,CPX(-1.0,0.0),sigma,1);
        delete[] BB;
        for (int IB=0;IB<4*bandwidth-1;IB++) {
            delete[] B[IB];
        }
        delete[] B;
        for (int ibw=0;ibw<nblocksband;ibw++) {
            delete H[ibw];
        }
        delete[] H;
        H = NULL;
        delete H0;
        H0 = NULL;
        delete H1;
        H1 = NULL;
        delete H1t;
        H1t = NULL;
if (!worldrank) cout << "TIME FOR SIGMA SOLVER " << get_time(sabtime) << endl;
// ONLY CHANGES THE RESULT SLIGHTLY BUT IS IMPORTANT FOR PEXSI
// /*
        if (complexenergypoint) {
            for (int i=0;i<ntriblock;i++) for (int j=0;j<i;j++) {
                sigma[i+ntriblock*j]+=sigma[j+ntriblock*i];
                sigma[i+ntriblock*j]*=0.5;
            }
            for (int i=0;i<ntriblock;i++) for (int j=0;j<i;j++) {
                sigma[j+ntriblock*i]=sigma[i+ntriblock*j];
            }
        }
// */
    }
    return 0;
}

int BoundarySelfEnergy::GetSigmaEig(MPI_Comm boundary_comm,transport_parameters transport_params)
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
    injection_methods::injection_method_type injection_method=transport_params.injection_method;
    double colzerothr=transport_params.colzero_threshold;
    double eps_limit=transport_params.eps_limit;
    if (complexenergypoint) eps_limit=transport_params.eps_limit_cc;
    int neigbeyn=max(1,int(transport_params.fac_neigbeyn*ndof));
    if (complexenergypoint) neigbeyn=max(1,int(transport_params.fac_neigbeyn_cc*ndof));
    double eps_decay=transport_params.eps_decay;
    double eps_eigval_degen=transport_params.eps_eigval_degen;
    double svd_fac=transport_params.svd_cutoff;
    double NQ_beyn=transport_params.n_points_beyn;
    double NCRC_beyn=transport_params.NCRC_beyn;
    int ntasks_beyn=transport_params.tasks_per_integration_point;
    int boundary_rank;
    MPI_Comm_rank(boundary_comm,&boundary_rank);
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
    CPX *Vtra=NULL;
    CPX *Vref=NULL;
    CPX *lambdatra=NULL;
    CPX *lambdaref=NULL;
    double *veltra=NULL;
    double *velref=NULL;
    CPX *lambdavec=NULL;
    CPX *eigvecc=NULL;
    int neigval=0;
    int ndectra=0;
    int ndecref=0;
    int nprotra=0;
    int nproref=0;
    if (!boundary_rank) {
        Vtra=new CPX[ntriblock*ntriblock];
        Vref=new CPX[ntriblock*ntriblock];
        lambdatra=new CPX[ntriblock];
        lambdaref=new CPX[ntriblock];
        veltra=new double[ntriblock];
        velref=new double[ntriblock];
        lambdavec=new CPX[2*bandwidth*ndof];
        eigvecc=new CPX[ndof*2*bandwidth*ndof];
    }
sabtime=get_time(d_zer);
    if (injection_method==injection_methods::BEYN) {
        if (complexenergypoint) {
            InjectionBeyn<CPX> *k_inj = new InjectionBeyn<CPX>(2*bandwidth,1.0/eps_limit,svd_fac,NQ_beyn,NCRC_beyn,ntasks_beyn);
            neigval = k_inj->execute(H,ndof,neigbeyn,lambdavec,eigvecc,inj_sign,boundary_comm);
            delete k_inj;
        } else {
            InjectionBeyn<double> *k_inj = new InjectionBeyn<double>(2*bandwidth,1.0/eps_limit,svd_fac,NQ_beyn,NCRC_beyn,ntasks_beyn);
            neigval = k_inj->execute(H,ndof,neigbeyn,lambdavec,eigvecc,inj_sign,boundary_comm);
            delete k_inj;
        }
    } else {
        CPX* KScpx=NULL;
        if (!boundary_rank) {
            KScpx = new CPX[ndofsq*nblocksband];
            for (int ibw=0;ibw<nblocksband;ibw++) {
                H[ibw]->sparse_to_full(&KScpx[ndofsq*ibw],ndof,ndof);
            }
        }
        InjectionIEV inj_iev;
        if (inj_iev.Compute(neigval,lambdavec,eigvecc,KScpx,ndof,bandwidth,inj_sign,complexenergypoint,colzerothr,boundary_comm)) return (LOGCERR, EXIT_FAILURE);
        if (!boundary_rank) {
            delete[] KScpx;
        }
    }
if (!worldrank) cout << "TIME FOR EIGENVALUE SOLVER " << get_time(sabtime) << endl;
sabtime=get_time(d_zer);
    if (!boundary_rank) {
// DETERMINE TYPE OF EIGENVALUE/VECTOR
        int *dectravec=new int[neigval];
        int *decrefvec=new int[neigval];
        int *protravec=new int[neigval];
        int *prorefvec=new int[neigval];
        CPX *matcdof=new CPX[ndofsq]();
        CPX *vecout=new CPX[ndof];
        for (int iindnzcoln=0;iindnzcoln<neigval;iindnzcoln++) {
            CPX lambda=lambdavec[iindnzcoln];
            if (abs(lambda)>eps_limit && abs(lambda)<1.0/eps_limit) {
                if ((abs(lambda)>d_one+eps_decay && inj_sign>0) || (abs(lambda)<d_one-eps_decay && inj_sign<0))
                    dectravec[ndectra++]=iindnzcoln;
                else if ((abs(lambda)<d_one-eps_decay && inj_sign>0) || (abs(lambda)>d_one+eps_decay && inj_sign<0))
                    decrefvec[ndecref++]=iindnzcoln;
                else {
                    c_zscal(ndofsq,z_zer,matcdof,1);
                    for (int ibandw=1;ibandw<=bandwidth;ibandw++)
                        H[bandwidth+ibandw]->add_sparse_to_full(matcdof,ndof,ndof,CPX(ibandw,d_zer)*pow(lambda,-(ibandw-1)));//ACHTUNG IST DAS AUCH FUER BEIDE INJ SIGNS SO RICHTIG MIT DEM VZ
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
            } // END IF k not infinite
        } // END FOR
        delete[] matcdof;
        delete[] vecout;
        CPX normy;
        for (int ipro=0;ipro<nprotra;ipro++) {
            lambdatra[ipro]=lambdavec[protravec[ipro]];
            c_zcopy(ndof,&eigvecc[ndof*protravec[ipro]],1,&Vtra[ndof*bandwidth*ipro+ndof*(bandwidth-1)],1);
            for (int ibandw=2;ibandw<=bandwidth;ibandw++) {
                c_zcopy(ndof,&Vtra[ndof*bandwidth*ipro+ndof*(bandwidth+1-ibandw)],1,&Vtra[ndof*bandwidth*ipro+ndof*(bandwidth-ibandw)],1);
                c_zscal(ndof,lambdavec[protravec[ipro]],&Vtra[ndof*bandwidth*ipro+ndof*(bandwidth-ibandw)],1);
            }
            normy=CPX(d_one/c_dznrm2(bandwidth*ndof,&Vtra[ipro*bandwidth*ndof],1),d_zer);
            c_zscal(bandwidth*ndof,normy,&Vtra[ipro*bandwidth*ndof],1);
            veltra[ipro]*=real(normy)*real(normy);
        }
        for (int idec=0;idec<ndectra;idec++) {
            lambdatra[idec+nprotra]=lambdavec[dectravec[idec]];
            c_zcopy(ndof,&eigvecc[ndof*dectravec[idec]],1,&Vtra[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth-1)],1);
            for (int ibandw=2;ibandw<=bandwidth;ibandw++) {
                c_zcopy(ndof,&Vtra[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth+1-ibandw)],1,&Vtra[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth-ibandw)],1);
                c_zscal(ndof,lambdavec[dectravec[idec]],&Vtra[ndof*bandwidth*(idec+nprotra)+ndof*(bandwidth-ibandw)],1);
            }
            normy=CPX(d_one/c_dznrm2(bandwidth*ndof,&Vtra[(idec+nprotra)*bandwidth*ndof],1),d_zer);
            c_zscal(bandwidth*ndof,normy,&Vtra[(idec+nprotra)*bandwidth*ndof],1);
        }
        for (int ipro=0;ipro<nproref;ipro++) {
            lambdaref[ipro]=lambdavec[prorefvec[ipro]];
            c_zcopy(ndof,&eigvecc[ndof*prorefvec[ipro]],1,&Vref[ndof*bandwidth*ipro+ndof*(bandwidth-1)],1);
            for (int ibandw=2;ibandw<=bandwidth;ibandw++) {
                c_zcopy(ndof,&Vref[ndof*bandwidth*ipro+ndof*(bandwidth+1-ibandw)],1,&Vref[ndof*bandwidth*ipro+ndof*(bandwidth-ibandw)],1);
                c_zscal(ndof,lambdavec[prorefvec[ipro]],&Vref[ndof*bandwidth*ipro+ndof*(bandwidth-ibandw)],1);
            }
            normy=CPX(d_one/c_dznrm2(bandwidth*ndof,&Vref[ipro*bandwidth*ndof],1),d_zer);
            c_zscal(bandwidth*ndof,normy,&Vref[ipro*bandwidth*ndof],1);
            velref[ipro]*=real(normy)*real(normy);
        }
        for (int idec=0;idec<ndecref;idec++) {
            lambdaref[idec+nproref]=lambdavec[decrefvec[idec]];
            c_zcopy(ndof,&eigvecc[ndof*decrefvec[idec]],1,&Vref[ndof*bandwidth*(idec+nproref)+ndof*(bandwidth-1)],1);
            for (int ibandw=2;ibandw<=bandwidth;ibandw++) {
                c_zcopy(ndof,&Vref[ndof*bandwidth*(idec+nproref)+ndof*(bandwidth+1-ibandw)],1,&Vref[ndof*bandwidth*(idec+nproref)+ndof*(bandwidth-ibandw)],1);
                c_zscal(ndof,lambdavec[decrefvec[idec]],&Vref[ndof*bandwidth*(idec+nproref)+ndof*(bandwidth-ibandw)],1);
            }
            normy=CPX(d_one/c_dznrm2(bandwidth*ndof,&Vref[(idec+nproref)*bandwidth*ndof],1),d_zer);
            c_zscal(bandwidth*ndof,normy,&Vref[(idec+nproref)*bandwidth*ndof],1);
        }
        delete[] dectravec;
        delete[] decrefvec;
        delete[] protravec;
        delete[] prorefvec;
        delete[] lambdavec;
        delete[] eigvecc;
if (!worldrank) cout << "TIME FOR EIGENVALUE VELOCITY AND NORM " << get_time(sabtime) << endl;
 /*
stringstream mysstream;
mysstream << "AllEigvals" << worldrank;
ofstream myfile(mysstream.str().c_str());
myfile.precision(15);
for (int iele=0;iele<nproref+ndecref;iele++) {
    CPX k_eigval=CPX(0.0,1.0)*log(lambdaref[iele]);
    myfile << real(k_eigval) << " " << imag(k_eigval) << endl;
}
myfile.close();
 */
        for (int ibw=0;ibw<nblocksband;ibw++) {
            delete H[ibw];
        }
        delete[] H;
        H = NULL;
        if (nprotra!=nproref) return (LOGCERR, EXIT_FAILURE);
        n_propagating=nprotra;
        n_dec=ndecref;
        int neigbas=nproref+ndecref;
        CPX *VT=new CPX[neigbas*ntriblock];
        CPX *matcpx=new CPX[ntriblock*neigbas];
        CPX *invgrs=new CPX[neigbas*neigbas];
        sabtime=get_time(d_zer);
        if (eps_limit<1.0E-4) {
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
        CPX *cworkcond=new CPX[2*neigbas];
        double *dworkcond=new double[2*neigbas];
        c_zgecon('1',neigbas,invgrs,neigbas,anorm,&rcond,cworkcond,dworkcond,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        delete[] cworkcond;
        delete[] dworkcond;
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
        int iter_max=1;
        if (inversion_with_sparse_mult_symm) {
            CPX *matctri  = new CPX[triblocksize];
            CPX *presigma  = new CPX[triblocksize];
            for (int iter=0;iter<iter_max;iter++){
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
            }
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
            for (int iter=0;iter<iter_max;iter++){
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

// ONLY CHANGES THE RESULT SLIGHTLY BUT IS IMPORTANT FOR PEXSI
// /*
        if (complexenergypoint) {
            for (int i=0;i<ntriblock;i++) for (int j=0;j<i;j++) {
                sigma[i+ntriblock*j]+=sigma[j+ntriblock*i];
                sigma[i+ntriblock*j]*=0.5;
            }
            for (int i=0;i<ntriblock;i++) for (int j=0;j<i;j++) {
                sigma[j+ntriblock*i]=sigma[i+ntriblock*j];
            }
        }
// */

        if (compute_inj) {
            sabtime=get_time(d_zer);
            c_dscal(nprotra*ntriblock,-d_one,((double*)Vref)+1,2);
//          swap(Vref,Vtra);
//          swap(lambdaref,lambdatra);
//          c_dscal(nprotra,-d_one,((double*)lambdaref)+1,2);
            for (int i=0;i<nprotra;i++) {
                int degeneracy=1;
                for (int j=0;j<i;j++) {
                    if (abs(lambdaref[i]-lambdaref[j])<eps_eigval_degen && abs(velref[i]-velref[j])<eps_eigval_degen) {
                        degeneracy++;
/*
                        CPX prod=c_zdotc(ntriblock,&Vref[i*ntriblock],1,&Vref[j*ntriblock],1);
                        c_zaxpy(ntriblock,-prod,&Vref[i*ntriblock],1,&Vref[j*ntriblock],1);
                        CPX norm=CPX(d_one/c_dznrm2(ntriblock,&Vref[j*ntriblock],1),d_zer);
                        c_zscal(ntriblock,norm,&Vref[j*ntriblock],1);
*/
                    }
                }
                if (degeneracy>2) eigval_degeneracy=i;
            }
            inj = new CPX[ntriblock*nprotra];
            lambdapro = new CPX[nprotra];
            for (int ipro=0;ipro<nprotra;ipro++) {
                lambdapro[ipro]=pow(lambdaref[ipro],+inj_sign);
            }
            H1t->mat_vec_mult(Vref,inj,nprotra);
            for (int ipro=0;ipro<nprotra;ipro++) {
                c_zscal(ntriblock,pow(lambdaref[ipro],+inj_sign*bandwidth),&Vref[ipro*ntriblock],1);
            }
            c_zgemm('N','N',ntriblock,nprotra,ntriblock,z_one,sigma,ntriblock,Vref,ntriblock,z_one,inj,ntriblock);
            for (int ipro=0;ipro<nprotra;ipro++) {
                c_zscal(ntriblock,CPX(d_one/sqrt(2*M_PI*velref[ipro]),d_zer),&inj[ipro*ntriblock],1);
            }
if (!worldrank) cout << "TIME FOR MATRIX MATRIX MULTIPLICATIONS FOR INJECTION " << get_time(sabtime) << endl;
        }
        delete[] Vtra;
        delete[] Vref;
        delete[] lambdatra;
        delete[] lambdaref;
        delete[] veltra;
        delete[] velref;
        delete H0;
        H0 = NULL;
        delete H1;
        H1 = NULL;
        delete H1t;
        H1t = NULL;
    }//end if !boundary_rank
    return 0;
}
