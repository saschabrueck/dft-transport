#include <mpi.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
using namespace std;
#include "ScaLapack.H"
#include "CSC.H"
#include "LinearSolver.H"
//#include "SuperLU.H"
#include "Umfpack.H"
//#include "MUMPS.H"
#include "Pardiso.H"
#include "GetSigma.H"

BoundarySelfEnergy::BoundarySelfEnergy()
{
    sigmal = NULL;
    sigmar = NULL;
    injl = NULL;
    injr = NULL;
    n_propagating = -1;
}

BoundarySelfEnergy::~BoundarySelfEnergy()
{
    delete[] sigmal;
    delete[] sigmar;
    delete[] injl;
    delete[] injr;
}

int BoundarySelfEnergy::GetSigma(TCSR<CPX> *SumHamC,CPX energy,int contact,transport_methods::transport_method method,c_transport_type parameter_sab)
{
    double d_one=1.0;
    double d_zer=0.0;
    CPX z_one=CPX(d_one,d_zer);
    CPX z_zer=CPX(d_zer,d_zer);
    double sabtime;
    int iinfo=0;
// get parameters
    int ncells=parameter_sab.n_cells;
    int bandwidth=parameter_sab.bandwidth;
    double colzerothr=parameter_sab.colzero_threshold;
    double eps_limit=parameter_sab.eps_limit;
    double eps_decay=parameter_sab.eps_decay;
// complex or real energy
    int complexenergypoint=0;
    if (imag(energy)) complexenergypoint=1;
    if (complexenergypoint && method==transport_methods::WF) return (LOGCERR, EXIT_FAILURE);
// inj_sign
    int inj_sign=-1;
// set parameters
    int ndof=SumHamC->size_tot/ncells;
    if (ndof*ncells!=SumHamC->size_tot) return (LOGCERR, EXIT_FAILURE);
    int ndofsq,ndofsqbandwidth;
    ndofsq=ndof*ndof;
    ndofsqbandwidth=ndofsq*(2*bandwidth+1);
    int ntriblock=bandwidth*ndof;
    int triblocksize=ntriblock*ntriblock;
    ndofsqbandwidth=ndofsq*(2*bandwidth+1);
// copy block
    CPX* KScpx=new CPX[ndofsqbandwidth];
    SumHamC->contactunitcell(KScpx,ndof,bandwidth,contact);
// assemble tridiagonalblocks
    CPX* H0cpx=new CPX[triblocksize];
    CPX* H1cpx=new CPX[triblocksize]();
    for (int ibandwidth=0;ibandwidth<bandwidth;ibandwidth++)
        for (int jbandwidth=0;jbandwidth<bandwidth;jbandwidth++)
            for (int jdof=0;jdof<ndof;jdof++)
                c_zcopy(ndof,&KScpx[(bandwidth-ibandwidth+jbandwidth)*ndofsq+jdof*ndof],1,&H0cpx[(bandwidth*(jdof+jbandwidth*ndof)+ibandwidth)*ndof],1);
    for (int ibandwidth=0;ibandwidth<bandwidth;ibandwidth++)
        for (int jbandwidth=0;jbandwidth<=ibandwidth;jbandwidth++)
            for (int jdof=0;jdof<ndof;jdof++)
                c_zcopy(ndof,&KScpx[(2*bandwidth-ibandwidth+jbandwidth)*ndofsq+jdof*ndof],1,&H1cpx[(bandwidth*(jdof+jbandwidth*ndof)+ibandwidth)*ndof],1);
// HERE STARTS THE NEW METHOD 
    CPX *mats=new CPX[ndofsq];
    c_zcopy(ndofsq,KScpx,1,mats,1);
    for (int ibandwidth=1;ibandwidth<2*bandwidth+1;ibandwidth++)
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
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
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
    }
    delete[] matmn;
    CPX *eigvecc=new CPX[ndof*nindnzcoln];
    if (blockbwpoe-blockbwpos==ndof-1) {
        c_zlacpy('A',ndof,nindnzcoln,&eigveccfull[blockbwpos],nindnzcoln,eigvecc,ndof);
    } else {
// RECONSTRUCT
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
    CPX *lambdavec=new CPX[nindnzcoln];
    int *dectravec=new int[nindnzcoln];
    int *decrefvec=new int[nindnzcoln];
    int *protravec=new int[nindnzcoln];
    int *prorefvec=new int[nindnzcoln];
    double *veltra=new double[nindnzcoln];
    double *velref=new double[nindnzcoln];
    int ndectra=0;
    int ndecref=0;
    int nprotra=0;
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
                } else return (LOGCERR, EXIT_FAILURE);
           } // END IF decaying or propagating
           lambdavec[iindnzcoln]=lambda;
        } // END IF k not infinite
    } // END FOR
    delete[] eigvalreal;
    delete[] eigvalimag;
    delete[] matcdof;
    delete[] vecout;
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
// KSeig is KScpx in V-base
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
        double* KSfull=new double[ndofsqbandwidth];
        c_dcopy(ndofsqbandwidth,(double*)KScpx,2,KSfull,1);
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
// inversion
    int *pivarrayg=new int[neigbas];
    sigmal=new CPX[triblocksize];
    full_conjugate_transpose(neigbas,ntriblock,Vtra,matcpx);
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
    full_conjugate_transpose(neigbas,ntriblock,Vref,matcpx);
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
    sabtime=get_time(d_zer);
    CPX *H1cpxt = new CPX[triblocksize];
    full_transpose(ntriblock,ntriblock,H1cpx,H1cpxt);
    CPX *matctri = new CPX[triblocksize];
    TCSR<CPX> *taut = new TCSR<CPX>(ntriblock,triblocksize,SumHamC->findx);
    taut->full_to_sparse(H1cpxt,ntriblock,ntriblock);
// i wouldnt need all the transpose if g00R, which is contained in sigma, was symmetric
    full_transpose(ntriblock,ntriblock,sigmal,matctri);
    taut->trans_mat_vec_mult(matctri,presigmal,ntriblock,1);
    taut->trans_mat_vec_mult(presigmal,matctri,ntriblock,1);
    full_transpose(ntriblock,ntriblock,matctri,sigmal);
    delete taut;
    TCSR<CPX> *tau = new TCSR<CPX>(ntriblock,triblocksize,SumHamC->findx);
    tau->full_to_sparse(H1cpx,ntriblock,ntriblock);
    full_transpose(ntriblock,ntriblock,sigmar,matctri);
    tau->trans_mat_vec_mult(matctri,presigmar,ntriblock,1);
    tau->trans_mat_vec_mult(presigmar,matctri,ntriblock,1);
    full_transpose(ntriblock,ntriblock,matctri,sigmar);
    delete tau;
    delete[] H1cpxt;
    delete[] matctri;
    cout << "MATRIX MATRIX MULTIPLICATIONS FOR SIGMA " << get_time(sabtime) << endl;
    if (method==transport_methods::WF) {
        sabtime=get_time(d_zer);
        n_propagating=nprotra;
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
        CPX *injl1=new CPX[ntriblock*nprotra];
        CPX *injr1=new CPX[ntriblock*nprotra];
        c_zcopy(ntriblock*nprotra,injl,1,injl1,1);
        c_zcopy(ntriblock*nprotra,injr,1,injr1,1);
        for (int ipro=0;ipro<nprotra;ipro++) {
            c_zscal(ntriblock,pow(lambdavec[protravec[ipro]],-bandwidth),&injl1[ipro*ntriblock],1);
            c_zscal(ntriblock,pow(lambdavec[prorefvec[ipro]],+bandwidth),&injr1[ipro*ntriblock],1);
        }
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,z_one,H0cpx,ntriblock,Vtracp,ntriblock,z_one,injl1,ntriblock);
        c_zgemm('N','N',ntriblock,nprotra,ntriblock,z_one,H0cpx,ntriblock,Vrefcp,ntriblock,z_one,injr1,ntriblock);
        delete[] Vtracp;
        delete[] Vrefcp;
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
    }
    delete[] H0cpx;
    delete[] H1cpx;
    delete[] presigmal;
    delete[] presigmar;
// maybe delete those earlier
    delete[] KScpx;
    delete[] Vtra;
    delete[] Vref;
// lil arrays that do not take much memory
    delete[] indnzcolvecn;
    delete[] indrzcolvecn;
    delete[] lambdavec;
    delete[] dectravec;
    delete[] decrefvec;
    delete[] protravec;
    delete[] prorefvec;
    delete[] veltra;
    delete[] velref;

    return 0;
}
