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
#include "Arpack.H"
#include "GetSigma.H"

class InjectionIEV {
public:
InjectionIEV();
~InjectionIEV();
int Compute(CPX*,CPX*,CPX*,CPX*,CPX*,int&,int&,int&,int&,double*,double*,int,int,int,int,double,double,double,MPI_Comm);
private:
int Pascal(int,int);
};

InjectionIEV::InjectionIEV()
{
}

InjectionIEV::~InjectionIEV()
{
}

int InjectionIEV::Pascal(int r,int n)
{
    if( n == 0 )
        return 1;
    if( r == 0 || r == n )
        return 1;
    return Pascal( r - 1, n - 1 ) + Pascal( r, n - 1 );
}

int InjectionIEV::Compute(CPX* KScpx,CPX* Vtra,CPX* Vref,CPX* lambdatra,CPX* lambdaref,int& ndectra,int& ndecref,int& nprotra,int& nproref,double* veltra,double* velref,int ndof,int bandwidth,int inj_sign,int complexenergypoint,double colzerothr,double eps_limit,double eps_decay,MPI_Comm boundary_comm)
{
    int boundary_rank,boundary_size;
    MPI_Comm_rank(boundary_comm,&boundary_rank);
    MPI_Comm_size(boundary_comm,&boundary_size);
if (!boundary_rank) {
    double d_one=1.0;
    double d_zer=0.0;
    CPX z_one=CPX(d_one,d_zer);
    CPX z_zer=CPX(d_zer,d_zer);
    double sabtime;
    int iinfo=0;
// set parameters
    int ndofsq=ndof*ndof;
    int nblocksband=2*bandwidth+1;
    int ntriblock=bandwidth*ndof;
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
//        cout << "TIME FOR CPX INVERSION " << get_time(sabtime) << endl;
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
//        cout << "TIME FOR REAL INVERSION " << get_time(sabtime) << endl;
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
    CPX *eigvalcpx=new CPX[nindnzcoln];
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
//        cout << "TIME FOR CPX DIAGONALIZATION " << get_time(sabtime) << endl;
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
//        cout << "TIME FOR REAL DIAGONALIZATION " << get_time(sabtime) << endl;
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
    CPX *eigvecc=new CPX[ndof*nindnzcoln];
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
//        cout << "TIME FOR RECONSTRUCTION " << get_time(sabtime) << endl;
    }
    delete[] indnzcolvecn;
    delete[] indrzcolvecn;
    delete[] eigveccfull;
    delete[] matmr;
    int neigval=nindnzcoln;
// DETERMINE TYPE OF EIGENVALUE/VECTOR
    CPX *lambdavec=new CPX[neigval];
    int *dectravec=new int[neigval];
    int *decrefvec=new int[neigval];
    int *protravec=new int[neigval];
    int *prorefvec=new int[neigval];
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
    delete[] lambdavec;
    delete[] dectravec;
    delete[] decrefvec;
    delete[] protravec;
    delete[] prorefvec;
    delete[] eigvecc;
}
    return 0;
}
