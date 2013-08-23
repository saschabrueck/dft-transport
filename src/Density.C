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
#include "tmprGF.H"
#include "GetSigma.H"
#include "Density.H"

extern "C" {
int ldlt_preprocess__(int *, int *, int *, int *, int *, int *, int *);
int ldlt_fact__(int *, int *, int *, CPX *);
int ldlt_free__(int *);
int ldlt_blkselinv__(int *, int*, int*, CPX *, int*);
}

int density(TCSR<double> *KohnSham,TCSR<double> *Overlap,TCSR<double> *Ps,CPX energy,CPX weight,double *muvec,int n_mu,transport_methods::transport_method method,c_transport_type parameter_sab)
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
    double evoltfactor=parameter_sab.evoltfactor;
// complex or real energy
    int complexenergypoint=0;
    if (imag(energy)) complexenergypoint=1;
    if (complexenergypoint && method==transport_methods::WF) return (LOGCERR, EXIT_FAILURE);
// H-E*S
    TCSR<CPX> *SumHamC;
    if (int addpotential=0) {
        cout << "ADDING POTENTIAL" << endl;
        double *Vvec = new double[Overlap->size_tot];
        int na=5;
        int nb=na+5;
        int nc=ncells;
        for (int ivvec=0;          ivvec<na*4*13*12; ivvec++) Vvec[ivvec] =  0.0;
        for (int ivvec=na*4*13*12; ivvec<nb*4*13*12; ivvec++) Vvec[ivvec] =  1.5;
        for (int ivvec=nb*4*13*12; ivvec<nc*4*13*12; ivvec++) Vvec[ivvec] =  0.5;
        TCSR<double> *Pot = new TCSR<double>(Overlap,Vvec);
        SumHamC = new TCSR<CPX>(CPX(evoltfactor,d_zer),KohnSham,-energy,Overlap,z_one,Pot);
        delete Pot;
    } else {
        SumHamC = new TCSR<CPX>(CPX(evoltfactor,d_zer),KohnSham,-energy,Overlap);
    }
// set parameters
    int ndof=SumHamC->size_tot/ncells;
    if (ndof*ncells!=SumHamC->size_tot) return (LOGCERR, EXIT_FAILURE);
    int ndofsq,ndofsqbandwidth;
    ndofsq=ndof*ndof;
    ndofsqbandwidth=ndofsq*(2*bandwidth+1);
    int ntriblock=bandwidth*ndof;
    int triblocksize=ntriblock*ntriblock;
    ndofsqbandwidth=ndofsq*(2*bandwidth+1);
// compute self energies
    BoundarySelfEnergy selfenergy;
    if ( selfenergy.GetSigma(SumHamC,energy,1,method,parameter_sab) ) return (LOGCERR, EXIT_FAILURE);
    BoundarySelfEnergy selfenergy2;
    if ( selfenergy2.GetSigma(SumHamC,energy,2,method,parameter_sab) ) return (LOGCERR, EXIT_FAILURE);
    CPX *sigmal=selfenergy.sigmal;
    CPX *sigmar=selfenergy2.sigmar;
    CPX *sigmar2=selfenergy.sigmar;
cout << "COMPARE SIGMAR " << sigmar[triblocksize-1] << " AND SIGMAR2 " << sigmar2[triblocksize-1] << endl;
// now add sigma to tridiagonalblocks and convert to sparse
    CPX *Hlcpx = new CPX[2*triblocksize]();
    c_zaxpy(triblocksize,-z_one,sigmal,1,Hlcpx,1);
    SumHamC->addcontact(Hlcpx,ndof,bandwidth,1);
    TCSR<CPX> *LeftCorner;
    LeftCorner = new TCSR<CPX>(ntriblock,2*triblocksize,SumHamC->findx);
    LeftCorner->full_to_sparse(Hlcpx,ntriblock,2*ntriblock);
    delete[] Hlcpx;
    CPX *Hrcpx = new CPX[2*triblocksize]();
    c_zaxpy(triblocksize,-z_one,sigmar,1,&Hrcpx[triblocksize],1);
    SumHamC->addcontact(Hrcpx,ndof,bandwidth,2);
    TCSR<CPX> *RightCorner;
    RightCorner = new TCSR<CPX>(ntriblock,2*triblocksize,SumHamC->findx);
    RightCorner->full_to_sparse(Hrcpx,ntriblock,2*ntriblock);
    delete[] Hrcpx;
// replace corners to complete sparse matrix
    TCSR<CPX> *HamSig;
    int nhaminterior=SumHamC->edge_i[SumHamC->size_tot-ntriblock]-SumHamC->edge_i[ntriblock];
    int nleftcorner=LeftCorner->edge_i[ntriblock]-LeftCorner->findx;
    int nrightcorner=RightCorner->edge_i[ntriblock]-RightCorner->findx;
    HamSig = new TCSR<CPX>(SumHamC->size_tot,nhaminterior+nleftcorner+nrightcorner,SumHamC->findx);
    HamSig->assembleshift(LeftCorner,SumHamC,RightCorner,ntriblock,SumHamC->size_tot-2*ntriblock);
    delete LeftCorner;
    delete RightCorner;
    CPX* KScpx=new CPX[ndofsqbandwidth];
    SumHamC->contactunitcell(KScpx,ndof,bandwidth,1);
    delete SumHamC;
    if (method==transport_methods::TEST) {
// this is only a lil debugging thing
        CPX* H0cpx=new CPX[triblocksize];
        for (int ibandwidth=0;ibandwidth<bandwidth;ibandwidth++)
            for (int jbandwidth=0;jbandwidth<bandwidth;jbandwidth++)
                for (int jdof=0;jdof<ndof;jdof++)
                    c_zcopy(ndof,&KScpx[(bandwidth-ibandwidth+jbandwidth)*ndofsq+jdof*ndof],1,&H0cpx[(bandwidth*(jdof+jbandwidth*ndof)+ibandwidth)*ndof],1);
        c_zaxpy(triblocksize,-z_one,sigmal,1,H0cpx,1);
        c_zaxpy(triblocksize,-z_one,sigmar,1,H0cpx,1);
        int *pivarrayd= new int[ntriblock];
        sabtime=get_time(d_zer);
        c_zgetrf(ntriblock,ntriblock,H0cpx,ntriblock,pivarrayd,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        CPX workytestc;
        c_zgetri(ntriblock,H0cpx,ntriblock,pivarrayd,&workytestc,-1,&iinfo);
        int lworkyd=int(real(workytestc));
        CPX *workyd=new CPX[lworkyd];
        c_zgetri(ntriblock,H0cpx,ntriblock,pivarrayd,workyd,lworkyd,&iinfo);
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        cout << "TIME FOR G INVERSION " << get_time(sabtime) << endl;
        delete[] pivarrayd;
        delete[] workyd;
//transm=-trace((sigmar-sigmar')*G*(sigmal-sigmal')*G'); i think its better not to use imag() because of symmetry error
        sabtime=get_time(d_zer);
        CPX *matctrit= new CPX[triblocksize];
        for (int jtri=0;jtri<ntriblock;jtri++)
            for (int itri=0;itri<ntriblock;itri++)
                matctrit[itri*ntriblock+jtri]=sigmal[itri*ntriblock+jtri]-conj(sigmal[jtri*ntriblock+itri]);
        c_zgemm('N','C',ntriblock,ntriblock,ntriblock,z_one,matctrit,ntriblock,H0cpx,ntriblock,z_zer,sigmal,ntriblock);
        for (int jtri=0;jtri<ntriblock;jtri++)
            for (int itri=0;itri<ntriblock;itri++)
                matctrit[itri*ntriblock+jtri]=sigmar[itri*ntriblock+jtri]-conj(sigmar[jtri*ntriblock+itri]);
        c_zgemm('N','N',ntriblock,ntriblock,ntriblock,z_one,matctrit,ntriblock,H0cpx,ntriblock,z_zer,sigmar,ntriblock);
        delete[] matctrit;
        delete[] H0cpx;
        CPX trace=0;
        for (int jtri=0;jtri<ntriblock;jtri++)
            for (int itri=0;itri<ntriblock;itri++)
                trace+=-sigmar[itri*ntriblock+jtri]*sigmal[jtri*ntriblock+itri];
        cout << "MATRIX MATRIX MULTIPLICATIONS AND TRACE FOR TRANSMISSION " << get_time(sabtime) << endl;
        cout << "Energy " << energy << " Transmission " << real(trace) << endl;
    } else if (method==transport_methods::NEGF) {
        sabtime=get_time(d_zer);
        int inversion_method=0;
        if (inversion_method==0) {
            if (ncells%bandwidth) return (LOGCERR, EXIT_FAILURE);
            std::vector<int> Bvec(ncells/bandwidth,0);
            for (int ii=0;ii<ncells/bandwidth;ii++) Bvec[ii]=ii*ntriblock;
            tmprGF::sparse_invert(HamSig,Bvec);
            Ps->add_imag(HamSig,-weight/M_PI);
            delete HamSig;
        } else if (inversion_method==1) {
#ifdef HAVE_PARDISO            
            Pardiso::sparse_invert(HamSig);
            Ps->add_imag(HamSig,-weight/M_PI);
            delete HamSig;
#endif
        } else if (inversion_method==2) {
#ifdef HAVE_SELINV            
            HamSig->change_findx(1);
            TCSR<CPX> * HamSigTri= new TCSR<CPX>(HamSig->size_tot,HamSig->n_nonzeros/2+HamSig->size_tot,HamSig->findx);
            HamSigTri->extract_lower_triangle(HamSig);
            delete HamSig;
            TCSC<CPX,int> *HamSigCSC=new TCSC<CPX,int>(HamSigTri,1);
            delete HamSigTri;
            int token=0;
            int *perm = NULL;
            int Lnnz;
            int order=-1;
            ldlt_preprocess__(&token, &HamSigCSC->size, HamSigCSC->edge_j, HamSigCSC->index_i, &Lnnz, &order, perm);   
            ldlt_fact__(&token, HamSigCSC->edge_j, HamSigCSC->index_i, HamSigCSC->nnz);
            delete HamSigCSC;
            TCSC<CPX,int> *InverseMat=new TCSC<CPX,int>(Overlap->size,Lnnz,1);
            int dumpL=0;
            ldlt_blkselinv__(&token, InverseMat->edge_j, InverseMat->index_i, InverseMat->nnz, &dumpL);
            ldlt_free__(&token);
// missing function to convert InverseMat from tridiag to complete and add to density mat Ps
            delete InverseMat;
#endif
        } else return (LOGCERR, EXIT_FAILURE);
        cout << "TIME FOR SPARSE INVERSION " << get_time(sabtime) << endl;
    } else if (method==transport_methods::WF) {
// sparse solver
        int nprol=selfenergy.n_propagating;
        int npror=selfenergy2.n_propagating;
        CPX *injl=selfenergy.injl;
        CPX *injr=selfenergy2.injr;
        CPX *injr2=selfenergy.injr;
        CPX *Sol=new CPX[HamSig->size_tot*(nprol+npror)];
        CPX *Inj=new CPX[HamSig->size_tot*(nprol+npror)]();
cout << "COMPARE INJR " << injr[ntriblock-1] << " AND INJR2 " << injr2[ntriblock-1] << endl;
        c_zlacpy('A',ntriblock,nprol,injl,ntriblock,Inj,HamSig->size_tot);
        c_zlacpy('A',ntriblock,npror,injr,ntriblock,&Inj[HamSig->size_tot*(nprol+1)-ntriblock],HamSig->size_tot);
        LinearSolver<CPX>* solver;
        solver = new Umfpack<CPX>(HamSig,MPI_COMM_WORLD);
        sabtime=get_time(d_zer);
        solver->prepare();
        solver->solve_equation(Sol,Inj,nprol+npror);
        cout << "TIME FOR WAVEFUNCTION SPARSE SOLVER WITH "<< ncells <<" UNIT CELLS " << get_time(sabtime) << endl;
        delete[] Inj;
        delete solver;
        delete HamSig;
        sabtime=get_time(d_zer);
        double fermil=0.0;
        if (real(energy)<=muvec[0]) fermil=1.0;
        double fermir=0.0;
        if (real(energy)<=muvec[1]) fermir=1.0;
        Ps->psipsidagger(Sol,nprol,weight*fermil);
        Ps->psipsidagger(&Sol[HamSig->size_tot*nprol],npror,weight*fermir);
        cout << "TIME FOR CONSTRUCTION OF S-PATTERN DENSITY MATRIX " << get_time(sabtime) << endl;
// transmission
        CPX *vecoutdof=new CPX[ndof];
        double transml=d_zer;
        for (int iband=1;iband<=bandwidth;iband++) {
            for (int ipro=0;ipro<nprol;ipro++) {
                c_zgemv('N',ndof,ndof,z_one,&KScpx[ndofsq*(bandwidth+iband)],ndof,&Sol[ndof*(iband+ncells*ipro)],1,z_zer,vecoutdof,1);
                transml+=iband*4*M_PI*imag(c_zdotc(ndof,&Sol[ndof*ncells*ipro],1,vecoutdof,1));
            }
        }
        delete[] Sol;
        delete[] vecoutdof;
        cout << "Energy " << energy << " Transmission " << transml << endl;
    } else return (LOGCERR, EXIT_FAILURE);
    delete[] KScpx;

    return 0;
}
