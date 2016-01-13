#include <mpi.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>
#include <algorithm>
using namespace std;
#include "ScaLapack.H"
#include "CSC.H"
#include "LinearSolver.H"
#include "tmprGF.H"
#ifdef HAVE_MUMPS
#include "MUMPS.H"
#endif
#include "SuperLU.H"
#include "c_pexsi_interface.h"
#include "GetSigma.H"
#include "Density.H"

int density(TCSR<double> *KohnSham,TCSR<double> *Overlap,TCSR<double> *Ps,CPX energy,CPX weight,transport_methods::transport_method method,std::vector<double> muvec,std::vector<contact_type> contactvec,double &current,double &transm,double &dos,std::vector<int> propnum,transport_parameters *parameter_sab,int evecpos,MPI_Comm matrix_comm)
{
    double d_one=1.0;
    double d_zer=0.0;
    CPX z_one=CPX(d_one,d_zer);
    CPX z_img=CPX(d_zer,d_one);
    double sabtime;
    TCSR<CPX> *HamSig;
    CPX* inj = NULL;
    int nprol, npror;
    CPX *lambda;
    int *dist_sol;
    int *displc_sol;
    int matrix_procs,matrix_rank;
    MPI_Comm_size(matrix_comm,&matrix_procs);
    MPI_Comm_rank(matrix_comm,&matrix_rank);
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
// get parameters always at the beginning of a subroutine (or maybe even better in the constructor) to make it independent of the structure containing the parameters
    int n_mu=muvec.size();
    int solver_method=parameter_sab->linear_solver;
    int cutout=parameter_sab->cutout;
int bandwidth=contactvec[0].bandwidth;
int ndof=contactvec[0].ndof;
int ncells=Overlap->size_tot/ndof;
int ntriblock=bandwidth*ndof;
int tra_block=0;
    std::vector<TCSR<CPX>*> Gamma(n_mu);
    if (method!=transport_methods::EQ) {
sabtime=get_time(d_zer);
        TCSR<CPX> *SumHamC = new TCSR<CPX>(Overlap->size,Overlap->n_nonzeros,Overlap->findx);
        SumHamC->copy_contain(Overlap,d_one);
        c_zscal(SumHamC->n_nonzeros,-energy,SumHamC->nnz,1);
        c_daxpy(SumHamC->n_nonzeros,d_one,KohnSham->nnz,1,(double*)SumHamC->nnz,2);
if (!worldrank) cout << "TIME FOR SumHamC " << get_time(sabtime) << endl;
// set pbc to zero 
        if (!cutout) SumHamC->settozeropbc(bandwidth,ndof);
//SumHamC->remove_thr(1.0E-6);
// compute self energies
        std::vector<BoundarySelfEnergy> selfenergies(n_mu);
        if (matrix_procs>1 && matrix_procs%n_mu) {
            if (!matrix_rank) cout << "Choose number of tasks per energy point as one or a multiple of number of contacts" << endl;
            return (LOGCERR, EXIT_FAILURE);
        }
        int boundary_id=matrix_rank*n_mu/matrix_procs;
        int n_bound_comm=min(matrix_procs,n_mu);
        MPI_Comm boundary_comm;
        MPI_Comm_split(matrix_comm,boundary_id,matrix_rank,&boundary_comm);
        for (int iseq=0;iseq<n_mu/n_bound_comm;iseq++) {
            int ipos=boundary_id+iseq*n_bound_comm;
            if (ipos<n_mu) {
                if ( selfenergies[ipos].Set_master(matrix_comm,boundary_comm) ) return (LOGCERR, EXIT_FAILURE);
            }
sabtime=get_time(d_zer);
            for (int i_bound_id=0;i_bound_id<n_bound_comm;i_bound_id++) {
                int ibpos=i_bound_id+iseq*n_bound_comm;
                if ( selfenergies[ibpos].Cutout(SumHamC,contactvec[ibpos],energy,method,matrix_comm) ) return (LOGCERR, EXIT_FAILURE);
            }
if (!worldrank) cout << "TIME FOR SIGMA CUTOUT " << get_time(sabtime) << endl;
sabtime=get_time(d_zer);
            if (ipos<n_mu) {
                if ( selfenergies[ipos].GetSigma(boundary_comm,evecpos,parameter_sab) ) return (LOGCERR, EXIT_FAILURE);
            }
if (!worldrank) cout << "TIME FOR SIGMA GETSIGMA " << get_time(sabtime) << endl;
sabtime=get_time(d_zer);
            for (int i_bound_id=0;i_bound_id<n_bound_comm;i_bound_id++) {
                int ibpos=i_bound_id+iseq*n_bound_comm;
                selfenergies[ibpos].Distribute(SumHamC,matrix_comm);
            }
if (!worldrank) cout << "TIME FOR SIGMA DISTRIBUTE " << get_time(sabtime) << endl;
        }
        MPI_Comm_free(&boundary_comm);
//add sigma to sumhamc
sabtime=get_time(d_zer);
        TCSR<CPX> **HamSigVec = new TCSR<CPX>*[n_mu+1];
        CPX *HamSigSigns = new CPX[n_mu+1];
        HamSigVec[0]=SumHamC;
        HamSigSigns[0]=z_one;
        for (int i_mu=0;i_mu<n_mu;i_mu++) {
            HamSigVec[i_mu+1]=selfenergies[i_mu].spsigmadist;
            HamSigSigns[i_mu+1]=-z_one;
        }
        HamSig = new TCSR<CPX>(n_mu+1,HamSigSigns,HamSigVec);
        delete[] HamSigVec;
        delete[] HamSigSigns;
        delete SumHamC;
if (!worldrank) cout << "TIME FOR ADDING SIGMA " << get_time(sabtime) << endl;
        if (method==transport_methods::WF) {
            nprol=selfenergies[0].n_propagating;
            npror=selfenergies[1].n_propagating;
            lambda = new CPX[nprol+npror];
            c_zcopy(nprol,selfenergies[0].lambdapro,1,lambda,1);
            c_zcopy(npror,selfenergies[1].lambdapro,1,&lambda[nprol],1);
if (nprol!=propnum[0] && propnum[0]>=0) if (!matrix_rank) cout << "WARNING: FOUND " << nprol << " OF " << propnum[0] << " MODES AT E=" << real(energy) << " POSITION " << evecpos << " LEFT" << endl;
if (npror!=propnum[1] && propnum[1]>=0) if (!matrix_rank) cout << "WARNING: FOUND " << npror << " OF " << propnum[1] << " MODES AT E=" << real(energy) << " POSITION " << evecpos << " RIGHT" << endl;
            dist_sol = new int[matrix_procs];
            MPI_Allgather(&HamSig->size,1,MPI_INT,dist_sol,1,MPI_INT,matrix_comm);
            displc_sol = new int[matrix_procs+1]();
            for (int iii=1;iii<matrix_procs+1;iii++) {
                displc_sol[iii]=displc_sol[iii-1]+dist_sol[iii-1];
            }
            int injsize_loc_max=*max_element(dist_sol,dist_sol+matrix_procs);
            inj = new CPX[injsize_loc_max*(nprol+npror)]();
            if (selfenergies[0].spainjdist->n_nonzeros) {
                selfenergies[0].spainjdist->sparse_to_full(inj,HamSig->size,nprol);
            }
            if (selfenergies[1].spainjdist->n_nonzeros) {
                selfenergies[1].spainjdist->sparse_to_full(&inj[HamSig->size*nprol],HamSig->size,npror);
            }
        }
        if (method==transport_methods::NEGF) {
            //Gamma
        }
    }
    if (method==transport_methods::EQ) {
        if (KohnSham->findx!=1 || Overlap->findx!=1) return (LOGCERR, EXIT_FAILURE);
 
        double *HS_nnz_inp = new double[2*Overlap->n_nonzeros]();
        double *HS_nnz_out = new double[2*Overlap->n_nonzeros]();
 
        c_dcopy(Overlap->n_nonzeros,Overlap->nnz,1,HS_nnz_inp,2);
        c_zscal(Overlap->n_nonzeros,-energy,(CPX*)HS_nnz_inp,1);
        c_daxpy(Overlap->n_nonzeros,1.0,KohnSham->nnz,1,HS_nnz_inp,2);
 
        int n_nonzeros_global;
        MPI_Allreduce(&Overlap->n_nonzeros,&n_nonzeros_global,1,MPI_INT,MPI_SUM,matrix_comm);
        int info;
 
        PPEXSIOptions options;
        PPEXSISetDefaultOptions(&options);
        options.npSymbFact = 1;
        options.ordering = 1;
  
        PPEXSIPlan plan;
        int nprowcol[2]={0,0};
        MPI_Dims_create(matrix_procs,2,nprowcol);
        plan = PPEXSIPlanInitialize(matrix_comm,nprowcol[0],nprowcol[1],-1,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
        PPEXSILoadRealSymmetricHSMatrix(plan,options,Overlap->size_tot,n_nonzeros_global,Overlap->n_nonzeros,Overlap->size,Overlap->edge_i,Overlap->index_j,HS_nnz_inp,1,NULL,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
        PPEXSISymbolicFactorizeComplexSymmetricMatrix(plan,options,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
        PPEXSISelInvComplexSymmetricMatrix(plan,options,HS_nnz_inp,HS_nnz_out,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
        PPEXSIPlanFinalize(plan,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
 
        delete[] HS_nnz_inp;
        c_zscal(Overlap->n_nonzeros,-weight/M_PI*CPX(0.0,1.0),(CPX*)HS_nnz_out,1);
        dos=c_ddot(Overlap->n_nonzeros,Overlap->nnz,1,HS_nnz_out,2);
        c_daxpy(Overlap->n_nonzeros,1.0,HS_nnz_out,2,Ps->nnz,1);
        delete[] HS_nnz_out;
    } else if (method==transport_methods::GF) {
        if (HamSig->findx!=1) return (LOGCERR, EXIT_FAILURE);
 
        double *HS_nnz_inp = new double[2*HamSig->n_nonzeros]();
        c_dcopy(2*HamSig->n_nonzeros,(double*)HamSig->nnz,1,HS_nnz_inp,1);
 
        int n_nonzeros_global;
        MPI_Allreduce(&HamSig->n_nonzeros,&n_nonzeros_global,1,MPI_INT,MPI_SUM,matrix_comm);
        int info;
 
        PPEXSIOptions options;
        PPEXSISetDefaultOptions(&options);
        options.npSymbFact = 1;
        options.ordering = 1;
  
        PPEXSIPlan plan;
        int nprowcol[2]={0,0};
        MPI_Dims_create(matrix_procs,2,nprowcol);
        plan = PPEXSIPlanInitialize(matrix_comm,nprowcol[0],nprowcol[1],-1,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
        PPEXSILoadRealSymmetricHSMatrix(plan,options,HamSig->size_tot,n_nonzeros_global,HamSig->n_nonzeros,HamSig->size,HamSig->edge_i,HamSig->index_j,HS_nnz_inp,1,NULL,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
        PPEXSISymbolicFactorizeComplexSymmetricMatrix(plan,options,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
        PPEXSISelInvComplexSymmetricMatrix(plan,options,HS_nnz_inp,(double*)HamSig->nnz,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
        PPEXSIPlanFinalize(plan,&info);
        if (info) return (LOGCERR, EXIT_FAILURE);
 
        delete[] HS_nnz_inp;
        Ps->add_real(HamSig,-weight/M_PI*CPX(0.0,1.0));
        delete HamSig;
    } else if (method==transport_methods::WF) {
        TCSR<CPX> *H1cut = new TCSR<CPX>(HamSig,tra_block*ntriblock,ntriblock,(tra_block+1)*ntriblock,ntriblock);
        TCSR<CPX> *H1 = new TCSR<CPX>(H1cut,0,matrix_comm);
        delete H1cut;
        sabtime=get_time(d_zer);
        LinearSolver<CPX>* solver;
        if (solver_method==12) {
            solver = new SuperLU<CPX>(HamSig,matrix_comm);
#ifdef HAVE_MUMPS
        } else if (solver_method==13) {
            solver = new MUMPS<CPX>(HamSig,matrix_comm);
#endif
        } else return (LOGCERR, EXIT_FAILURE);
        solver->prepare();
        CPX* sol = new CPX[dist_sol[matrix_rank]*(nprol+npror)]();
        solver->solve_equation(sol, inj, nprol+npror);
        delete[] inj;
        delete solver;
        delete HamSig;
if (!worldrank) cout << "TIME FOR WAVEFUNCTION SPARSE SOLVER " << get_time(sabtime) << endl;
        int solsize=displc_sol[matrix_procs];
        CPX* Sol = new CPX[solsize*(nprol+npror)];
        for (int icol=0;icol<nprol+npror;icol++) {
            MPI_Allgatherv(&sol[dist_sol[matrix_rank]*icol],dist_sol[matrix_rank],MPI_DOUBLE_COMPLEX,&Sol[solsize*icol],dist_sol,displc_sol,MPI_DOUBLE_COMPLEX,matrix_comm);
        }
        delete[] sol;
        delete[] dist_sol;
        delete[] displc_sol;
        if (parameter_sab->method==3) {
sabtime=get_time(d_zer);
            if (muvec[0]>muvec[1]) {
                double fermil = fermi(real(energy),muvec[0],parameter_sab->temperature,0)-fermi(real(energy),muvec[1],parameter_sab->temperature,0);
                CPX* SolT = new CPX[solsize*nprol];
                full_transpose(nprol,solsize,Sol,SolT);
                Ps->psipsidagger_transpose(SolT,nprol,+weight*fermil);
                delete[] SolT;
            } else {
                double fermir = fermi(real(energy),muvec[1],parameter_sab->temperature,0)-fermi(real(energy),muvec[0],parameter_sab->temperature,0);
                CPX* SolT = new CPX[solsize*npror];
                full_transpose(npror,solsize,&Sol[solsize*nprol],SolT);
                Ps->psipsidagger_transpose(SolT,npror,+weight*fermir);
                delete[] SolT;
            }
if (!worldrank) cout << "TIME FOR CONSTRUCTION OF S-PATTERN DENSITY MATRIX " << get_time(sabtime) << endl;
        }
        if (parameter_sab->method==2) {
            CPX* Soll = new CPX[solsize*(nprol+npror)]();
            CPX* Solr = new CPX[solsize*(nprol+npror)]();
            for (int ibw=bandwidth+1;ibw<ncells;ibw++) {
//            for (int ibw=0;ibw<ncells;ibw++) {
                c_zlacpy('A',ndof,nprol+npror,Sol,solsize,&Soll[ndof*ibw],solsize);
                c_zlacpy('A',ndof,nprol+npror,&Sol[solsize-ndof],solsize,&Solr[solsize-ndof*(ibw+1)],solsize);
                for (int ipro=0;ipro<nprol+npror;ipro++) {
                    c_zscal(solsize,pow(lambda[ipro],-1),&Soll[ipro*solsize],1);
                    c_zscal(solsize,pow(lambda[ipro],+1),&Solr[ipro*solsize],1);
                }
            }
            double fermil=fermi(real(energy),muvec[0],parameter_sab->temperature,0);
            double fermir=fermi(real(energy),muvec[1],parameter_sab->temperature,0);
sabtime=get_time(d_zer);
            CPX* SolT = new CPX[solsize*max(nprol,npror)];
            CPX* SollT = new CPX[solsize*max(nprol,npror)];
            CPX* SolrT = new CPX[solsize*max(nprol,npror)];
            full_transpose(nprol,solsize,Sol,SolT);
            full_transpose(nprol,solsize,Soll,SollT);
            full_transpose(nprol,solsize,Solr,SolrT);
            double dosl=Ps->psipsidagger_transpose(Overlap,SolT,SollT,SolrT,nprol,ndof,bandwidth,+weight*fermil);
            full_transpose(npror,solsize,&Sol[solsize*nprol],SolT);
            full_transpose(npror,solsize,&Soll[solsize*nprol],SollT);
            full_transpose(npror,solsize,&Solr[solsize*nprol],SolrT);
            double dosr=Ps->psipsidagger_transpose(Overlap,SolT,SollT,SolrT,npror,ndof,bandwidth,+weight*fermir);
            delete[] SolT;
            delete[] SollT;
            delete[] SolrT;
/*
            double dosl=Ps->psipsidagger(Overlap,Sol,Soll,Solr,nprol,ndof,bandwidth,+weight*fermil);
            double dosr=Ps->psipsidagger(Overlap,&Sol[Ps->size_tot*nprol],&Soll[Ps->size_tot*nprol],&Solr[Ps->size_tot*nprol],npror,ndof,bandwidth,+weight*fermir);
*/
            dos=dosl+dosr;
if (!worldrank) cout << "TIME FOR CONSTRUCTION OF S-PATTERN DENSITY MATRIX " << get_time(sabtime) << endl;
            delete[] Soll;
            delete[] Solr;
        }
        delete[] lambda;
// transmission
        if (!matrix_rank) {
            CPX *vecoutdof=new CPX[ntriblock];
            double transml=d_zer;
            for (int ipro=0;ipro<nprol;ipro++) {
                H1->mat_vec_mult(&Sol[Ps->size_tot*ipro+(tra_block+1)*ntriblock],vecoutdof,1);
                transml+=4*M_PI*imag(c_zdotc(ntriblock,&Sol[Ps->size_tot*ipro+tra_block*ntriblock],1,vecoutdof,1));
            }
            double transmr=d_zer;
            for (int ipro=nprol;ipro<nprol+npror;ipro++) {
                H1->mat_vec_mult(&Sol[Ps->size_tot*ipro+(tra_block+1)*ntriblock],vecoutdof,1);
                transmr+=4*M_PI*imag(c_zdotc(ntriblock,&Sol[Ps->size_tot*ipro+tra_block*ntriblock],1,vecoutdof,1));
            }
            delete[] vecoutdof;
if (abs(abs(transml)-abs(transmr))/max(1.0,min(abs(transml),abs(transmr)))>0.1) cout << "CAUTION: TRANSMISSION " << transml << " " << transmr << " AT ENERGY " << real(energy) << " AT POSITION " << evecpos << endl;
//if (abs(transml)<1.0E-2 && nprol>0 && npror>0) cout << "RED ALERT: ZERO TRANSMISSION IN SPITE OF " << nprol << " AND " << npror << " AT ENERGY " << real(energy) << " AT POSITION " << evecpos << endl;
            transm=transml;
            double diff_fermi=fermi(real(energy),muvec[0],parameter_sab->temperature,0)-fermi(real(energy),muvec[1],parameter_sab->temperature,0);
            current=2.0*E_ELECTRON*E_ELECTRON/(2.0*M_PI*H_BAR)*diff_fermi*transm;
        }
        delete[] Sol;
        delete H1;
    } else return (LOGCERR, EXIT_FAILURE);

    return 0;
}
