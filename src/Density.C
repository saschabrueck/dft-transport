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
#include "Umfpack.H"
#include "MUMPS.H"
#include "SuperLU.H"
#include "tmprGF.H"
#include "c_pexsi_interface.h"
#include "GetSigma.H"
#ifdef HAVE_PARDISO
#include "Pardiso.H"
#include "SpikeSolver.H"
#endif
#include "Density.H"

extern "C" {
int ldlt_preprocess__(int *, int *, int *, int *, int *, int *, int *);
int ldlt_fact__(int *, int *, int *, CPX *);
int ldlt_free__(int *);
int ldlt_blkselinv__(int *, int*, int*, CPX *, int*);
}

int density(TCSR<double> *KohnSham,TCSR<double> *Overlap,TCSR<double> *OverlapPBC,TCSR<double> *Ps,CPX energy,CPX weight,transport_methods::transport_method method,int n_mu,double *muvec,int *contactvec,double &current,double &transm,double &dos,std::vector<int> propnum,int *atom_of_bf,double *erhoperatom,double *drhoperatom,double *Vatom,c_transport_type parameter_sab,int distribute_pmat,int evecpos,MPI_Comm matrix_comm)
{
    double d_one=1.0;
    double d_zer=0.0;
    CPX z_one=CPX(d_one,d_zer);
    CPX z_zer=CPX(d_zer,d_zer);
    CPX z_img=CPX(d_zer,d_one);
    double sabtime;
    int iinfo;
    TCSR<CPX> *HamSig;
    CPX* inj = NULL;
    int nprol, npror;
    CPX *lambda;
    int *dist_sol;
    int *displc_sol;
    vector<TCSR<CPX>*> Gamma(n_mu);
    int matrix_procs,matrix_rank;
    MPI_Comm_size(matrix_comm,&matrix_procs);
    MPI_Comm_rank(matrix_comm,&matrix_rank);
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
// get parameters always at the beginning of a subroutine (or maybe even better in the constructor) to make it independent of the structure containing the parameters
// but ndof bandwidth etc will be included in contactvec, which will become a more complicated structure
    int ncells=parameter_sab.n_cells;
    int bandwidth=parameter_sab.bandwidth;
    int ndof=Overlap->size_tot/parameter_sab.n_cells;
    int ntriblock=bandwidth*ndof;
    {
        TCSR<CPX> *SumHamC = new TCSR<CPX>(z_one,KohnSham,-energy,Overlap);
// set pbc to zero 
        SumHamC->settozeropbc(bandwidth,ndof);
// compute self energies
        vector<BoundarySelfEnergy> selfenergies(n_mu);
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
            for (int i_bound_id=0;i_bound_id<n_bound_comm;i_bound_id++) {
                int ibpos=i_bound_id+iseq*n_bound_comm;
                if ( selfenergies[ibpos].Cutout(SumHamC,contactvec[ibpos],energy,method,parameter_sab,matrix_comm) ) return (LOGCERR, EXIT_FAILURE);
            }
            if (ipos<n_mu) {
                if ( selfenergies[ipos].GetSigma(boundary_comm,evecpos) ) return (LOGCERR, EXIT_FAILURE);
            }
            for (int i_bound_id=0;i_bound_id<n_bound_comm;i_bound_id++) {
                int ibpos=i_bound_id+iseq*n_bound_comm;
                selfenergies[ibpos].Distribute(SumHamC,matrix_comm);
            }
        }
        MPI_Comm_free(&boundary_comm);
//add sigma to sumhamc
        TCSR<CPX> *SumHamCplusSigma = new TCSR<CPX>(z_one,SumHamC,-z_one,selfenergies[0].spsigmadist);
        delete SumHamC;
        HamSig = new TCSR<CPX>(z_one,SumHamCplusSigma,-z_one,selfenergies[1].spsigmadist);
        delete SumHamCplusSigma;
        if (method==transport_methods::WF) {
            nprol=selfenergies[0].n_propagating;
            npror=selfenergies[1].n_propagating;
            lambda = new CPX[nprol+npror];
            c_zcopy(nprol,selfenergies[0].lambdapro,1,lambda,1);
            c_zcopy(npror,selfenergies[1].lambdapro,1,&lambda[nprol],1);
if (nprol!=propnum[0]) cout << "WARNING: FOUND " << nprol << " OF " << propnum[0] << " LEFT MODES AT " << real(energy) << endl;
if (npror!=propnum[1]) cout << "WARNING: FOUND " << npror << " OF " << propnum[1] << " RIGHT MODES AT " << real(energy) << endl;
            inj = new CPX[HamSig->size*(nprol+npror)]();
            if (selfenergies[0].spainjdist->n_nonzeros) {
                selfenergies[0].spainjdist->sparse_to_full(inj,HamSig->size,nprol);
            }
            if (selfenergies[1].spainjdist->n_nonzeros) {
                selfenergies[1].spainjdist->sparse_to_full(&inj[HamSig->size*nprol],HamSig->size,npror);
            }
            dist_sol = new int[matrix_procs];
            MPI_Allgather(&HamSig->size,1,MPI_INT,dist_sol,1,MPI_INT,matrix_comm);
            displc_sol = new int[matrix_procs+1]();
            for (int iii=1;iii<matrix_procs+1;iii++) {
                displc_sol[iii]=displc_sol[iii-1]+dist_sol[iii-1];
            }
        }
        if (method==transport_methods::NEGF) {
            if (matrix_procs>1) return (LOGCERR, EXIT_FAILURE);
            for (int i_mu=0;i_mu<n_mu;i_mu++) {
                TCSR<CPX> *Sigma_H = new TCSR<CPX>(selfenergies[i_mu].spsigmadist);
                Sigma_H->sparse_transpose(selfenergies[i_mu].spsigmadist);
                c_dscal(Sigma_H->n_nonzeros,-d_one,((double*)Sigma_H->nnz)+1,2);
                Gamma[i_mu] = new TCSR<CPX>(z_one,selfenergies[i_mu].spsigmadist,-z_one,Sigma_H);
                delete Sigma_H;
//OR JUST USE IMAGINARY PART OF SIGMA INSTEAD
            }
        }
    }
    if (method==transport_methods::NEGF) {
        int HamSigsize_tot=HamSig->size_tot;
        int HamSigfindx=HamSig->findx;
        int *pivarrays = new int[HamSigsize_tot];
        CPX *HamBig0 = new CPX[HamSigsize_tot*HamSigsize_tot]();
        HamSig->sparse_to_full(HamBig0,HamSigsize_tot,HamSigsize_tot);
        delete HamSig;
        sabtime=get_time(d_zer);
        c_zgetrf(HamSigsize_tot,HamSigsize_tot,HamBig0,HamSigsize_tot,pivarrays,&iinfo);
        if (!worldrank) cout << "TIME FOR LU " << get_time(sabtime) << endl;
        if (iinfo) return (LOGCERR, EXIT_FAILURE);
        CPX nworks;
        c_zgetri(HamSigsize_tot,HamBig0,HamSigsize_tot,pivarrays,&nworks,-1,&iinfo);
        int lworks=int(real(nworks));
        CPX* works=new CPX[lworks];
        sabtime=get_time(d_zer);
        c_zgetri(HamSigsize_tot,HamBig0,HamSigsize_tot,pivarrays,works,lworks,&iinfo);
        if (!worldrank) cout << "TIME FOR INV " << get_time(sabtime) << endl;
        delete[] works;
        delete[] pivarrays;
        CPX *HamBig1 = new CPX[HamSigsize_tot*HamSigsize_tot]();
        CPX *HamBig2 = new CPX[HamSigsize_tot*HamSigsize_tot]();
        c_zcopy(HamSigsize_tot*HamSigsize_tot,HamBig0,1,HamBig1,1);
        c_dscal(HamSigsize_tot*HamSigsize_tot,-d_one,((double*)HamBig1)+1,2);
        sabtime=get_time(d_zer);
        Gamma[1]->trans_mat_vec_mult(HamBig1,HamBig2,HamSigsize_tot,1);
        if (!worldrank) cout << "TIME FOR MULT A " << get_time(sabtime) << endl;
        full_transpose(HamSigsize_tot,HamSigsize_tot,HamBig2,HamBig1);
        TCSR<CPX> *Gamma_G = new TCSR<CPX>(HamSigsize_tot,HamSigsize_tot*ntriblock,HamSigfindx);
        Gamma_G->full_to_sparse(HamBig1,HamSigsize_tot,HamSigsize_tot);
        c_zcopy(HamSigsize_tot*HamSigsize_tot,HamBig0,1,HamBig1,1);
        sabtime=get_time(d_zer);
        Gamma_G->trans_mat_vec_mult(HamBig1,HamBig2,HamSigsize_tot,1);
        if (!worldrank) cout << "TIME FOR MULT B " << get_time(sabtime) << endl;
        delete Gamma_G;
        sabtime=get_time(d_zer);
        Gamma[0]->trans_mat_vec_mult(HamBig2,HamBig1,HamSigsize_tot,1);
        if (!worldrank) cout << "TIME FOR MULT C " << get_time(sabtime) << endl;
        for (int i=0;i<HamSigsize_tot;i++) transm+=-real(HamBig1[(HamSigsize_tot+1)*i]);
        full_transpose(HamSigsize_tot,HamSigsize_tot,HamBig2,HamBig1);
        delete[] HamBig2;
        TCSR<CPX> *HamBigSparse = new TCSR<CPX>(HamSigsize_tot,HamSigsize_tot*HamSigsize_tot,HamSigfindx);
        HamBigSparse->full_to_sparse(HamBig1,HamSigsize_tot,HamSigsize_tot);
        delete[] HamBig1;
//we need the real part of hambig but my routines give me the imaginary one
//        c_zscal(HamBigSparse->n_nonzeros,CPX(0.0,1.0),HamBigSparse->nnz,1);
        Ps->add_real(HamBigSparse,-weight/M_PI*z_img);
        Overlap->atomdensity(HamBigSparse,-2.0*weight/M_PI*z_img,atom_of_bf,erhoperatom);
        delete HamBigSparse;
        delete[] HamBig0;
    } else if (method==transport_methods::GF) {
        sabtime=get_time(d_zer);
        int inversion_method=0;
        if (inversion_method==0) {
            HamSig->change_findx(1);
            int n_nonzeros_global;
            MPI_Allreduce(&HamSig->n_nonzeros,&n_nonzeros_global,1,MPI_INT,MPI_SUM,matrix_comm);
            double *HamSig_nnz = new double[2*HamSig->n_nonzeros];
            c_dcopy(2*HamSig->n_nonzeros,(double*)HamSig->nnz,1,HamSig_nnz,1);
            PSelInvComplexSymmetricInterface(HamSig->size_tot,n_nonzeros_global,HamSig->n_nonzeros,HamSig->size,HamSig->edge_i,HamSig->index_j,HamSig_nnz,2,1,matrix_comm,matrix_procs,1,(double*)HamSig->nnz,&iinfo);
            Ps->add_real(HamSig,-weight/M_PI*z_img);
            Overlap->atomdensity(HamSig,-2.0*weight/M_PI*z_img,atom_of_bf,erhoperatom);
            delete HamSig;
        } else if (inversion_method==3) {
            TCSR<CPX> *HamSigG = new TCSR<CPX>(HamSig,0,matrix_comm);
            if (ncells%bandwidth) return (LOGCERR, EXIT_FAILURE);
            if (distribute_pmat && matrix_procs>1) {
                if (!matrix_rank) {
                    std::vector<int> Bvec(ncells/bandwidth,0);
                    for (int ii=0;ii<ncells/bandwidth;ii++) Bvec[ii]=(ii+1)*ntriblock-1;
                    omp_set_num_threads(omp_get_max_threads()*matrix_procs);
                    tmprGF::sparse_invert(HamSigG,Bvec);
                    omp_set_num_threads(omp_get_max_threads()/matrix_procs);
                }
                HamSigG->scatter(HamSig,0,matrix_comm);
                delete HamSigG;
                Ps->add_real(HamSig,-weight/M_PI*z_img);
                Overlap->atomdensity(HamSig,-2.0*weight/M_PI*z_img,atom_of_bf,erhoperatom);
                delete HamSig;
            } else {
                delete HamSig;
                if (!matrix_rank) {
                    std::vector<int> Bvec(ncells/bandwidth,0);
                    for (int ii=0;ii<ncells/bandwidth;ii++) Bvec[ii]=(ii+1)*ntriblock-1;
                    omp_set_num_threads(omp_get_max_threads()*matrix_procs);
                    tmprGF::sparse_invert(HamSigG,Bvec);
                    omp_set_num_threads(omp_get_max_threads()/matrix_procs);
//                    Ps->add_real(HamSigG,-weight/M_PI*z_img);
                }
/*
                TCSR<CPX> *P_PBC = new TCSR<CPX>(OverlapPBC->size,OverlapPBC->n_nonzeros,OverlapPBC->findx);
                P_PBC->copy_index(OverlapPBC);
                for (int i_mu=0;i_mu<n_mu;i_mu++) {
                    P_PBC->fill_pbc_block(HamSigG,ndof,bandwidth,contactvec[i_mu],matrix_comm);
                }
                Ps->add_real(P_PBC,-weight/M_PI*z_img);
                OverlapPBC->atomdensity(P_PBC,-2.0*weight/M_PI*z_img,atom_of_bf,erhoperatom);
                delete P_PBC;
*/
// /*
                TCSR<CPX> *Green = new TCSR<CPX>(KohnSham->size,KohnSham->n_nonzeros,KohnSham->findx);
                Green->copy_index(KohnSham);
                Green->add(HamSigG,CPX(1.0,0.0));
                for (int i_mu=0;i_mu<n_mu;i_mu++) {
                    Green->replicate_cell(ndof,bandwidth,contactvec[i_mu],bandwidth,matrix_comm);
                }
                Ps->add_real(Green,-weight/M_PI*z_img);
                Overlap->atomdensity(Green,-2.0*weight/M_PI*z_img,atom_of_bf,erhoperatom);
                OverlapPBC->atomdensity(Green,-2.0*weight/M_PI*z_img,atom_of_bf,erhoperatom);
                for (int i_diag=0;i_diag<Green->size;i_diag++) dos+=imag(Green->nnz[Green->diag_pos[i_diag]]);
// */ 
//                if (matrix_procs==1) Overlap->atomdensity(HamSigG,-2.0*weight/M_PI*z_img,atom_of_bf,erhoperatom);
                delete HamSigG;
            }
#ifdef HAVE_PARDISO            
        } else if (inversion_method==1) {
            TCSR<CPX> *HamSigG = new TCSR<CPX>(HamSig,0,matrix_comm);
            if (distribute_pmat && matrix_procs>1) {
                if (!matrix_rank) {
                    Pardiso::sparse_invert(HamSigG);
                }
                HamSigG->scatter(HamSig,0,matrix_comm);
                delete HamSigG;
                Ps->add_real(HamSig,-weight/M_PI*z_img);
                Overlap->atomdensity(HamSig,-2.0*weight/M_PI*z_img,atom_of_bf,erhoperatom);
                delete HamSig;
            } else {
                delete HamSig;
                if (!matrix_rank) {
                    Pardiso::sparse_invert(HamSigG);
                    Ps->add_real(HamSigG,-weight/M_PI*z_img);
                }
                delete HamSigG;
            }
#endif
#ifdef HAVE_SELINV            
        } else if (inversion_method==2) {
            TCSR<CPX> *HamSigG = new TCSR<CPX>(HamSig,0,matrix_comm);
            if (distribute_pmat && matrix_procs>1) return (LOGCERR, EXIT_FAILURE);
            delete HamSig;
            if (!matrix_rank) {
                HamSigG->change_findx(1);
                TCSR<CPX> *HamSigTri= new TCSR<CPX>(HamSigG->size,HamSigG->n_nonzeros/2+HamSigG->size,HamSigG->findx);
                HamSigTri->extract_lower_triangle(HamSigG);
                delete HamSigG;
                TCSC<CPX,int> *HamSigCSC=new TCSC<CPX,int>(HamSigTri,1);
                delete HamSigTri;
                int token=0;
                int *perm = NULL;
                int Lnnz;
                int order=-1;
                ldlt_preprocess__(&token, &HamSigCSC->size, HamSigCSC->edge_j, HamSigCSC->index_i, &Lnnz, &order, perm);   
                ldlt_fact__(&token, HamSigCSC->edge_j, HamSigCSC->index_i, HamSigCSC->nnz);
                delete HamSigCSC;
                int dumpL=0;
                TCSR<CPX> *InverseMatTrans=new TCSR<CPX>(Ps->size_tot,Lnnz,1);
                ldlt_blkselinv__(&token, InverseMatTrans->edge_i, InverseMatTrans->index_j, InverseMatTrans->nnz, &dumpL);
                ldlt_free__(&token);
                TCSR<CPX> *InverseMat=new TCSR<CPX>(Ps->size_tot,Lnnz,1);
                InverseMat->sparse_transpose(InverseMatTrans);
                TCSR<CPX> *InverseMatTotal=new TCSR<CPX>(z_one,InverseMat,z_one,InverseMatTrans);
                for (int idiag=0;idiag<InverseMatTotal->size;idiag++)
                    InverseMatTotal->nnz[InverseMatTotal->diag_pos[idiag]]*=CPX(0.5,0.0);
                delete InverseMat;
                delete InverseMatTrans;
                Ps->add_real(InverseMatTotal,-weight/M_PI*z_img);
                if (matrix_procs==1) Overlap->atomdensity(InverseMatTotal,-2.0*weight/M_PI*z_img,atom_of_bf,erhoperatom);
                delete InverseMatTotal;
            } else {
                delete HamSigG;
            }
#endif
        } else return (LOGCERR, EXIT_FAILURE);
if (!worldrank) cout << "TIME FOR SPARSE INVERSION " << get_time(sabtime) << endl;
    } else if (method==transport_methods::WF) {
        int tra_block=0;
        TCSR<CPX> *H1cut = new TCSR<CPX>(HamSig,tra_block*ntriblock,ntriblock,(tra_block+1)*ntriblock,ntriblock);
        TCSR<CPX> *H1 = new TCSR<CPX>(H1cut,0,matrix_comm);
        delete H1cut;
        sabtime=get_time(d_zer);
        LinearSolver<CPX>* solver;
        int solver_method=0;
        if (solver_method==0) {
            solver = new SuperLU<CPX>(HamSig,matrix_comm);
        } else if (solver_method==1) {
            solver = new MUMPS<CPX>(HamSig,matrix_comm);
#ifdef HAVE_PARDISO            
        } else if (solver_method==2) {
            solver = new SpikeSolver<CPX>(HamSig,matrix_comm);
#endif
        } else if (solver_method==3 && matrix_procs==1) {
            solver = new Umfpack<CPX>(HamSig,matrix_comm);
        } else return (LOGCERR, EXIT_FAILURE);
        solver->prepare();
        CPX* sol = new CPX[dist_sol[matrix_rank]*(nprol+npror)]();
        solver->solve_equation(sol, inj, nprol+npror);
        delete[] inj;
        delete solver;
        delete HamSig;
if (!worldrank) cout << "TIME FOR WAVEFUNCTION SPARSE SOLVER WITH "<< ncells <<" UNIT CELLS " << get_time(sabtime) << endl;
        int solsize=displc_sol[matrix_procs];
        CPX* Sol = new CPX[solsize*(nprol+npror)];
        for (int icol=0;icol<nprol+npror;icol++) {
            MPI_Allgatherv(&sol[dist_sol[matrix_rank]*icol],dist_sol[matrix_rank],MPI_DOUBLE_COMPLEX,&Sol[solsize*icol],dist_sol,displc_sol,MPI_DOUBLE_COMPLEX,matrix_comm);
        }
        delete[] sol;
        delete[] dist_sol;
        delete[] displc_sol;
        CPX* Soll = new CPX[solsize*(nprol+npror)]();
        CPX* Solr = new CPX[solsize*(nprol+npror)]();
        for (int ibw=bandwidth+1;ibw<ncells;ibw++) {
//        for (int ibw=0;ibw<ncells;ibw++) {
            c_zlacpy('A',ndof,nprol+npror,Sol,solsize,&Soll[ndof*ibw],solsize);
            c_zlacpy('A',ndof,nprol+npror,&Sol[solsize-ndof],solsize,&Solr[solsize-ndof*(ibw+1)],solsize);
            for (int ipro=0;ipro<nprol+npror;ipro++) {
                c_zscal(solsize,pow(lambda[ipro],-1),&Soll[ipro*solsize],1);
                c_zscal(solsize,pow(lambda[ipro],+1),&Solr[ipro*solsize],1);
            }
        }
        delete[] lambda;
        double Temp=parameter_sab.temperature;
        double fermil=0.0;
        double fermir=0.0;
        if (!parameter_sab.n_abscissae) {
            fermil=fermi(real(energy),muvec[0],Temp,0);
            fermir=fermi(real(energy),muvec[1],Temp,0);
        } else if (muvec[0]>muvec[1]) {
            fermil=fermi(real(energy),muvec[0],Temp,0)-fermi(real(energy),muvec[1],Temp,0);
        } else {
            fermir=fermi(real(energy),muvec[1],Temp,0)-fermi(real(energy),muvec[0],Temp,0);
        }
//        double dfermil=fermi(real(energy),muvec[0],Temp,2);
//        double dfermir=fermi(real(energy),muvec[1],Temp,2);
        sabtime=get_time(d_zer);
        if (distribute_pmat || matrix_procs==1) {
            Ps->psipsidagger(Sol,Soll,Solr,nprol,ndof,bandwidth,+weight*fermil);
            Ps->psipsidagger(&Sol[Ps->size_tot*nprol],&Soll[Ps->size_tot*nprol],&Solr[Ps->size_tot*nprol],npror,ndof,bandwidth,+weight*fermir);
            dos=Overlap->psipsidagger(Sol,nprol+npror);
//            dos=Overlap->psipsidagger(Sol,nprol);
//            dos=Overlap->psipsidagger(&Sol[Ps->size_tot*nprol],npror);
/*
if (parameter_sab.n_abscissae) {
*/
            Overlap->psipsidagger(Sol,nprol,+2.0*weight*fermil,atom_of_bf,erhoperatom);
            OverlapPBC->psipsidagger(Sol,Soll,nprol,+2.0*weight*fermil,atom_of_bf,erhoperatom);
            OverlapPBC->psipsidagger(Sol,Solr,nprol,+2.0*weight*fermil,atom_of_bf,erhoperatom);
            Overlap->psipsidagger(&Sol[Ps->size_tot*nprol],npror,+2.0*weight*fermir,atom_of_bf,erhoperatom);
            OverlapPBC->psipsidagger(&Sol[Ps->size_tot*nprol],&Soll[Ps->size_tot*nprol],npror,+2.0*weight*fermir,atom_of_bf,erhoperatom);
            OverlapPBC->psipsidagger(&Sol[Ps->size_tot*nprol],&Solr[Ps->size_tot*nprol],npror,+2.0*weight*fermir,atom_of_bf,erhoperatom);
/*
} else {
            Overlap->psipsidagger(Sol,nprol,+2.0*weight*fermil,atom_of_bf,erhoperatom,Vatom,real(energy));
            Overlap->psipsidagger(&Sol[Ps->size_tot*nprol],npror,+2.0*weight*fermir,atom_of_bf,erhoperatom,Vatom,real(energy));
}
*/
/*
            Overlap->psipsidagger(Sol,nprol,-2.0*weight*dfermil,atom_of_bf,drhoperatom);
            OverlapPBC->psipsidagger(Sol,Soll,nprol,-2.0*weight*dfermil,atom_of_bf,drhoperatom);
            OverlapPBC->psipsidagger(Sol,Solr,nprol,-2.0*weight*dfermil,atom_of_bf,drhoperatom);
            Overlap->psipsidagger(&Sol[Ps->size_tot*nprol],npror,-2.0*weight*dfermir,atom_of_bf,drhoperatom);
            OverlapPBC->psipsidagger(&Sol[Ps->size_tot*nprol],&Soll[Ps->size_tot*nprol],npror,-2.0*weight*dfermir,atom_of_bf,drhoperatom);
            OverlapPBC->psipsidagger(&Sol[Ps->size_tot*nprol],&Solr[Ps->size_tot*nprol],npror,-2.0*weight*dfermir,atom_of_bf,drhoperatom);
*/
        } else {
            if (!matrix_rank) {
                Ps->psipsidagger(Sol,Soll,Solr,nprol,ndof,bandwidth,+weight*fermil);
                Ps->psipsidagger(&Sol[Ps->size_tot*nprol],&Soll[Ps->size_tot*nprol],&Solr[Ps->size_tot*nprol],npror,ndof,bandwidth,+weight*fermir);
                dos=Overlap->psipsidagger(Sol,nprol+npror);
            }
        }
        delete[] Soll;
        delete[] Solr;
if (!worldrank) cout << "TIME FOR CONSTRUCTION OF S-PATTERN DENSITY MATRIX " << get_time(sabtime) << endl;
// transmission
        if (!matrix_rank) {
            CPX *vecoutdof=new CPX[ntriblock];
            H1->shift_resize(tra_block*ntriblock,ntriblock,(tra_block+1)*ntriblock,ntriblock);
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
            if (abs(abs(transml)-abs(transmr))>0.1) return (LOGCERR, EXIT_FAILURE);
            transm=transml;
            double diff_fermi=fermi(real(energy),muvec[0],Temp,0)-fermi(real(energy),muvec[1],Temp,0);
            current=2.0*E_ELECTRON*E_ELECTRON/(2.0*M_PI*H_BAR)*diff_fermi*transm;
        }
        delete[] Sol;
        delete H1;
    } else return (LOGCERR, EXIT_FAILURE);

    return 0;
}
