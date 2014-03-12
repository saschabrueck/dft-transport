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
#include "Pardiso.H"
#include "tmprGF.H"
#include "GetSigma.H"
#include "SpikeSolver.H"
#include "Density.H"

extern "C" {
int ldlt_preprocess__(int *, int *, int *, int *, int *, int *, int *);
int ldlt_fact__(int *, int *, int *, CPX *);
int ldlt_free__(int *);
int ldlt_blkselinv__(int *, int*, int*, CPX *, int*);
}

int density(TCSR<double> *KohnSham,TCSR<double> *Overlap,TCSR<double> *Ps,CPX energy,CPX weight,transport_methods::transport_method method,int n_mu,double *muvec,int *contactvec,double &current,int propnum,int *atom_of_bf,double *erhoperatom,double *drhoperatom,c_transport_type parameter_sab,MPI_Comm matrix_comm)
{
    double d_one=1.0;
    double d_zer=0.0;
    CPX z_one=CPX(d_one,d_zer);
    CPX z_zer=CPX(d_zer,d_zer);
    double sabtime;
    int iinfo=0;
    TCSR<CPX> *HamSig;
    CPX* inj = NULL;
    int nprol, npror;
    int *dist_sol;
    int *displc_sol;
    int matrix_procs,matrix_rank;
    MPI_Comm_size(matrix_comm,&matrix_procs);
    MPI_Comm_rank(matrix_comm,&matrix_rank);
// get parameters
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
        int matrix_procs,matrix_rank;
        MPI_Comm_size(matrix_comm,&matrix_procs);
        MPI_Comm_rank(matrix_comm,&matrix_rank);
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
                if ( selfenergies[ipos].GetSigma(boundary_comm) ) return (LOGCERR, EXIT_FAILURE);
            }
            for (int i_bound_id=0;i_bound_id<n_bound_comm;i_bound_id++) {
                int ibpos=i_bound_id+iseq*n_bound_comm;
                selfenergies[ibpos].Distribute(SumHamC,matrix_comm);
            }
        }
        MPI_Comm_free(&boundary_comm);
//add sigma to sumhamc
        TCSR<CPX> *SumHamCplusSigmal = new TCSR<CPX>(z_one,SumHamC,-z_one,selfenergies[0].spsigmaldist);
        delete SumHamC;
        HamSig = new TCSR<CPX>(z_one,SumHamCplusSigmal,-z_one,selfenergies[1].spsigmardist);
        delete SumHamCplusSigmal;
        if (method==transport_methods::WF) {
            nprol=selfenergies[0].n_propagating;
            npror=selfenergies[1].n_propagating;
if (npror!=propnum) cout << "WARNING: FOUND " << npror << " OF " << propnum << " MODES AT " << real(energy) << endl;
            if (selfenergies[0].spainjldist->n_nonzeros || selfenergies[1].spainjrdist->n_nonzeros) {
                inj = new CPX[HamSig->size*(nprol+npror)]();
            }
            if (selfenergies[0].spainjldist->n_nonzeros) {
                selfenergies[0].spainjldist->sparse_to_full(inj,HamSig->size,nprol);
            }
            if (selfenergies[1].spainjrdist->n_nonzeros) {
                selfenergies[1].spainjrdist->sparse_to_full(&inj[HamSig->size*nprol],HamSig->size,npror);
            }
            dist_sol=new int[matrix_procs];
            MPI_Allgather(&HamSig->size,1,MPI_INT,dist_sol,1,MPI_INT,matrix_comm);
            displc_sol=new int[matrix_procs+1]();
            for (int iii=1;iii<matrix_procs+1;iii++) {
                displc_sol[iii]=displc_sol[iii-1]+dist_sol[iii-1];
            }
        }
    }
    if (method==transport_methods::TEST) {
// this is only a lil debugging thing
        int triblocksize=ntriblock*ntriblock;
        CPX* H0cpx=new CPX[triblocksize];
        CPX* sigmal=new CPX[triblocksize];
        CPX* sigmar=new CPX[triblocksize];
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
        delete[] sigmal;
        delete[] sigmar;
        CPX trace=0;
        for (int jtri=0;jtri<ntriblock;jtri++)
            for (int itri=0;itri<ntriblock;itri++)
                trace+=-sigmar[itri*ntriblock+jtri]*sigmal[jtri*ntriblock+itri];
        cout << "MATRIX MATRIX MULTIPLICATIONS AND TRACE FOR TRANSMISSION " << get_time(sabtime) << endl;
        cout << "Energy " << energy << " Transmission " << real(trace) << endl;
    } else if (method==transport_methods::NEGF) {
TCSR<CPX> *HamSigG = new TCSR<CPX>(HamSig,0,matrix_comm);
delete HamSig;
HamSig = HamSigG;
        sabtime=get_time(d_zer);
        int inversion_method=0;
        if (inversion_method==0) {
            if (ncells%bandwidth) return (LOGCERR, EXIT_FAILURE);
if (!matrix_rank) {
            std::vector<int> Bvec(ncells/bandwidth,0);
            for (int ii=0;ii<ncells/bandwidth;ii++) Bvec[ii]=(ii+1)*ntriblock-1;
            tmprGF::sparse_invert(HamSig,Bvec);
            Ps->add_imag(HamSig,+weight/M_PI);
            Overlap->atomdensity(HamSig,+2.0*weight/M_PI,atom_of_bf,erhoperatom);
            Ps->contactdensity(ndof,bandwidth);
}
            delete HamSig;
#ifdef HAVE_PARDISO            
        } else if (inversion_method==1) {
if (!matrix_rank) {
            Pardiso::sparse_invert(HamSig);
            Ps->add_imag(HamSig,+weight/M_PI);
            Ps->contactdensity(ndof,bandwidth);
}
            delete HamSig;
#endif
#ifdef HAVE_SELINV            
        } else if (inversion_method==2) {
if (!matrix_rank) {
            HamSig->change_findx(1);
            TCSR<CPX> *HamSigTri= new TCSR<CPX>(HamSig->size,HamSig->n_nonzeros/2+HamSig->size,HamSig->findx);
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
            Ps->add_imag(InverseMatTotal,+weight/M_PI);
            delete InverseMatTotal;
            Ps->contactdensity(ndof,bandwidth);
} else {
            delete HamSig;
}
#endif
        } else return (LOGCERR, EXIT_FAILURE);
        cout << "TIME FOR SPARSE INVERSION " << get_time(sabtime) << endl;
    } else if (method==transport_methods::WF) {
// sparse solver
        TCSR<CPX> *H1cut = new TCSR<CPX>(HamSig,0,ntriblock,ntriblock,ntriblock);
        TCSR<CPX> *H1 = new TCSR<CPX>(H1cut,0,matrix_comm);
        delete H1cut;
        sabtime=get_time(d_zer);
        CPX *Sol = NULL;
        int solver_method=0;
        if (solver_method==0) { 
            CPX *Inj = NULL;
            if (!matrix_rank) {
                Inj = new CPX[displc_sol[matrix_procs]*(nprol+npror)];
            }
            if (inj==NULL) inj = new CPX[dist_sol[matrix_rank]*(nprol+npror)](); 
            for (int icol=0;icol<nprol+npror;icol++) {
                MPI_Gatherv(&inj[dist_sol[matrix_rank]*icol],dist_sol[matrix_rank],MPI_DOUBLE_COMPLEX,&Inj[displc_sol[matrix_procs]*icol],dist_sol,displc_sol,MPI_DOUBLE_COMPLEX,0,matrix_comm);
            }
            delete[] inj;
            LinearSolver<CPX>* solver;
            int do_umfpack=0;
            if (matrix_procs==1 && do_umfpack==1) {
                solver = new Umfpack<CPX>(HamSig,matrix_comm);
            } else {
                solver = new MUMPS<CPX>(HamSig,matrix_comm,2);
            }
            delete HamSig;
            solver->prepare();
if (!matrix_rank) {
            Sol = new CPX[displc_sol[matrix_procs]*(nprol+npror)];
}
            solver->solve_equation(Sol,Inj,nprol+npror);
            delete solver;
            if (!matrix_rank) {
                delete[] Inj;
            }
            sabtime=get_time(d_zer);
        } else if (solver_method==1) {
            if (matrix_procs != 2) {
              printf("WARNING: SpikeSolver needs to be tested with more than two partitions\n");
            }
            LinearSolver<CPX>* solver;
            // SPIKE solver
            solver = new SpikeSolver<CPX>(HamSig, matrix_comm);
            int bandwidth = 1092;     // TODO: write routine to find bandwidth
            std::vector<int> p_lines(matrix_procs+1);
            for (int i = 0; i < matrix_procs + 1; ++i) {
              p_lines[i] = displc_sol[i];
            }
            solver->prepare(bandwidth, p_lines);
            CPX* sol = new CPX[dist_sol[matrix_rank]*(nprol+npror)]();
            solver->solve_equation(sol, inj, nprol+npror);

            delete[] inj;
            //delete solver; // TODO: probably
            delete HamSig;
if (!matrix_rank) {
            Sol = new CPX[displc_sol[matrix_procs]*(nprol+npror)];
}
            for (int icol=0;icol<nprol+npror;icol++) {
                MPI_Gatherv(&sol[dist_sol[matrix_rank]*icol],dist_sol[matrix_rank],MPI_DOUBLE_COMPLEX,&Sol[displc_sol[matrix_procs]*icol],dist_sol,displc_sol,MPI_DOUBLE_COMPLEX,0,matrix_comm);
            }
            delete[] sol;
        } else return (LOGCERR, EXIT_FAILURE);
        cout << "TIME FOR WAVEFUNCTION SPARSE SOLVER WITH "<< ncells <<" UNIT CELLS " << get_time(sabtime) << endl;
        delete[] dist_sol;
        delete[] displc_sol;
        double Temp=parameter_sab.extra_param3;
        double fermil=fermi(real(energy),muvec[0],Temp,0);
        double fermir=fermi(real(energy),muvec[1],Temp,0);
if (!matrix_rank) {
// WARNING THOSE ONLY WORK FOR COMPLETE MATRIX ON THIS RANK
        Ps->psipsidagger(Sol,nprol,-weight*fermil);
        Ps->psipsidagger(&Sol[Ps->size_tot*nprol],npror,-weight*fermir);
        Ps->contactdensity(ndof,bandwidth);//do this later at the end of Energyvector in parallel
        current=Overlap->psipsidaggerdosdebug(Sol,nprol+npror);
        Overlap->psipsidagger(Sol,nprol,-2.0*weight*fermil,atom_of_bf,erhoperatom);
        Overlap->psipsidagger(&Sol[Ps->size_tot*nprol],npror,-2.0*weight*fermir,atom_of_bf,erhoperatom);
        double dfermil=fermi(real(energy),muvec[0],Temp,2);
        double dfermir=fermi(real(energy),muvec[1],Temp,2);
        Overlap->psipsidagger(Sol,nprol,-2.0*weight*dfermil,atom_of_bf,drhoperatom);
        Overlap->psipsidagger(&Sol[Ps->size_tot*nprol],npror,-2.0*weight*dfermir,atom_of_bf,drhoperatom);
        cout << "TIME FOR CONSTRUCTION OF S-PATTERN DENSITY MATRIX " << get_time(sabtime) << endl;
}//end if !matrix_rank
// transmission
        double transml=d_zer;
        CPX *vecoutdof=new CPX[ntriblock];
        if (!matrix_rank) {
            H1->shift_resize(0,ntriblock,ntriblock,ntriblock);
            for (int ipro=0;ipro<nprol;ipro++) {
                H1->mat_vec_mult(&Sol[Ps->size_tot*ipro+ntriblock],vecoutdof,1);
                transml+=4*M_PI*imag(c_zdotc(ntriblock,&Sol[Ps->size_tot*ipro],1,vecoutdof,1));
            }
            delete[] Sol;
        }
        delete H1;
        delete[] vecoutdof;
int worldrank; MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
        cout << worldrank << " Energy " << energy << " Transmission " << transml << endl;
    } else return (LOGCERR, EXIT_FAILURE);

    return 0;
}
