#include "Types.H"
#include "libcp2k.h"
#include "DiagScaLapack.H"
#include "SemiSelfConsistent.H"
#include "EnergyVector.H"
#include "WriteMatrix.H"

/**  
 *   \brief Takes the overlap (S) and Kohn-Sham (KS) matrices as input and returns a density matrix (P).
 *          This function acts as the gate to the CP2K's world. 
 *          Here only call functions that evaluate a P matrix.   
 */

extern int scf_iter_sab;

void add_to_cp2k_csr(cp2k_csr_interop_type& cp2kCSRmat,double factor,int start_i_fr,int start_j_fr,int start_i_to,int start_j_to,int length_i,int length_j,MPI_Comm cp2k_comm)
{
    TCSR<double> *Pmc = new TCSR<double>(cp2kCSRmat,start_i_fr,length_i,start_j_fr,length_j);
    TCSR<double> *Pm  = new TCSR<double>(Pmc,cp2k_comm);
    delete Pmc;
    double *Pf = new double[length_i*length_j];
    Pm->sparse_to_full(Pf,length_i,length_j);
    delete Pm;
    for(int r=0;r<cp2kCSRmat.nrows_local;r++){
        int i=r+cp2kCSRmat.first_row-start_i_to;
        if(i>=0 && i<length_i){
            for(int e=cp2kCSRmat.rowptr_local[r]-1;e<cp2kCSRmat.rowptr_local[r+1]-1;e++){
                int j=cp2kCSRmat.colind_local[e]-1-start_j_to;
                if(j>=0 && j<length_j){
                    cp2kCSRmat.nzvals_local[e]+=factor*Pf[i+j*length_i];
                }
            }
        }
    }
    delete[] Pf;
}

void copy_to_cp2k_csr(cp2k_csr_interop_type& cp2kCSRmat,double factor,int start_i_fr,int start_j_fr,int start_i_to,int start_j_to,int length_i,int length_j,MPI_Comm cp2k_comm)
{
    TCSR<double> *Pmc = new TCSR<double>(cp2kCSRmat,start_i_fr,length_i,start_j_fr,length_j);
    TCSR<double> *Pm  = new TCSR<double>(Pmc,cp2k_comm);
    delete Pmc;
    double *Pf = new double[length_i*length_j];
    Pm->sparse_to_full(Pf,length_i,length_j);
    delete Pm;
    for(int r=0;r<cp2kCSRmat.nrows_local;r++){
        int i=r+cp2kCSRmat.first_row-start_i_to;
        if(i>=0 && i<length_i){
            for(int e=cp2kCSRmat.rowptr_local[r]-1;e<cp2kCSRmat.rowptr_local[r+1]-1;e++){
                int j=cp2kCSRmat.colind_local[e]-1-start_j_to;
                if(j>=0 && j<length_j){
                    cp2kCSRmat.nzvals_local[e]=factor*Pf[i+j*length_i];
                }
            }
        }
    }
    delete[] Pf;
}

#ifdef HAVE_PIMAG
void c_scf_method(cp2k_transport_parameters cp2k_transport_params, cp2k_csr_interop_type S, cp2k_csr_interop_type KS, cp2k_csr_interop_type * P, cp2k_csr_interop_type * PImag)
{
#else
void c_scf_method(cp2k_transport_parameters cp2k_transport_params, cp2k_csr_interop_type S, cp2k_csr_interop_type KS, cp2k_csr_interop_type * P)
{
    cp2k_csr_interop_type * PImag = NULL;
#endif
    MPI::Intercomm Comm;
    Comm = MPI::COMM_WORLD;
    int rank = Comm.Get_rank();

    c_dscal(P->nze_local,0.5,P->nzvals_local,1);
 
    transport_parameters* transport_params = new transport_parameters();
    transport_params->n_occ                       = cp2k_transport_params.n_occ;
    transport_params->method                      = cp2k_transport_params.method;
    transport_params->injection_method            = cp2k_transport_params.injection_method;
    transport_params->linear_solver               = cp2k_transport_params.linear_solver;
    transport_params->rlaxis_integration_method   = cp2k_transport_params.rlaxis_integration_method;
    transport_params->n_abscissae                 = cp2k_transport_params.n_abscissae;
    transport_params->n_kpoint                    = cp2k_transport_params.n_kpoint;
    transport_params->num_interval                = cp2k_transport_params.num_interval;
    transport_params->num_contacts                = cp2k_transport_params.num_contacts;
    transport_params->tasks_per_point             = cp2k_transport_params.tasks_per_energy_point;
    transport_params->gpus_per_point              = cp2k_transport_params.gpus_per_point;
    transport_params->colzero_threshold           = cp2k_transport_params.colzero_threshold;
    transport_params->eps_limit                   = cp2k_transport_params.eps_limit;
    transport_params->eps_limit_cc                = cp2k_transport_params.eps_limit_cc;
    transport_params->eps_decay                   = cp2k_transport_params.eps_decay;
    transport_params->eps_singularity_curvatures  = cp2k_transport_params.eps_singularity_curvatures;
    transport_params->eps_mu                      = cp2k_transport_params.eps_mu;
    transport_params->eps_eigval_degen            = cp2k_transport_params.eps_eigval_degen;
    transport_params->energy_interval             = cp2k_transport_params.energy_interval;
    transport_params->min_interval                = cp2k_transport_params.min_interval;
    transport_params->temperature                 = cp2k_transport_params.temperature*cp2k_transport_params.evoltfactor/K_BOLTZMANN;
    transport_params->svd_cutoff                  = cp2k_transport_params.svd_cutoff;
    transport_params->n_points_beyn               = cp2k_transport_params.n_points_beyn;
    transport_params->tasks_per_integration_point = cp2k_transport_params.tasks_per_integration_point;
    transport_params->extra_scf                   = cp2k_transport_params.extra_scf;
    transport_params->n_atoms                     = cp2k_transport_params.n_atoms-cp2k_transport_params.cutout[0]-cp2k_transport_params.cutout[1];
    transport_params->evoltfactor                 = cp2k_transport_params.evoltfactor;
    transport_params->NCRC_beyn                   = cp2k_transport_params.ncrc_beyn;
    transport_params->obc                         = cp2k_transport_params.cutout[0] || cp2k_transport_params.cutout[1];
    if (cp2k_transport_params.eps_fermi<=(numeric_limits<double>::epsilon)()) {
        transport_params->eps_fermi               = (numeric_limits<double>::epsilon)();
    } else {
        transport_params->eps_fermi               = cp2k_transport_params.eps_fermi;
    }
    if (cp2k_transport_params.n_rand_beyn>1.0 || cp2k_transport_params.n_rand_beyn<=0.0) {
        transport_params->fac_neigbeyn            = 1.0;
    } else {
        transport_params->fac_neigbeyn            = cp2k_transport_params.n_rand_beyn;
    }
    if (cp2k_transport_params.n_rand_cc_beyn>1.0 || cp2k_transport_params.n_rand_cc_beyn<=0.0) {
        transport_params->fac_neigbeyn_cc         = 1.0;
    } else {
        transport_params->fac_neigbeyn_cc         = cp2k_transport_params.n_rand_cc_beyn;
    }

    std::vector<int> Bsizes;
    for (int i=0;i<cp2k_transport_params.n_blocks;i++) {
        Bsizes.push_back(cp2k_transport_params.tridiag_blocks[i]);
    }
    Bsizes[0]-=cp2k_transport_params.cutout[0];
    Bsizes[Bsizes.size()-1]-=cp2k_transport_params.cutout[1];
    std::vector<int> orb_per_atom(1,0);
    for (int iii=0;iii<transport_params->n_atoms;iii++) {
        orb_per_atom.push_back(orb_per_atom[iii]+cp2k_transport_params.nsgf[cp2k_transport_params.cutout[0]+iii]);//partial_sum
    }
 
    std::vector<contact_type> contactvec(transport_params->num_contacts);
    int size_muvec=0;
    int stride = 6;
    for (int i_c=0;i_c<transport_params->num_contacts;i_c++) if (cp2k_transport_params.contacts_data[4+stride*i_c]) size_muvec++;
    std::vector<double> muvec(size_muvec);
    int i_m=0;
    int i_n=0;
    for (int i_c=0;i_c<transport_params->num_contacts;i_c++) {
        int i_t;
        if (cp2k_transport_params.contacts_data[4+stride*i_c]) {
            i_t=i_m++;
        } else {
            i_t=muvec.size()+i_n++;
        }
        int bandwidth  = cp2k_transport_params.contacts_data[0+stride*i_c];
        int inj_sign   = cp2k_transport_params.contacts_data[3+stride*i_c];
        int sigma_atom = cp2k_transport_params.contacts_data[5+stride*i_c];
        int atom_start = cp2k_transport_params.contacts_data[1+stride*i_c];
        int n_atoms_c  = cp2k_transport_params.contacts_data[2+stride*i_c];
        int atom_stop  = atom_start + n_atoms_c;
        contactvec[i_t].bandwidth = bandwidth;
        contactvec[i_t].inj_sign  = inj_sign;
        contactvec[i_t].natoms    = n_atoms_c;
        contactvec[i_t].atomstart = atom_start-cp2k_transport_params.cutout[0];
        contactvec[i_t].start     = std::accumulate(&cp2k_transport_params.nsgf[cp2k_transport_params.cutout[0]],&cp2k_transport_params.nsgf[atom_start],0);
        contactvec[i_t].start_bs  = std::accumulate(&cp2k_transport_params.nsgf[0],&cp2k_transport_params.nsgf[atom_start],0);
        contactvec[i_t].ndof      = std::accumulate(&cp2k_transport_params.nsgf[atom_start],&cp2k_transport_params.nsgf[atom_stop],0);
        contactvec[i_t].n_ele     = std::accumulate(&cp2k_transport_params.zeff[atom_start],&cp2k_transport_params.zeff[atom_stop],0.0);
    }

    transport_params->cutl = std::accumulate(&cp2k_transport_params.nsgf[0],&cp2k_transport_params.nsgf[cp2k_transport_params.cutout[0]],0);
    transport_params->cutr = std::accumulate(&cp2k_transport_params.nsgf[cp2k_transport_params.n_atoms-cp2k_transport_params.cutout[1]],&cp2k_transport_params.nsgf[cp2k_transport_params.n_atoms],0);

//    transport_params->obc = 1;
    transport_params->n_points_inv = 256;

    transport_params->cp2k_scf_iter = ++scf_iter_sab;

    if (!rank) for (uint i=0;i<contactvec.size();i++) cout << "Contact BW " << contactvec[i].bandwidth << endl;

    transport_params->negf_solver=0;
    if (transport_params->eps_decay<0.0) {
        transport_params->negf_solver=1;
        transport_params->eps_decay=-transport_params->eps_decay;
    }

    int do_double=0;
    if (transport_params->method==3) do_double=1;
    for (int i=0;i<=do_double;i++) {
        if (do_double) {
            if (!i) {
                transport_params->cutl=0;
                transport_params->cutr=P->nrows_total/2;
                transport_params->obc=1;
                contactvec[1].start+=-P->nrows_total/2;
                contactvec[1].start_bs+=-P->nrows_total/2;
                transport_params->cp2k_scf_iter=scf_iter_sab;
            } else {
                transport_params->cutl=P->nrows_total/2;
                transport_params->cutr=0;
                transport_params->obc=1;
                contactvec[0].start_bs+=P->nrows_total/2;
                contactvec[1].start_bs+=P->nrows_total/2;
                transport_params->cp2k_scf_iter=-scf_iter_sab;
            }
        }

        int wr_cutblocksize = 0;
        int wr_bw           = 0;
        int wr_ndof         = 0;
        if (!transport_params->obc && transport_params->method==0 && contactvec.size()>0) {
            wr_bw   = contactvec[0].bandwidth;
            wr_ndof = contactvec[0].ndof;
        }
        switch (transport_params->method) {
            case 0:
                if (!rank) cout << "Writing Matrices" << endl;
                //write_matrix(Overlap,KohnSham,wr_cutblocksize,wr_bw,wr_ndof);
                break;
            case 1:
                if (!rank) cout << "Starting ScaLaPackDiag" << endl;
                //if (diagscalapack(Overlap,KohnSham,transport_params)) throw std::exception();
                break;
            case 2:
                if (!rank) cout << "Starting CP2K core/valence Hamiltonian + OMEN Poisson local self consistent code" << endl;
                if (semiselfconsistent(S,KS,P,PImag,muvec,contactvec,Bsizes,orb_per_atom,transport_params)) throw std::exception();
                break;
            case 3:
            case 4:
            default:
                if (!rank) cout << "Starting Transport " << transport_params->method << endl;
                Energyvector energyvector;
                if (energyvector.Execute(S,KS,P,PImag,muvec,contactvec,Bsizes,orb_per_atom,NULL,NULL,transport_params)) throw std::exception();
        }

    }

    if (!transport_params->extra_scf && transport_params->obc) {
        for (int i=0;i<contactvec[0].bandwidth;i++) {
            for (int ii=0;ii<=i;ii++) {
                int size_tot=P->nrows_total;
                int offset=0;
                copy_to_cp2k_csr(*P,0.5,\
                        (i+1+offset)*contactvec[0].ndof,\
                        offset*contactvec[0].ndof,\
                        ii*contactvec[0].ndof,\
                        size_tot+(ii-(i+1))*contactvec[0].ndof,\
                        contactvec[0].ndof,contactvec[0].ndof,MPI_COMM_WORLD);
                add_to_cp2k_csr(*P,0.5,\
                        size_tot-(offset+1)*contactvec[0].ndof,\
                        size_tot-(i+2+offset)*contactvec[0].ndof,\
                        ii*contactvec[0].ndof,\
                        size_tot+(ii-(i+1))*contactvec[0].ndof,\
                        contactvec[0].ndof,contactvec[0].ndof,MPI_COMM_WORLD);
                copy_to_cp2k_csr(*P,0.5,\
                        offset*contactvec[0].ndof,\
                        (i+1+offset)*contactvec[0].ndof,\
                        size_tot+(ii-(i+1))*contactvec[0].ndof,\
                        ii*contactvec[0].ndof,\
                        contactvec[0].ndof,contactvec[0].ndof,MPI_COMM_WORLD);
                add_to_cp2k_csr(*P,0.5,\
                        size_tot-(i+2+offset)*contactvec[0].ndof,\
                        size_tot-(offset+1)*contactvec[0].ndof,\
                        size_tot+(ii-(i+1))*contactvec[0].ndof,\
                        ii*contactvec[0].ndof,\
                        contactvec[0].ndof,contactvec[0].ndof,MPI_COMM_WORLD);
            }
        }
    }
    if (do_double) {
        for (int i=0;i<contactvec[0].bandwidth;i++) {
            for (int ii=0;ii<=i;ii++) {
                int size_tot=P->nrows_total;
                int mid=size_tot/2;
                copy_to_cp2k_csr(*P,0.5,\
                        mid+(-i-2)*contactvec[0].ndof,\
                        mid+(-1)*contactvec[0].ndof,\
                        mid+(-i-1+ii)*contactvec[0].ndof,\
                        mid+(0+ii)*contactvec[0].ndof,\
                        contactvec[0].ndof,contactvec[0].ndof,MPI_COMM_WORLD);
                add_to_cp2k_csr(*P,0.5,\
                        mid+0*contactvec[0].ndof,\
                        mid+(i+1)*contactvec[0].ndof,\
                        mid+(-i-1+ii)*contactvec[0].ndof,\
                        mid+(0+ii)*contactvec[0].ndof,\
                        contactvec[0].ndof,contactvec[0].ndof,MPI_COMM_WORLD);
                copy_to_cp2k_csr(*P,0.5,\
                        mid+(-1)*contactvec[0].ndof,\
                        mid+(-i-2)*contactvec[0].ndof,\
                        mid+(0+ii)*contactvec[0].ndof,\
                        mid+(-i-1+ii)*contactvec[0].ndof,\
                        contactvec[0].ndof,contactvec[0].ndof,MPI_COMM_WORLD);
                add_to_cp2k_csr(*P,0.5,\
                        mid+(i+1)*contactvec[0].ndof,\
                        mid+0*contactvec[0].ndof,\
                        mid+(0+ii)*contactvec[0].ndof,\
                        mid+(-i-1+ii)*contactvec[0].ndof,\
                        contactvec[0].ndof,contactvec[0].ndof,MPI_COMM_WORLD);
            }
        }
    }

    if (!transport_params->extra_scf) {
        double *mulli=new double[S.nrows_total]();
        for (int i=0;i<S.nrows_local;i++) {
            for (int e=S.rowptr_local[i]-1;e<S.rowptr_local[i+1]-1;e++) {
                mulli[i+S.first_row]+=S.nzvals_local[e]*P->nzvals_local[e];
            }
        }
        MPI_Allreduce(MPI_IN_PLACE,mulli,S.nrows_total,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if (!rank) {
            stringstream mysstream;
            mysstream << "Density_BF_" << scf_iter_sab;
            ofstream myfile(mysstream.str().c_str());
            myfile.precision(8);
            for (int i=0;i<S.nrows_total;i++) myfile<<mulli[i]<<"\n";
            myfile.close();
        }
        delete[] mulli;
    }
 
    delete transport_params;
}

