#include "c_scf.H"

/**  
 *   \brief Takes the overlap (S) and Kohn-Sham (KS) matrices as input and returns a density matrix (P).
 *          This function acts as the gate to the CP2K's world. 
 *          Here only call functions that evaluate a P matrix.   
 */
void c_scf_method(cp2k_transport_parameters cp2k_transport_params, cp2k_csr_interop_type S, cp2k_csr_interop_type KS, cp2k_csr_interop_type * P)
{
    MPI::Intercomm Comm;
    Comm = MPI::COMM_WORLD;
    int rank = Comm.Get_rank();
 
    transport_parameters* transport_params = new transport_parameters();
    transport_params->n_occ                      = cp2k_transport_params.n_occ;
    transport_params->n_atoms                    = cp2k_transport_params.n_atoms;
    transport_params->method                     = cp2k_transport_params.method;
    transport_params->injection_method           = cp2k_transport_params.injection_method;
    transport_params->linear_solver              = cp2k_transport_params.linear_solver;
    transport_params->rlaxis_integration_method  = cp2k_transport_params.rlaxis_integration_method;
    transport_params->n_abscissae                = cp2k_transport_params.n_abscissae;
    transport_params->n_kpoint                   = cp2k_transport_params.n_kpoint;
    transport_params->num_interval               = cp2k_transport_params.num_interval;
    transport_params->num_contacts               = cp2k_transport_params.num_contacts;
    transport_params->tasks_per_point            = cp2k_transport_params.tasks_per_point;
    transport_params->colzero_threshold          = cp2k_transport_params.colzero_threshold;
    transport_params->eps_limit                  = cp2k_transport_params.eps_limit;
    transport_params->eps_limit_cc               = cp2k_transport_params.eps_limit_cc;
    transport_params->eps_decay                  = cp2k_transport_params.eps_decay;
    transport_params->eps_singularity_curvatures = cp2k_transport_params.eps_singularity_curvatures;
    transport_params->eps_mu                     = cp2k_transport_params.eps_mu;
    transport_params->eps_eigval_degen           = cp2k_transport_params.eps_eigval_degen;
    transport_params->energy_interval            = cp2k_transport_params.energy_interval;
    transport_params->min_interval               = cp2k_transport_params.min_interval;
    transport_params->temperature                = cp2k_transport_params.temperature;
    transport_params->svd_cutoff                 = cp2k_transport_params.svd_cutoff;
    transport_params->n_points_beyn              = cp2k_transport_params.n_points_beyn;
    transport_params->extra_scf                  = cp2k_transport_params.extra_scf;
    transport_params->cutout                     = cp2k_transport_params.cutout[0] || cp2k_transport_params.cutout[1];
    transport_params->NCRC_beyn                  = cp2k_transport_params.ncrc_beyn;
    if (cp2k_transport_params.eps_fermi<=(numeric_limits<double>::epsilon)()) {
        transport_params->eps_fermi              = (numeric_limits<double>::epsilon)();
    } else {
        transport_params->eps_fermi              = cp2k_transport_params.eps_fermi;
    }
    if (cp2k_transport_params.n_rand_beyn>1.0 || cp2k_transport_params.n_rand_beyn<=0.0) {
        transport_params->fac_neigbeyn           = 1.0;
    } else {
        transport_params->fac_neigbeyn           = cp2k_transport_params.n_rand_beyn;
    }
    if (cp2k_transport_params.n_rand_cc_beyn>1.0 || cp2k_transport_params.n_rand_cc_beyn<=0.0) {
        transport_params->fac_neigbeyn_cc        = 1.0;
    } else {
        transport_params->fac_neigbeyn_cc        = cp2k_transport_params.n_rand_cc_beyn;
    }
 
    std::vector<contact_type> contactvec(transport_params->num_contacts);
    int size_muvec=0;
    for (int i_c=0;i_c<transport_params->num_contacts;i_c++) if (cp2k_transport_params.contacts_data[4+5*i_c]) size_muvec++;
    std::vector<double> muvec(size_muvec);
    int i_m=0;
    int i_n=0;
    for (int i_c=0;i_c<transport_params->num_contacts;i_c++) {
        int i_t;
        if (cp2k_transport_params.contacts_data[4+5*i_c]) {
            i_t=i_m++;
        } else {
            i_t=muvec.size()+i_n++;
        }
        int bandwidth  = cp2k_transport_params.contacts_data[0+5*i_c];
        int inj_sign   = cp2k_transport_params.contacts_data[3+5*i_c];
        int atom_start = cp2k_transport_params.contacts_data[1+5*i_c];
        int n_atoms_c  = cp2k_transport_params.contacts_data[2+5*i_c];
        if (atom_start==-1) {
            if (!cp2k_transport_params.contacts_data[4+5*i_c]) {
                atom_start = cp2k_transport_params.n_atoms/2;
            } else if (inj_sign==-1) {
                atom_start = cp2k_transport_params.n_atoms-n_atoms_c*bandwidth-cp2k_transport_params.cutout[1];
            } else {
                atom_start = cp2k_transport_params.cutout[0];
            }
        }
        int atom_stop  = atom_start + n_atoms_c;
        contactvec[i_t].bandwidth = bandwidth;
        contactvec[i_t].inj_sign  = inj_sign;
        contactvec[i_t].start     = std::accumulate(&cp2k_transport_params.nsgf[cp2k_transport_params.cutout[0]],&cp2k_transport_params.nsgf[atom_start],0);
        contactvec[i_t].ndof      = std::accumulate(&cp2k_transport_params.nsgf[atom_start],&cp2k_transport_params.nsgf[atom_stop],0);
        contactvec[i_t].n_ele     = std::accumulate(&cp2k_transport_params.zeff[atom_start],&cp2k_transport_params.zeff[atom_stop],0.0);
    }

    int cut_l = std::accumulate(&cp2k_transport_params.nsgf[0],&cp2k_transport_params.nsgf[cp2k_transport_params.cutout[0]],0);
    int cut_r = std::accumulate(&cp2k_transport_params.nsgf[cp2k_transport_params.n_atoms-cp2k_transport_params.cutout[1]],&cp2k_transport_params.nsgf[cp2k_transport_params.n_atoms],0);

    TCSR<double> *Overlap, *KohnSham;
    if (transport_params->cutout) {
        Overlap  = new TCSR<double>( S,cut_l,S.nrows_total-cut_l-cut_r,cut_l,S.nrows_total-cut_l-cut_r);
        KohnSham = new TCSR<double>(KS,cut_l,S.nrows_total-cut_l-cut_r,cut_l,S.nrows_total-cut_l-cut_r);
    } else {
        Overlap  = new TCSR<double>( S.nrows_local, S.nze_local,1);
        KohnSham = new TCSR<double>(KS.nrows_local,KS.nze_local,1);
        cp2kCSR_to_CSR( S, Overlap);
        cp2kCSR_to_CSR(KS,KohnSham);
    }
    c_dscal(KohnSham->n_nonzeros,cp2k_transport_params.evoltfactor,KohnSham->nnz,1);

    int wr_cutblocksize = 0;
    int wr_bw           = 0;
    int wr_ndof         = 0;
    if (!transport_params->cutout) {
        wr_bw   = contactvec[0].bandwidth;
        wr_ndof = contactvec[0].ndof;
    }
    switch (transport_params->method) {
        case 0:
            if (!rank) cout << "Writing Matrices" << endl;
            write_matrix(Overlap,KohnSham,wr_cutblocksize,wr_bw,wr_ndof);
            break;
        case 1:
            if (!rank) cout << "Starting ScaLaPackDiag" << endl;
            if (diagscalapack(Overlap,KohnSham,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
            break;
        case 2:
            if (!rank) cout << "Starting CP2K core/valence Hamiltonian + OMEN Poisson local self consistent code" << endl;
            if (semiselfconsistent(Overlap,KohnSham,muvec,contactvec,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
            break;
        case 3:
        case 4:
        default:
            if (!rank) cout << "Starting Transport " << transport_params->method << endl;
            Energyvector energyvector;
            if (energyvector.Execute(Overlap,KohnSham,muvec,contactvec,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
    }
 
    if (cp2k_transport_params.extra_scf) {
        c_dscal(P->nze_local,0.5,P->nzvals_local,1);
    } else {
        if (transport_params->cutout) {
            c_dscal(P->nze_local,0.5,P->nzvals_local,1);
            Overlap->copy_shifted(*P,cut_l,S.nrows_total-cut_l-cut_r,cut_l,S.nrows_total-cut_l-cut_r);
        } else {
            CSR_to_cp2kCSR(Overlap,*P);
        }
    }

    delete Overlap;
    delete KohnSham;
    delete transport_params;
}

