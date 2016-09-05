#include "libcp2k.h"
#include "CSR.H"
#ifdef HAVE_OMEN_POISSON
#include "SemiSelfConsistent.H"
#endif
#include "EnergyVector.H"
//#include "DiagScaLapack.H"
//#include "WriteMatrix.H"

/**  
 *   \brief Takes the overlap (S) and Kohn-Sham (KS) matrices as input and returns a density matrix (P).
 *          This function acts as the gate to the CP2K's world. 
 *          Here only call functions that evaluate a P matrix.   
 */

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
    int rank,mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (!rank) cout << "Starting Transport" << endl;

    c_dscal(P->nze_local,0.5,P->nzvals_local,1);
 
    std::vector<int> atom_of_bf;
    for (int a=0;a<cp2k_transport_params.n_atoms;a++) {
        for (int i=0;i<cp2k_transport_params.nsgf[a];i++) {
            atom_of_bf.push_back(a);
        }
    }
 
    int cndof_pbc=0;
    int cndof_mid=0;
    int cbw_pbc=0;
    int cbw_mid=0;
    for (int system=0;system<=int(static_cast<cp2k_methods::cp2k_method_type>(cp2k_transport_params.method)==cp2k_methods::TRANSPORT);system++) {
        transport_parameters transport_params;
        transport_params.cp2k_scf_iter               = cp2k_transport_params.iscf;
        transport_params.n_occ                       = cp2k_transport_params.n_occ;
        transport_params.n_abscissae                 = cp2k_transport_params.n_abscissae;
        transport_params.n_kpoint                    = cp2k_transport_params.n_kpoint;
        transport_params.num_interval                = cp2k_transport_params.num_interval;
        transport_params.tasks_per_point             = cp2k_transport_params.tasks_per_energy_point;
        transport_params.gpus_per_point              = cp2k_transport_params.gpus_per_point;
        transport_params.colzero_threshold           = cp2k_transport_params.colzero_threshold;
        transport_params.eps_limit                   = cp2k_transport_params.eps_limit;
        transport_params.eps_limit_cc                = cp2k_transport_params.eps_limit_cc;
        transport_params.eps_decay                   = cp2k_transport_params.eps_decay;
        transport_params.eps_singularity_curvatures  = cp2k_transport_params.eps_singularity_curvatures;
        transport_params.eps_mu                      = cp2k_transport_params.eps_mu;
        transport_params.eps_eigval_degen            = cp2k_transport_params.eps_eigval_degen;
        transport_params.energy_interval             = cp2k_transport_params.energy_interval;
        transport_params.min_interval                = cp2k_transport_params.min_interval;
        transport_params.svd_cutoff                  = cp2k_transport_params.svd_cutoff;
        transport_params.n_points_beyn               = cp2k_transport_params.n_points_beyn;
        transport_params.tasks_per_integration_point = cp2k_transport_params.tasks_per_integration_point;
        transport_params.evoltfactor                 = cp2k_transport_params.evoltfactor;
        transport_params.NCRC_beyn                   = cp2k_transport_params.ncrc_beyn;
        transport_params.n_points_inv                = cp2k_transport_params.n_points_inv;
        transport_params.n_atoms                     = cp2k_transport_params.n_atoms;
        transport_params.pexsi_ordering              = cp2k_transport_params.ordering;
        transport_params.pexsi_row_ordering          = cp2k_transport_params.row_ordering;
        transport_params.pexsi_verbosity             = cp2k_transport_params.verbosity;
        transport_params.pexsi_np_symb_fact          = cp2k_transport_params.pexsi_np_symb_fact;
        transport_params.temperature                 = cp2k_transport_params.temperature*cp2k_transport_params.evoltfactor;
        transport_params.boltzmann_ev                = cp2k_transport_params.boltzmann/cp2k_transport_params.e_charge;
        transport_params.conduct_quant               = 2.0*cp2k_transport_params.e_charge*cp2k_transport_params.e_charge/(2.0*M_PI*cp2k_transport_params.h_bar);
        transport_params.extra_scf                   = cp2k_transport_params.extra_scf;
        transport_params.obc                         = cp2k_transport_params.obc_equilibrium || cp2k_transport_params.cutout[0] || cp2k_transport_params.cutout[1];
        transport_params.injection_method            = static_cast<injection_methods::injection_method_type>(cp2k_transport_params.injection_method);
        transport_params.lin_solver_method           = static_cast<lin_solver_methods::lin_solver_method_type>(cp2k_transport_params.linear_solver);
        transport_params.real_int_method             = static_cast<real_int_methods::real_int_method_type>(cp2k_transport_params.rlaxis_integration_method);
        transport_params.cp2k_method                 = static_cast<cp2k_methods::cp2k_method_type>(cp2k_transport_params.method);
        if (cp2k_transport_params.eps_fermi<=(numeric_limits<double>::epsilon)()) {
            transport_params.eps_fermi               = (numeric_limits<double>::epsilon)();
        } else {
            transport_params.eps_fermi               = cp2k_transport_params.eps_fermi;
        }
        if (cp2k_transport_params.n_rand_beyn>1.0 || cp2k_transport_params.n_rand_beyn<=0.0) {
            transport_params.fac_neigbeyn            = 1.0;
        } else {
            transport_params.fac_neigbeyn            = cp2k_transport_params.n_rand_beyn;
        }
        if (cp2k_transport_params.n_rand_cc_beyn>1.0 || cp2k_transport_params.n_rand_cc_beyn<=0.0) {
            transport_params.fac_neigbeyn_cc         = 1.0;
        } else {
            transport_params.fac_neigbeyn_cc         = cp2k_transport_params.n_rand_cc_beyn;
        }
        transport_params.negf_solver                 = false;
        if (cp2k_transport_params.qt_formalism==41) {
            transport_params.negf_solver             = true;
        }
 
        std::vector<contact_type> contactvec(cp2k_transport_params.num_contacts);
        int size_muvec=0;
        int stride = cp2k_transport_params.stride_contacts;
        for (uint i_c=0;i_c<contactvec.size();i_c++) if (cp2k_transport_params.contacts_data[4+stride*i_c]) size_muvec++;
        std::vector<double> muvec(size_muvec);
        int i_m=0;
        int i_n=0;
        for (uint i_c=0;i_c<contactvec.size();i_c++) {
            int i_t;
            if (cp2k_transport_params.contacts_data[4+stride*i_c]) {
                i_t=i_m++;
            } else {
                i_t=muvec.size()+i_n++;
            }
            contactvec[i_t].bandwidth = cp2k_transport_params.contacts_data[0+stride*i_c];
            contactvec[i_t].inj_sign  = cp2k_transport_params.contacts_data[3+stride*i_c];
            contactvec[i_t].natoms    = cp2k_transport_params.contacts_data[2+stride*i_c];
            contactvec[i_t].atomstart = cp2k_transport_params.contacts_data[1+stride*i_c];
            if (contactvec[i_t].atomstart<0) {
                if (contactvec[i_t].inj_sign==1) {
                    contactvec[i_t].atomstart = cp2k_transport_params.cutout[0];
                } else if (contactvec[i_t].inj_sign==-1) {
                    if (transport_params.cp2k_method==cp2k_methods::TRANSPORT) {
                        contactvec[i_t].atomstart = transport_params.n_atoms/2-contactvec[i_t].natoms;
                    } else {
                        contactvec[i_t].atomstart = transport_params.n_atoms-contactvec[i_t].natoms-cp2k_transport_params.cutout[1];
                    }
                }
            }
        }

        std::vector<int> natoms_start(mpi_size);
        MPI_Allgather(&atom_of_bf[S.first_row],1,MPI_INT,&natoms_start[0],1,MPI_INT,MPI_COMM_WORLD);
        natoms_start.push_back(transport_params.n_atoms);
        std::vector<int> natoms_local;
        for (int i=0;i<mpi_size;i++) {
            natoms_local.push_back(natoms_start[i+1]-natoms_start[i]);
        }
 
        std::vector<std::vector<int>> pairlist(natoms_local[rank]);
        for (int i=0;i<S.nrows_local;i++) {
            int a=atom_of_bf[S.first_row+i]-atom_of_bf[S.first_row];
            for (int e=S.rowptr_local[i]-1;e<S.rowptr_local[i+1]-1;e++) {
                pairlist[a].push_back(atom_of_bf[S.colind_local[e]-1]+1);
            }
        }
        int maxpair=0;
        for (uint a=0;a<pairlist.size();a++) {
            sort(pairlist[a].begin(),pairlist[a].end());
            pairlist[a].erase(unique(pairlist[a].begin(),pairlist[a].end()),pairlist[a].end());
            maxpair=max(maxpair,int(pairlist[a].size()));
        }
        MPI_Allreduce(MPI_IN_PLACE,&maxpair,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
 
        int *pairmatrix_local = new int[maxpair*pairlist.size()]();
        for (uint a=0;a<pairlist.size();a++) {
            for (uint i=0;i<pairlist[a].size();i++) {
                pairmatrix_local[a*maxpair+i]=pairlist[a][i];
            }
        }
        pairlist=std::vector<std::vector<int>>();
        for (int i=0;i<mpi_size;i++) {
            natoms_local[i]*=maxpair;
            natoms_start[i]*=maxpair;
        }
        int *pairmatrix_global=NULL;
        if (!rank) {
            pairmatrix_global = new int[maxpair*transport_params.n_atoms];
        }
        MPI_Gatherv(pairmatrix_local,natoms_local[rank],MPI_INT,pairmatrix_global,&natoms_local[0],&natoms_start[0],MPI_INT,0,MPI_COMM_WORLD);
        delete[] pairmatrix_local;
 
        std::vector<int> Bsizes;
        int num_tridiag_blocks = 1;
        if (!rank) {
            std::vector<bool> adjmat(transport_params.n_atoms*transport_params.n_atoms);
            for (int a=0;a<transport_params.n_atoms;a++) {
                for (int i=0;i<maxpair;i++) {
                    int b=pairmatrix_global[a*maxpair+i]-1;
                    if (b+1) {
                        adjmat[b*transport_params.n_atoms+a]=true;
                        adjmat[a*transport_params.n_atoms+b]=true;
                    }
                }
            }
            delete[] pairmatrix_global;
            pairmatrix_global = NULL;
 
            for (uint i=0;i<contactvec.size();i++) {
                if (contactvec[i].bandwidth<=0) {
                    contactvec[i].bandwidth=0;
                    int block_exists = 1;
                    while (block_exists) {
                        block_exists = 0;
                        for (int irow=contactvec[i].atomstart;irow<contactvec[i].atomstart+contactvec[i].natoms;irow++) {
                            int icol_start = contactvec[i].atomstart+contactvec[i].inj_sign*(contactvec[i].bandwidth+1)*contactvec[i].natoms;
                            int icol_end = contactvec[i].atomstart+contactvec[i].natoms+contactvec[i].inj_sign*(contactvec[i].bandwidth+1)*contactvec[i].natoms;
                            for (int icol=icol_start;icol<icol_end;icol++) {
                                if (adjmat[icol*transport_params.n_atoms+irow]) {
                                    block_exists = 1;
                                    contactvec[i].bandwidth++;
                                    break;
                                }
                            }
                            if (block_exists) break;
                        }
                    }
                }
                contactvec[i].sigma_natoms = contactvec[i].natoms*contactvec[i].bandwidth;
 
                int optimize_sigma_size=0;
                if (optimize_sigma_size) {
                    int tridiag_start = contactvec[i].atomstart;
                    if (contactvec[i].inj_sign==-1) tridiag_start = contactvec[i].atomstart-(contactvec[i].bandwidth-1)*contactvec[i].natoms;
                    int mincol = transport_params.n_atoms;
                    int maxcol = 0;
                    for (int irow=tridiag_start;irow<tridiag_start+contactvec[i].bandwidth*contactvec[i].natoms;irow++) {
                        int icol_start = tridiag_start+contactvec[i].inj_sign*contactvec[i].bandwidth*contactvec[i].natoms;
                        int icol_end = tridiag_start+contactvec[i].bandwidth*contactvec[i].natoms+contactvec[i].inj_sign*contactvec[i].bandwidth*contactvec[i].natoms;
                        for (int icol=icol_start;icol<icol_end;icol++) {
                            if (adjmat[icol*transport_params.n_atoms+irow]) {
                                if (icol < mincol) mincol = icol;
                                if (icol > maxcol) maxcol = icol;
                            }
                        }
                    }
                    contactvec[i].sigma_natoms = maxcol+1-(tridiag_start+contactvec[i].bandwidth*contactvec[i].natoms);
                    if (contactvec[i].inj_sign==-1) contactvec[i].sigma_natoms = tridiag_start-mincol;
                }
            }
 
            std::vector<int> tridiag_blocks_start;
            tridiag_blocks_start.push_back(cp2k_transport_params.cutout[0]);
            tridiag_blocks_start.push_back(contactvec[0].atomstart+contactvec[0].sigma_natoms);
            int tridiag_end = contactvec[1].atomstart+contactvec[1].natoms-contactvec[1].sigma_natoms;
            int maxcol = 0;
            while (tridiag_blocks_start[num_tridiag_blocks] < tridiag_end) {
                for (int irow=tridiag_blocks_start[num_tridiag_blocks-1];irow<tridiag_blocks_start[num_tridiag_blocks];irow++) {
                    int colend = transport_params.n_atoms-cp2k_transport_params.cutout[1];
                    if (irow<contactvec[0].bandwidth*contactvec[0].natoms) {
                        colend = contactvec[0].atomstart+2*contactvec[0].bandwidth*contactvec[0].natoms;
                    }
                    for (int icol=irow;icol<colend;icol++) {
                        if (adjmat[icol*transport_params.n_atoms+irow]) {
                            if (icol > maxcol) {
                                maxcol = icol;
                            }
                        }
                    }
                }
                num_tridiag_blocks++;
                tridiag_blocks_start.push_back(maxcol+1);
            }
            tridiag_blocks_start[0] = 0;
            tridiag_blocks_start[num_tridiag_blocks] = transport_params.n_atoms;
            for (int i=0;i<num_tridiag_blocks;i++) {
                Bsizes.push_back(tridiag_blocks_start[i+1]-tridiag_blocks_start[i]);
            }
        }
 
        for (uint i=0;i<contactvec.size();i++) {
             if (!rank) cout << "Contact BW " << contactvec[i].bandwidth << endl;
             MPI_Bcast(&contactvec[i].bandwidth,1,MPI_INT,0,MPI_COMM_WORLD);
             MPI_Bcast(&contactvec[i].sigma_natoms,1,MPI_INT,0,MPI_COMM_WORLD);
        }
 
        MPI_Bcast(&num_tridiag_blocks,1,MPI_INT,0,MPI_COMM_WORLD);
        Bsizes.resize(num_tridiag_blocks);
        MPI_Bcast(&Bsizes[0],num_tridiag_blocks,MPI_INT,0,MPI_COMM_WORLD);
 
        Bsizes[0]-=cp2k_transport_params.cutout[0];
        Bsizes[Bsizes.size()-1]-=cp2k_transport_params.cutout[1];
 
        for (uint i=0;i<contactvec.size();i++) {
            int atom_start = contactvec[i].atomstart;
            int atom_stop  = atom_start + contactvec[i].natoms;
            contactvec[i].atomstart = atom_start-cp2k_transport_params.cutout[0];
            contactvec[i].start     = std::accumulate(&cp2k_transport_params.nsgf[cp2k_transport_params.cutout[0]],&cp2k_transport_params.nsgf[atom_start],0);
            contactvec[i].start_bs  = std::accumulate(&cp2k_transport_params.nsgf[0],&cp2k_transport_params.nsgf[atom_start],0);
            contactvec[i].ndof      = std::accumulate(&cp2k_transport_params.nsgf[atom_start],&cp2k_transport_params.nsgf[atom_stop],0);
            contactvec[i].n_ele     = std::accumulate(&cp2k_transport_params.zeff[atom_start],&cp2k_transport_params.zeff[atom_stop],0.0);
        }
 
        transport_params.cutl = std::accumulate(&cp2k_transport_params.nsgf[0],&cp2k_transport_params.nsgf[cp2k_transport_params.cutout[0]],0);
        transport_params.cutr = std::accumulate(&cp2k_transport_params.nsgf[cp2k_transport_params.n_atoms-cp2k_transport_params.cutout[1]],&cp2k_transport_params.nsgf[cp2k_transport_params.n_atoms],0);
 
        transport_params.n_atoms = transport_params.n_atoms-cp2k_transport_params.cutout[0]-cp2k_transport_params.cutout[1];
        std::vector<int> orb_per_atom(1,0);
        for (int i=0;i<transport_params.n_atoms;i++) {
            orb_per_atom.push_back(orb_per_atom[i]+cp2k_transport_params.nsgf[cp2k_transport_params.cutout[0]+i]);
        }

        if (transport_params.cp2k_method==cp2k_methods::TRANSPORT) {
            transport_params.obc=true;
            cndof_mid=contactvec[1].ndof;
            cbw_mid=contactvec[1].bandwidth;
            if (!system) {
                transport_params.cutl=0;
                transport_params.cutr=S.nrows_total/2;
            } else {
                transport_params.cutl=S.nrows_total/2;
                transport_params.cutr=0;
                transport_params.cp2k_scf_iter*=-1;
                if (contactvec[0].inj_sign==+1) {
                    contactvec[0].inj_sign=-1;
                    contactvec[0].start   =S.nrows_total/2-contactvec[0].ndof;
                    contactvec[0].start_bs=S.nrows_total-contactvec[0].ndof;
                    contactvec[1].inj_sign=+1;
                    contactvec[1].start   =0;
                    contactvec[1].start_bs=S.nrows_total/2;
                } else if (contactvec[0].inj_sign==-1) {
                    contactvec[0].inj_sign=+1;
                    contactvec[0].start   =0;
                    contactvec[0].start_bs=S.nrows_total/2;
                    contactvec[1].inj_sign=-1;
                    contactvec[1].start   =S.nrows_total/2-contactvec[1].ndof;
                    contactvec[1].start_bs=S.nrows_total-contactvec[1].ndof;
                } else throw std::exception();
                if (contactvec.size()>muvec.size()) {
                    contactvec[2].start_bs+=S.nrows_total/2;
                }
            }
        }
        if (!transport_params.extra_scf && transport_params.obc) {
            cndof_pbc=contactvec[0].ndof;
            cbw_pbc=contactvec[0].bandwidth;
        }

        //write_matrix(Overlap,KohnSham,wr_cutblocksize,wr_bw,wr_ndof);
        //if (diagscalapack(Overlap,KohnSham,transport_params)) throw std::exception();
        if (transport_params.cp2k_method==cp2k_methods::LOCAL_SCF) {
#ifdef HAVE_OMEN_POISSON
            if (semiselfconsistent(S,KS,P,PImag,muvec,contactvec,Bsizes,orb_per_atom,transport_params)) throw std::exception();
#endif
        } else {
            Energyvector energyvector;
            if (energyvector.Execute(S,KS,P,PImag,muvec,contactvec,Bsizes,orb_per_atom,NULL,NULL,transport_params)) throw std::exception();
        }

    }

     if (cndof_pbc) {
        int size_tot=S.nrows_total;
        int offset=0;
        int cndof=cndof_pbc;
        for (int i=0;i<cbw_pbc;i++) {
            for (int ii=0;ii<=i;ii++) {
                copy_to_cp2k_csr(*P,0.5,\
                        (i+1+offset)*cndof,\
                        offset*cndof,\
                        ii*cndof,\
                        size_tot+(ii-(i+1))*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_to_cp2k_csr(*P,0.5,\
                        size_tot-(offset+1)*cndof,\
                        size_tot-(i+2+offset)*cndof,\
                        ii*cndof,\
                        size_tot+(ii-(i+1))*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                copy_to_cp2k_csr(*P,0.5,\
                        offset*cndof,\
                        (i+1+offset)*cndof,\
                        size_tot+(ii-(i+1))*cndof,\
                        ii*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_to_cp2k_csr(*P,0.5,\
                        size_tot-(i+2+offset)*cndof,\
                        size_tot-(offset+1)*cndof,\
                        size_tot+(ii-(i+1))*cndof,\
                        ii*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
            }
        }
    }
    if (cndof_mid) {
        int size_tot=S.nrows_total;
        int mid=size_tot/2;
        int cndof=cndof_mid;
        for (int i=0;i<cbw_mid;i++) {
            for (int ii=0;ii<=i;ii++) {
                copy_to_cp2k_csr(*P,0.5,\
                        mid+(-i-2)*cndof,\
                        mid+(-1)*cndof,\
                        mid+(-i-1+ii)*cndof,\
                        mid+(0+ii)*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_to_cp2k_csr(*P,0.5,\
                        mid+0*cndof,\
                        mid+(i+1)*cndof,\
                        mid+(-i-1+ii)*cndof,\
                        mid+(0+ii)*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                copy_to_cp2k_csr(*P,0.5,\
                        mid+(-1)*cndof,\
                        mid+(-i-2)*cndof,\
                        mid+(0+ii)*cndof,\
                        mid+(-i-1+ii)*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_to_cp2k_csr(*P,0.5,\
                        mid+(i+1)*cndof,\
                        mid+0*cndof,\
                        mid+(0+ii)*cndof,\
                        mid+(-i-1+ii)*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
            }
        }
    }

    if (!cp2k_transport_params.extra_scf) {
        std::vector<double> mulli(cp2k_transport_params.n_atoms,0.0);
        for (int i=0;i<S.nrows_local;i++) {
            for (int e=S.rowptr_local[i]-1;e<S.rowptr_local[i+1]-1;e++) {
                mulli[atom_of_bf[i+S.first_row]]+=S.nzvals_local[e]*P->nzvals_local[e];
            }
        }
        MPI_Allreduce(MPI_IN_PLACE,&mulli[0],cp2k_transport_params.n_atoms,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if (!rank) {
            stringstream mysstream;
            mysstream << "Mulliken_" << cp2k_transport_params.iscf;
            ofstream myfile(mysstream.str().c_str());
            myfile.precision(8);
            for (int i=0;i<cp2k_transport_params.n_atoms;i++) myfile<<mulli[i]<<"\n";
            myfile.close();
        }
    }
}

