#include "libcp2k.h"
#include "CSR.H"
#ifdef HAVE_OMEN_POISSON
#include "SemiSelfConsistent.H"
#endif
#include "EnergyVector.H"
//#include "DiagScaLapack.H"
//#include "WriteMatrix.H"

void write_cp2k_csr(cp2k_csr_interop_type& cp2kCSRmat,const char* filename)
{
    ofstream myfile;
    myfile.open(filename);
    myfile.precision(14);
    for(int i=0;i<cp2kCSRmat.nrows_local;i++){
        for(int e=cp2kCSRmat.rowptr_local[i]-1;e<cp2kCSRmat.rowptr_local[i+1]-1;e++){
            myfile<<i+1+cp2kCSRmat.first_row<<" "<<cp2kCSRmat.colind_local[e]<<" "<<cp2kCSRmat.nzvals_local[e]<<"\n";
        }
    }
    myfile.close();                                                   
}

void read_nzvals_bin(cp2k_csr_interop_type& cp2kCSRmat,const char* filename,MPI_Comm cp2k_comm)
{
    MPI_File file;
    MPI_Status status;
    MPI_File_open(cp2k_comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&file);
    int rank,mpi_size;
    MPI_Comm_size(cp2k_comm,&mpi_size);
    MPI_Comm_rank(cp2k_comm,&rank);

    int *dist = new int[mpi_size];
    int *disp = new int[mpi_size];
    MPI_Allgather(&cp2kCSRmat.nze_local,1,MPI_INT,dist,1,MPI_INT,cp2k_comm);
    disp[0]=0;
    for (int i=0;i<mpi_size-1;i++) {
        disp[i+1]=disp[i]+dist[i];
    }

    MPI_Offset offset = disp[rank]*sizeof(double);

    delete[] dist;
    delete[] disp;

    MPI_File_seek(file,offset,MPI_SEEK_SET);
    MPI_File_read(file,cp2kCSRmat.nzvals_local,cp2kCSRmat.nze_local,MPI_DOUBLE,&status);

    MPI_File_close(&file);
}

void add_full_to_scaled_cp2k_csr(cp2k_csr_interop_type& cp2kCSRmat,double* Pf,int start_i_to,int start_j_to,int length_i,int length_j,double a, double b)
{
    for(int r=0;r<cp2kCSRmat.nrows_local;r++){
        int i=r+cp2kCSRmat.first_row-start_i_to;
        if(i>=0 && i<length_i){
            for(int e=cp2kCSRmat.rowptr_local[r]-1;e<cp2kCSRmat.rowptr_local[r+1]-1;e++){
                int j=cp2kCSRmat.colind_local[e]-1-start_j_to;
                if(j>=0 && j<length_j){
                    cp2kCSRmat.nzvals_local[e]=a*Pf[i+j*length_i]+b*cp2kCSRmat.nzvals_local[e];
                }
            }
        }
    }
}

void add_from_to_scaled_cp2k_csr(cp2k_csr_interop_type& cp2kCSRmat,double a, double b,int start_i_fr,int start_j_fr,int start_i_to,int start_j_to,int length_i,int length_j,MPI_Comm cp2k_comm)
{
    TCSR<double> *Pmc = new TCSR<double>(cp2kCSRmat,start_i_fr,length_i,start_j_fr,length_j);
    TCSR<double> *Pm  = new TCSR<double>(Pmc,cp2k_comm);
    delete Pmc;
    double *Pf = new double[length_i*length_j];
    Pm->sparse_to_full(Pf,length_i,length_j);
    delete Pm;
    add_full_to_scaled_cp2k_csr(cp2kCSRmat,Pf,start_i_to,start_j_to,length_i,length_j,a,b);
    delete[] Pf;
}

void write_scaled_cp2k_csr_bin(cp2k_csr_interop_type& cp2kCSRmat,const char* filename,double factor,MPI_Comm cp2k_comm)
{
    MPI_File file;
    MPI_Status status;
    MPI_File_open(cp2k_comm,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&file);
    int rank,mpi_size;
    MPI_Comm_size(cp2k_comm,&mpi_size);
    MPI_Comm_rank(cp2k_comm,&rank);

    int *dist = new int[mpi_size];
    int *disp = new int[mpi_size];
    MPI_Allgather(&cp2kCSRmat.nze_local,1,MPI_INT,dist,1,MPI_INT,cp2k_comm);
    disp[0]=0;
    for (int i=0;i<mpi_size-1;i++) {
        disp[i+1]=disp[i]+dist[i];
    }

    MPI_Offset offset = 0;
    if (rank) offset=(4*disp[rank]+3)*sizeof(double);
    delete[] dist;
    delete[] disp;

    MPI_File_seek(file,offset,MPI_SEEK_SET);

    if (!rank) {
        double head_1=(double) cp2kCSRmat.nrows_total;
        double head_2=(double) cp2kCSRmat.nze_total;
        double head_3=(double) 1;
        MPI_File_write(file,&head_1,1,MPI_DOUBLE,&status);
        MPI_File_write(file,&head_2,1,MPI_DOUBLE,&status);
        MPI_File_write(file,&head_3,1,MPI_DOUBLE,&status);
    }

    for (int i=0;i<cp2kCSRmat.nrows_local;i++) {
        for (int e=cp2kCSRmat.rowptr_local[i]-1;e<cp2kCSRmat.rowptr_local[i+1]-1;e++) {
            double i_val=(double) i+cp2kCSRmat.first_row+1;
            double j_val=(double) cp2kCSRmat.colind_local[e];
            double r_val=factor*cp2kCSRmat.nzvals_local[e];
            double m_val=0.0;
            MPI_File_write(file,&i_val,1,MPI_DOUBLE,&status);
            MPI_File_write(file,&j_val,1,MPI_DOUBLE,&status);
            MPI_File_write(file,&r_val,1,MPI_DOUBLE,&status);
            MPI_File_write(file,&m_val,1,MPI_DOUBLE,&status);
        }
    }

    MPI_File_close(&file);
}

void write_scaled_cp2k_csr_bin_remove_pbc(cp2k_csr_interop_type& cp2kCSRmat,const char* filename,double factor,int bw,int ndof,MPI_Comm cp2k_comm)
{
    MPI_File file;
    MPI_Status status;
    MPI_File_open(cp2k_comm,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&file);
    int rank,mpi_size;
    MPI_Comm_size(cp2k_comm,&mpi_size);
    MPI_Comm_rank(cp2k_comm,&rank);

    int *remove_element = new int[max(cp2kCSRmat.nze_local,1)]();
    int n_removed_elements = 0;

    for (int i=0;i<max(0,min(cp2kCSRmat.nrows_local,bw*ndof-cp2kCSRmat.first_row));i++) {
        for (int e=cp2kCSRmat.rowptr_local[i]-1;e<cp2kCSRmat.rowptr_local[i+1]-1;e++) {
            if ( (cp2kCSRmat.colind_local[e]-1)/ndof - (i+cp2kCSRmat.first_row)/ndof > bw ) {
                remove_element[e] = 1;
                n_removed_elements++;
            }
        }
    }

    for (int i=max(0,min(cp2kCSRmat.nrows_local,cp2kCSRmat.nrows_total-bw*ndof-cp2kCSRmat.first_row));i<cp2kCSRmat.nrows_local;i++) {
        for (int e=cp2kCSRmat.rowptr_local[i]-1;e<cp2kCSRmat.rowptr_local[i+1]-1;e++) {
            if ( (cp2kCSRmat.colind_local[e]-1)/ndof - (i+cp2kCSRmat.first_row)/ndof < -bw ) {
                remove_element[e] = 1;
                n_removed_elements++;
            }
        }
    }

    int remaining_elements = cp2kCSRmat.nze_local - n_removed_elements;

    int *dist = new int[mpi_size];
    int *disp = new int[mpi_size+1];
    MPI_Allgather(&remaining_elements,1,MPI_INT,dist,1,MPI_INT,cp2k_comm);
    disp[0]=0;
    for (int i=0;i<mpi_size;i++) {
        disp[i+1]=disp[i]+dist[i];
    }
    delete[] dist;

    MPI_Offset offset = 0;
    if (rank) offset=(4*disp[rank]+3)*sizeof(double);

    MPI_File_seek(file,offset,MPI_SEEK_SET);

    if (!rank) {
        double head_1=(double) cp2kCSRmat.nrows_total;
        double head_2=(double) disp[mpi_size];
        double head_3=(double) 1;
        MPI_File_write(file,&head_1,1,MPI_DOUBLE,&status);
        MPI_File_write(file,&head_2,1,MPI_DOUBLE,&status);
        MPI_File_write(file,&head_3,1,MPI_DOUBLE,&status);
    }

    delete[] disp;

    for (int i=0;i<cp2kCSRmat.nrows_local;i++) {
        for (int e=cp2kCSRmat.rowptr_local[i]-1;e<cp2kCSRmat.rowptr_local[i+1]-1;e++) {
            if (!remove_element[e]) {
                double i_val=(double) i+cp2kCSRmat.first_row+1;
                double j_val=(double) cp2kCSRmat.colind_local[e];
                double r_val=factor*cp2kCSRmat.nzvals_local[e];
                double m_val=0.0;
                MPI_File_write(file,&i_val,1,MPI_DOUBLE,&status);
                MPI_File_write(file,&j_val,1,MPI_DOUBLE,&status);
                MPI_File_write(file,&r_val,1,MPI_DOUBLE,&status);
                MPI_File_write(file,&m_val,1,MPI_DOUBLE,&status);
            }
        }
    }

    delete[] remove_element;
    MPI_File_close(&file);
}

/*!  
 *   \brief Takes the overlap (S) and Kohn-Sham (KS) matrices as input and returns a density matrix (P).
 *          This function acts as the gate to the CP2K's world. 
 */
#ifdef HAVE_PIMAG
void c_scf_method(cp2k_transport_parameters cp2k_transport_params, cp2k_csr_interop_type S, cp2k_csr_interop_type KS, cp2k_csr_interop_type * P, cp2k_csr_interop_type * PImag)
{
#else
void c_scf_method(cp2k_transport_parameters cp2k_transport_params, cp2k_csr_interop_type S, cp2k_csr_interop_type KS, cp2k_csr_interop_type * P)
{
    cp2k_csr_interop_type * PImag = NULL;
#endif
#ifdef HAVE_SPLITSOLVE
    char gpu_string[255];
    set_gpu(0,gpu_string);
#endif

    int rank,mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (!rank) cout << "Starting Transport" << endl;

    c_dscal(P->nze_local,0.5,P->nzvals_local,1);

/*
    write_scaled_cp2k_csr_bin( S,"S_4.bin",1.0,MPI_COMM_WORLD);
    write_scaled_cp2k_csr_bin(KS,"H_4.bin",cp2k_transport_params.evoltfactor,MPI_COMM_WORLD);
    return;
*/

    double dens_mixing=cp2k_transport_params.dens_mixing;
    double *P_save=NULL;
    if (dens_mixing<1.0 || dens_mixing>0.0) {
        P_save = new double[P->nze_local];
        c_dcopy(P->nze_local,P->nzvals_local,1,P_save,1);
    }

    int cndof_pbc=0;
    int cndof_mid=0;
    int cbw_pbc=0;
    int cbw_mid=0;
    for (int system=0;system<=int(static_cast<cp2k_methods::cp2k_method_type>(cp2k_transport_params.method)==cp2k_methods::TRANSPORT);system++) {
        transport_parameters transport_params;
        transport_params.cp2k_scf_iter               = (1-2*system)*cp2k_transport_params.iscf;
        transport_params.n_occ                       = cp2k_transport_params.n_occ;
        transport_params.n_abscissae                 = cp2k_transport_params.num_pole;
        transport_params.n_kpoint                    = cp2k_transport_params.n_kpoint;
        transport_params.num_interval                = cp2k_transport_params.num_interval;
        transport_params.tasks_per_point             = cp2k_transport_params.tasks_per_energy_point;
        transport_params.tasks_per_point_cc          = cp2k_transport_params.tasks_per_pole;
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
        transport_params.pexsi_ordering              = cp2k_transport_params.ordering;
        transport_params.pexsi_row_ordering          = cp2k_transport_params.row_ordering;
        transport_params.pexsi_verbosity             = cp2k_transport_params.verbosity;
        transport_params.pexsi_np_symb_fact          = cp2k_transport_params.pexsi_np_symb_fact;
        transport_params.temperature                 = cp2k_transport_params.temperature*cp2k_transport_params.evoltfactor;
        transport_params.boltzmann_ev                = cp2k_transport_params.boltzmann/cp2k_transport_params.e_charge;
        transport_params.conduct_quant               = 2.0*cp2k_transport_params.e_charge*cp2k_transport_params.e_charge/(2.0*M_PI*cp2k_transport_params.h_bar);
        transport_params.extra_scf                   = cp2k_transport_params.extra_scf;
        transport_params.injection_method            = static_cast<injection_methods::injection_method_type>(cp2k_transport_params.injection_method);
        transport_params.lin_solver_method           = static_cast<lin_solver_methods::lin_solver_method_type>(cp2k_transport_params.linear_solver);
        transport_params.inv_solver_method           = static_cast<inv_solver_methods::inv_solver_method_type>(cp2k_transport_params.matrixinv_method);
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
        transport_params.update_fermi                = true;
        transport_params.get_fermi_neutral           = false;
        if (cp2k_transport_params.transport_neutral==52) {
            transport_params.get_fermi_neutral       = true;
        }

        int cutout[2]={0,0};
        cutout[0]=cp2k_transport_params.cutout[0];
        cutout[1]=cp2k_transport_params.cutout[1];
        if (transport_params.cp2k_method==cp2k_methods::TRANSPORT) {
            transport_params.extra_scf = false;
            cutout[system]=0;
            cutout[1-system]=cp2k_transport_params.n_atoms/2;
        }
        transport_params.obc                         = cp2k_transport_params.obc_equilibrium || cutout[0] || cutout[1];

        std::vector<contact_type> contactvec(cp2k_transport_params.num_contacts);
        int num_p=0;
        int num_m=0;
        int stride = cp2k_transport_params.stride_contacts;
        for (uint i_c=0;i_c<contactvec.size();i_c++) {
            if (cp2k_transport_params.contacts_data[4+stride*i_c]) {
                if (cp2k_transport_params.contacts_data[3+stride*i_c]==+1) {
                    num_p++;
                } else {
                    num_m++;
                }
            }
        }
        std::vector<double> muvec(num_p+num_m);
        int i_p=0;
        int i_m=0;
        int i_n=0;
        for (uint i_c=0;i_c<contactvec.size();i_c++) {
            int i_t;
            if (cp2k_transport_params.contacts_data[4+stride*i_c]) {
                if (cp2k_transport_params.contacts_data[3+stride*i_c]==+1) {
                    i_t=i_p++;
                } else {
                    i_t=num_p+i_m++;
                }
            } else {
                i_t=muvec.size()+i_n++;
            }
            contactvec[i_t].bandwidth = cp2k_transport_params.contacts_data[0+stride*i_c];
            contactvec[i_t].inj_sign  = cp2k_transport_params.contacts_data[3+stride*i_c];
            contactvec[i_t].natoms    = cp2k_transport_params.contacts_data[2+stride*i_c];
            contactvec[i_t].atomstart = cp2k_transport_params.contacts_data[1+stride*i_c];
            if (contactvec[i_t].atomstart<0) {
                if (contactvec[i_t].inj_sign==+1) {
                    contactvec[i_t].atomstart = cutout[0];
                } else if (contactvec[i_t].inj_sign==-1) {
                    contactvec[i_t].atomstart = cp2k_transport_params.n_atoms-contactvec[i_t].natoms-cutout[1];
                }
            }
        }
        if (transport_params.cp2k_method==cp2k_methods::TRANSPORT) {
            if (system) {
                swap(contactvec[0].natoms,contactvec[1].natoms);
                swap(contactvec[0].bandwidth,contactvec[1].bandwidth);
            }
            contactvec[0].atomstart = cutout[0];
            contactvec[1].atomstart = cp2k_transport_params.n_atoms-contactvec[1].natoms-cutout[1];
            for (uint i=muvec.size();i<contactvec.size();i++) {
                contactvec[i].atomstart+=cutout[0];
            }
        }

        if (transport_params.cp2k_method==cp2k_methods::TRANSMISSION) {
            if (!transport_params.extra_scf && !transport_params.obc) {
                contactvec.resize(1);
                muvec.resize(1);
                contactvec[0].bandwidth = cp2k_transport_params.contacts_data[0];
                contactvec[0].inj_sign  = cp2k_transport_params.contacts_data[3];
                contactvec[0].natoms    = cp2k_transport_params.contacts_data[2];
                contactvec[0].atomstart = cp2k_transport_params.contacts_data[1];
                if (contactvec[0].atomstart<0) {
                    if (contactvec[0].inj_sign==+1) {
                        contactvec[0].atomstart = cutout[0];
                    } else if (contactvec[0].inj_sign==-1) {
                        contactvec[0].atomstart = cp2k_transport_params.n_atoms-contactvec[0].natoms-cutout[1];
                    }
                }
            } else {
                cutout[0] = contactvec[0].atomstart;
                cutout[1] = cp2k_transport_params.n_atoms-contactvec[1].atomstart-contactvec[1].natoms;
            }
        }

        std::vector<int> atom_of_bf;
        for (int a=0;a<cp2k_transport_params.n_atoms;a++) {
            for (int i=0;i<cp2k_transport_params.nsgf[a];i++) {
                atom_of_bf.push_back(a);
            }
        }
        atom_of_bf.push_back(cp2k_transport_params.n_atoms);
 
        std::vector<int> natoms_start(mpi_size);
        MPI_Allgather(&atom_of_bf[S.first_row],1,MPI_INT,&natoms_start[0],1,MPI_INT,MPI_COMM_WORLD);
        natoms_start.push_back(cp2k_transport_params.n_atoms);
        std::vector<int> natoms_local;
        for (int i=0;i<mpi_size;i++) {
            natoms_local.push_back(natoms_start[i+1]-natoms_start[i]);
        }
 
        std::vector<std::vector<int>> pairlist(natoms_local[rank]);
        int i_bf=0;
        int maxpair=0;
        for (uint a=0;a<pairlist.size();a++) {
            for (int i=0;i<cp2k_transport_params.nsgf[natoms_start[rank]+a];i++) {
                for (int e=S.rowptr_local[i_bf]-1;e<S.rowptr_local[i_bf+1]-1;e++) {
                    pairlist[a].push_back(atom_of_bf[S.colind_local[e]-1]+1);
                }
                i_bf++;
                pairlist[a].erase(unique(pairlist[a].begin(),pairlist[a].end()),pairlist[a].end());
            }
            sort(pairlist[a].begin(),pairlist[a].end());
            pairlist[a].erase(unique(pairlist[a].begin(),pairlist[a].end()),pairlist[a].end());
            maxpair=max(maxpair,int(pairlist[a].size()));
        }
        if (i_bf!=S.nrows_local) throw std::exception();
        MPI_Allreduce(MPI_IN_PLACE,&maxpair,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
        atom_of_bf=std::vector<int>();
 
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
            pairmatrix_global = new int[maxpair*cp2k_transport_params.n_atoms];
        }
        MPI_Gatherv(pairmatrix_local,natoms_local[rank],MPI_INT,pairmatrix_global,&natoms_local[0],&natoms_start[0],MPI_INT,0,MPI_COMM_WORLD);
        delete[] pairmatrix_local;
 
        std::vector<int> Bsizes;
        int num_tridiag_blocks = 1;
        if (!rank) {
            std::vector<bool> adjmat(cp2k_transport_params.n_atoms*cp2k_transport_params.n_atoms);
            for (int a=0;a<cp2k_transport_params.n_atoms;a++) {
                for (int i=0;i<maxpair;i++) {
                    int b=pairmatrix_global[a*maxpair+i]-1;
                    if (b+1) {
                        adjmat[b*cp2k_transport_params.n_atoms+a]=true;
                        adjmat[a*cp2k_transport_params.n_atoms+b]=true;
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
                                if (adjmat[icol*cp2k_transport_params.n_atoms+irow]) {
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
                    int mincol = cp2k_transport_params.n_atoms;
                    int maxcol = 0;
                    for (int irow=tridiag_start;irow<tridiag_start+contactvec[i].bandwidth*contactvec[i].natoms;irow++) {
                        int icol_start = tridiag_start+contactvec[i].inj_sign*contactvec[i].bandwidth*contactvec[i].natoms;
                        int icol_end = tridiag_start+contactvec[i].bandwidth*contactvec[i].natoms+contactvec[i].inj_sign*contactvec[i].bandwidth*contactvec[i].natoms;
                        for (int icol=icol_start;icol<icol_end;icol++) {
                            if (adjmat[icol*cp2k_transport_params.n_atoms+irow]) {
                                if (icol < mincol) mincol = icol;
                                if (icol > maxcol) maxcol = icol;
                            }
                        }
                    }
                    contactvec[i].sigma_natoms = maxcol+1-(tridiag_start+contactvec[i].bandwidth*contactvec[i].natoms);
                    if (contactvec[i].inj_sign==-1) contactvec[i].sigma_natoms = tridiag_start-mincol;
                }
            }

            if (!(transport_params.cp2k_method==cp2k_methods::TRANSMISSION && !transport_params.extra_scf && !transport_params.obc)) {
                std::vector<int> tridiag_blocks_start;
                tridiag_blocks_start.push_back(cutout[0]);
                tridiag_blocks_start.push_back(contactvec[0].atomstart+contactvec[0].sigma_natoms);
                int tridiag_end = contactvec[1].atomstart+contactvec[1].natoms-contactvec[1].sigma_natoms;
                int maxcol = 0;
                while (tridiag_blocks_start[num_tridiag_blocks] < tridiag_end) {
                    for (int irow=tridiag_blocks_start[num_tridiag_blocks-1];irow<tridiag_blocks_start[num_tridiag_blocks];irow++) {
                        int colend = cp2k_transport_params.n_atoms-cutout[1];
                        if (irow<contactvec[0].bandwidth*contactvec[0].natoms) {
                            colend = contactvec[0].atomstart+2*contactvec[0].bandwidth*contactvec[0].natoms;
                        }
                        for (int icol=irow;icol<colend;icol++) {
                            if (adjmat[icol*cp2k_transport_params.n_atoms+irow]) {
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
                tridiag_blocks_start[num_tridiag_blocks] = cp2k_transport_params.n_atoms;
                for (int i=0;i<num_tridiag_blocks;i++) {
                    Bsizes.push_back(tridiag_blocks_start[i+1]-tridiag_blocks_start[i]);
                }
 
                Bsizes[0]-=cutout[0];
                Bsizes[Bsizes.size()-1]-=cutout[1];
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
 
        for (uint i=0;i<contactvec.size();i++) {
            int atom_start = contactvec[i].atomstart;
            int atom_stop  = atom_start + contactvec[i].natoms;
            contactvec[i].atomstart = atom_start-cutout[0];
            contactvec[i].start     = std::accumulate(cp2k_transport_params.nsgf+cutout[0], cp2k_transport_params.nsgf+atom_start,0);
            contactvec[i].start_bs  = std::accumulate(cp2k_transport_params.nsgf,           cp2k_transport_params.nsgf+atom_start,0);
            contactvec[i].ndof      = std::accumulate(cp2k_transport_params.nsgf+atom_start,cp2k_transport_params.nsgf+atom_stop ,0);
            contactvec[i].n_ele     = std::accumulate(cp2k_transport_params.zeff+atom_start,cp2k_transport_params.zeff+atom_stop ,0.0);
        }
 
        transport_params.cutl = std::accumulate(cp2k_transport_params.nsgf,                                        cp2k_transport_params.nsgf+cutout[0]                    ,0);
        transport_params.cutr = std::accumulate(cp2k_transport_params.nsgf+cp2k_transport_params.n_atoms-cutout[1],cp2k_transport_params.nsgf+cp2k_transport_params.n_atoms,0);
 
        if (!transport_params.extra_scf && transport_params.obc && !system) {
            cndof_pbc=contactvec[0].ndof;
            cbw_pbc=contactvec[0].bandwidth;
        }
        if (system) {
            cndof_mid=contactvec[0].ndof;
            cbw_mid=contactvec[0].bandwidth;
        }

        std::vector<int> orb_per_atom(1,0);
        for (int i=0;i<cp2k_transport_params.n_atoms-cutout[0]-cutout[1];i++) {
            orb_per_atom.push_back(orb_per_atom[i]+cp2k_transport_params.nsgf[cutout[0]+i]);
        }

        if (transport_params.cp2k_method==cp2k_methods::WRITE_OUT) {
            write_scaled_cp2k_csr_bin_remove_pbc( S,"S_4.bin",1.0,contactvec[0].bandwidth,contactvec[0].ndof,MPI_COMM_WORLD);
            write_scaled_cp2k_csr_bin_remove_pbc(KS,"H_4.bin",cp2k_transport_params.evoltfactor,contactvec[0].bandwidth,contactvec[0].ndof,MPI_COMM_WORLD);
            return;
        }
        if (transport_params.cp2k_method==cp2k_methods::LOCAL_SCF) {
#ifdef HAVE_OMEN_POISSON
            transport_params.update_fermi=false;
            if (semiselfconsistent(S,KS,P,PImag,muvec,contactvec,Bsizes,orb_per_atom,dens_mixing,transport_params)) throw std::exception();
#endif
        } else {
            if (transport_params.get_fermi_neutral) {
                Energyvector energyvector;
                if (energyvector.Execute(S,KS,P,PImag,muvec,contactvec,Bsizes,orb_per_atom,NULL,transport_params)) throw std::exception();
                transport_params.update_fermi=false;
                transport_params.get_fermi_neutral=false;
            }
            Energyvector energyvector;
            if (energyvector.Execute(S,KS,P,PImag,muvec,contactvec,Bsizes,orb_per_atom,NULL,transport_params)) throw std::exception();
        }

    }

    int overwrite_first_last_diag_block=0;
    if (cndof_pbc) {
        int size_tot=S.nrows_total;
        int offset=0;
        int cndof=cndof_pbc;
        for (int i=0;i<cbw_pbc;i++) {
            for (int ii=0;ii<=i;ii++) {
                add_from_to_scaled_cp2k_csr(*P,0.5,0.0,\
                        (i+1+offset)*cndof,\
                        offset*cndof,\
                        ii*cndof,\
                        size_tot+(ii-(i+1))*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_from_to_scaled_cp2k_csr(*P,0.5,1.0,\
                        size_tot-(offset+1)*cndof,\
                        size_tot-(i+2+offset)*cndof,\
                        ii*cndof,\
                        size_tot+(ii-(i+1))*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_from_to_scaled_cp2k_csr(*P,0.5,0.0,\
                        offset*cndof,\
                        (i+1+offset)*cndof,\
                        size_tot+(ii-(i+1))*cndof,\
                        ii*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_from_to_scaled_cp2k_csr(*P,0.5,1.0,\
                        size_tot-(i+2+offset)*cndof,\
                        size_tot-(offset+1)*cndof,\
                        size_tot+(ii-(i+1))*cndof,\
                        ii*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
            }
        }
        if (overwrite_first_last_diag_block) {
            add_from_to_scaled_cp2k_csr(*P,1.0,0.0,\
                    size_tot-(offset+2)*cndof,\
                    size_tot-(offset+2)*cndof,\
                    size_tot-(offset+1)*cndof,\
                    size_tot-(offset+1)*cndof,\
                    cndof,cndof,MPI_COMM_WORLD);
            add_from_to_scaled_cp2k_csr(*P,1.0,0.0,\
                    (offset+1)*cndof,\
                    (offset+1)*cndof,\
                    offset*cndof,\
                    offset*cndof,\
                    cndof,cndof,MPI_COMM_WORLD);
        }
    }
    if (cndof_mid) {
        int size_tot=S.nrows_total;
        int mid=size_tot/2;
        int cndof=cndof_mid;
        for (int i=0;i<cbw_mid;i++) {
            for (int ii=0;ii<=i;ii++) {
                add_from_to_scaled_cp2k_csr(*P,0.5,0.0,\
                        mid+(-i-2)*cndof,\
                        mid+(-1)*cndof,\
                        mid+(-i-1+ii)*cndof,\
                        mid+(0+ii)*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_from_to_scaled_cp2k_csr(*P,0.5,1.0,\
                        mid+0*cndof,\
                        mid+(i+1)*cndof,\
                        mid+(-i-1+ii)*cndof,\
                        mid+(0+ii)*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_from_to_scaled_cp2k_csr(*P,0.5,0.0,\
                        mid+(-1)*cndof,\
                        mid+(-i-2)*cndof,\
                        mid+(0+ii)*cndof,\
                        mid+(-i-1+ii)*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
                add_from_to_scaled_cp2k_csr(*P,0.5,1.0,\
                        mid+(i+1)*cndof,\
                        mid+0*cndof,\
                        mid+(0+ii)*cndof,\
                        mid+(-i-1+ii)*cndof,\
                        cndof,cndof,MPI_COMM_WORLD);
            }
        }
        if (overwrite_first_last_diag_block) {
            add_from_to_scaled_cp2k_csr(*P,1.0,0.0,\
                    mid-2*cndof,\
                    mid-2*cndof,\
                    mid-1*cndof,\
                    mid-1*cndof,\
                    cndof,cndof,MPI_COMM_WORLD);
            add_from_to_scaled_cp2k_csr(*P,1.0,0.0,\
                    mid+1*cndof,\
                    mid+1*cndof,\
                    mid+0*cndof,\
                    mid+0*cndof,\
                    cndof,cndof,MPI_COMM_WORLD);
        }
    }

    std::vector<int> atom_of_bf;
    for (int a=0;a<cp2k_transport_params.n_atoms;a++) {
        for (int i=0;i<cp2k_transport_params.nsgf[a];i++) {
            atom_of_bf.push_back(a);
        }
    }
    atom_of_bf.push_back(cp2k_transport_params.n_atoms);

    std::vector<double> mulli(cp2k_transport_params.n_atoms,0.0);
    for (int i=0;i<S.nrows_local;i++) {
        for (int e=S.rowptr_local[i]-1;e<S.rowptr_local[i+1]-1;e++) {
            mulli[atom_of_bf[i+S.first_row]]+=2.0*S.nzvals_local[e]*P->nzvals_local[e];
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
    if (!rank) cout << "CHARGE " << cp2k_transport_params.iscf << " IS " << std::accumulate(cp2k_transport_params.zeff,cp2k_transport_params.zeff+cp2k_transport_params.n_atoms,-std::accumulate(mulli.begin(),mulli.end(),0.0)) << endl;

    mulli.assign(cp2k_transport_params.n_atoms,0.0);
    for (int i=0;i<S.nrows_local;i++) {
        for (int e=S.rowptr_local[i]-1;e<S.rowptr_local[i+1]-1;e++) {
            mulli[atom_of_bf[i+S.first_row]]+=2.0*S.nzvals_local[e]*KS.nzvals_local[e];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,&mulli[0],cp2k_transport_params.n_atoms,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (!rank) {
        stringstream mysstream;
        mysstream << "Mullipot_" << cp2k_transport_params.iscf;
        ofstream myfile(mysstream.str().c_str());
        myfile.precision(8);
        for (int i=0;i<cp2k_transport_params.n_atoms;i++) myfile<<mulli[i]<<"\n";
        myfile.close();
    }

    if (dens_mixing<1.0 || dens_mixing>0.0) {
        c_dscal(P->nze_local,dens_mixing,P->nzvals_local,1);
        c_daxpy(P->nze_local,1.0-dens_mixing,P_save,1,P->nzvals_local,1);
        delete[] P_save;
    }

    if (!rank) cout << "Transport iteration " << cp2k_transport_params.iscf << " finished" << endl;
}

