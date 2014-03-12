#include "Utilities.H"
#include "CSR.H"
#include "Density.H"
#include "GetSingularities.H"
#include "Quadrature.H"

int energyvector(TCSR<double> *Overlap,TCSR<double> *KohnSham,int n_mu,double* muvec, int* contactvec,double* dopingvec,double* electronchargeperatom,double* derivativechargeperatom,int n_atoms,int* atom_of_bf,c_transport_type transport_params)
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    MPI_Comm matrix_comm;
    MPI_Comm eq_rank_matrix_comm;
    double sabtime;
// check parameters
    if ( Overlap->size_tot%transport_params.n_cells || transport_params.bandwidth<1 ) return (LOGCERR, EXIT_FAILURE);
    int tasks_per_point=1;
    if (tasks_per_point>1) if (!iam) cout<<"Distributing matrix over "<<tasks_per_point<<" tasks"<<endl;
    if ( nprocs%tasks_per_point ) {
        if (!iam) cout << "Choose number of tasks per energy point as a divider of total number of tasks" << endl;
        return (LOGCERR, EXIT_FAILURE);
    }
// additional parameters not needed in the future
    int n_cells=transport_params.n_cells;
// only to write out unit cell Hamiltonian and Overlap matrix
#ifdef _CONTACT_WRITEOUT
    int ndof=Overlap->size_tot/n_cells;
    TCSR<double> *Hcut = new TCSR<double>(KohnSham,0,ndof,0,(transport_params.bandwidth+1)*ndof);
    TCSR<double> *Hsp = new TCSR<double>(Hcut,0,MPI_COMM_WORLD);
    delete Hcut;
    if (!iam) {
        Hsp->shift_resize(0,ndof,0,(transport_params.bandwidth+1)*ndof);
        Hsp->change_findx(1);
        Hsp->write("H012");
    }
    delete Hsp;
    TCSR<double> *Scut = new TCSR<double>(Overlap,0,ndof,0,(transport_params.bandwidth+1)*ndof);
    TCSR<double> *Ssp = new TCSR<double>(Scut,MPI_COMM_WORLD);
    delete Scut;
    if (!iam) {
        Ssp->shift_resize(0,ndof,0,(transport_params.bandwidth+1)*ndof);
        Ssp->change_findx(1);
        Ssp->write("S012");
    }
    delete Ssp;
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
#endif
    c_dscal(KohnSham->n_nonzeros,transport_params.evoltfactor,KohnSham->nnz,1);
// allocate matrices to gather on every node
    sabtime=get_time(0.0);
    TCSR<double> *KohnShamCollect;
    TCSR<double> *OverlapCollect;
    TCSR<double> *Ps=NULL;
    if (tasks_per_point > 1) {
        MPI_Comm KS_matrix_comm;
        MPI_Comm KS_eq_rank_matrix_comm;
        KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD,tasks_per_point,&KS_matrix_comm,&KS_eq_rank_matrix_comm);
        OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD,tasks_per_point,&matrix_comm,&eq_rank_matrix_comm);
        int comm_result, eq_comm_result;
        MPI_Comm_compare(KS_matrix_comm,matrix_comm,&comm_result);
        MPI_Comm_compare(KS_eq_rank_matrix_comm,eq_rank_matrix_comm,&eq_comm_result);
        if (comm_result!=MPI_CONGRUENT || eq_comm_result!=MPI_CONGRUENT) return (LOGCERR, EXIT_FAILURE);
        MPI_Comm_free(&KS_matrix_comm);
        MPI_Comm_free(&KS_eq_rank_matrix_comm);
        int matrix_rank, n_mat_comm;
        MPI_Comm_size(eq_rank_matrix_comm,&n_mat_comm);
        MPI_Comm_rank(matrix_comm,&matrix_rank);
        int *master_ranks = new int[n_mat_comm];
        if (matrix_rank==0) {
            MPI_Gather(&iam,1,MPI_INT,master_ranks,1,MPI_INT,0,eq_rank_matrix_comm);
        }
        MPI_Bcast(master_ranks,n_mat_comm,MPI_INT,0,MPI_COMM_WORLD);
        Ps = new TCSR<double>(Overlap,master_ranks,n_mat_comm,MPI_COMM_WORLD);
        Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
        delete[] master_ranks;
    } else {
        KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
        OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
        MPI_Comm_split(MPI_COMM_WORLD,iam,iam,&matrix_comm);
        MPI_Comm_dup(MPI_COMM_WORLD,&eq_rank_matrix_comm);
        Ps = new TCSR<double>(OverlapCollect->size,OverlapCollect->n_nonzeros,OverlapCollect->findx);
        Ps->copy_index(OverlapCollect);
        Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
#ifdef _CP2K_WRITEOUT
if (!iam) {
c_dscal(KohnShamCollect->n_nonzeros,1.0/transport_params.evoltfactor,KohnShamCollect->nnz,1);
KohnShamCollect->write_CSR("KohnSham");
OverlapCollect->write_CSR("Overlap");
ofstream paramoutfile("TransportParams");
paramoutfile << transport_params.method << endl;
paramoutfile << transport_params.bandwidth << endl;
paramoutfile << transport_params.n_cells << endl;
paramoutfile << transport_params.n_occ << endl;
paramoutfile << transport_params.n_abscissae << endl;
paramoutfile << transport_params.n_kpoint << endl;
paramoutfile << transport_params.extra_int_param1 << endl;
paramoutfile << transport_params.extra_int_param2 << endl;
paramoutfile << transport_params.extra_int_param3 << endl;
paramoutfile << transport_params.evoltfactor << endl;
paramoutfile << transport_params.colzero_threshold << endl;
paramoutfile << transport_params.eps_limit << endl;
paramoutfile << transport_params.eps_decay << endl;
paramoutfile << transport_params.eps_singularities << endl;
paramoutfile << transport_params.extra_param1 << endl;
paramoutfile << transport_params.extra_param2 << endl;
paramoutfile << transport_params.extra_param3 << endl;
paramoutfile.close();
}
MPI_Barrier(MPI_COMM_WORLD);
exit(0);
#endif
#ifdef _OMEN_WRITEOUT
if (!iam) {
int ndof=Overlap->size_tot/n_cells;
KohnShamCollect->removepbc(transport_params.bandwidth,ndof);
KohnShamCollect->change_findx(1);
KohnShamCollect->write_CSR_bin("H_4.bin");
OverlapCollect->removepbc(transport_params.bandwidth,ndof);
OverlapCollect->change_findx(1);
OverlapCollect->write_CSR_bin("S_4.bin");
}
MPI_Barrier(MPI_COMM_WORLD);
exit(0);
#endif
    }
    if (!iam) cout << "TIME FOR DISTRIBUTING MATRICES " << get_time(sabtime) << endl;
// determine singularity stuff
    double Temp=transport_params.extra_param3;
    sabtime=get_time(0.0);
    Singularities singularities(transport_params);
    if ( singularities.Execute(KohnSham,Overlap,n_mu,muvec,dopingvec,contactvec) ) return (LOGCERR, EXIT_FAILURE);
muvec[0]=(muvec[0]+muvec[1])/2.0;muvec[1]=muvec[0];
    if (!iam) cout << "TIME FOR SINGULARITIES " << get_time(sabtime) << endl;
// determine elements in nonequilibrium range
    double k_b=K_BOLTZMANN;
    double nonequi_start=muvec[0]-35.0*k_b*Temp;
int integral_over_real_axis=transport_params.extra_int_param2;
if (!iam) if (integral_over_real_axis) cout << "Integral over real axis" << endl;
if (integral_over_real_axis) nonequi_start=singularities.energy_gs;
    double nonequi_end=muvec[n_mu-1]+35.0*k_b*Temp;
    std::vector<double> energylist;
    energylist.push_back(nonequi_start);
    for (uint i_energies=0;i_energies<singularities.energies.size();i_energies++)
        if (singularities.energies[i_energies]>nonequi_start && singularities.energies[i_energies]<nonequi_end)
            energylist.push_back(singularities.energies[i_energies]);
    energylist.push_back(nonequi_end);
    if (!iam) cout << energylist.size() << " energies in nonequilibrium range from a total of " << singularities.energies.size() << " singularity points" << endl;
    if (!iam) {
        ofstream myfile;
        myfile.open("SingularityList");
        myfile.precision(15);
        for (uint iele=0;iele<energylist.size();iele++)
            myfile << energylist[iele] << endl;
        myfile.close();
    }
// get integration points along complex contour with gauss legendre
    int num_points_on_contour=transport_params.n_abscissae;
    Quadrature gausslegendre(quadrature_types::CCGL,singularities.energy_gs,nonequi_start,Temp,muvec[0],num_points_on_contour);
    std::vector<CPX> energyvector(gausslegendre.abscissae);
    std::vector<CPX> stepvector(gausslegendre.weights);
    std::vector<transport_methods::transport_method> methodvector(gausslegendre.abscissae.size(),transport_methods::NEGF);
// get integration points along real axis with gauss cheby
    double smallest_energy_distance=transport_params.extra_param2;
if (!iam) cout<<"Smallest enery distance "<<smallest_energy_distance<<endl;
if (!iam) cout<<"Max number of points per small interval "<<transport_params.extra_int_param1<<endl;
if (!iam) cout<<"Average distance for big intervals "<<transport_params.extra_param1<<endl;
    for (uint i_energies=1;i_energies<energylist.size();i_energies++) {
        int num_points_per_interval=max(transport_params.extra_int_param1,int(ceil(abs(energylist[i_energies]-energylist[i_energies-1])/transport_params.extra_param1)));
        while ((energylist[i_energies]-energylist[i_energies-1])*(1.0-cos(M_PI/(2.0*num_points_per_interval)))/2.0<smallest_energy_distance && num_points_per_interval>1)
            --num_points_per_interval;
        if (num_points_per_interval>1) {
            Quadrature gausscheby(quadrature_types::GC,energylist[i_energies-1],energylist[i_energies],0.0,0.0,num_points_per_interval);
            energyvector.insert(energyvector.end(),gausscheby.abscissae.begin(),gausscheby.abscissae.end());
            stepvector.insert(stepvector.end(),gausscheby.weights.begin(),gausscheby.weights.end());
            std::vector<transport_methods::transport_method> methodblock(gausscheby.abscissae.size(),transport_methods::WF);
            methodvector.insert(methodvector.end(),methodblock.begin(),methodblock.end());
        }
    }
/*
int ntrapez=2000;
Quadrature trapez(quadrature_types::TR,nonequi_start,nonequi_end,Temp,muvec[0],ntrapez);
energyvector.resize(ntrapez);
stepvector.resize(ntrapez);
methodvector.resize(ntrapez);
copy(trapez.abscissae.begin(),trapez.abscissae.end(),energyvector.begin());
copy(trapez.weights.begin(),trapez.weights.end(),stepvector.begin());
fill(methodvector.begin(),methodvector.end(),transport_methods::WF);
*/
/*
energyvector.resize(1);
energyvector[0]=CPX(-9.12376,0.0);
stepvector.resize(1);
stepvector[0]=CPX(1.0,0.0);
methodvector.resize(1);
methodvector[0]=transport_methods::WF;
*/
    if (!iam) cout << "Size of Energyvector " << energyvector.size() << endl;
// get propagating modes from bandstructure OF RIGHT CONTACT ONLY -> ONLY FOR EQUILIBRIUM OR I HAVE TO IMPLEMENT LEFT CONTACT AS WELL
    std::vector< std::vector<double> > propagating;
    if (!iam) singularities.get_propagating(propagating,energyvector);
    singularities.delete_matrices();
    int *propagating_sizes = new int[energyvector.size()];
    if (!iam) {
        for (int ie=0;ie<energyvector.size();ie++) {
            propagating_sizes[ie]=propagating[ie].size();
        }
    }
    MPI_Bcast(propagating_sizes,energyvector.size(),MPI_INT,0,MPI_COMM_WORLD);
if (!iam) {
ofstream myfile("Propagating");
for (int ii=0;ii<energyvector.size();ii++) myfile << real(energyvector[ii]) << " " << propagating[ii].size() << endl;
myfile.close();
}
// run distributed
    std::vector<double> currentvector(energyvector.size(),0.0);
    double *eperatom = new double[n_atoms]();
    double *dperatom = new double[n_atoms]();
    int matrix_id, n_mat_comm;
    MPI_Comm_size(eq_rank_matrix_comm,&n_mat_comm);
    MPI_Comm_rank(eq_rank_matrix_comm,&matrix_id);
    uint jpos;
    for (int iseq=0;iseq<int(ceil(double(energyvector.size())/n_mat_comm));iseq++)
//      for (int iseq=0;iseq<1;iseq++)
        if ( (jpos=matrix_id+iseq*n_mat_comm)<energyvector.size())
            if (abs(stepvector[jpos])>0.0)
                if (density(KohnShamCollect,OverlapCollect,Ps,energyvector[jpos],stepvector[jpos],methodvector[jpos],n_mu,muvec,contactvec,currentvector[jpos],propagating_sizes[jpos],atom_of_bf,eperatom,dperatom,transport_params,matrix_comm))
                    return (LOGCERR, EXIT_FAILURE);
    delete KohnShamCollect;
    delete OverlapCollect;
    delete[] propagating_sizes;
    MPI_Allreduce(eperatom,electronchargeperatom,n_atoms,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(dperatom,derivativechargeperatom,n_atoms,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    delete[] eperatom;
    delete[] dperatom;
    for (int ibw=0;ibw<transport_params.bandwidth;ibw++) {
        c_dcopy(n_atoms/n_cells,&electronchargeperatom[transport_params.bandwidth*n_atoms/n_cells],1,&electronchargeperatom[ibw*n_atoms/n_cells],1);
        c_dcopy(n_atoms/n_cells,&electronchargeperatom[(n_cells-transport_params.bandwidth-1)*n_atoms/n_cells],1,&electronchargeperatom[(n_cells-ibw-1)*n_atoms/n_cells],1);
        c_dcopy(n_atoms/n_cells,&derivativechargeperatom[transport_params.bandwidth*n_atoms/n_cells],1,&derivativechargeperatom[ibw*n_atoms/n_cells],1);
        c_dcopy(n_atoms/n_cells,&derivativechargeperatom[(n_cells-transport_params.bandwidth-1)*n_atoms/n_cells],1,&derivativechargeperatom[(n_cells-ibw-1)*n_atoms/n_cells],1);
    }
    if (!iam) {
        cout << "Number of electrons per unit cell NEW ";
        double e_total=0.0;
        for (int icell=0;icell<n_cells;icell++) {
            double e_per_unit_cell=0.0;
            for (int iatom=icell*n_atoms/n_cells;iatom<(icell+1)*n_atoms/n_cells;iatom++) {
                e_per_unit_cell+=electronchargeperatom[iatom];
                e_total+=electronchargeperatom[iatom];
            }
            cout << " " << e_per_unit_cell;
        }
        cout << endl;
        cout << "Total number of electrons NEW " << e_total << endl;
    }
// trPS per energy point
    ofstream myfile;
    std::vector<double> currentvector2(energyvector.size(),0.0);
    MPI_Allreduce(&currentvector[0],&currentvector2[0],energyvector.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) {
        myfile.open("DOS_Profile");
        myfile.precision(15);
        for (uint iele=0;iele<energyvector.size();iele++)
            if (abs(stepvector[iele])>0.0)
                myfile << real(energyvector[iele]) << " " << real(stepvector[iele]) << " " << currentvector2[iele] << endl;
        myfile.close();
    }

    if (tasks_per_point > 1) {
        int matrix_rank, matrix_size;
        MPI_Comm_size(matrix_comm,&matrix_size);
        MPI_Comm_rank(matrix_comm,&matrix_rank);
        if (matrix_rank==0) {
            Ps->reduce(0,eq_rank_matrix_comm);
        }
        if (transport_params.method!=3) {
            Ps->scatter(Overlap,0,MPI_COMM_WORLD);
        }
    } else {
        if (transport_params.method!=3) {
            Ps->reducescatter(Overlap,MPI_COMM_WORLD);
        }
    }
    delete Ps;
    MPI_Comm_free(&matrix_comm);
    MPI_Comm_free(&eq_rank_matrix_comm);

    return 0;
}
