#include "Utilities.H"
#include "Density.H"
#include "GetSingularities.H"
#include "Quadrature.H"
#include "EnergyVector.H"
#include "pole.hpp"
#include <iterator>

Energyvector::Energyvector()
{
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);
}

Energyvector::~Energyvector()
{
}

int Energyvector::Execute(TCSR<double> *Overlap,TCSR<double> *KohnSham,int n_mu,double* muvec, int* contactvec,double* electronchargeperatom,double* derivativechargeperatom,double *Vatom,int n_atoms,int* atom_of_bf,c_transport_type transport_params)
{
    double sabtime;
    if ( Overlap->size_tot%transport_params.n_cells || transport_params.bandwidth<1 ) return (LOGCERR, EXIT_FAILURE);
    int tasks_per_point=transport_params.extra_int_param2;
    if (!iam) cout << "Distributing matrix over " << tasks_per_point << " tasks" << endl;
    if ( nprocs%tasks_per_point ) {
        if (!iam) cout << "Choose number of tasks per energy point as a divider of total number of tasks" << endl;
        return (LOGCERR, EXIT_FAILURE);
    }
    int distribute_pmat=1;
// allocate matrices to gather on every node
    sabtime=get_time(0.0);
    MPI_Comm matrix_comm;
    MPI_Comm eq_rank_matrix_comm;
    TCSR<double> *KohnShamCollect;
    TCSR<double> *OverlapCollect;
    TCSR<double> *Ps = NULL;
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
        if (distribute_pmat) {
            Ps = new TCSR<double>(OverlapCollect);
        } else {
            int *master_ranks = new int[n_mat_comm];
            if (matrix_rank==0) {
                MPI_Gather(&iam,1,MPI_INT,master_ranks,1,MPI_INT,0,eq_rank_matrix_comm);
            }
            MPI_Bcast(master_ranks,n_mat_comm,MPI_INT,0,MPI_COMM_WORLD);
            Ps = new TCSR<double>(Overlap,master_ranks,n_mat_comm,MPI_COMM_WORLD);
            delete[] master_ranks;
        }
        Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
    } else {
        KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
        OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
        MPI_Comm_split(MPI_COMM_WORLD,iam,iam,&matrix_comm);
        MPI_Comm_dup(MPI_COMM_WORLD,&eq_rank_matrix_comm);
        Ps = new TCSR<double>(OverlapCollect->size,OverlapCollect->n_nonzeros,OverlapCollect->findx);
        Ps->copy_index(OverlapCollect);
        Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
    }
    if (!iam) cout << "TIME FOR DISTRIBUTING MATRICES " << get_time(sabtime) << endl;
    std::vector<CPX> energyvector;
    std::vector<CPX> stepvector;
    std::vector<transport_methods::transport_method> methodvector;
    std::vector< std::vector<int> > propagating_sizes;
    if (determine_energyvector(energyvector,stepvector,methodvector,propagating_sizes,KohnSham,Overlap,muvec,contactvec,transport_params,n_mu)) return (LOGCERR, EXIT_FAILURE);
    std::vector<double> currentvector(energyvector.size(),0.0);
    std::vector<double> transmission(energyvector.size(),0.0);
    std::vector<double> dos_profile(energyvector.size(),0.0);
    double *eperatom = new double[n_atoms]();
    double *dperatom = new double[n_atoms]();
    int n_cells=transport_params.n_cells;
    int ndof=Overlap->size_tot/n_cells;
    TCSR<double> *OverlapCollectSave = new TCSR<double>(OverlapCollect);
    KohnShamCollect->settozeropbc(transport_params.bandwidth,ndof);
    OverlapCollect->settozeropbc(transport_params.bandwidth,ndof);
    TCSR<double> *OverlapCollectPBC = new TCSR<double>(1.0,OverlapCollectSave,-1.0,OverlapCollect);
    delete OverlapCollectSave;
    int matrix_id, n_mat_comm;
    MPI_Comm_size(eq_rank_matrix_comm,&n_mat_comm);
    MPI_Comm_rank(eq_rank_matrix_comm,&matrix_id);
    sabtime=get_time(0.0);
    unsigned int jpos;
    for (int iseq=0;iseq<int(ceil(double(energyvector.size())/n_mat_comm));iseq++)
        if ( (jpos=matrix_id+iseq*n_mat_comm)<energyvector.size())
            if (abs(stepvector[jpos])>0.0)
                if (density(KohnShamCollect,OverlapCollect,OverlapCollectPBC,Ps,energyvector[jpos],stepvector[jpos],methodvector[jpos],n_mu,muvec,contactvec,currentvector[jpos],transmission[jpos],dos_profile[jpos],propagating_sizes[jpos],atom_of_bf,eperatom,dperatom,Vatom,transport_params,distribute_pmat,jpos,matrix_comm))
                    return (LOGCERR, EXIT_FAILURE);
    if (!iam) cout << "TIME FOR DENSITY " << get_time(sabtime) << endl;
    delete KohnShamCollect;
    delete OverlapCollect;
    delete OverlapCollectPBC;
    MPI_Allreduce(eperatom,electronchargeperatom,n_atoms,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//    MPI_Allreduce(dperatom,derivativechargeperatom,n_atoms,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    delete[] eperatom;
    delete[] dperatom;
    if (!iam) {
        cout << "Number of electrons per unit cell";
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
        cout << "Total number of electrons " << e_total << endl;
    }
// trPS per energy point
    ofstream myfile;
    std::vector<double> currentvector2(energyvector.size(),0.0);
    std::vector<double> transmission2(energyvector.size(),0.0);
    std::vector<double> dos_profile2(energyvector.size(),0.0);
    MPI_Allreduce(&currentvector[0],&currentvector2[0],energyvector.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&transmission[0],&transmission2[0],energyvector.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&dos_profile[0],&dos_profile2[0],energyvector.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) {
        myfile.open("DOS_Profile");
        myfile.precision(15);
        for (uint iele=0;iele<energyvector.size();iele++)
            myfile << real(energyvector[iele]) << " " << imag(energyvector[iele]) << " " << real(stepvector[iele]) << " " << imag(stepvector[iele]) << " " << dos_profile2[iele] << endl;
        myfile.close();
    }
    if (!iam) {
        myfile.open("Transmission");
        myfile.precision(15);
        for (uint iele=0;iele<energyvector.size();iele++)
            myfile << real(energyvector[iele]) << " " << imag(energyvector[iele]) << " " << real(stepvector[iele]) << " " << imag(stepvector[iele]) << " " << transmission2[iele] << endl;
        myfile.close();
    }
    if (!iam) cout << "CURRENT IS " << c_ddot(energyvector.size(),&currentvector2[0],1,(double*)&stepvector[0],2) << endl;

    if (tasks_per_point > 1) {
        int matrix_rank, matrix_size;
        MPI_Comm_size(matrix_comm,&matrix_size);
        MPI_Comm_rank(matrix_comm,&matrix_rank);
        Ps->reduce(0,eq_rank_matrix_comm);
        for (int i_rank=0;i_rank<matrix_size;i_rank++) {
            Ps->scatter(Overlap,i_rank,MPI_COMM_WORLD);
        }
    } else {
        Ps->reducescatter(Overlap,MPI_COMM_WORLD);
    }
    delete Ps;
    MPI_Comm_free(&matrix_comm);
    MPI_Comm_free(&eq_rank_matrix_comm);

    return 0;
}

int Energyvector::determine_energyvector(std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,std::vector<transport_methods::transport_method> &methodvector,std::vector< std::vector<int> > &propagating_sizes,TCSR<double> *KohnSham,TCSR<double> *Overlap,double *muvec,int *contactvec,c_transport_type transport_params,int n_mu)
{
    double sabtime=get_time(0.0);
    double Temp=transport_params.extra_param3;
    Singularities singularities(transport_params,n_mu);
    if ( singularities.Execute(KohnSham,Overlap,contactvec) ) return (LOGCERR, EXIT_FAILURE);
    if (!iam) cout << "TIME FOR SINGULARITIES " << get_time(sabtime) << endl;
    for (int i_mu=0;i_mu<n_mu;i_mu++) singularities.write_bandstructure(i_mu);

    double delta_eps_fermi=-log((numeric_limits<double>::epsilon)())*K_BOLTZMANN*Temp;
    double muvec_min=*min_element(muvec,muvec+n_mu);
    double muvec_max=*max_element(muvec,muvec+n_mu);
    double nonequi_start=muvec_min-delta_eps_fermi;
    double nonequi_end=muvec_max+delta_eps_fermi;
    double energy_vb=*max_element(singularities.energies_vb.begin(),singularities.energies_vb.end());
    double energy_cb=*min_element(singularities.energies_cb.begin(),singularities.energies_cb.end());

    if (!transport_params.n_abscissae) {
        add_real_axis_energies(energy_cb,nonequi_end,energyvector,stepvector,methodvector,singularities.energies_extremum,transport_params,n_mu);
    } else if (transport_params.method==2) {
        add_cmpx_cont_energies(singularities.energy_gs,singularities.energy_gs,muvec_min,energyvector,stepvector,methodvector,transport_params,n_mu); //FOR PEX ONLY MU AND EM
    } else {
// all localized states with lowest fermi level corresponding to occupation of localized states in bandgap
        add_cmpx_cont_energies(singularities.energy_gs,singularities.energy_gs,muvec_min,energyvector,stepvector,methodvector,transport_params,n_mu);
        add_real_axis_energies(nonequi_start,nonequi_end,energyvector,stepvector,methodvector,singularities.energies_extremum,transport_params,n_mu);
    }

    ifstream evecfile("OMEN_E");
    if (evecfile) {
        energyvector.clear();
        istream_iterator<double> start_evec(evecfile), end_evec;
        energyvector.assign(start_evec,end_evec);
        methodvector.clear();
        methodvector.assign(energyvector.size(),transport_methods::WF);
        stepvector.clear();
        if (energyvector.size()==1) {
            stepvector.assign(1,CPX(1.0,0.0));
        } else {
            stepvector.assign(1,(energyvector[1]-energyvector[0])/2.0);
            for (uint istep=1;istep<energyvector.size()-1;istep++) {
                stepvector.push_back((energyvector[istep+1]-energyvector[istep-1])/2.0);
            }
            stepvector.push_back((energyvector[energyvector.size()-1]-energyvector[energyvector.size()-2])/2.0);
        }
    }
    evecfile.close();
    if (!iam) cout << "Size of Energyvector " << energyvector.size() << endl;
    if (!iam) {
        ofstream myfile;
        myfile.open("E_dat");
        myfile.precision(15);
        myfile << energyvector.size() << endl;
        for (uint iele=0;iele<energyvector.size();iele++)
            myfile << real(energyvector[iele]) << endl;
        myfile.close();
    }
// get propagating modes from bandstructure
    propagating_sizes.resize(energyvector.size());
    for (uint ie=0;ie<energyvector.size();ie++) propagating_sizes[ie].resize(n_mu);
    if (!iam) {
        std::vector< std::vector< std::vector<double> > > propagating = singularities.get_propagating(energyvector);
        for (uint ie=0;ie<energyvector.size();ie++) {
            for (int i_mu=0;i_mu<n_mu;i_mu++) {
                propagating_sizes[ie][i_mu]=propagating[i_mu][ie].size();
            }
        }
    }
    for (uint ie=0;ie<energyvector.size();ie++) {
        MPI_Bcast(&propagating_sizes[ie][0],n_mu,MPI_INT,0,MPI_COMM_WORLD);
    }
    return 0;
}

void Energyvector::add_real_axis_energies(double nonequi_start,double nonequi_end,std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,std::vector<transport_methods::transport_method> &methodvector,const std::vector< std::vector<double> > &energies_extremum,c_transport_type transport_params,int n_mu)
{
    std::vector<double> energylist;
    int n_energies;
    if (!iam) {
        energylist.push_back(nonequi_start);
        for (int i_mu=0;i_mu<n_mu;i_mu++)
            for (uint i_energies=0;i_energies<energies_extremum[i_mu].size();i_energies++)
                if (energies_extremum[i_mu][i_energies]>nonequi_start && energies_extremum[i_mu][i_energies]<nonequi_end)
                    energylist.push_back(energies_extremum[i_mu][i_energies]);
        energylist.push_back(nonequi_end);
        std::sort(energylist.begin(),energylist.end());
        n_energies=energylist.size();
    }
    MPI_Bcast(&n_energies,1,MPI_INT,0,MPI_COMM_WORLD);
    energylist.resize(n_energies);
    MPI_Bcast(&energylist[0],n_energies,MPI_DOUBLE,0,MPI_COMM_WORLD);
    double smallest_energy_distance=transport_params.extra_param2;
    if (!iam) cout<<"Smallest enery distance "<<smallest_energy_distance<<endl;
    if (!iam) cout<<"Max number of points per small interval "<<transport_params.extra_int_param1<<endl;
    if (!iam) cout<<"Average distance for big intervals "<<transport_params.extra_param1<<endl;
    if (!iam) cout<<"Singularities in range "<< n_energies-2 << endl;
    for (uint i_energies=1;i_energies<energylist.size();i_energies++) {
        int num_points_per_interval=max(transport_params.extra_int_param1,int(ceil(abs(energylist[i_energies]-energylist[i_energies-1])/transport_params.extra_param1)));
        while ((energylist[i_energies]-energylist[i_energies-1])/2.0*(1.0-cos(M_PI/(2.0*num_points_per_interval)))<smallest_energy_distance && num_points_per_interval>1)
//        while ((energylist[i_energies]-energylist[i_energies-1])/2.0*(1.0-tanh(M_PI/2.0*sinh(3.0)))<smallest_energy_distance && num_points_per_interval>1)
            --num_points_per_interval;
        if (num_points_per_interval>1) {
            Quadrature quadrature(quadrature_types::GC,energylist[i_energies-1],energylist[i_energies],num_points_per_interval);
            energyvector.insert(energyvector.end(),quadrature.abscissae.begin(),quadrature.abscissae.end());
            stepvector.insert(stepvector.end(),quadrature.weights.begin(),quadrature.weights.end());
            methodvector.resize(energyvector.size(),transport_methods::WF);
        }
    }
}

void Energyvector::add_cmpx_cont_energies(double start,double end,double mu,std::vector<CPX> &energyvector,std::vector<CPX> &stepvector,std::vector<transport_methods::transport_method> &methodvector,c_transport_type transport_params,int n_mu)
{
    double Temp=transport_params.extra_param3;
    enum choose_method {do_pexsi,do_pole_summation,do_contour,do_line};
    choose_method method=do_pexsi;
    if (method==do_pexsi) {
        int num_points_on_contour=transport_params.n_abscissae;
        energyvector.resize(num_points_on_contour);
        stepvector.resize(num_points_on_contour);
        double EM=abs(mu-start); // Max|E-mu| for all Eigenvalues E of 
        if (PEXSI::GetPoleDensity(&energyvector[0],&stepvector[0],num_points_on_contour,K_BOLTZMANN*Temp,0.0,EM,mu)) LOGCERR;
        c_zscal(num_points_on_contour,CPX(M_PI/2.0,0.0),&stepvector[0],1);
    } else if (method==do_pole_summation) {
        double Temp_r=1.0*Temp;
        double Temp_i=1.0*Temp;
        double eps=-log((numeric_limits<double>::epsilon)())*K_BOLTZMANN*Temp;
        double eps_r=-log((numeric_limits<double>::epsilon)())*K_BOLTZMANN*Temp_r;
        double mu_r=start-eps_r;
        double eps_i=-log((numeric_limits<double>::epsilon)())*K_BOLTZMANN*Temp_i;
        CPX    mu_i=CPX(mu_r-eps_r,eps_i);

        int nu;
        CPX zval;
        nu=0;
        while (imag(zval=CPX(mu,(2*nu+++1)*M_PI*K_BOLTZMANN*Temp))<2.0*eps_i) {
            energyvector.push_back(zval);
            stepvector.push_back(-CPX(0.0,2.0*M_PI*K_BOLTZMANN*Temp)*fermi(CPX(0.0,-1.0)*zval,CPX(0.0,-1.0)*mu_i,Temp_i,0));
        }

        nu=0;
        while (imag(zval=CPX(mu_r,(2*nu+++1)*M_PI*K_BOLTZMANN*Temp_r))<2.0*eps_i) {
            energyvector.push_back(zval);
            stepvector.push_back(+CPX(0.0,2.0*M_PI*K_BOLTZMANN*Temp_r)*fermi(CPX(0.0,-1.0)*zval,CPX(0.0,-1.0)*mu_i,Temp_i,0));
        }

        nu=0;
        while (real(zval=mu_i+CPX((2*nu+++1)*M_PI*K_BOLTZMANN*Temp_i,0.0))<mu+eps) {
            energyvector.push_back(zval);
            stepvector.push_back(CPX(2.0*M_PI*K_BOLTZMANN*Temp_i,0.0)*(fermi(zval,CPX(1.0,0.0)*mu,Temp,0)-fermi(zval,CPX(1.0,0.0)*mu_r,Temp_r,0)));
        }
    } else if (method==do_contour) {
        Quadrature quadrature(quadrature_types::CCGL,start,end,transport_params.n_abscissae);
        energyvector.insert(energyvector.end(),quadrature.abscissae.begin(),quadrature.abscissae.end());
        stepvector.insert(stepvector.end(),quadrature.weights.begin(),quadrature.weights.end());
    } else if (method==do_line) {
        Quadrature quadrature(quadrature_types::MR,start,end,transport_params.n_abscissae);
        energyvector.insert(energyvector.end(),quadrature.abscissae.begin(),quadrature.abscissae.end());
        stepvector.insert(stepvector.end(),quadrature.weights.begin(),quadrature.weights.end());
    }
    methodvector.resize(energyvector.size(),transport_methods::GF);
}
