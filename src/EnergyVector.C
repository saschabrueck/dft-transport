#include "Utilities.H"
#include "CSR.H"
#include "Density.H"
#include "GetSingularities.H"
#include "Quadrature.H"

int energyvector(TCSR<double> *Overlap,TCSR<double> *KohnSham,TCSR<double> *P_Matrix,c_transport_type transport_params)
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    double sabtime;
// check parameters
    if ( Overlap->size_tot%transport_params.n_cells || transport_params.bandwidth<1 ) return (LOGCERR, EXIT_FAILURE);
// additional parameters not needed in the future
    int ncells=transport_params.n_cells;
    int size=Overlap->size_tot;
    int ndof=size/ncells;
    int begin;
    int number;
// allocate matrices to gather on every node
    TCSR<double> *KohnShamCollectPrimal = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
    TCSR<double> *KohnShamCollect;
    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
// add the potential 
    int addpotential=0;
    if (addpotential) {
        cout << "ADDING POTENTIAL" << endl;
        double *Vvec = new double[Overlap->size_tot];
        int na=5;
        int nb=na+5;
        int nc=ncells;
        for (int ivvec=0;          ivvec<na*4*13*12; ivvec++) Vvec[ivvec] =  0.0;
        for (int ivvec=na*4*13*12; ivvec<nb*4*13*12; ivvec++) Vvec[ivvec] =  1.5;
        for (int ivvec=nb*4*13*12; ivvec<nc*4*13*12; ivvec++) Vvec[ivvec] =  0.5;
        TCSR<double> *Pot = new TCSR<double>(Overlap,Vvec);
        delete[] Vvec;
        KohnShamCollect = new TCSR<double>(transport_params.evoltfactor,KohnShamCollectPrimal,1.0,Pot);
        delete Pot;
        delete KohnShamCollectPrimal;
    } else {
        KohnShamCollect = KohnShamCollectPrimal;
        KohnShamCollectPrimal = NULL;
        c_dscal(KohnShamCollect->n_nonzeros,transport_params.evoltfactor,KohnShamCollect->nnz,1);
    }
    TCSR<double> *Ps = new TCSR<double>(OverlapCollect->size,OverlapCollect->n_nonzeros,OverlapCollect->findx);
    Ps->copy_index(OverlapCollect);
    Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
// determine singularity stuff
    sabtime=get_time(0.0);
    Singularities *singularities = new Singularities(KohnShamCollect,OverlapCollect,transport_params);
    cout << "TIME FOR SINGULARITIES " << get_time(sabtime) << endl;
    double Temp=0.0;
    int n_mu=2;
    double *muvec = new double[n_mu];
    muvec[0]=singularities->energy_vbe;
//    muvec[1]=singularities->energy_vbe+0.5;
    muvec[1]=singularities->energy_vbe;
// determine elements in nonequilibrium range
    double k_b=K_BOLTZMANN;
    double nonequi_start=muvec[0]-35.0*k_b*Temp;
    double nonequi_end=muvec[n_mu-1]+35.0*k_b*Temp;
    std::vector<double> energylist;
    energylist.push_back(nonequi_start);
    for (int i_energies=0;i_energies<singularities->n_energies;i_energies++)
        if (singularities->energies[i_energies]>nonequi_start && singularities->energies[i_energies]<nonequi_end)
            energylist.push_back(singularities->energies[i_energies]);
    energylist.push_back(nonequi_end);
    if (!iam) cout << energylist.size() << " energies in nonequilibrium range from a total of " << singularities->n_energies << " singularity points" << endl;
// get integration points along complex contour with gauss legendre
    int num_points_on_contour=transport_params.n_abscissae;
    Quadrature gausslegendre(quadrature_types::CCGL,singularities->energy_gs,nonequi_start,Temp,muvec[0],num_points_on_contour);
    std::vector<CPX> energyvector(gausslegendre.abscissae);
    std::vector<CPX> stepvector(gausslegendre.weights);
    std::vector<transport_methods::transport_method> methodvector(num_points_on_contour,transport_methods::NEGF);
// get integration points along real axis with gauss cheby
    int num_points_per_interval=transport_params.n_abscissae;
    std::vector<transport_methods::transport_method> methodblock(num_points_per_interval,transport_methods::WF);
    Quadrature *gausscheby;
    double smallest_energy_distance=transport_params.extra_param1;
    for (uint i_energies=1;i_energies<energylist.size();i_energies++) {
        if (energylist[i_energies-1]<energylist[i_energies]-smallest_energy_distance) {
            gausscheby = new Quadrature(quadrature_types::GC,energylist[i_energies-1],energylist[i_energies],0.0,0.0,num_points_per_interval);
            energyvector.insert(energyvector.end(),gausscheby->abscissae.begin(),gausscheby->abscissae.end());
            stepvector.insert(stepvector.end(),gausscheby->weights.begin(),gausscheby->weights.end());
            methodvector.insert(methodvector.end(),methodblock.begin(),methodblock.end());
            delete gausscheby;
        }
    }
    delete singularities;
    cout << "Size of Energyvector " << energyvector.size() << endl;
// run distributed
    double skip_point_weight_thr=0.0;
    uint jpos;
    for (int iseq=0;iseq<int(ceil(double(energyvector.size())/nprocs));iseq++)
        if ( (jpos=iam+iseq*nprocs)<energyvector.size())
            if (abs(stepvector[jpos])>skip_point_weight_thr)
                if (density(KohnShamCollect,OverlapCollect,Ps,energyvector[jpos],stepvector[jpos],muvec,n_mu,methodvector[jpos],transport_params)) return (LOGCERR, EXIT_FAILURE);
// copy density to corners
//    Ps->contactdensity(ndof,bw);
// trPS
    double trPS;
    double trPS_Sum;
// trPS per unit cell
    ofstream myfile;
    if (!iam) myfile.open("NumberOfElectrons");
    if (!iam) myfile << "Number of electrons per unit cell";
    for (int idens=0;idens<ncells;idens++) {
        begin=Ps->edge_i[idens*ndof]-Ps->findx;
        number=Ps->edge_i[idens*ndof+ndof]-Ps->findx-begin;
        trPS=c_ddot(number,&Ps->nnz[begin],1,&OverlapCollect->nnz[begin],1);
        MPI_Allreduce(&trPS,&trPS_Sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if (!iam) myfile << " " << trPS_Sum;
    }
    if (!iam) myfile << endl;
    if (!iam) myfile.close();
// trPS total
    trPS=c_ddot(Ps->n_nonzeros,Ps->nnz,1,OverlapCollect->nnz,1);
    MPI_Allreduce(&trPS,&trPS_Sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) cout << "Total number of electrons " << trPS_Sum << endl;
// sum and scatter to distributed p matrix
    Ps->reducescatter(P_Matrix,MPI_COMM_WORLD);
    delete Ps;
    delete KohnShamCollect;
    delete OverlapCollect;
    delete[] muvec;

    return 0;
}
