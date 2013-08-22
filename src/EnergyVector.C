#include "Utilities.H"
#include "CSR.H"
#include "Density.H"
#include "GetSingularities.H"
#include "Quadrature.H"

void cgqf ( int nt, int kind, double alpha, double beta, double a, double b,
  double t[], double wts[] );

int energyvector(TCSR<double> *Overlap,TCSR<double> *KohnSham,TCSR<double> *P_Matrix,c_transport_type transport_params)
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    double sabtime;
// check parameters
    if ( Overlap->size_tot%transport_params.n_cells || transport_params.bandwidth<1 ) return (LOGCERR, EXIT_FAILURE);
// allocate matrices to gather on every node
    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
    TCSR<double> *Ps = new TCSR<double>(OverlapCollect->size,OverlapCollect->n_nonzeros,OverlapCollect->findx);
    Ps->copy_index(OverlapCollect);
    Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
// additional parameters not needed in the future
    int bw=transport_params.bandwidth;
    int ncells=transport_params.n_cells;
    int size=OverlapCollect->size_tot;
    int ndof=size/ncells;
    int begin;
    int number;
// determine singularity stuff
    sabtime=get_time(0.0);
    Singularities *singularities = new Singularities(KohnShamCollect,OverlapCollect,transport_params);
    cout << "TIME FOR SINGULARITIES " << get_time(sabtime) << endl;
    int n_energies_below_vbe=0;
    int start_energy_wf_method=0;
    double muvec[2];
    muvec[0]=singularities->energy_vbe;
//    muvec[1]=singularities->energy_vbe+0.5;
    muvec[1]=singularities->energy_vbe;
    while(singularities->energies[n_energies_below_vbe++]<=muvec[1]);
//    while(singularities->energies[start_energy_wf_method++]<=muvec[0]); // I HAVE TO CHECK IF THIS IS RIGHT OR IF IT IS TOO HIGH
    if (!iam) cout << n_energies_below_vbe-start_energy_wf_method << "computed, " << n_energies_below_vbe << " below, " << singularities->n_energies << " total" << endl;
    if (!iam) for (int iii=0;iii<n_energies_below_vbe;iii++) cout << singularities->energies[iii] << endl;
// get integration points along real axis with simple uniform trapezoidal rule
/*
    int size_energyvector=1;
//    int size_energyvector=nprocs;
    CPX *energyvector = new CPX[size_energyvector];
    CPX *stepvector = new CPX[size_energyvector];
    transport_methods::transport_method *methodvector = new transport_methods::transport_method[size_energyvector];
    double skip_point_weight_thr=0.0;
    if (size_energyvector==1) {
        stepvector[0]=CPX(1.0,0.0);
        energyvector[0]=CPX(1.2,0.0);
    } else {
        double step=(singularities->energy_vbe-singularities->energy_gs)/(size_energyvector-1);
        for (int i_step=0;i_step<size_energyvector;i_step++) {
            methodvector[i_step]=transport_methods::WF;
            energyvector[i_step]=CPX(singularities->energy_gs+i_step*step,0.0);
            if (!i_step || i_step==size_energyvector-1)
                stepvector[i_step]=CPX(step/2.0,0.0);
            else
                stepvector[i_step]=CPX(step,0.0);
        }
    }
*/
// get integration points along real axis with gauss cheby

    Quadrature *gausscheby;
    int num_points_per_interval=transport_params.n_abscissae;
    double skip_point_weight_thr=10E-12;
    int size_energyvector;
    if (singularities->energy_gs<singularities->energies[0])
        size_energyvector=num_points_per_interval*singularities->n_energies;
    else
        size_energyvector=num_points_per_interval*(singularities->n_energies-1);
    CPX *energyvector = new CPX[size_energyvector];
    CPX *stepvector = new CPX[size_energyvector];
    transport_methods::transport_method *methodvector = new transport_methods::transport_method[size_energyvector];
    double *energyvectord = new double[size_energyvector]();
    double *stepvectord = new double[size_energyvector]();
    int i_start=0;
    if (singularities->energy_gs<singularities->energies[0]) {
        cgqf(num_points_per_interval,2,0.0,0.0,singularities->energy_gs,singularities->energies[0],energyvectord,stepvectord);
        for (int ipoint=0;ipoint<num_points_per_interval;ipoint++)
            stepvectord [ipoint] *= sqrt ( ( energyvectord[ipoint] - singularities->energy_gs ) * ( singularities->energies[0] - energyvectord[ipoint] ) );
/*
        gausscheby = new Quadrature(quadrature_types::GC,singularities->energy_gs,singularities->energies[0],0.0,singularities->energy_vbe,num_points_per_interval);
        std::copy(gausscheby->abscissae.begin(),gausscheby->abscissae.end(),energyvector);
        std::copy(gausscheby->weights.begin(),gausscheby->weights.end(),stepvector);
        delete gausscheby;
*/
        i_start=num_points_per_interval;
    }
    for (int i_energies=0;i_energies<n_energies_below_vbe-1;i_energies++) {
        cgqf(num_points_per_interval,2,0.0,0.0,singularities->energies[i_energies],singularities->energies[i_energies+1],&energyvectord[i_start+i_energies*num_points_per_interval],&stepvectord[i_start+i_energies*num_points_per_interval]);
        for (int ipoint=0;ipoint<num_points_per_interval;ipoint++) {
            stepvectord [ipoint+i_start+i_energies*num_points_per_interval] *= sqrt ( ( energyvectord[ipoint+i_start+i_energies*num_points_per_interval] - singularities->energies[i_energies] ) * ( singularities->energies[i_energies+1] - energyvectord[ipoint+i_start+i_energies*num_points_per_interval] ) );
        if (singularities->energies[i_energies+1]-singularities->energies[i_energies]<0.001) stepvectord [ipoint+i_start+i_energies*num_points_per_interval] = 0;
        }
/*
        gausscheby = new Quadrature(quadrature_types::GC,singularities->energies[i_energies],singularities->energies[i_energies+1],0.0,singularities->energy_vbe,num_points_per_interval);
        std::copy(gausscheby->abscissae.begin(),gausscheby->abscissae.end(),&energyvector[i_start+i_energies*num_points_per_interval]);
        std::copy(gausscheby->weights.begin(),gausscheby->weights.end(),&stepvector[i_start+i_energies*num_points_per_interval]);
        delete gausscheby;
*/
    }
    for (int i_evec=0;i_evec<size_energyvector;i_evec++) {
        energyvector[i_evec]=CPX(energyvectord[i_evec],0.0);
        stepvector[i_evec]=CPX(stepvectord[i_evec],0.0);
        methodvector[i_evec]=transport_methods::WF;
    }

// get integration points along complex contour with gauss legendre
/*
    int num_points_on_contour=20;
    Quadrature gausslegendre(quadrature_types::CCGL,singularities->energy_gs,singularities->energy_vbe,0.0,singularities->energy_vbe,num_points_on_contour);
    int size_energyvector=num_points_on_contour;
    transport_methods::transport_method *methodvector = new transport_methods::transport_method[size_energyvector];
    CPX *energyvector = new CPX[size_energyvector];
    CPX *stepvector = new CPX[size_energyvector];
    std::fill_n (methodvector,size_energyvector,transport_methods::NEGF);
    std::copy(gausslegendre.abscissae.begin(),gausslegendre.abscissae.end(),energyvector);
    std::copy(gausslegendre.weights.begin(),gausslegendre.weights.end(),stepvector);
    double skip_point_weight_thr=0.0;
*/
// delete singularities
    delete singularities;
// run distributed
    int seq_per_cpu=int(ceil(double(size_energyvector)/nprocs));
    int jpos;
    for (int iseq=0;iseq<seq_per_cpu;iseq++)
        if ( (jpos=iam+iseq*nprocs)<size_energyvector)
            if (abs(stepvector[jpos])>skip_point_weight_thr)
                if (density(KohnShamCollect,OverlapCollect,Ps,energyvector[jpos],stepvector[jpos],muvec,2,methodvector[jpos],transport_params)) return (LOGCERR, EXIT_FAILURE);
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

    return 0;
}
