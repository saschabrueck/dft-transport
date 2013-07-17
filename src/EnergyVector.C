#include "CSR.H"
#include "Density.H"
#include "GetSingularities.H"
#include "Quadrature.H"

int energyvector(TCSR<double> *Overlap,TCSR<double> *KohnSham,TCSR<double> *P_Matrix,c_transport_type transport_params)
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);
// check parameters
    if ( Overlap->size_tot%transport_params.n_cells || transport_params.bandwidth<1 ) return (cerr<<__LINE__<<endl, EXIT_FAILURE);
// allocate matrices to gather on every node
    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
    TCSR<CPX> *Ps = new TCSR<CPX>(OverlapCollect->size,OverlapCollect->n_nonzeros,OverlapCollect->findx);
    Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
// determine singularity stuff
    Singularities *singularities = new Singularities(KohnShamCollect,OverlapCollect,transport_params);
// get integration points along real axis with simple uniform trapezoidal rule
/*    double step;
    if (nprocs==1)
        step=2.0;
    else
        step=(singularities->energy_vbe-singularities->energy_gs)/(nprocs-1);
    CPX energy=CPX(singularities->energy_gs+iam*step,0.0);
    if (!iam || iam==nprocs-1)
        step/=2;*/
// here is the real quadrature
    int num_points_per_interval=transport_params.n_abscissae;
    int size_energyvector=num_points_per_interval*singularities->n_energies;
    CPX *energyvector = new CPX[size_energyvector];
    CPX *stepvector = new CPX[size_energyvector];
    Quadrature *gausscheby = new Quadrature(quadrature_types::GC,singularities->energy_gs,singularities->energies[0],0.0,singularities->energy_vbe,num_points_per_interval);
    std::copy(gausscheby->abscissae.begin(),gausscheby->abscissae.end(),energyvector);
    std::copy(gausscheby->weights.begin(),gausscheby->weights.end(),stepvector);
    delete gausscheby;
    for (int i_energies=1;i_energies<singularities->n_energies;i_energies++) {
        gausscheby = new Quadrature(quadrature_types::GC,singularities->energies[i_energies-1],singularities->energies[i_energies],0.0,singularities->energy_vbe,num_points_per_interval);
        std::copy(gausscheby->abscissae.begin(),gausscheby->abscissae.end(),&energyvector[i_energies*num_points_per_interval]);
        std::copy(gausscheby->weights.begin(),gausscheby->weights.end(),&stepvector[i_energies*num_points_per_interval]);
        delete gausscheby;
    }
    delete singularities;
// run distributed
    int method=2;
    int seq_per_cpu=int(ceil(double(size_energyvector)/nprocs));
    int jpos;
    for (int iseq=0;iseq<seq_per_cpu;iseq++)
        if ( (jpos=iam+iseq*nprocs)<size_energyvector)
            if (density(KohnShamCollect,OverlapCollect,Ps,energyvector[jpos],stepvector[jpos],method,transport_params)) return (cerr<<__LINE__<<endl, EXIT_FAILURE);
// trPS
    CPX trPScpx=c_ddot(Ps->n_nonzeros,(double*) Ps->nnz,2,OverlapCollect->nnz,1);
    CPX trPScpxSum=CPX(0.0,0.0);
    MPI_Allreduce(&trPScpx,&trPScpxSum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) cout << "Number of electrons from density matrix tr(PS) " << real(trPScpxSum) << endl;
// sum and scatter to distributed p matrix
    Ps->reducescatterconvert(P_Matrix,MPI_COMM_WORLD);

    return 0;
}
