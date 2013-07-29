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
    if ( Overlap->size_tot%transport_params.n_cells || transport_params.bandwidth<1 ) return (LOGCERR, EXIT_FAILURE);
// allocate matrices to gather on every node
    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
    TCSR<double> *Ps = new TCSR<double>(OverlapCollect->size,OverlapCollect->n_nonzeros,OverlapCollect->findx);
    Ps->copy_index(OverlapCollect);
    Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
// determine singularity stuff
    Singularities *singularities = new Singularities(KohnShamCollect,OverlapCollect,transport_params);
// get integration points along real axis with simple uniform trapezoidal rule
/*
    int size_energyvector=1;
//    int size_energyvector=nprocs;
    CPX *energyvector = new CPX[size_energyvector];
    CPX *stepvector = new CPX[size_energyvector];
    double skip_point_weight_thr=0.0;
    if (size_energyvector==1) {
        stepvector[0]=CPX(1.0,0.0);
        energyvector[0]=CPX(1.2,0.0);
    } else {
        double step=(singularities->energy_vbe-singularities->energy_gs)/(size_energyvector-1);
        for (int i_step=0;i_step<size_energyvector;i_step++) {
            energyvector[i_step]=CPX(singularities->energy_gs+i_step*step,0.0);
            if (!i_step || i_step==size_energyvector-1)
                stepvector[i_step]=CPX(step/2.0,0.0);
            else
                stepvector[i_step]=CPX(step,0.0);
        }
    }
*/
// get integration points along real axis with gauss cheby
/*
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
    int i_start=0;
    if (singularities->energy_gs<singularities->energies[0]) {
        gausscheby = new Quadrature(quadrature_types::GC,singularities->energy_gs,singularities->energies[0],0.0,singularities->energy_vbe,num_points_per_interval);
        std::copy(gausscheby->abscissae.begin(),gausscheby->abscissae.end(),energyvector);
        std::copy(gausscheby->weights.begin(),gausscheby->weights.end(),stepvector);
        delete gausscheby;
        i_start=num_points_per_interval;
    }
    for (int i_energies=0;i_energies<singularities->n_energies-1;i_energies++) {
        gausscheby = new Quadrature(quadrature_types::GC,singularities->energies[i_energies],singularities->energies[i_energies+1],0.0,singularities->energy_vbe,num_points_per_interval);
        std::copy(gausscheby->abscissae.begin(),gausscheby->abscissae.end(),&energyvector[i_start+i_energies*num_points_per_interval]);
        std::copy(gausscheby->weights.begin(),gausscheby->weights.end(),&stepvector[i_start+i_energies*num_points_per_interval]);
        delete gausscheby;
    }
*/
// get integration points along complex contour with gauss legendre

    int num_points_on_contour=48;
    Quadrature gausslegendre(quadrature_types::CCGL,singularities->energy_gs,singularities->energy_vbe,0.0,singularities->energy_vbe,num_points_on_contour);
    int size_energyvector=num_points_on_contour;
    CPX *energyvector = new CPX[size_energyvector];
    CPX *stepvector = new CPX[size_energyvector];
    std::copy(gausslegendre.abscissae.begin(),gausslegendre.abscissae.end(),energyvector);
    std::copy(gausslegendre.weights.begin(),gausslegendre.weights.end(),stepvector);
    double skip_point_weight_thr=0.0;

// delete singularities
    delete singularities;
// run distributed
    int seq_per_cpu=int(ceil(double(size_energyvector)/nprocs));
    int jpos;
    for (int iseq=0;iseq<seq_per_cpu;iseq++)
        if ( (jpos=iam+iseq*nprocs)<size_energyvector)
            if (abs(stepvector[jpos])>skip_point_weight_thr)
                if (density(KohnShamCollect,OverlapCollect,Ps,energyvector[jpos],stepvector[jpos],transport_methods::NEGF,transport_params)) return (LOGCERR, EXIT_FAILURE);
// trPS
// ich hab ja noch den Ps hier, der generell und nicht tri ist aber nullen hat, also kann ich das normale overlapcollect nehmen
// aber ich muss halt mal zwei nehmen und einmal die diagonalelementproduktsumme wieder abziehen
    int bw=transport_params.bandwidth;
    int ncells=transport_params.n_cells;
    int size=OverlapCollect->size_tot;
    int ndof=size/ncells;
    int middlebegin=Ps->edge_i[bw*ndof]-Ps->findx;
    int middlenumber=Ps->edge_i[bw*ndof+ndof]-Ps->findx-middlebegin;
    double trPS_1=2*c_ddot(middlenumber,&Ps->nnz[middlebegin],1,&OverlapCollect->nnz[middlebegin],1);
    double trPS_2=0.0;
    for (int ii=bw*ndof;ii<bw*ndof+ndof;ii++) trPS_2+=Ps->nnz[Ps->diag_pos[ii]]*OverlapCollect->nnz[OverlapCollect->diag_pos[ii]];
    double trPS=trPS_1-trPS_2;
//das ist das alte
//    double trPS=c_ddot(Ps->n_nonzeros,Ps->nnz,1,OverlapCollect->nnz,1);
    double trPS_Sum;
    MPI_Allreduce(&trPS,&trPS_Sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) cout << "Number of electrons from density matrix tr(PS) " << trPS_Sum*ncells << endl;
// sum and scatter to distributed p matrix
    Ps->reducescatter(P_Matrix,MPI_COMM_WORLD);
    delete Ps;
    delete KohnShamCollect;
    delete OverlapCollect;

    return 0;
}
