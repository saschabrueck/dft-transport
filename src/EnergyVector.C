#include "CSR.H"
#include "Types.H"
#include "Density.H"

int energyvector(TCSR<double> *Overlap,TCSR<double> *KohnSham,TCSR<double> *P_Matrix)
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);
// get parameters
    ParameterStruct parameter_sab;
    CPX energy;
    double energystart;
    double energyend;
    ifstream paramfile("Parameters");
    paramfile >> parameter_sab.bandwidth >> parameter_sab.ndof >> parameter_sab.evoltfactor >> parameter_sab.colzerothr >> parameter_sab.eps_limit >> parameter_sab.eps_decay >> energystart >> energyend;
    paramfile.close();
    if (parameter_sab.bandwidth<1) return (cerr<<__LINE__<<endl, EXIT_FAILURE);
// get integration points along real axis with simple uniform trapezoidal rule
    double energystep;
    if (nprocs==1)
        energystep=2.0;
    else
        energystep=(energyend-energystart)/(nprocs-1);
    energy=CPX(energystart+iam*energystep,0.0);
    double step;
    if (!iam || iam==nprocs-1)
        step=energystep/2;
    else
        step=energystep;
// allocate matrices to gather on every node
    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
//    TCSR<double> *KohnShamCollect = new TCSR<double>("CNT5_h.mtx");
//    TCSR<double> *OverlapCollect  = new TCSR<double>("CNT5_s.mtx");
    TCSR<CPX> *Ps = new TCSR<CPX>(OverlapCollect->size,OverlapCollect->n_nonzeros,OverlapCollect->findx);
// run 
    int method=2;
    if (density(KohnShamCollect,OverlapCollect,Ps,energy,method,parameter_sab)) return (cerr<<__LINE__<<endl, EXIT_FAILURE);
// sum and scatter to distributed p matrix
    Ps->reducescatterfactorconvert(P_Matrix,MPI_COMM_WORLD,step);

    return 0;
}
