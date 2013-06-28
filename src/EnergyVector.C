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
    paramfile >> parameter_sab.bandwidth >> parameter_sab.ncells >> parameter_sab.evoltfactor >> parameter_sab.colzerothr >> parameter_sab.eps_limit >> parameter_sab.eps_decay >> energystart >> energyend;
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
    Ps->init_variable(Ps->nnz,Ps->n_nonzeros);
// run 
    int method=2;
    if (density(KohnShamCollect,OverlapCollect,Ps,energy,step,method,parameter_sab)) return (cerr<<__LINE__<<endl, EXIT_FAILURE);
// trPS
    CPX trPScpx=c_ddot(Ps->n_nonzeros,(double*) Ps->nnz,2,OverlapCollect->nnz,1);
    CPX trPScpxSum=CPX(0.0,0.0);
    MPI_Allreduce(&trPScpx,&trPScpxSum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) cout << "Number of electrons from density matrix tr(PS) " << real(trPScpxSum) << endl;
// sum and scatter to distributed p matrix
    Ps->reducescatterconvert(P_Matrix,MPI_COMM_WORLD);

    return 0;
}
