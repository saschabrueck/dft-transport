#include "DiagScaLapack.H"
#include "p_eig.H"
#include "array_tools.H"

int diagscalapack(TCSR<double> *Overlap,TCSR<double> *KohnSham,transport_parameters *parameters_transport)
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);

    int nocc=parameters_transport->n_occ;

    double sabtime;
    double *OVfull, *KSfull;

    int nvec=KohnSham->size_tot;
    if (nvec!=Overlap->size_tot) return (LOGCERR, EXIT_FAILURE);

    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,0,MPI_COMM_WORLD);
    if (!iam) {
        KSfull = new double[nvec*nvec];
        KohnShamCollect->sparse_to_full(KSfull,nvec,nvec);
    } // end if
    delete KohnShamCollect;

    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);
    if (!iam) {
        OVfull = new double[nvec*nvec];
        OverlapCollect ->sparse_to_full(OVfull,nvec,nvec);
    } // end if

    double *eigval = new double [nvec];

    if (!iam) sabtime=get_time(0.0);
    if (p_eig(KSfull,OVfull,eigval,nvec,MPI_COMM_WORLD)) return (LOGCERR, EXIT_FAILURE);
    if (!iam) cout << "Time for p_eig " << get_time(sabtime) << endl;

    if (!iam) full_transpose(nvec,nvec,KSfull,OVfull);

// THIS IS OF COURSE UNNECESSARILY COMPLICATED BUT IT IS TO TRY MPI ROUTINES
    if (!iam) sabtime=get_time(0.0);
    TCSR<double> *Density;
    if (!iam)
        Density = new TCSR<double>(OverlapCollect,OVfull,1.0/nprocs,nocc);
    else
        Density = new TCSR<double>(OverlapCollect->size,OverlapCollect->n_nonzeros,OverlapCollect->findx);
    if (!iam) cout << "Time for Dens " << get_time(sabtime) << endl;
    delete OverlapCollect;

    if (!iam) {
        delete[] KSfull;
        delete[] OVfull;
    }
    if (!iam) for (int i_out=0;i_out<nocc;i_out++) cout << eigval[i_out] << "\n";
    delete[] eigval;

    MPI_Bcast(Density->nnz,Density->n_nonzeros,MPI_DOUBLE,0,MPI_COMM_WORLD);
    Density->reducescatter(Overlap,MPI_COMM_WORLD);

    delete Density;

    return 0;
}

