#include "DiagScaLapack.H"
#include "p_eig.H"

int diagscalapack(TCSR<double> *Overlap,TCSR<double> *KohnSham,transport_parameters *parameters_transport)
{
    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);

    int nocc=parameters_transport->n_occ;

    double sabtime;
    double *OVfull=NULL;
    double *KSfull=NULL;

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

sabtime=get_time(0.0);
    if (p_eig(KSfull,OVfull,eigval,nvec,MPI_COMM_WORLD)) return (LOGCERR, EXIT_FAILURE);
if (!iam) cout << "Time for p_eig " << get_time(sabtime) << endl;

    c_dscal(OverlapCollect->n_nonzeros,0.0,OverlapCollect->nnz,1);
if (!iam) sabtime=get_time(0.0);
//    if (!iam) OverlapCollect->psipsidagger(KSfull,nocc,1.0); THIS TAKES REALLY A LOT OF TIME
    if (!iam) full_transpose(nocc,nvec,KSfull,OVfull);
    if (!iam) OverlapCollect->psipsidagger_transpose(OVfull,nocc,1.0);

if (!iam) cout << "Time for Dens " << get_time(sabtime) << endl;

    if (!iam) {
        delete[] KSfull;
        delete[] OVfull;
    }
//if (!iam) for (int i_out=0;i_out<nocc;i_out++) cout << eigval[i_out] << "\n";
    delete[] eigval;

    OverlapCollect->reducescatter(Overlap,MPI_COMM_WORLD);
    delete OverlapCollect;

    return 0;
}

