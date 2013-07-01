#include "ScaLapack.H"

int kpointintegration(TCSR<double> *Overlap,TCSR<double> *KohnSham,TCSR<double> *P_Matrix,int nocc)
{
    int bw=2;
    int ncells=5;

    int size=Overlap->size_tot;
    int ndof=size/ncells;
    if (ndof*ncells!=size) return (cerr<<__LINE__<<endl, EXIT_FAILURE);
    int sizesq=size*size;
    CPX iunit=CPX(0.0,1.0);

    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);

    double step,value;
    if (nprocs==1)
        value=0;
    else {
        step=2*M_PI/(nprocs-1);
        value=-M_PI+iam*step;
    }

    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);

    double *ham=new double[sizesq];
    double *tau=new double[sizesq];
    double *taut=new double[sizesq];
    KohnShamCollect->moveawaypbc(taut,tau,bw,ndof);
    KohnShamCollect->sparse_to_full(ham,size,size);
    delete KohnShamCollect;

    CPX *A=new CPX[sizesq];
    for (int i=0;i<size;i++)
        for (int j=0;j<size;j++)
            A[j*size+i]=ham[j*size+i]+tau[j*size+i]*exp(iunit*value)+tau[i*size+j]*exp(-iunit*value);

    delete[] ham;
    delete[] tau;
    delete[] taut;

    double *sam=new double[sizesq];
    double *sau=new double[sizesq];
    double *saut=new double[sizesq];
    OverlapCollect->moveawaypbc(saut,sau,bw,ndof);
    OverlapCollect->sparse_to_full(sam,size,size);

    CPX *B=new CPX[sizesq];
    for (int i=0;i<size;i++)
        for (int j=0;j<size;j++)
            B[j*size+i]=sam[j*size+i]+sau[j*size+i]*exp(iunit*value)+sau[i*size+j]*exp(-iunit*value);

    delete[] sam;
    delete[] sau;
    delete[] saut;

    double *eigval=new double[size];
    double rwork[3*size-2];
    int info;
    CPX twork;
    c_zhegv(1,'V','U',size,A,size,B,size,eigval,&twork,-1,rwork,&info);
    int lwork=int(real(twork));
    CPX *work=new CPX[lwork];
    c_zhegv(1,'V','U',size,A,size,B,size,eigval,work,lwork,rwork,&info);
    if (info) { cout<<info<<endl; return info; }

    double alphastep;
    if (nprocs==1)
        alphastep=1.0;
    else if (!iam || iam==nprocs-1)
        alphastep=step/2;
    else
        alphastep=step;

    full_transpose(size,size,A,B);
    delete[] A;
    TCSR<CPX> *Ps = new TCSR<CPX>(OverlapCollect,B,CPX(alphastep,0.0),nocc);
    delete[] B;
    delete OverlapCollect;

    TCSR<CPX> *P_Matrix_cpx = new TCSR<CPX>(P_Matrix->size,P_Matrix->n_nonzeros,P_Matrix->findx);
    P_Matrix_cpx->copy_index(P_Matrix);
    Ps->reducescatter(P_Matrix_cpx,MPI_COMM_WORLD);

    c_dcopy(P_Matrix->n_nonzeros,(double*)P_Matrix_cpx->nnz,2,P_Matrix->nnz,1);

    double trPS=c_ddot(P_Matrix->n_nonzeros,P_Matrix->nnz,1,Overlap->nnz,1);
    double trPS_Sum=0;
    MPI_Allreduce(&trPS,&trPS_Sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) cout << "Number of electrons from density matrix tr(PS) " << trPS_Sum << endl;

    return 0;

}
