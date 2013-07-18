#include "ScaLapack.H"

int kpointintegration(TCSR<double> *Overlap,TCSR<double> *KohnSham,TCSR<double> *P_Matrix,c_transport_type trans_params)
{
    int bw=trans_params.bandwidth;
    int ncells=trans_params.n_cells;
    int nocc=trans_params.n_occ;

    int size=Overlap->size_tot;
    int ndof=size/ncells;
    if (ndof*ncells!=size) return (LOGCERR, EXIT_FAILURE);
    int sizesq=size*size;
    CPX iunit=CPX(0.0,1.0);

    int iam, nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&iam);

    double step,value;
    if (nprocs==1) {
        step=1.0;
        value=0.0;
    } else {
        step=1.0/(nprocs-1);
        value=-M_PI+iam*step*2*M_PI;
    }

    double alphastep;
    if (nprocs==1)
        alphastep=1.0;
    else if (!iam || iam==nprocs-1)
        alphastep=1.0/2;
    else
        alphastep=1.0;

    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,MPI_COMM_WORLD);
    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,MPI_COMM_WORLD);

    double *ham=new double[sizesq];
    double *tau=new double[sizesq];
    double *taut=new double[sizesq];
    KohnShamCollect->moveawaypbc(taut,tau,bw,ndof);
    KohnShamCollect->sparse_to_full(ham,size,size);
    delete KohnShamCollect;

    CPX *overlapnnzsave=new CPX[OverlapCollect->n_nonzeros];
    c_dcopy(OverlapCollect->n_nonzeros,OverlapCollect->nnz,1,(double*)overlapnnzsave,2);

    CPX *A=new CPX[sizesq];
    for (int i=0;i<size;i++)
        for (int j=0;j<size;j++)
            A[j*size+i]=ham[j*size+i]+tau[j*size+i]*exp(iunit*value)+taut[j*size+i]*exp(-iunit*value);

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
            B[j*size+i]=sam[j*size+i]+sau[j*size+i]*exp(iunit*value)+saut[j*size+i]*exp(-iunit*value);

    delete[] sam;
    delete[] sau;
    delete[] saut;
    CPX *B2=new CPX[sizesq];
    c_zcopy(sizesq,B,1,B2,1);

    double *eigval=new double[size];
    double *rwork = new double[3*size-2];
    int info;
    CPX twork;
    c_zhegv(1,'V','U',size,A,size,B,size,eigval,&twork,-1,rwork,&info);
    int lwork=int(real(twork));
    CPX *work=new CPX[lwork];
    c_zhegv(1,'V','U',size,A,size,B,size,eigval,work,lwork,rwork,&info);
    if (info) { cout<<info<<endl; return info; }
    delete[] work;
    delete[] rwork;

    double fermi;
    if (!iam) fermi=(eigval[nocc-1]+eigval[nocc])/2;
    MPI_Bcast(&fermi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    int nocclocal=-1;
    while (eigval[++nocclocal]<fermi);
    int nocctry=0;
    MPI_Allreduce(&nocclocal,&nocctry,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) cout << "New Fermi Level " << fermi << endl;
    cout << iam << " contributes " << nocclocal << endl;
    if (!iam) cout << "New Number of Electrons " << (double) nocctry/nprocs << endl;
    while (nocctry!=nocc*nprocs) {
        if (nocctry<nocc*nprocs) fermi=(fermi+eigval[nocc])/2;
        if (nocctry>nocc*nprocs) fermi=(fermi+eigval[nocc-1])/2;
        MPI_Bcast(&fermi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if (!iam) cout << "New Fermi Level " << fermi << endl;
        nocclocal=-1;
        while (eigval[++nocclocal]<fermi);
        cout << iam << " contributes " << nocclocal << endl;
        MPI_Allreduce(&nocclocal,&nocctry,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        if (!iam) cout << "New Number of Electrons " << (double) nocctry/nprocs << endl;
    }

    full_transpose(size,size,A,B);
    delete[] A;
    TCSR<CPX> *Ps = new TCSR<CPX>(OverlapCollect,B,CPX(alphastep*step,0.0),nocclocal);
    delete OverlapCollect;
    Ps->sparse_to_full(B,size,size);
    CPX trPS_k=c_zdotc(sizesq,B,1,B2,1);
    CPX trPS_kSum=CPX(0.0,0.0);
    delete[] B;
    delete[] B2;
    MPI_Allreduce(&trPS_k,&trPS_kSum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) cout << "Number of electrons from Sum_k(tr(P[k]S[k])) " << real(trPS_kSum) << endl;

    int middlebegin=Ps->edge_i[bw*ndof]-Ps->findx;
    int middlenumber=Ps->edge_i[bw*ndof+ndof]-Ps->findx-middlebegin;
    CPX trPScpx=c_zdotc(middlenumber,&Ps->nnz[middlebegin],1,&overlapnnzsave[middlebegin],1);
    CPX trPScpxSum=CPX(0.0,0.0);
    MPI_Allreduce(&trPScpx,&trPScpxSum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
    if (!iam) cout << "Number of electrons from density matrix strip " << real(trPScpxSum)*ncells << endl;

    TCSR<CPX> *P_Matrix_cpx = new TCSR<CPX>(P_Matrix->size,P_Matrix->n_nonzeros,P_Matrix->findx);
    P_Matrix_cpx->copy_index(P_Matrix);
    Ps->reducescatter(P_Matrix_cpx,MPI_COMM_WORLD);

    c_dcopy(P_Matrix->n_nonzeros,(double*)P_Matrix_cpx->nnz,2,P_Matrix->nnz,1);

    return 0;

}
