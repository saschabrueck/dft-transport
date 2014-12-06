#include "ScaLapack.H"

void pdgeev(int n,double *a,int *desca,double *wr,double *wi,double *vr,int *descvr,\
            double *work,int lwork,MPI_Comm comm,int *info)
{
    int ictxt,NPROW,NPCOL,mycol,myrow,nout,ROWMAX,COLMAX;
    int MB,NB,itau,descq[9],descqz[9],descvt[9],nr_loc,nc_loc;
    double *tau,*q,*qz,*vt,*aglobal,*qzglobal,*vtglobal;

    ictxt    = desca[1];
    MB       = desca[4];
    NB       = desca[5];
    
    Cblacs_gridinfo(ictxt,&NPROW,&NPCOL,&myrow,&mycol);
    itau     = c_numroc(n-1,NB,mycol,0,NPCOL);
    tau      = new double[itau];
    
    nr_loc   = c_numroc(n,MB,myrow,0,NPROW);
    nc_loc   = c_numroc(n,NB,mycol,0,NPCOL);
    ROWMAX   = (n-1)/(NPROW*MB)+1;
    COLMAX   = (n-1)/(NPCOL*NB)+1;
    
    c_descinit(descq,n,n,MB,NB,0,0,ictxt,nr_loc,info);
    c_descinit(descqz,n,n,MB,NB,0,0,ictxt,nr_loc,info);
    c_descinit(descvt,n,n,MB,NB,0,0,ictxt,nr_loc,info);

    q        = new double[nr_loc*nc_loc];
    qz       = new double[nr_loc*nc_loc];
    vt       = new double[nr_loc*nc_loc];
    aglobal  = new double[n*n];
    qzglobal = new double[n*n];
    vtglobal = new double[n*n];
    
    c_pdgehrd(n,1,n,a,1,1,desca,tau,work,lwork,info);
    c_pdlacpy('L',n,n,a,1,1,desca,q,1,1,descq);
    get_hess(myrow,mycol,n,n,nr_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,a);
    c_pdlaset('A',n,n,0.0,1.0,qz,1,1,descqz);
    c_pdormhr('L','N',n,n,1,n,q,1,1,descq,tau,qz,1,1,descqz,work,lwork,info);
    c_pdlahqr(true,true,n,1,n,a,desca,wr,wi,1,n,qz,descqz,work,lwork,0,0,info);
    
    receive(myrow,mycol,n,n,nr_loc,nc_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,aglobal,\
            a,comm);
    receive(myrow,mycol,n,n,nr_loc,nc_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,qzglobal,\
            qz,comm);
    
    c_dhseqr('S','V',n,1,n,aglobal,n,wr,wi,qzglobal,n,work,lwork,info);
    c_dtrevc('R','A',false,n,aglobal,n,NULL,1,vtglobal,n,n,&nout,work,info);

    distribute(myrow,mycol,n,n,nr_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,qzglobal,qz);
    distribute(myrow,mycol,n,n,nr_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,vtglobal,vt);
    
    c_pdgemm('N','N',n,n,n,1.0,qz,1,1,descqz,vt,1,1,descvt,0.0,vr,1,1,descvr);

    delete[] tau;
    delete[] q;
    delete[] qz;
    delete[] vt;
    delete[] aglobal;
    delete[] qzglobal;
    delete[] vtglobal;

}

/************************************************************************************************/

void distribute(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
                int ROWMAX,int COLMAX,double *A,double *Aloc)
{
    int IC,IR,ic0,icN,ir0,irM,sizeM,ic;

    if((myrow>=0)&&(mycol>=0)){
        for(IC=0;IC<COLMAX;IC++){    
            for(IR=0;IR<ROWMAX;IR++){    
                ic0   = mycol*NB+IC*NPCOL*NB;
                icN   = min(ic0+NB,N);
                ir0   = myrow*MB+IR*NPROW*MB;
                irM   = min(ir0+MB,M);
                sizeM = irM-ir0;
                if(sizeM){
                    for(ic=ic0;ic<icN;ic++){
                        c_dcopy(sizeM,&A[ic*M+ir0],1,&Aloc[(IC*NB+ic-ic0)*NR_loc+IR*MB],1);
                    }
                }
                else{
                    break;
                }
            }
        }
    }
            
}

/************************************************************************************************/

void zdistribute(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
                 int ROWMAX,int COLMAX,CPX *A,CPX *Aloc)
{
    int IC,IR,ic0,icN,ir0,irM,sizeM,ic;

    if((myrow>=0)&&(mycol>=0)){
        for(IC=0;IC<COLMAX;IC++){    
            for(IR=0;IR<ROWMAX;IR++){    
                ic0   = mycol*NB+IC*NPCOL*NB;
                icN   = min(ic0+NB,N);
                ir0   = myrow*MB+IR*NPROW*MB;
                irM   = min(ir0+MB,M);
                sizeM = irM-ir0;
                if(sizeM){
                    for(ic=ic0;ic<icN;ic++){
                        c_zcopy(sizeM,&A[ic*M+ir0],1,&Aloc[(IC*NB+ic-ic0)*NR_loc+IR*MB],1);
                    }
                }
                else{
                    break;
                }
            }
        }
    }
            
}

/************************************************************************************************/

void get_hess(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
              int ROWMAX,int COLMAX,double *Aloc)
{
    int IC,IR,ic0,icN,ir0,irM,ic,ir;

    if((myrow>=0)&&(mycol>=0)){
        for(IC=0;IC<COLMAX;IC++){    
            for(IR=0;IR<ROWMAX;IR++){    
                ic0   = mycol*NB+IC*NPCOL*NB;
                icN   = min(ic0+NB,N);
                ir0   = myrow*MB+IR*NPROW*MB;
                irM   = min(ir0+MB,M);
                for(ir=ir0;ir<irM;ir++){
                    for(ic=ic0;ic<icN;ic++){
                        if(ic+1<ir){
                            Aloc[(IC*NB+ic-ic0)*NR_loc+IR*MB+ir-ir0]=0;
                        }
                    }
                }
            }
        }
    }
            
}

/************************************************************************************************/

void get_zhess(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
              int ROWMAX,int COLMAX,CPX *Aloc)
{
    int IC,IR,ic0,icN,ir0,irM,ic,ir;

    if((myrow>=0)&&(mycol>=0)){
        for(IC=0;IC<COLMAX;IC++){    
            for(IR=0;IR<ROWMAX;IR++){    
                ic0   = mycol*NB+IC*NPCOL*NB;
                icN   = min(ic0+NB,N);
                ir0   = myrow*MB+IR*NPROW*MB;
                irM   = min(ir0+MB,M);
                for(ir=ir0;ir<irM;ir++){
                    for(ic=ic0;ic<icN;ic++){
                        if(ic+1<ir){
                            Aloc[(IC*NB+ic-ic0)*NR_loc+IR*MB+ir-ir0]=CPX(0.0,0.0);
                        }
                    }
                }
            }
        }
    }
            
}

/************************************************************************************************/

void receive(int myrow,int mycol,int M, int N,int NR_loc,int NC_loc,int MB,int NB,int NPROW,\
             int NPCOL,int ROWMAX,int COLMAX,double *A,double *Aloc,MPI_Comm comm)
{
    MPI_Status status;
    int mr,mc;

    if((myrow>=0)&&(mycol>=0)){
        if((myrow==0)&&(mycol==0)){
            insert_block(myrow,mycol,M,N,NR_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,A,Aloc);
            for(mr=0;mr<NPROW;mr++){
                for(mc=0;mc<NPCOL;mc++){
                    if((mr+mc)>0){
                        MPI_Recv(&NR_loc,1,MPI_INT,mr*NPCOL+mc,0,comm,&status);
                        MPI_Recv(&NC_loc,1,MPI_INT,mr*NPCOL+mc,1,comm,&status);
                        MPI_Recv(Aloc,NR_loc*NC_loc,MPI_DOUBLE,mr*NPCOL+mc,2,comm,&status);
                        insert_block(mr,mc,M,N,NR_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,A,Aloc);
                    }
                }
            }
        }else{
            MPI_Send(&NR_loc,1,MPI_INT,0,0,comm);
            MPI_Send(&NC_loc,1,MPI_INT,0,1,comm);
            MPI_Send(Aloc,NR_loc*NC_loc,MPI_DOUBLE,0,2,comm);
        }
    }
    MPI_Bcast(A,M*N,MPI_DOUBLE,0,comm);
}

/************************************************************************************************/

void zreceive(int myrow,int mycol,int M, int N,int NR_loc,int NC_loc,int MB,int NB,int NPROW,\
              int NPCOL,int ROWMAX,int COLMAX,CPX *A,CPX *Aloc,MPI_Comm comm)
{
    MPI_Status status;
    int mr,mc;

    if((myrow>=0)&&(mycol>=0)){
        if((myrow==0)&&(mycol==0)){
            zinsert_block(myrow,mycol,M,N,NR_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,A,Aloc);
            for(mr=0;mr<NPROW;mr++){
                for(mc=0;mc<NPCOL;mc++){
                    if((mr+mc)>0){
                        MPI_Recv(&NR_loc,1,MPI_INT,mr*NPCOL+mc,0,comm,&status);
                        MPI_Recv(&NC_loc,1,MPI_INT,mr*NPCOL+mc,1,comm,&status);
                        MPI_Recv(Aloc,NR_loc*NC_loc,MPI_DOUBLE_COMPLEX,mr*NPCOL+mc,2,\
                                 comm,&status);
                        zinsert_block(mr,mc,M,N,NR_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,A,Aloc);
                    }
                }
            }
        }else{
            MPI_Send(&NR_loc,1,MPI_INT,0,0,comm);
            MPI_Send(&NC_loc,1,MPI_INT,0,1,comm);
            MPI_Send(Aloc,NR_loc*NC_loc,MPI_DOUBLE_COMPLEX,0,2,comm);
        }
    }
    MPI_Bcast(A,M*N,MPI_DOUBLE_COMPLEX,0,comm);
}

/************************************************************************************************/

void insert_block(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
                  int ROWMAX,int COLMAX,double *A,double *Aloc)
{
    int IC,IR,ic0,icN,ir0,irM,sizeM,ic;

    for(IC=0;IC<=COLMAX;IC++){    
        for(IR=0;IR<=ROWMAX;IR++){    
            ic0   = mycol*NB+IC*NPCOL*NB;
            icN   = min(ic0+NB,N);
            ir0   = myrow*MB+IR*NPROW*MB;
            irM   = min(ir0+MB,M);
            sizeM = irM-ir0;
            if(sizeM){
                for(ic=ic0;ic<icN;ic++){
                    c_dcopy(sizeM,&Aloc[(IC*NB+ic-ic0)*NR_loc+IR*MB],1,&A[ic*M+ir0],1);
                }
            }
            else{
                break;
            }
        }
    }
}

/************************************************************************************************/

void zinsert_block(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
                   int ROWMAX,int COLMAX,CPX *A,CPX *Aloc)
{
    int IC,IR,ic0,icN,ir0,irM,sizeM,ic;

    for(IC=0;IC<=COLMAX;IC++){    
        for(IR=0;IR<=ROWMAX;IR++){    
            ic0   = mycol*NB+IC*NPCOL*NB;
            icN   = min(ic0+NB,N);
            ir0   = myrow*MB+IR*NPROW*MB;
            irM   = min(ir0+MB,M);
            sizeM = irM-ir0;
            if(sizeM){
                for(ic=ic0;ic<icN;ic++){
                    c_zcopy(sizeM,&Aloc[(IC*NB+ic-ic0)*NR_loc+IR*MB],1,&A[ic*M+ir0],1);
                }
            }
            else{
                break;
            }
        }
    }
}

/************************************************************************************************/

void pdgeev_driver(int n, double *a, double *wr, double *wi,double *vr,int NPROW,int NPCOL,\
                   MPI_Comm comm,int *info)
{
    int ictxt,mycol,myrow,nr_loc,nc_loc;
    int ROWMAX,COLMAX,desca[9],descvr[9];
    int MB,NB;
    double *aloc,*vrloc;
    
    ictxt        = comm;
    get_blocking_factors(NPROW,NPCOL,n,n,&MB,&NB);
    Cblacs_gridinit(&ictxt,"Row-major",NPROW,NPCOL);
    Cblacs_gridinfo(ictxt,&NPROW,&NPCOL,&myrow,&mycol);

    nr_loc       = c_numroc(n,MB,myrow,0,NPROW);
    nc_loc       = c_numroc(n,NB,mycol,0,NPCOL);

    ROWMAX       = (n-1)/(NPROW*MB)+1;
    COLMAX       = (n-1)/(NPCOL*NB)+1;

    aloc         = new double[nr_loc*nc_loc];
    vrloc        = new double[nr_loc*nc_loc];

    distribute(myrow,mycol,n,n,nr_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,a,aloc);

    c_descinit(desca,n,n,MB,NB,0,0,ictxt,nr_loc,info);
    c_descinit(descvr,n,n,MB,NB,0,0,ictxt,nr_loc,info);

    int ihip     = c_numroc(n,MB,myrow,0,NPROW);
    int ilrow    = c_indxg2p(1,MB,myrow,0,NPROW);
    int ilcol    = c_indxg2p(1,NB,mycol,0,NPCOL);
    int ihlp     = c_numroc(n,MB,myrow,ilrow,NPROW);
    int inlq     = c_numroc(n,NB,mycol,ilcol,NPCOL);
    int lwork    = NB*(MB+max(ihip+1,ihlp+inlq));
    double *work = new double[lwork];
    
    pdgeev(n,aloc,desca,wr,wi,vrloc,descvr,work,lwork,comm,info);

    receive(myrow,mycol,n,n,nr_loc,nc_loc,MB,NB,NPROW,NPCOL,ROWMAX,COLMAX,vr,vrloc,comm);

    Cblacs_gridexit(ictxt);

    delete[] aloc;
    delete[] vrloc;
    delete[] work;

}

/************************************************************************************************/

void pdgemm_driver(char transa,char transb,int m,int n,int k,double alpha,double *a,\
                   int lda,double *b,int ldb,double beta,double *c,int ldc,int NPROW,\
                   int NPCOL,MPI_Comm comm)
{
    int ictxt,MB,NB,myrow,mycol,info;
    int anr_loc,anc_loc,aROWMAX,aCOLMAX,na;
    int bnr_loc,bnc_loc,bROWMAX,bCOLMAX,nb;
    int cnr_loc,cnc_loc,cROWMAX,cCOLMAX,nc;
    int desca[9],descb[9],descc[9];
    double *aloc,*bloc,*cloc;

    if(transa == 'N'){
        na       = k;
    }else{
        na       = m;
    }

    if(transb == 'N'){
        nb       = n;
        nc       = n;
    }else{
        nb       = k;
        nc       = k;
    }
    
    ictxt        = comm;
    Cblacs_gridinit(&ictxt,"Row-major",NPROW,NPCOL);
    Cblacs_gridinfo(ictxt,&NPROW,&NPCOL,&myrow,&mycol);
    get_blocking_factors(NPROW,NPCOL,max(m,k),max(k,n),&MB,&NB);
    
    anr_loc      = c_numroc(lda,MB,myrow,0,NPROW);
    anc_loc      = c_numroc(na,NB,mycol,0,NPCOL);
    aROWMAX      = (lda-1)/(NPROW*MB)+1;
    aCOLMAX      = (na-1)/(NPCOL*NB)+1;
    aloc         = new double[anr_loc*anc_loc];

    distribute(myrow,mycol,lda,na,anr_loc,MB,NB,NPROW,NPCOL,aROWMAX,aCOLMAX,a,aloc);
    c_descinit(desca,lda,na,MB,NB,0,0,ictxt,anr_loc,&info);

    bnr_loc      = c_numroc(ldb,MB,myrow,0,NPROW);
    bnc_loc      = c_numroc(nb,NB,mycol,0,NPCOL);
    bROWMAX      = (ldb-1)/(NPROW*MB)+1;
    bCOLMAX      = (nb-1)/(NPCOL*NB)+1;
    bloc         = new double[bnr_loc*bnc_loc];

    distribute(myrow,mycol,ldb,nb,bnr_loc,MB,NB,NPROW,NPCOL,bROWMAX,bCOLMAX,b,bloc);
    c_descinit(descb,ldb,nb,MB,NB,0,0,ictxt,bnr_loc,&info);

    cnr_loc      = c_numroc(ldc,MB,myrow,0,NPROW);
    cnc_loc      = c_numroc(nc,NB,mycol,0,NPCOL);
    cROWMAX      = (ldc-1)/(NPROW*MB)+1;
    cCOLMAX      = (nc-1)/(NPCOL*NB)+1;
    cloc         = new double[cnr_loc*cnc_loc];

    distribute(myrow,mycol,ldc,nc,cnr_loc,MB,NB,NPROW,NPCOL,cROWMAX,cCOLMAX,c,cloc);
    c_descinit(descc,ldc,nc,MB,NB,0,0,ictxt,cnr_loc,&info);
    
    c_pdgemm(transa,transb,m,n,k,alpha,aloc,1,1,desca,bloc,1,1,descb,beta,cloc,1,1,descc);

    receive(myrow,mycol,ldc,nc,cnr_loc,cnc_loc,MB,NB,NPROW,NPCOL,cROWMAX,cCOLMAX,c,cloc,comm);

    Cblacs_gridexit(ictxt);

    delete[] aloc;
    delete[] bloc;
    delete[] cloc;

}

/************************************************************************************************/

void pzgemm_driver(char transa,char transb,int m,int n,int k,CPX alpha,CPX *a,\
                   int lda,CPX *b,int ldb,CPX beta,CPX *c,int ldc,int NPROW,\
                   int NPCOL,MPI_Comm comm)
{
    int ictxt,MB,NB,myrow,mycol,info;
    int anr_loc,anc_loc,aROWMAX,aCOLMAX,na;
    int bnr_loc,bnc_loc,bROWMAX,bCOLMAX,nb;
    int cnr_loc,cnc_loc,cROWMAX,cCOLMAX,nc;
    int desca[9],descb[9],descc[9];
    CPX *aloc,*bloc,*cloc;

    if(transa == 'N'){
        na       = k;
    }else{
        na       = m;
    }

    if(transb == 'N'){
        nb       = n;
        nc       = n;
    }else{
        nb       = k;
        nc       = k;
    }
    
    ictxt        = comm;
    Cblacs_gridinit(&ictxt,"Row-major",NPROW,NPCOL);
    Cblacs_gridinfo(ictxt,&NPROW,&NPCOL,&myrow,&mycol);
    get_blocking_factors(NPROW,NPCOL,max(m,k),max(k,n),&MB,&NB);
    
    anr_loc      = c_numroc(lda,MB,myrow,0,NPROW);
    anc_loc      = c_numroc(na,NB,mycol,0,NPCOL);
    aROWMAX      = (lda-1)/(NPROW*MB)+1;
    aCOLMAX      = (na-1)/(NPCOL*NB)+1;
    aloc         = new CPX[anr_loc*anc_loc];

    zdistribute(myrow,mycol,lda,na,anr_loc,MB,NB,NPROW,NPCOL,aROWMAX,aCOLMAX,a,aloc);
    c_descinit(desca,lda,na,MB,NB,0,0,ictxt,anr_loc,&info);

    bnr_loc      = c_numroc(ldb,MB,myrow,0,NPROW);
    bnc_loc      = c_numroc(nb,NB,mycol,0,NPCOL);
    bROWMAX      = (ldb-1)/(NPROW*MB)+1;
    bCOLMAX      = (nb-1)/(NPCOL*NB)+1;
    bloc         = new CPX[bnr_loc*bnc_loc];

    zdistribute(myrow,mycol,ldb,nb,bnr_loc,MB,NB,NPROW,NPCOL,bROWMAX,bCOLMAX,b,bloc);
    c_descinit(descb,ldb,nb,MB,NB,0,0,ictxt,bnr_loc,&info);

    cnr_loc      = c_numroc(ldc,MB,myrow,0,NPROW);
    cnc_loc      = c_numroc(nc,NB,mycol,0,NPCOL);
    cROWMAX      = (ldc-1)/(NPROW*MB)+1;
    cCOLMAX      = (nc-1)/(NPCOL*NB)+1;
    cloc         = new CPX[cnr_loc*cnc_loc];

    zdistribute(myrow,mycol,ldc,nc,cnr_loc,MB,NB,NPROW,NPCOL,cROWMAX,cCOLMAX,c,cloc);
    c_descinit(descc,ldc,nc,MB,NB,0,0,ictxt,cnr_loc,&info);
    
    c_pzgemm(transa,transb,m,n,k,alpha,aloc,1,1,desca,bloc,1,1,descb,beta,cloc,1,1,descc);

    zreceive(myrow,mycol,ldc,nc,cnr_loc,cnc_loc,MB,NB,NPROW,NPCOL,cROWMAX,cCOLMAX,c,cloc,comm);

    Cblacs_gridexit(ictxt);

    delete[] aloc;
    delete[] bloc;
    delete[] cloc;

}

/************************************************************************************************/

void pdsyev_driver(char jobz,char uplo,int n,double *a, double *lambda,double *z,int NPROW,\
                   int NPCOL,MPI_Comm comm)
{

    int ictxt,MB,NB,myrow,mycol,info;
    int anr_loc,anc_loc,aROWMAX,aCOLMAX,lda;
    int znr_loc,znc_loc,zROWMAX,zCOLMAX,ldz;
    int desca[9],descz[9];
    int lwork;
    double *aloc,*zloc;
    double *work;

    ictxt        = comm;
    Cblacs_gridinit(&ictxt,"Row-major",NPROW,NPCOL);
    Cblacs_gridinfo(ictxt,&NPROW,&NPCOL,&myrow,&mycol);
    get_blocking_factors(NPROW,NPCOL,n,n,&MB,&NB);
    
    lda          = n;
    anr_loc      = c_numroc(lda,MB,myrow,0,NPROW);
    anc_loc      = c_numroc(n,NB,mycol,0,NPCOL);
    aROWMAX      = (lda-1)/(NPROW*MB)+1;
    aCOLMAX      = (n-1)/(NPCOL*NB)+1;
    aloc         = new double[anr_loc*anc_loc];

    distribute(myrow,mycol,lda,n,anr_loc,MB,NB,NPROW,NPCOL,aROWMAX,aCOLMAX,a,aloc);
    c_descinit(desca,lda,n,MB,NB,0,0,ictxt,anr_loc,&info);

    ldz          = n;
    znr_loc      = c_numroc(ldz,MB,myrow,0,NPROW);
    znc_loc      = c_numroc(n,NB,mycol,0,NPCOL);
    zROWMAX      = (ldz-1)/(NPROW*MB)+1;
    zCOLMAX      = (n-1)/(NPCOL*NB)+1;
    zloc         = new double[znr_loc*znc_loc];

    distribute(myrow,mycol,ldz,n,znr_loc,MB,NB,NPROW,NPCOL,zROWMAX,zCOLMAX,z,zloc);
    c_descinit(descz,ldz,n,MB,NB,0,0,ictxt,znr_loc,&info);

    lwork        = 20*n;
    work         = new double[lwork];

    c_pdsyev(jobz,uplo,n,aloc,1,1,desca,lambda,z,1,1,descz,work,lwork,&info);

    receive(myrow,mycol,ldz,n,znr_loc,znc_loc,MB,NB,NPROW,NPCOL,zROWMAX,zCOLMAX,z,zloc,comm);

    Cblacs_gridexit(ictxt);

    delete[] aloc;
    delete[] zloc;
    delete[] work;

}

/************************************************************************************************/

void get_blocking_factors(int NPROW,int NPCOL,int m,int n,int* MB,int* NB)
{
    int row_per_processor    = m/NPROW;
    int block_per_rprocessor = row_per_processor/ROWBLOCK;
    int mb                   = row_per_processor/block_per_rprocessor;

    int col_per_processor    = n/NPCOL;
    int block_per_cprocessor = col_per_processor/COLBLOCK;
    int nb                   = col_per_processor/block_per_cprocessor;

    *MB                      = (mb+nb)/2;
    *NB                      = (mb+nb)/2;
}

/************************************************************************************************/
