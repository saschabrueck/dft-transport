#include "Arpack.H"


extern "C" {

    int fortran_name(dnaupd,DNAUPD)(int *ido,char *bmat,int *n,char *which,int *nev,double *tol,\
                                    double *resid,int *ncv,double *v,int *ldv,int iparam[11],\
                                    int ipntr[14],double *workd,double *workl,int *lworkl,\
                                    int *info);

    int fortran_name(dneupd,DNEUPD)(int *rvec,char *howmny,int *select,double *dr,double* di,\
                                    double *Z,int *ldz,double *sigmar,double *sigmai,double *workev,\
                                    char *bmat,int *n,char *which,int *nev,double *tol,\
                                    double *resid,int *ncv,double *v,int *ldv,int iparam[11],\
                                    int ipntr[14],double *workd,double *workl,int *lworkl,\
                                    double *rwork,int *info);

    int fortran_name(znaupd,ZNAUPD)(int *ido,char *bmat,int *n,char *which,int *nev,double *tol,\
                                    CPX *resid,int *ncv,CPX *v,int *ldv,int iparam[11],\
                                    int ipntr[14],CPX *workd,CPX *workl,int *lworkl,\
                                    double *rwork,int *info);

    int fortran_name(zneupd,ZNEUPD)(int *rvec,char *howmny,int *select,CPX *evals,CPX *evecs,\
                                    int *ldz,CPX *sigma,CPX *workev,char *bmat,int *n,\
                                    char *which,int *nev,double *tol,CPX *resid,int *ncv,\
                                    CPX *v,int *ldv,int iparam[11],int ipntr[14],CPX *workd,\
                                    CPX *workl,int *lworkl,double *rwork,int *info);

}



Arpack::Arpack()
{
}

/************************************************************************************************/

Arpack::~Arpack()
{
}

/************************************************************************************************/

bool Arpack::eigs(CPX *Zcpx,CPX *Dcpx,double *A,int ndof,int n,int nev,double &inptime)
{
    double sabtime;
    
    int  ido        = 0;
    char bmat[1]    = {'I'};
    char which[2]   = {'L','M'};
    double tol      = 1.0e-15;
    double* resid   = new double[n];
    int  ncv        = min(8*nev,n-1);
    double* v       = new double[n*ncv];
    int  ldv        = n;
    int iparam[11];
    iparam[0]       = 1;
    iparam[2]       = 100;
//    iparam[3]       = 1;
    iparam[6]       = 1;
    int ipntr[14];
    double* workd   = new double[3*n];
    int lworkl      = 3*ncv*ncv + 6 * ncv;
    double* workl   = new double[lworkl];
    double* rwork   = new double[ncv];
    int info        = 1;
    
    srand(0);
    for (int ir=0;ir<n;ir++) {
        resid[ir]=(double)rand()/RAND_MAX;
    }

    int iternum=0;

    int worldrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);

    while(ido != 99) {

        ++iternum;
        
        fortran_name(dnaupd,DNAUPD)(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, \
                                    v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

        if (info) cout << "WARNING INFO " << info << " AT ITERATION " << iternum << endl;
        
        if(ido == -1 || ido == 1) {
            sabtime=get_time(0.0);

if (n/ndof==1) {
            c_dgemv('N',n,n,-1.0,A,n,&workd[ipntr[0] - 1],1,0.0,&workd[ipntr[1] - 1],1);
} else if (n/ndof==2) {
            c_dcopy(ndof,&workd[ipntr[0] - 1],1,&workd[ipntr[1] - 1],1);
            c_dgemv('N',ndof,n,-1.0,A,ndof,&workd[ipntr[0] - 1],1,0.0,&workd[ipntr[1] - 1 + ndof],1);
            c_daxpy(ndof,-1.0,&workd[ipntr[1] - 1 + ndof],1,&workd[ipntr[1] - 1],1);
} else {

            c_dcopy(n - ndof,&workd[ipntr[0] - 1 + ndof],1,&workd[ipntr[1] - 1],1);
//            c_dgemv('N',ndof,n,1.0,A,ndof,&workd[ipntr[0] - 1],1,0.0,&workd[ipntr[1] - 1 + n - ndof],1);
            c_dgemv('T',n,ndof,-1.0,A,n,&workd[ipntr[0] - 1],1,0.0,&workd[ipntr[1] - 1 + n - ndof],1);

}
            
            inptime+=get_time(sabtime);
        }else if(ido==99){
            cout << worldrank << " NEEDED " << iparam[2] << " ITERATIONS AND " << iternum << " MATRIX VECTOR MULTIPLICATION" << endl;
            break;
        }
    }

    int nev2       = max(iparam[4],nev);
    double sigmar;
    double sigmai;
    int   rvec     = 0;
    char  howmny[] = {'A'};
    int*  select   = new int[ncv];
    int   ldz      = n;
    double* workev = new double[3*ncv];
    double *Di     = new double[nev2+4];
    double *D      = new double[nev2+4];
    double *Z;

    if(Zcpx!=NULL) {
        rvec=1;
        Z = new double[n*(nev2+4)];
    }

    switch(info) {
    case 0:
        
        if(iparam[4]<nev){
            cout << "the eigensolver found " << iparam[4] << " of " << nev << " eigenvalues\n";
        }

        fortran_name(dneupd,DNEUPD)(&rvec, howmny, select, D, Di, Z, &ldz, &sigmar, &sigmai, workev,\
                                    bmat, &n, which, &nev2, &tol, resid, &ncv, v, &ldv, \
                                    iparam, ipntr, workd, workl, &lworkl, rwork, &info);
        break;
    case 1:
        printf("eigensolver quit after reaching max iterations\n");
        break;
    case 3:
        printf("eigensolver quit because he could not apply shifts\n");
        break;
    }

    delete[] select;
    delete[] workev;

    delete[] resid;
    delete[] v;
    delete[] workd;
    delete[] workl;
    delete[] rwork;

    for (int i=0;i<nev;i++) {
        Dcpx[i]=CPX(D[i],Di[i]);
    }

    if(Zcpx!=NULL) {
        int jev=0;
        while (jev<nev) {
            if (!Di[jev]) {
                for (int i=0;i<n;i++) {
                    Zcpx[jev*n+i]=CPX(Z[jev*n+i],0.0);
                }
                jev++;
            }
            else if (D[jev]==D[jev+1]) {
                for (int i=0;i<n;i++) {
                    Zcpx[ jev   *n+i]=CPX(Z[jev*n+i], Z[(jev+1)*n+i]);
                    Zcpx[(jev+1)*n+i]=CPX(Z[jev*n+i],-Z[(jev+1)*n+i]);
                }
                jev+=2;
            }
            else {
                printf("%i HAS ERROR IN RESULT OF ARPACK AT %i\nIMAGEIG IS %f AND ONE MORE %f\nAND FOUND ARE %i\n",worldrank,jev,Di[jev],Di[jev+1],iparam[4]);
                jev++;
            }
        }
        delete[] Z;
    }
    delete[] Di;
    delete[] D;

    return true;
}

/************************************************************************************************/

bool Arpack::eigs(CPX *Z,CPX *D,CPX *A,int ndof,int n,int nev,double &inptime)
{
    double sabtime;

    int  ido        = 0;
    char bmat[1]    = {'I'};
    char which[2]   = {'L', 'M'};
    double tol      = 1.0e-15;
    int  ncv        = min(4*nev,n-1);
    CPX* resid      = new CPX[n];
    CPX* v          = new CPX[n*ncv];
    int  ldv        = n;
    int iparam[11];
    iparam[0]       = 1;
    iparam[2]       = 100;
//    iparam[3]       = 1;
    iparam[6]       = 1;
    int ipntr[14];
    int lworkl      = 3*ncv*ncv + 6 * ncv;
    CPX* workd      = new CPX[3*n];
    CPX* workl      = new CPX[lworkl];
    double* rwork   = new double[ncv];
    int info        = 1;

    srand(0);
    for (int ir=0;ir<n;ir++) {
        resid[ir]=CPX((double)rand()/RAND_MAX,(double)rand()/RAND_MAX);
    }

    int iternum=0;

    int worldrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);

    while(ido != 99) {

        ++iternum;

#ifdef _AIX

        znaupd(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, \
                                    v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
#else

        fortran_name(znaupd,ZNAUPD)(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, \
                                    v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

#endif
        if (info) cout << "WARNING INFO " << info << " AT ITERATION " << iternum << endl;

        if(ido == -1 || ido == 1) {
            sabtime=get_time(0.0);

            c_zcopy(n - ndof,&workd[ipntr[0] - 1 + ndof],1,&workd[ipntr[1] - 1],1);
//            c_zgemv('N',ndof,n,CPX(1.0,0.0),A,ndof,&workd[ipntr[0] - 1],1,CPX(0.0,0.0),&workd[ipntr[1] - 1 + n - ndof],1);
            c_zgemv('T',n,ndof,CPX(-1.0,0.0),A,n,&workd[ipntr[0] - 1],1,CPX(0.0,0.0),&workd[ipntr[1] - 1 + n - ndof],1);

            inptime+=get_time(sabtime);
        }else if(ido==99){
            cout << worldrank << " NEEDED " << iparam[2] << " ITERATIONS AND " << iternum << " MATRIX VECTOR MULTIPLICATION" << endl;
            break;
        }
        
    }

    CPX sigma;
    int   rvec     = 0;
    char  howmny[] = {'A'};
    int*  select   = new int[ncv];
    CPX* workev    = new CPX[2*ncv];
    int   ldz      = n;

    if(Z!=NULL) rvec = 1;

    switch(info) {
    case 0:
        
        if(iparam[4]<nev){
            cout << "the eigensolver found " << iparam[4] << " of " << nev << " eigenvalues\n";
        }

#ifdef _AIX

        zneupd(&rvec, howmny, select, D, Z, &ldz, &sigma, workev,\
                                    bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, \
                                    iparam, ipntr, workd, workl, &lworkl, rwork, &info);
#else

        fortran_name(zneupd,ZNEUPD)(&rvec, howmny, select, D, Z, &ldz, &sigma, workev,\
                                    bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, \
                                    iparam, ipntr, workd, workl, &lworkl, rwork, &info);

#endif
        if (info) printf("problem in postproc\n");
        break;
    case 1:
        printf("eigensolver quit after reaching max iterations\n");
        break;
    case 3:
        printf("eigensolver quit because he could not apply shifts\n");
        break;
    }

    delete[] select;
    delete[] workev;
    delete[] resid;
    delete[] v;
    delete[] workd;
    delete[] workl;
    delete[] rwork;

    return true;
}

/************************************************************************************************/
