#include "Arpack.H"


extern "C" {

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

bool Arpack::eigs(CPX *Z,CPX *D,CPX *A,int n,int nev)
{

    int  ido        = 0;
    char bmat[1]    = {'I'};
    char which[2]   = {'L', 'M'};
    double tol      = 1.0e-10;
    int  ncv        = min(2*nev,n-1);
    CPX* resid      = new CPX[n];
    CPX* v          = new CPX[n*ncv];
    int  ldv        = n;
    int iparam[11];
    iparam[0]       = 1;
    iparam[2]       = 500;
    iparam[3]       = 1;
    iparam[6]       = 1;
    int ipntr[14];
    int lworkl      = 3*ncv*ncv + 6 * ncv;
    CPX* workd      = new CPX[3*n];
    CPX* workl      = new CPX[lworkl];
    double* rwork   = new double[ncv];
    int info        = 0;

    while(ido != 99) {

#ifdef _AIX

        znaupd(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, \
                                    v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
#else

        fortran_name(znaupd,ZNAUPD)(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, \
                                    v, &ldv, iparam, ipntr, workd, workl, &lworkl, rwork, &info);

#endif

        if(ido == -1 || ido == 1) {

            c_zgemv('N',n,n,CPX(1.0,0.0),A,n,&workd[ipntr[0] - 1],1,CPX(0.0,0.0),&workd[ipntr[1] - 1],1);
            
        }else if(ido==99){
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
