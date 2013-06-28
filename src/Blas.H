#ifndef __BLAS
#define __BLAS

#include "mpi.h"
#include <iostream>
#include "Types.H"

extern "C" {
//Blas
    void fortran_name(dcopy,DCOPY)(int *n,double *dx,int *incx,double *dy,int *incy);
    void fortran_name(daxpy,DAXPY)(int *n,double *da,double *dx,int *incx,double *dy,\
                                      int *incy);
    double fortran_name(dnrm2,DNRM2)(int *n,double *dx,int *incx);
    double fortran_name(dznrm2,DZNRM2)(int *n,CPX *dx,int *incx);
    void fortran_name(dscal,DSCAL)(int *n,double *da,double *dx,int *incx);
    void fortran_name(zscal,ZSCAL)(int *n,CPX *da,CPX *dx,int *incx);
    void fortran_name(dgemm,DGEMM)(char *transa, char *transb, int *m, int *n, int *k, \
                                   double *alpha, double *a, int *lda, double *b, int *ldb, \
                                   double *beta, double *c, int *ldc);
    void fortran_name(dsymm,DSYMM)(char *side, char *uplo, int *m, int *n, double *alpha, \
                                   double *a, int *lda, double *b, int *ldb, double *beta, \
                                   double *c, int *ldc);
    void fortran_name(zgemm,ZGEMM)(char *transa, char *transb, int *m, int *n, int *k, \
                                   CPX *alpha, CPX *a, int *lda, CPX *b, int *ldb, \
                                   CPX *beta, CPX *c, int *ldc);
    void fortran_name(zhemm,ZHEMM)(char *side, char *uplo, int *m, int *n, CPX *alpha, \
                                   CPX *a, int *lda, CPX *b, int *ldb, CPX *beta, \
                                   CPX *c, int *ldc);
    void fortran_name(dgemv,DGEMV)(char *trans, int *m, int *n, double *alpha, double *a, \
                                   int *lda, double *x, int *incx, double *beta, double *y, \
                                   int *incy);
    
    void fortran_name(zgemv,ZGEMV)(char *trans, int *m, int *n, CPX *alpha, CPX *a, \
                                   int *lda, CPX *x, int *incx, CPX *beta, CPX *y, \
                                   int *incy);

    double fortran_name(ddot,DDOT)(int *n,double *x,int *incx,double *y,int *incy);
    CPX fortran_name(zdotc,ZDOTC)(int *n,CPX *x,int *incx,CPX *y,int *incy);
    void fortran_name(zcopy,ZCOPY)(int *n,CPX *dx,int *incx,CPX *dy,int *incy);
    void fortran_name(zaxpy,ZAXPY)(int *n, CPX *alpha, CPX *x, int *incx, CPX *y, int *incy);
    double fortran_name(dasum,DASUM)(int *n,double *dx,int *incx);
// SAB
    double fortran_name(dzasum,DZASUM)(int *n,CPX *dx,int *incx);

//Lapack
    void fortran_name(dgetrf,DGETRF)(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
    void fortran_name(dgetrs,DGETRS)(char *trans, int *n, int *nrhs, double *a, int *lda, \
                                     int *ipiv, double *b, int *ldb, int *info);
    void fortran_name(zgetrf,ZGETRF)(int *m, int *n, CPX *a, int *lda, int *ipiv, int *info);
    void fortran_name(zgetrs,ZGETRS)(char *trans, int *n, int *nrhs, CPX *a, int *lda, \
                                     int *ipiv, CPX *b, int *ldb, int *info);
    void fortran_name(zgetri,ZGETRI)(int *n,CPX *a,int *lda,int *ipiv,CPX *work,int *lwork,\
                                     int *info);
    void fortran_name(dgeev,DGEEV)(char *jobvl, char *jobvr, int *n, double *a, int *lda, \
                                   double *wr, double *wi, double *vl, int *ldvl, double *vr, \
                                   int *ldvr, double *work, int *lwork, int *info);
    void fortran_name(dsyev,DSYEV)(char *JOBZ,char *UPLO,int *N,double *A,int *LDA,double *W,\
				   double *WORK,int *LWORK,int *INFO);
    void fortran_name(dggev,DGGEV)(char* jobvl, char *jobvr, int *n, double *a, int *lda, \
                                   double *b, int *ldb, double *alphar, double *alphai, \
                                   double *beta, double *vl, int *ldvl, double *vr, int *ldvr, \
                                   double *work, int *lwork, int *info);
    void fortran_name(zgeev,ZGEEV)(char *jobvl, char *jobvr, int *n, CPX *a, int *lda, \
                                   CPX *w, CPX *vl, int *ldvl, CPX *vr, int *ldvr, CPX *work, \
                                   int *lwork, double *rwork, int *info);
    void fortran_name(zheev,ZHEEV)(char *jobvl, char *uplo, int *n, CPX *a, int *lda, \
                                   double *w, CPX *work, int *lwork, double *rwork, int *info);
    void fortran_name(dgetri,DGETRI)(int *n,double *a,int *lda,int *ipiv,double *work,int *lwork,\
                                     int *info);
    void fortran_name(dsytri,DSYTRI)(char *uplo,int *n,double *a,int *lda,int *ipiv,double *work,\
                                     int *info);
    void fortran_name(zhetrf,ZHETRF)(char *uplo,int *n,CPX *a,int *lda,int *ipiv,CPX *work,\
                                     int *lwork,int *info);
    void fortran_name(zhetri,ZHETRI)(char *uplo,int *n,CPX *a,int *lda,int *ipiv,CPX *work,\
                                     int *info);
    void fortran_name(zhetrs,ZHETRS)(char *uplo,int *n,int *nrhs,CPX *a,int *lda,int *ipiv,\
                                     CPX *b,int *ldb,int *info);
    void fortran_name(dsysv,DSYSV)(char *uplo,int *n,int *nrhs,double *a,int *lda,int *ipiv,\
                                   double *b,int *ldb,double *work,int *lwork,int *info);
    void fortran_name(dsytrf,DSYTRF)(char *uplo,int *n,double *a,int *lda,int *ipiv,double *work,\
                                     int *lwork,int *info);
    void fortran_name(dsytrs,DSYTRS)(char *uplo,int *n,int *nrhs,double *a,int *lda,int *ipiv,\
                                     double *b,int *ldb,int *info);
    void fortran_name(dstebz,DSTEBZ)(char *range,char *order,int *iter,double *vl,double *vu,int *il,int *iu,\
                                     double *abstol,double *diag,double *offd,int *neval,int *nsplit,\
                                     double *eval,int *iblock,int *isplit,double *work,int *iwork,int *info);
    void fortran_name(zlarnv,ZLARNV)(int*,int*,int*,CPX*);
// SAB
    void fortran_name(dsygv,DSYGV)(int *itype,char *jobz,char *uplo,int *n,double *a,int *lda,double *b,int *ldb,\
                    double *w,double *work,int *lwork,int *info);
    void fortran_name(zhegv,ZHEGV)(int *itype,char *jobz,char *uplo,int *n,CPX *a,int *lda,CPX *b,int *ldb,\
                    double *w,CPX *work,int *lwork,double *rwork,int *info);
}

/*My Blas*******************************************************************************************/

inline void c_icopy(int n,int *dx,int incx,int *dy,int incy)
{
    int i;
    
    for(i=0;i<n;i++) dy[i*incy] = dx[i*incx];
}

/*Blas*******************************************************************************************/
    
inline void c_dcopy(int n,double *dx,int incx,double *dy,int incy)
{
    fortran_name(dcopy,DCOPY)(&n,dx,&incx,dy,&incy);
}

/************************************************************************************************/

inline void c_daxpy(int n,double da,double *dx,int incx,double *dy,int incy)
{
    fortran_name(daxpy,DAXPY)(&n,&da,dx,&incx,dy,&incy);
}

/************************************************************************************************/

inline double c_dnrm2(int n,double* dx,int incx)
{
    return fortran_name(dnrm2,DNRM2)(&n,dx,&incx);
}

/************************************************************************************************/

inline double c_dznrm2(int n,CPX* dx,int incx)
{
    return fortran_name(dznrm2,DZNRM2)(&n,dx,&incx);
}

/************************************************************************************************/

inline void c_dscal(int n,double da,double *dx,int incx)
{
    fortran_name(dscal,DSCAL)(&n,&da,dx,&incx);
}

/************************************************************************************************/

inline void c_zscal(int n,CPX da,CPX *dx,int incx)
{
    fortran_name(zscal,ZSCAL)(&n,&da,dx,&incx);
}

/************************************************************************************************/

inline void c_dgemm(char transa, char transb, int m, int n, int k, double alpha, double *a, \
                    int lda, double *b, int ldb, double beta, double *c, int ldc)
{
    fortran_name(dgemm,DGEMM)(&transa,&transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
}

/************************************************************************************************/

inline void c_dsymm(char side,char uplo, int m, int n, double alpha, double *a, int lda, \
                    double *b, int ldb, double beta, double *c, int ldc)
{
    fortran_name(dsymm,DSYMM)(&side,&uplo,&m,&n,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
}

/************************************************************************************************/

inline void c_zgemm(char transa, char transb, int m, int n, int k, CPX alpha, CPX *a, \
                    int lda, CPX *b, int ldb, CPX beta, CPX *c, int ldc)
{
    fortran_name(zgemm,ZGEMM)(&transa,&transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
}

/************************************************************************************************/

inline void c_zhemm(char side, char uplo, int m, int n, CPX alpha, CPX *a, int lda, \
                    CPX *b, int ldb, CPX beta, CPX *c, int ldc)
{
    fortran_name(zhemm,ZHEMM)(&side,&uplo,&m,&n,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
}

/************************************************************************************************/

inline void c_dgemv(char transa, int m, int n, double alpha, double *a, int lda, double *x, \
                    int incx, double beta, double *y, int incy)
{
    fortran_name(dgemv,DGEMV)(&transa,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}

/************************************************************************************************/

inline void c_zgemv(char transa, int m, int n, CPX alpha, CPX *a, int lda, CPX *x, \
                    int incx, CPX beta, CPX *y, int incy)
{
    fortran_name(zgemv,ZGEMV)(&transa,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);

}

/***********************************************************************************************/
inline double c_ddot(int n, double *x, int incx, double *y, int incy)
{
    return fortran_name(ddot,DDOT)(&n,x,&incx,y,&incy);
}

/************************************************************************************************/

inline CPX c_zdotc(int n, CPX *x, int incx, CPX *y, int incy)
{
        double real,imag;

        real = c_ddot(n,(double*)x,2*incx,(double*)y,2*incy)+c_ddot(n,(double*)&x[0]+1,2*incx,(double*)&y[0]+1,2*incy);
        imag = -c_ddot(n,(double*)&x[0]+1,2*incx,(double*)y,2*incy)+c_ddot(n,(double*)x,2*incx,(double*)&y[0]+1,2*incy);

        return CPX(real,imag);

// FOR UNKNOWN REASON INTEL COMPILER SEGFAULTS WITH THAT
//        return fortran_name(zdotc,ZDOTC)(&n,x,&incx,y,&incy);
}

/************************************************************************************************/

inline void c_zcopy(int n,CPX *dx,int incx,CPX *dy,int incy)
{
    fortran_name(zcopy,ZCOPY)(&n,dx,&incx,dy,&incy);
}

/************************************************************************************************/

inline void c_zaxpy(int n, CPX alpha, CPX *x, int incx, CPX *y, int incy)
{
    fortran_name(zaxpy,ZAXPY)(&n,&alpha,x,&incx,y,&incy);
}

/************************************************************************************************/

inline double c_dasum(int n,double *dx,int incx)
{
    return fortran_name(dasum,DASUM)(&n,dx,&incx);
}

/*SAB********************************************************************************************/

inline double c_dzasum(int n,CPX *dx,int incx)
{
    return fortran_name(dzasum,DZASUM)(&n,dx,&incx);
}

/*SAB*END****************************************************************************************/

template <typename T,typename W>
inline void c_tcopy(int n,T *dx,int incx,W *dy,int incy);

template<>
inline void c_tcopy(int n,double *dx,int incx,double *dy,int incy)
{
    c_dcopy(n,dx,incx,dy,incy);
}

template<>
inline void c_tcopy(int n,CPX *dx,int incx,CPX *dy,int incy)
{
    c_zcopy(n,dx,incx,dy,incy);
}

template<>
inline void c_tcopy(int n,double *dx,int incx,CPX *dy,int incy)
{
    c_dcopy(n,dx,incx,(double*)dy,2*incy);
}

/************************************************************************************************/

template <typename T,typename W>
inline void c_taxpy(int n, T alpha, T *x, int incx, W *y, int incy);

template <>
inline void c_taxpy(int n, double alpha, double *x, int incx, double *y, int incy)
{
    c_daxpy(n,alpha,x,incx,y,incy);
}

template <>
inline void c_taxpy(int n, CPX alpha, CPX *x, int incx, CPX *y, int incy)
{
    c_zaxpy(n,alpha,x,incx,y,incy);
}

template <>
inline void c_taxpy(int n, double alpha, double *x, int incx, CPX *y, int incy)
{
    c_daxpy(n,alpha,x,incx,(double*)y,2*incy);
}

/************************************************************************************************/

template <typename T>
inline void c_tscal(int n,T da,T *dx,int incx);

template <>
inline void c_tscal(int n,double da,double *dx,int incx)
{
    c_dscal(n,da,dx,incx);
}

template <>
inline void c_tscal(int n,CPX da,CPX *dx,int incx)
{
    c_zscal(n,da,dx,incx);
}

/************************************************************************************************/

template <typename T>
inline void c_tgemm(char transa, char transb, int m, int n, int k, T alpha, T *a, \
                    int lda, T *b, int ldb, T beta, T *c, int ldc);

template <>
inline void c_tgemm(char transa, char transb, int m, int n, int k, double alpha, double *a, \
                    int lda, double *b, int ldb, double beta, double *c, int ldc)
{
    c_dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}

template <>
inline void c_tgemm(char transa, char transb, int m, int n, int k, CPX alpha, CPX *a, \
                    int lda, CPX *b, int ldb, CPX beta, CPX *c, int ldc)
{
    c_zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}

/*Lapack*****************************************************************************************/

inline void c_dgetrf(int m, int n, double *a, int lda, int *ipiv, int *info)
{
    fortran_name(dgetrf,DGETRF)(&m,&n,a,&lda,ipiv,info);
}

/************************************************************************************************/

inline void c_dgetrs(char transa, int n, int nrhs, double *a, int lda, int *ipiv, double *b, \
                   int ldb, int *info)
{
    fortran_name(dgetrs,DGETRS)(&transa,&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}

/************************************************************************************************/

inline void c_zgetrf(int m, int n, CPX *a, int lda, int *ipiv, int *info)
{
    fortran_name(zgetrf,ZGETRF)(&m,&n,a,&lda,ipiv,info);
}

/************************************************************************************************/

inline void c_zgetrs(char transa, int n, int nrhs, CPX *a, int lda, int *ipiv, CPX *b, \
                   int ldb, int *info)
{
    fortran_name(zgetrs,ZGETRS)(&transa,&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}

/************************************************************************************************/

inline void c_zgetri(int n,CPX *a,int lda,int *ipiv,CPX *work,int lwork,int *info)
{
    fortran_name(zgetri,ZGETRI)(&n,a,&lda,ipiv,work,&lwork,info);
}

/************************************************************************************************/

inline void c_dgeev(char jobvl, char jobvr, int n, double *a, int lda, double *wr, double *wi, \
                    double *vl, int ldvl, double *vr, int ldvr, double *work, int lwork, \
                    int *info)
{
    fortran_name(dgeev,DGEEV)(&jobvl,&jobvr,&n,a,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,info);
}

/************************************************************************************************/

inline void c_dsyev(char JOBZ,char UPLO,int N,double *A,int LDA,double *W,double *WORK,int LWORK,\
                    int *INFO)
{
  fortran_name(dsyev,DSYEV)(&JOBZ,&UPLO,&N,A,&LDA,W,WORK,&LWORK,INFO);
}

/************************************************************************************************/

inline void c_dggev(char jobvl, char jobvr, int n, double *a, int lda, double *b, int ldb, \
                    double *alphar, double *alphai, double *beta, double *vl, int ldvl, \
                    double *vr, int ldvr, double *work, int lwork, int *info)
{
    fortran_name(dggev,DGGEV)(&jobvl,&jobvr,&n,a,&lda,b,&ldb,alphar,alphai,beta,vl,&ldvl,vr,\
                              &ldvr,work,&lwork,info);
}

/************************************************************************************************/

inline void c_zgeev(char jobvl, char jobvr, int n, CPX *a, int lda, CPX *w, CPX *vl, int ldvl, \
                    CPX *vr, int ldvr, CPX *work, int lwork, double *rwork, int *info)
{
    fortran_name(zgeev,ZGEEV)(&jobvl,&jobvr,&n,a,&lda,w,vl,&ldvl,vr,&ldvr,work,&lwork,rwork,info);
}

/************************************************************************************************/

inline void c_zheev(char jobvl, char uplo, int n, CPX *a, int lda, double *w, CPX *work, \
		    int lwork, double *rwork, int *info)
{
    fortran_name(zheev,ZHEEV)(&jobvl,&uplo,&n,a,&lda,w,work,&lwork,rwork,info);
}

/************************************************************************************************/

inline void c_dgetri(int n,double *a,int lda,int *ipiv,double *work,int lwork,int *info)
{
    fortran_name(dgetri,DGETRI)(&n,a,&lda,ipiv,work,&lwork,info);
}

/************************************************************************************************/

inline void c_dsytri(char uplo,int n,double *a,int lda,int *ipiv,double *work,int *info)
{
    fortran_name(dsytri,DSYTRI)(&uplo,&n,a,&lda,ipiv,work,info);
}

/************************************************************************************************/

inline void c_zhetrf(char uplo,int n,CPX *a,int lda,int *ipiv,CPX *work,int lwork,int *info)
{
    fortran_name(zhetrf,ZHETRF)(&uplo,&n,a,&lda,ipiv,work,&lwork,info);
}

/************************************************************************************************/

inline void c_zhetri(char uplo,int n,CPX *a,int lda,int *ipiv,CPX *work,int *info)
{
    fortran_name(zhetri,ZHETRI)(&uplo,&n,a,&lda,ipiv,work,info);
}

/************************************************************************************************/

inline void c_zhetrs(char uplo,int n,int nrhs,CPX *a,int lda,int *ipiv,CPX *b,int ldb,\
                     int *info)
{
    fortran_name(zhetrs,ZHETRS)(&uplo,&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}

/************************************************************************************************/

inline void c_dsysv(char uplo,int n,int nrhs,double *a,int lda,int *ipiv,double *b,int ldb,\
                    double *work,int lwork,int *info)
{
    fortran_name(dsysv,DSYSV)(&uplo,&n,&nrhs,a,&lda,ipiv,b,&ldb,work,&lwork,info);
}

/************************************************************************************************/

inline void c_dsytrf(char uplo,int n,double *a,int lda,int *ipiv,double *work,int lwork,int *info)
{
    fortran_name(dsytrf,DSYTRF)(&uplo,&n,a,&lda,ipiv,work,&lwork,info);
}

/************************************************************************************************/

inline void c_dsytrs(char uplo,int n,int nrhs,double *a,int lda,int *ipiv,double *b,int ldb,\
                     int *info)
{
    fortran_name(dsytrs,DSYTRS)(&uplo,&n,&nrhs,a,&lda,ipiv,b,&ldb,info);
}

/************************************************************************************************/

inline void c_dstebz(char *range,char *order,int *iter,double *vl,double *vu,int *il,int *iu,\
                     double *abstol,double *diag,double *offd,int *neval,int *nsplit,\
                     double *eval,int *iblock,int *isplit,double *work,int *iwork,int *info)
{
    fortran_name(dstebz,DSTEBZ)(range,order,iter,vl,vu,il,iu,abstol,diag,offd,neval,nsplit,eval,\
	                            iblock,isplit,work,iwork,info);
}

/************************************************************************************************/

inline void c_zlarnv(int *idist,int *iseed,int n,CPX *x)
{
	fortran_name(zlarnv,ZLARNV)(idist,iseed,&n,x);
}

/*SAB********************************************************************************************/

inline void c_dsygv(int itype,char jobz,char uplo,int n,double *a,int lda,double *b,int ldb,\
                    double *w,double *work,int lwork,int *info)
{
    fortran_name(dsygv,DSYGV)(&itype,&jobz,&uplo,&n,a,&lda,b,&ldb,w,work,&lwork,info);
}

/************************************************************************************************/

inline void c_zhegv(int itype,char jobz,char uplo,int n,CPX *a,int lda,CPX *b,int ldb,\
                    double *w,CPX *work,int lwork,double *rwork,int *info)
{
    fortran_name(zhegv,ZHEGV)(&itype,&jobz,&uplo,&n,a,&lda,b,&ldb,w,work,&lwork,rwork,info);
}

/*SAB*END****************************************************************************************/

#endif




