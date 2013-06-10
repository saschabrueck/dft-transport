#ifndef __SCALAPACK
#define __SCALAPACK

#ifdef _AIX

#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <time.h>
#include "Blas.H"
#include "Utilities.H"

#else

#include <mpi.h>
#include <iostream.h>
#include <stdio.h>
#include <fstream.h>
#include <time.h>
#include "Blas.H"
#include "Utilities.H"

#endif

extern "C" {
    
    void fortran_name(dgehrd,DGEHRD)(int*,int*,int*,double*,int*,double*,double*,int*,int*);
    void fortran_name(dlacpy,DLACPY)(char*,int*,int*,double*,int*,double*,int*);
    void fortran_name(dorghr,DORGHR)(int*,int*,int*,double*,int*,double*,double*,int*,int*);
    void fortran_name(dhseqr,DHSEQR)(char*,char*,int*,int*,int*,double*,int*,double*,double*,\
                                     double*,int*,double*,int*,int*);
    void fortran_name(dtrevc,DTREVC)(char*,char*,bool*,int*,double*,int*,double*,int*,double*,\
                                     int*,int*,int*,double*,int*);
    void fortran_name(ztrevc,ZTREVC)(char*,char*,bool*,int*,CPX*,int*,CPX*,int*,CPX*,\
                                     int*,int*,int*,CPX*,double*,int*);
    void fortran_name(dgebal,DGEBAL)(char*,int*,double*,int*,int*,int*,double*,int*);

    int fortran_name(numroc,NUMROC)(int*,int*,int*,int*,int*);
    int fortran_name(indxg2p,INDXG2P)(int*,int*,int*,int*,int*);
    void fortran_name(descinit,DESCINIT)(int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
    void fortran_name(pdsyev,PDSYEV)(char*,char*,int*,double*,int*,int*,int*,double*,\
				     double*,int*,int*,int*,double*,int*,int*);
    void fortran_name(pdgehrd,PDGEHRD)(int*,int*,int*,double*,int*,int*,int*,double*,\
                                       double*,int*,int*);
    void fortran_name(pzgehrd,PZGEHRD)(int*,int*,int*,CPX*,int*,int*,int*,CPX*,\
                                       CPX*,int*,int*);
    void fortran_name(pdlacpy,PDLACPY)(char*,int*,int*,double*,int*,int*,int*,double*,\
                                       int*,int*,int*);
    void fortran_name(pzlacpy,PZLACPY)(char*,int*,int*,CPX*,int*,int*,int*,CPX*,\
                                       int*,int*,int*);
    void fortran_name(pdlaset,PDLASET)(char*,int*,int*,double*,double*,double*,int*,int*,int*);
    void fortran_name(pzlaset,PZLASET)(char*,int*,int*,CPX*,CPX*,CPX*,int*,int*,int*);
    void fortran_name(pdlahqr,PDLAHQR)(bool*,bool*,int*,int*,int*,double*,int*,double*,double*,\
                                       int*,int*,double*,int*,double*,int*,int*,int*,int*);
    void fortran_name(pzlahqr,PZLAHQR)(bool*,bool*,int*,int*,int*,CPX*,int*,CPX*,\
                                       int*,int*,CPX*,int*,CPX*,int*,int*,int*,int*);
    void fortran_name(pdormhr,PDORMHR)(char*,char*,int*,int*,int*,int*,double*,int*,int*,int*,\
                                       double*,double*,int*,int*,int*,double*,int*,int*);
    void fortran_name(pzunmhr,PZUNMHR)(char*,char*,int*,int*,int*,int*,CPX*,int*,int*,int*,\
                                       CPX*,CPX*,int*,int*,int*,CPX*,int*,int*);
    void fortran_name(pztrevc,PZTREVC)(char*,char*,bool*,int*,CPX*,int*,CPX*,int*,CPX*,int*,\
                                       int*,int*,CPX*,double*,int*);
    void fortran_name(pdgemm,PDGEMM)(char*,char*,int*,int*,int*,double*,double*,int*,int*,int*,\
                                     double*,int*,int*,int*,double*,double*,int*,int*,int*);
    void fortran_name(pzgemm,PZGEMM)(char*,char*,int*,int*,int*,CPX*,CPX*,int*,int*,int*,\
                                     CPX*,int*,int*,int*,CPX*,CPX*,int*,int*,int*);
    void fortran_name(pzgetrf,PZGETRF)(int*,int*,CPX*,int*,int*,int*,int*,int*);
    void fortran_name(pzgetrs,PZGETRS)(char*,int*,int*,CPX*,int*,int*,int*,int*,CPX*,int*,\
                                       int*,int*,int*);

    void fortran_name(blacs_pinfo,BLACS_pinfo)(int*,int*);
    void fortran_name(blacs_get,BLACS_GET)(int*,int*,int*);
    void fortran_name(blacs_gridinit,BLACS_GRIDINIT)(int*,char*,int*,int*);
    void fortran_name(blacs_gridinfo,BLACS_GRIDINFO)(int*,int*,int*,int*,int*);
    void fortran_name(blacs_gridexit,BLACS_GRIDEXIT)(int*);
}

/*Lapack*****************************************************************************************/

inline void c_dgehrd(int N,int ILO,int IHI,double *A,int LDA,double *TAU,double *WORK,\
                     int LWORK,int *INFO)
{
    fortran_name(dgehrd,DGEHRD)(&N,&ILO,&IHI,A,&LDA,TAU,WORK,&LWORK,INFO);
}

/************************************************************************************************/

inline void c_dlacpy(char UPLO,int M,int N,double *A,int LDA,double *B,int LDB)
{
    fortran_name(dlacpy,DLACPY)(&UPLO,&M,&N,A,&LDA,B,&LDB);
}

/************************************************************************************************/

inline void c_dorghr(int N,int ILO,int IHI,double *A,int LDA,double *TAU,double *WORK,\
                     int LWORK,int *INFO)
{
    fortran_name(dorghr,DORGHR)(&N,&ILO,&IHI,A,&LDA,TAU,WORK,&LWORK,INFO);
}

/************************************************************************************************/

inline void c_dhseqr(char JOB,char COMPZ,int N,int ILO,int IHI,double *H,int LDH,double *WR,\
                     double *WI,double *Z,int LDZ,double *WORK,int LWORK,int *INFO)
{
    fortran_name(dhseqr,DHSEQR)(&JOB,&COMPZ,&N,&ILO,&IHI,H,&LDH,WR,WI,Z,&LDZ,WORK,&LWORK,INFO);
}

/************************************************************************************************/

inline void c_dtrevc(char SIDE,char HOWMNY,bool SELECT,int N,double *T,int LDT,double *VL,\
                     int LDVL,double *VR,int LDVR,int MM,int *M,double *WORK,int *INFO)
{
    fortran_name(dtrevc,DTREVC)(&SIDE,&HOWMNY,&SELECT,&N,T,&LDT,VL,&LDVL,VR,&LDVR,&MM,M,\
                                WORK,INFO);
}

/************************************************************************************************/

inline void c_ztrevc(char SIDE,char HOWMNY,bool SELECT,int N,CPX *T,int LDT,CPX *VL,\
                     int LDVL,CPX *VR,int LDVR,int MM,int *M,CPX *WORK,double *RWORK,int *INFO)
{
    fortran_name(ztrevc,ZTREVC)(&SIDE,&HOWMNY,&SELECT,&N,T,&LDT,VL,&LDVL,VR,&LDVR,&MM,M,\
                                WORK,RWORK,INFO);
}

/************************************************************************************************/

inline void c_dgebal(char JOB,int N,double *A,int LDA,int ILO,int IHI,double *SCALE,int *INFO )
{
    fortran_name(dgebal,DGEBAL)(&JOB,&N,A,&LDA,&ILO,&IHI,SCALE,INFO);
}

/*ScaLapack**************************************************************************************/

inline int c_numroc(int N,int NB,int IPROC,int ISRCPROC,int NPROCS)
{
    return fortran_name(numroc,NUMROC)(&N,&NB,&IPROC,&ISRCPROC,&NPROCS);
}

/************************************************************************************************/

inline int c_indxg2p(int INDXGLOB,int NB,int IPROC,int ISRCPROC,int NPROCS)
{
    return fortran_name(indxg2p,INDXG2P)(&INDXGLOB,&NB,&IPROC,&ISRCPROC,&NPROCS);
}

/************************************************************************************************/

inline void c_descinit(int *DESC,int M,int N,int MB,int NB,int IRSRC,int ICSRC,int ICTXT,\
                       int LLD, int *INFO)
{
    fortran_name(descinit,DESCINIT)(DESC,&M,&N,&MB,&NB,&IRSRC,&ICSRC,&ICTXT,&LLD,INFO);
}

/************************************************************************************************/

inline void c_pdsyev(char JOBZ,char UPLO,int N,double *A,int IA,int JA,int *DESCA,double *W,
		     double *Z,int IZ,int JZ,int *DESCZ,double *WORK,int LWORK,int *INFO)
{

    fortran_name(pdsyev,PDSYEV)(&JOBZ,&UPLO,&N,A,&IA,&JA,DESCA,W,Z,&IZ,&JZ,DESCZ,WORK,\
				&LWORK,INFO);

}

/************************************************************************************************/

inline void c_pdgehrd(int N,int ILO,int IHI,double *A,int IA,int JA,int *DESCA,double* TAU,\
                      double *WORK,int LWORK,int *INFO)
{
    fortran_name(pdgehrd,PDGEHRD)(&N,&ILO,&IHI,A,&IA,&JA,DESCA,TAU,WORK,&LWORK,INFO);
}

/************************************************************************************************/

inline void c_pzgehrd(int N,int ILO,int IHI,CPX *A,int IA,int JA,int *DESCA,CPX* TAU,\
                      CPX *WORK,int LWORK,int *INFO)
{
    fortran_name(pzgehrd,PZGEHRD)(&N,&ILO,&IHI,A,&IA,&JA,DESCA,TAU,WORK,&LWORK,INFO);
}

/************************************************************************************************/

inline void c_pdlacpy(char UPLO,int M,int N,double *A,int IA,int JA,int *DESCA,double *B,\
                      int IB,int JB, int *DESCB)
{
    fortran_name(pdlacpy,PDLACPY)(&UPLO,&M,&N,A,&IA,&JA,DESCA,B,&IB,&JB,DESCB);
}

/************************************************************************************************/

inline void c_pzlacpy(char UPLO,int M,int N,CPX *A,int IA,int JA,int *DESCA,CPX *B,\
                      int IB,int JB, int *DESCB)
{
    fortran_name(pzlacpy,PZLACPY)(&UPLO,&M,&N,A,&IA,&JA,DESCA,B,&IB,&JB,DESCB);
}

/************************************************************************************************/

inline void c_pdlaset(char UPLO,int M,int N,double ALPHA,double BETA,double *A,int IA,\
                      int JA,int *DESCA)
{
    fortran_name(pdlaset,PDLASET)(&UPLO,&M,&N,&ALPHA,&BETA,A,&IA,&JA,DESCA);
}

/************************************************************************************************/

inline void c_pzlaset(char UPLO,int M,int N,CPX ALPHA,CPX BETA,CPX *A,int IA,\
                      int JA,int *DESCA)
{
    fortran_name(pzlaset,PZLASET)(&UPLO,&M,&N,&ALPHA,&BETA,A,&IA,&JA,DESCA);
}

/************************************************************************************************/

inline void c_pdlahqr(bool WANTT,bool WANTZ,int N,int ILO,int IHI,double *A,int *DESCA,\
                      double *WR,double *WI,int ILOZ,int IHIZ,double *Z,int *DESCZ,double *WORK,\
                      int LWORK,int *IWORK,int ILWORK,int *INFO)
{
    fortran_name(pdlahqr,PDLAHQR)(&WANTT,&WANTZ,&N,&ILO,&IHI,A,DESCA,WR,WI,&ILOZ,&IHIZ,Z,DESCZ,\
                                  WORK,&LWORK,IWORK,&ILWORK,INFO);
}

/************************************************************************************************/

inline void c_pzlahqr(bool WANTT,bool WANTZ,int N,int ILO,int IHI,CPX *A,int *DESCA,\
                      CPX *W,int ILOZ,int IHIZ,CPX *Z,int *DESCZ,CPX *WORK,\
                      int LWORK,int *IWORK,int ILWORK,int *INFO)
{
    fortran_name(pzlahqr,PZLAHQR)(&WANTT,&WANTZ,&N,&ILO,&IHI,A,DESCA,W,&ILOZ,&IHIZ,Z,DESCZ,\
                                  WORK,&LWORK,IWORK,&ILWORK,INFO);
}

/************************************************************************************************/

inline void c_pdormhr(char SIDE,char TRANS,int M,int N,int ILO,int IHI,double *A,int IA,int JA,\
                      int *DESCA,double *TAU,double *C,int IC,int JC,int *DESCC,double *WORK,\
                      int LWORK,int *INFO)
{
    fortran_name(pdormhr,PDORMHR)(&SIDE,&TRANS,&M,&N,&ILO,&IHI,A,&IA,&JA,DESCA,\
                                  TAU,C,&IC,&JC,DESCC,WORK,&LWORK,INFO);
}

/************************************************************************************************/

inline void c_pzunmhr(char SIDE,char TRANS,int M,int N,int ILO,int IHI,CPX *A,int IA,int JA,\
                      int *DESCA,CPX *TAU,CPX *C,int IC,int JC,int *DESCC,CPX *WORK,\
                      int LWORK,int *INFO)
{
    fortran_name(pzunmhr,PZUNMHR)(&SIDE,&TRANS,&M,&N,&ILO,&IHI,A,&IA,&JA,DESCA,\
                                  TAU,C,&IC,&JC,DESCC,WORK,&LWORK,INFO);
}

/************************************************************************************************/

inline void c_pztrevc(char SIDE,char HOWMNY,bool SELECT,int N,CPX *T,int *DESCT,CPX *VL,\
                      int *DESCVL,CPX *VR,int *DESCVR,int MM,int *M,CPX *WORK,double *RWORK,\
                      int *INFO)
{
    fortran_name(pztrevc,PZTREVC)(&SIDE,&HOWMNY,&SELECT,&N,T,DESCT,VL,DESCVL,VR,DESCVR,\
                                  &MM,M,WORK,RWORK,INFO);
}

/************************************************************************************************/

inline void c_pdgemm(char transa,char transb,int m,int n,int k,double alpha,double *a,int ia,\
                     int ja,int *desc_a,double *b,int ib,int jb,int *desc_b,double beta,\
                     double *c,int ic,int jc,int *desc_c)
{
    fortran_name(pdgemm,PDGEMM)(&transa,&transb,&m,&n,&k,&alpha,a,&ia,&ja,desc_a,b,&ib,&jb,\
                                desc_b,&beta,c,&ic,&jc,desc_c);
}

/************************************************************************************************/

inline void c_pzgemm(char transa,char transb,int m,int n,int k,CPX alpha,CPX *a,int ia,\
                     int ja,int *desc_a,CPX *b,int ib,int jb,int *desc_b,CPX beta,\
                     CPX *c,int ic,int jc,int *desc_c)
{
    fortran_name(pzgemm,PZGEMM)(&transa,&transb,&m,&n,&k,&alpha,a,&ia,&ja,desc_a,b,&ib,&jb,\
                                desc_b,&beta,c,&ic,&jc,desc_c);
}

/************************************************************************************************/

inline void c_pzgetrf(int M,int N,CPX *A,int IA,int JA,int *DESCA,int *IPIV,int *INFO)
{
    fortran_name(pzgetrf,PZGETRF)(&M,&N,A,&IA,&JA,DESCA,IPIV,INFO);
}

/************************************************************************************************/

inline void c_pzgetrs(char TRANS,int N,int NRHS,CPX *A,int IA,int JA,int *DESCA,int *IPIV,\
                      CPX *B,int IB,int JB,int *DESCB,int *INFO)
{
    fortran_name(pzgetrs,PZGETRS)(&TRANS,&N,&NRHS,A,&IA,&JA,DESCA,IPIV,B,&IB,&JB,DESCB,INFO);
}

/*Blacs******************************************************************************************/

inline void Cblacs_pinfo(int *MYPNUM,int *NPROCS)
{
    fortran_name(blacs_pinfo,BLACS_pinfo)(MYPNUM,NPROCS);
}

/************************************************************************************************/

inline void Cblacs_get(int ICONTXT,int WHAT,int *VAL)
{
    fortran_name(blacs_get,BLASCS_GET)(&ICONTXT,&WHAT,VAL);
}

/************************************************************************************************/

inline void Cblacs_gridinit(int *ICONTXT,char *ORDER,int NPROW,int NPCOL)
{
    fortran_name(blacs_gridinit,BLACS_GRIDINIT)(ICONTXT,ORDER,&NPROW,&NPCOL);
}

/************************************************************************************************/

inline void Cblacs_gridinfo(int ICONTXT,int *NPROW,int *NPCOL,int *MYPROW,int *MYPCOL)
{
    fortran_name(blacs_gridinfo,BLACS_GRIDINFO)(&ICONTXT,NPROW,NPCOL,MYPROW,MYPCOL);
}

/************************************************************************************************/

inline void Cblacs_gridexit(int ICONTXT)
{
    fortran_name(blacs_gridexit,BLACS_GRIDEXIT)(&ICONTXT);
}

/*Own********************************************************************************************/

void my_dgeev(int n,double *a,int lda,double *wr,double *wi,double *vr,int ldvr,double *work,\
              int lwork,int *info);
void pdgeev(int n,double *a,int *desca,double *wr,double *wi,double *vr,int *descvr,\
            double *work,int lwork,MPI_Comm comm,int *info);
void pzgeev(int n,CPX *a,int *desca,CPX *w,CPX *vr,int *descvr,CPX *work,int lwork,\
            MPI_Comm comm,int *info);
void pdgeev_driver(int n, double *a, double *wr, double *wi,double *vr,int NPROW,int NPCOL,\
                   MPI_Comm comm,int *info);
void pzgeev_driver(int n, CPX *a, CPX *w,CPX *vr,int NPROW,int NPCOL,MPI_Comm comm,int *info);
void pdgemm_driver(char transa,char transb,int m,int n,int k,double alpha,double *a,\
                   int lda,double *b,int ldb,double beta,double *c,int ldc,int NPROW,\
                   int NPCOL,MPI_Comm comm);
void pzgemm_driver(char transa,char transb,int m,int n,int k,CPX alpha,CPX *a,\
                   int lda,CPX *b,int ldb,CPX beta,CPX *c,int ldc,int NPROW,\
                   int NPCOL,MPI_Comm comm);
void distribute(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
                int ROWMAX,int COLMAX,double *A,double *Aloc);
void zdistribute(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
                 int ROWMAX,int COLMAX,CPX *A,CPX *Aloc);
void receive(int myrow,int mycol,int M, int N,int NR_loc,int NC_loc,int MB,int NB,int NPROW,\
             int NPCOL,int ROWMAX,int COLMAX,double *A,double *Aloc,MPI_Comm comm);
void zreceive(int myrow,int mycol,int M, int N,int NR_loc,int NC_loc,int MB,int NB,int NPROW,\
              int NPCOL,int ROWMAX,int COLMAX,CPX *A,CPX *Aloc,MPI_Comm comm);
void insert_block(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
                  int ROWMAX,int COLMAX,double *A,double *Aloc);
void zinsert_block(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
                   int ROWMAX,int COLMAX,CPX *A,CPX *Aloc);
void get_hess(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
              int ROWMAX,int COLMAX,double *Aloc);
void get_zhess(int myrow,int mycol,int M, int N,int NR_loc,int MB,int NB,int NPROW,int NPCOL,\
              int ROWMAX,int COLMAX,CPX *Aloc);
void get_blocking_factors(int NPROW,int NPCOL,int m,int n,int* MB,int* NB);

#endif


