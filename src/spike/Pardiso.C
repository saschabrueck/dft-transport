#include "Pardiso.H"
#include <cstddef>
#include <omp.h>


extern "C" {

    // Idea: use the complex interface as an interpretation of a 2xdouble
    // interface. So we cast to double if we have a CPX matrix and set
    // the corresponding parameter mtype to the right value.
    //int fortran_name(pardiso,PARDISO)(void*,int*,int*, int*,int*,int*,
    //    double*,int*,int*,int*,int*,int*,int*,double*,double*,int*,double*);

    // This was correct for the newest pardiso:
    //void fortran_name(pardisoinit,PARDISOINIT)(void*,int*,int*,int*,
    //    double*,int*);
    //int fortran_name(pardiso,PARDISO)(void*,int*,int*, int*,int*,int*,
    //    CPX*,int*,int*,int*,int*,int*,int*,CPX*,CPX*,int*,double*);



//    // But intel doesn't ship that
//    void fortran_name(pardisoinit,PARDISOINIT)(void*,int*,int*);
//    int fortran_name(pardiso,PARDISO)(void*,int*,int*, int*,int*,int*,
//        CPX*,int*,int*,int*,int*,int*,int*,CPX*,CPX*,int*);
    
void pardiso( void*, int *, int *, int *,
                   int *, int *,       void *, int *,
                   int *, int *,    int *, int *,
                   int *,    void *,       void *, int * );

void PARDISO( void*, int *, int *, int *,
                   int *,    int *,    void *, int *,
                   int *,    int *, int *, int *,
                   int *,    void *,       void *, int * );

void pardisoinit( void*, int *, int * );

void PARDISOINIT( void*, int *, int * );

/*
*  Note: The pardiso_64 interface is not supported on IA-32 architecture.
*        If called on IA-32, error = -12 is returned.
*/

void pardiso_64( void*, long long int *, long long int *, long long int *,
                   long long int *, long long int *,          void *, long long int *,
                   long long int *, long long int *, long long int *, long long int *,
                   long long int *,          void *,          void *, long long int * );

void PARDISO_64( void*, long long int *, long long int *, long long int *,
                   long long int *, long long int *,          void *, long long int *,
                   long long int *, long long int *, long long int *, long long int *,
                   long long int *,          void *,          void *, long long int * );
    
    
    
}

namespace Pardiso {

/** \brief Function to solve a sparse linear system
 * 
 *  \param[in]   matrix   The coefficient matrix.
 *
 *  \param[in]   RHS      The right hand side.
 *
 *  \param[in]   RHS_col  The number of columns in the right hand side.
 *
 *  \param[in|out] X      The array to contain the solution of the system.
 *
 *  \pre   The given matrix is a non-singular matrix (complex or real) 
 *         in CSR format, RHS is an array containing the right hand side to
 *         solve for.
 *
 *  \post  X contains the solution of the system.
 */
void sparse_solve(TCSR<CPX> *matrix, CPX *RHS, int RHS_col, CPX *X) {
    int numthreads = omp_get_max_threads();
     /* Matrix data. */
    int n = matrix->size;
    int * ia = matrix->edge_i;
    int * ja = matrix->index_j;
    CPX * a = matrix->nnz;
    int mtype = 13;       // complex unsymetric
    /* RHS and solution vectors. */
    CPX * b = RHS;
    CPX * x = X;
    int nrhs = RHS_col;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    int iparm[64];
    int maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    int i;
    double ddum;          /* Double dummy */
    int idum;         /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[2] = numthreads;
    iparm[3] = 0; // 0         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[17] = 0;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = 0;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    iparm[34] = 1;        /* PARDISO use C-style indexing for ia and ja arrays */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }
//    printf ("\nReordering completed ... ");
//    printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
//    printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during numerical factorization: %d", error);
        exit (2);
    }
//    printf ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
    /* Set right hand side to one. */
//    for ( i = 0; i < n; i++ )
//    {
//        b[i] = 1;
//    }
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
  
//    for ( i = 0; i < n; i++ )
//    {
//        printf ("\n x [%d] = % f", i, x[i]);
//    }
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
    
}

void sparse_solve(TCSR<double> *matrix, double *RHS, int RHS_col, double *X) {
     int numthreads = omp_get_max_threads();
     /* Matrix data. */
    int n = matrix->size;
    int * ia = matrix->edge_i;
    int * ja = matrix->index_j;
    double * a = matrix->nnz;
    int mtype =      11;       // real unsymetric
    /* RHS and solution vectors. */
    double * b = RHS;
    double * x = X;
    int nrhs = RHS_col;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    int iparm[64];
    int maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    int i;
    double ddum;          /* Double dummy */
    int idum;         /* Integer dummy. */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[2] = numthreads;
    iparm[3] = 0; // 0         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[17] = 0;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = 0;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    iparm[34] = 1;        /* PARDISO use C-style indexing for ia and ja arrays */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }
//    printf ("\nReordering completed ... ");
//    printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
//    printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during numerical factorization: %d", error);
        exit (2);
    }
//    printf ("\nFactorization completed ... ");
/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement. */
/* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2;         /* Max numbers of iterative refinement steps. */
    /* Set right hand side to one. */
//    for ( i = 0; i < n; i++ )
//    {
//        b[i] = 1;
//    }
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
  
//    for ( i = 0; i < n; i++ )
//    {
//        printf ("\n x [%d] = % f", i, x[i]);
//    }
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
}

} // namespace
