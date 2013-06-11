// contains the parts common to the sbrueck/mcalderara codes later to be merged

#ifdef __cplusplus
extern "C" {
#endif

   typedef void(*ptr2function) (c_DBCSR, c_DBCSR, c_DBCSR *);

   void c_scf_method(c_DBCSR S, c_DBCSR KS, c_DBCSR * P);

   void yyrestart(FILE *);
   void yyparse();
   extern FILE *yyin;

#ifdef __cplusplus
}
#endif

#include <mpi.h>
#include "cDBCSR.H"
#include "cDBCSR_tools.H"
#include "conversion.H"
#include "array_tools.H"

using namespace std;

#include <math.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>

#include "Types.H"
#include "Utilities.H"
#include "InputParameter.H"
#include "MPI_Timing.H"
#include "ScaLapack.H"
#include "Blas.H"
#include "CSR.H"
