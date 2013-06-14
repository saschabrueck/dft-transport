#ifndef c_scf_H
#define c_scf_H

#include <mpi.h>
#include "cDBCSR.H"
#include "cDBCSR_tools.H"
#include "conversion.H"
#include "array_tools.H"
//#include "GreenSolver.H"

int greensolver(TCSR<double>*,TCSR<double>*,TCSR<double>*);
int diagscalapack(TCSR<double>*,TCSR<double>*,TCSR<double>*,int);
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

#endif
