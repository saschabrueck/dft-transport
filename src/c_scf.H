#ifndef c_scf_H
#define c_scf_H

#include <mpi.h>
#include "cDBCSR.H"
#include "cDBCSR_tools.H"
#include "conversion.H"
#include "array_tools.H"

#ifdef __cplusplus
extern "C" {
#endif

   typedef void(*ptr2function) (int, c_DBCSR, c_DBCSR, c_DBCSR *);

   void c_scf_method(int nelectron, c_DBCSR S, c_DBCSR KS, c_DBCSR * P);

#ifdef __cplusplus
}
#endif

#endif
