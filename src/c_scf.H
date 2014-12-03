#ifndef __c_scf
#define __c_scf


#include <mpi.h>
#include <iostream>
#include "cDBCSR.H"
#include "cDBCSR_tools.H"
#include "conversion.H"
#include "array_tools.H"
#include "Types.H"

#include "KPointIntegration.H"
#include "DiagScaLapack.H"
#include "SemiSelfConsistent.H"

struct SCF_Exception{
    SCF_Exception(const int line,const char* file) {std::cerr<<"Error in line "<<line<<" of file "<<file<<std::endl;}
};

#ifdef __cplusplus
extern "C" {
#endif

   typedef void(*ptr2function) (c_transport_type, c_DBCSR, c_DBCSR, c_DBCSR *);
   void c_scf_method(c_transport_type trasnport_env_params, c_DBCSR S, c_DBCSR KS, c_DBCSR * P);

#ifdef __cplusplus
}
#endif

#endif
