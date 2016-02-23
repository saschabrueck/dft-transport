/**  
 *   \brief auxiliary types and functions for coupling CP2K and OMEN 
 *   \author Mohammad Hossein Bani-Hashemian
 */
#ifndef __c_scf
#define __c_scf

/* header files required by template functions defined in this file */
#include <mpi.h>
#include "CSR.H"
#include "libcp2k.H"
#include "Types.H"

/* here include header files of OMEN source files that evaluate a density matrix    *
 * the density matrix will be passed to CP2K ...                                    */
#include "DiagScaLapack.H"
#include "SemiSelfConsistent.H"
#include "EnergyVector.H"
#include "WriteMatrix.H"


/* types -------------------------------------------------------------------------- */
struct SCF_Exception{
    SCF_Exception(const int line,const char* file) {std::cerr<<"Error in line "<<line<<" of file "<<file<<std::endl;}
};

/* functions ---------------------------------------------------------------------- */
/**  
 *   \brief initializes an OMEN CSR matrix using a CP2K C-interoperable CSR matrix
 *   \author Mohammad Hossein Bani-Hashemian
 */
template <class T>
void cp2kCSR_to_CSR (cp2k_csr_interop_type& cp2kCSRmat, TCSR<T>* CSRmat)
{
    int N, n_nnz;
  
    N = cp2kCSRmat.nrows_local;
    n_nnz = cp2kCSRmat.nze_local;
  
    CSRmat->size = N;
    CSRmat->size_tot = cp2kCSRmat.nrows_total;
    CSRmat->n_nonzeros = n_nnz;
    CSRmat->type = 0;
    CSRmat->findx = 1;
    CSRmat->first_row = cp2kCSRmat.first_row;
  
    CSRmat->set_nnz(n_nnz, cp2kCSRmat.nzvals_local);
    CSRmat->set_index_i(N, cp2kCSRmat.nzerow_local);
    CSRmat->set_index_j(n_nnz, cp2kCSRmat.colind_local);
    CSRmat->set_edge_i(N, cp2kCSRmat.rowptr_local);
}

/**  
 *   \brief initializes a CP2K C-interoperable matrix using an OMEN CSR matrix
 *   \author Mohammad Hossein Bani-Hashemian
 */
template <class T>
void CSR_to_cp2kCSR (TCSR<T>* CSRmat, cp2k_csr_interop_type& cp2kCSRmat)
{
    MPI::Intercomm Comm;
    Comm = MPI::COMM_WORLD;
 
    c_icopy(CSRmat->size+1, CSRmat->edge_i, 1, cp2kCSRmat.rowptr_local, 1);
    c_icopy(CSRmat->n_nonzeros, CSRmat->index_j, 1, cp2kCSRmat.colind_local, 1);
    c_icopy(CSRmat->size, CSRmat->index_i, 1, cp2kCSRmat.nzerow_local, 1);
    c_tcopy(CSRmat->n_nonzeros, CSRmat->nnz, 1, cp2kCSRmat.nzvals_local, 1);
 
    cp2kCSRmat.nrows_total = CSRmat->size_tot;
    cp2kCSRmat.ncols_total = CSRmat->size_tot;
    cp2kCSRmat.nze_local = CSRmat->n_nonzeros;
    Comm.Allreduce(&cp2kCSRmat.nze_local, &cp2kCSRmat.nze_total, 1, MPI::INT, MPI::SUM);
    cp2kCSRmat.nrows_local = CSRmat->size;
    cp2kCSRmat.data_type = 3; /* dbcsr_types : dbcsr_type_real_8 = 3 */
    cp2kCSRmat.first_row = CSRmat->first_row;
}

#endif
