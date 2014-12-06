#include "c_scf.H"

/**  
 *   \brief Takes the overlap (S) and Kohn-Sham (KS) matrices as input and returns a density matrix (P).
 *          This function acts as the gate to the CP2K's world. 
 *          Here only call functions that evaluate a P matrix.   
 *   \author Mohammad Hossein Bani-Hashemian
 */
void c_scf_method(cp2k_transport_parameters cp2k_transport_params, cp2k_csr_interop_type S, cp2k_csr_interop_type KS, cp2k_csr_interop_type * P)
{
   MPI::Intercomm Comm;
   Comm = MPI::COMM_WORLD;
   int rank = Comm.Get_rank();

   TCSR<double> *Overlap, *KohnSham;
   Overlap = new TCSR<double>(S.nrows_local, S.nze_local, 1);
   KohnSham = new TCSR<double>(KS.nrows_local, KS.nze_local, 1);
   cp2kCSR_to_CSR(S, Overlap);
   cp2kCSR_to_CSR(KS, KohnSham);

   transport_parameters* transport_params = new transport_parameters();
   intit_transport_parameters_from_cp2k(cp2k_transport_params, transport_params);

   switch (transport_params->method) {
      case 0:
         if (!rank) cout << "Writing Matrices" << endl;
         if (semiselfconsistent(Overlap,KohnSham,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
      case 1:
         if (!rank) cout << "Starting ScaLaPackDiag" << endl;
         if (diagscalapack(Overlap,KohnSham,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
     case 2:
         if (!rank) cout << "Starting CP2K core/valence Hamiltonian + OMEN Poisson semi self consistent code" << endl;
         if (semiselfconsistent(Overlap,KohnSham,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
     case 3:
         if (!rank) cout << "Starting Transport" << endl;
         if (semiselfconsistent(Overlap,KohnSham,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
     default:
         if (!rank) cout << "Starting Transport" << endl;
         if (semiselfconsistent(Overlap,KohnSham,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
   }

   CSR_to_cp2kCSR(Overlap, *P);

   delete Overlap;
   delete KohnSham;
   delete transport_params;
}

