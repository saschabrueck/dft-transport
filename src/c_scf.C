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
   TCSR<double> *OverlapCut, *KohnShamCut;
   Overlap = new TCSR<double>(S.nrows_local, S.nze_local, 1);
   KohnSham = new TCSR<double>(KS.nrows_local, KS.nze_local, 1);
   cp2kCSR_to_CSR(S, Overlap);
   cp2kCSR_to_CSR(KS, KohnSham);

   transport_parameters* transport_params = new transport_parameters();
   intit_transport_parameters_from_cp2k(cp2k_transport_params, transport_params);

   c_dscal(KohnSham->n_nonzeros,transport_params->evoltfactor,KohnSham->nnz,1);

   switch (transport_params->method) {
      case 0:
         if (!rank) cout << "Writing Matrices" << endl;
         write_matrix(Overlap,KohnSham,transport_params);
      case 1:
         if (!rank) cout << "Starting ScaLaPackDiag" << endl;
         if (diagscalapack(Overlap,KohnSham,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
     case 2:
         if (!rank) cout << "Starting CP2K core/valence Hamiltonian + OMEN Poisson semi self consistent code" << endl;
         if (semiselfconsistent(Overlap,KohnSham,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
     case 3:
     default:
         if (!rank) cout << "Starting Transport" << endl;
         int cut_l=0;
         int cut_r=0;
         if (cut_l+cut_r) {
             OverlapCut  = new TCSR<double>(Overlap, cut_l,Overlap->size_tot-cut_l-cut_r,cut_l,Overlap->size_tot-cut_l-cut_r);
             KohnShamCut = new TCSR<double>(KohnSham,cut_l,Overlap->size_tot-cut_l-cut_r,cut_l,Overlap->size_tot-cut_l-cut_r);
         } else {
             OverlapCut  = Overlap;
             KohnShamCut = KohnSham;
         }
         if ( OverlapCut->size_tot%transport_params->n_cells || transport_params->bandwidth<1 ) throw SCF_Exception(__LINE__,__FILE__);
         std::vector<double> muvec(transport_params->num_contacts);
         std::vector<contact_type> contactvec(transport_params->num_contacts+1);
         for (uint i_mu=0;i_mu<contactvec.size();i_mu++) {
             contactvec[i_mu].bandwidth=transport_params->bandwidth;
             contactvec[i_mu].ndof=OverlapCut->size_tot/transport_params->n_cells; // ONLY IF ALL CELLS EQUAL
             contactvec[i_mu].n_occ=transport_params->n_occ/transport_params->n_cells; // THIS IS AN INTEGER DIVISION, IN GENERAL THE RESULT IS NOT CORRECT AND FOR CUT IT IS NOT CORRECT
         }
         contactvec[0].start=0;
         contactvec[0].inj_sign=+1;
         contactvec[1].start=OverlapCut->size_tot-contactvec[1].ndof*contactvec[1].bandwidth;
         contactvec[1].inj_sign=-1;
         contactvec[2].start=OverlapCut->size_tot/2;
         contactvec[2].inj_sign=+1;

         Energyvector energyvector;
         if (energyvector.Execute(OverlapCut,KohnShamCut,muvec,contactvec,transport_params)) throw SCF_Exception(__LINE__,__FILE__);

         if (cut_l+cut_r) {
             cp2kCSR_to_CSR(*P, Overlap);
             c_dscal(Overlap->n_nonzeros,0.5,Overlap->nnz,1);
             Overlap->copy_shifted(OverlapCut,cut_l,Overlap->size_tot-cut_r,cut_l,Overlap->size_tot-cut_r);
             delete OverlapCut;
             delete KohnShamCut;
         }
   }

   CSR_to_cp2kCSR(Overlap, *P);

   delete Overlap;
   delete KohnSham;
   delete transport_params;
}

