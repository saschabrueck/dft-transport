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

   c_dscal(P->nze_local,0.5,P->nzvals_local,1);

   c_dscal(KohnSham->n_nonzeros,cp2k_transport_params.evoltfactor,KohnSham->nnz,1);

   int cut_l=cp2k_transport_params.cutout[0];
   int cut_r=cp2k_transport_params.cutout[1];
   if (cut_l+cut_r) {
       OverlapCut  = new TCSR<double>(Overlap, cut_l,Overlap->size_tot-cut_l-cut_r,cut_l,Overlap->size_tot-cut_l-cut_r);
       KohnShamCut = new TCSR<double>(KohnSham,cut_l,Overlap->size_tot-cut_l-cut_r,cut_l,Overlap->size_tot-cut_l-cut_r);
   } else {
       OverlapCut  = Overlap;
       KohnShamCut = KohnSham;
   }

   transport_parameters* transport_params = new transport_parameters();
   transport_params->n_occ                      = cp2k_transport_params.n_occ;
   transport_params->n_atoms                    = cp2k_transport_params.n_atoms;
   transport_params->method                     = cp2k_transport_params.method;
   transport_params->injection_method           = cp2k_transport_params.injection_method;
   transport_params->linear_solver              = cp2k_transport_params.linear_solver;
   transport_params->n_abscissae                = cp2k_transport_params.n_abscissae;
   transport_params->n_kpoint                   = cp2k_transport_params.n_kpoint;
   transport_params->num_interval               = cp2k_transport_params.num_interval;
   transport_params->num_contacts               = cp2k_transport_params.num_contacts;
   transport_params->tasks_per_point            = cp2k_transport_params.tasks_per_point;
   transport_params->colzero_threshold          = cp2k_transport_params.colzero_threshold;
   transport_params->eps_limit                  = cp2k_transport_params.eps_limit;
   transport_params->eps_decay                  = cp2k_transport_params.eps_decay;
   transport_params->eps_singularity_curvatures = cp2k_transport_params.eps_singularity_curvatures;
   transport_params->eps_mu                     = cp2k_transport_params.eps_mu;
   transport_params->eps_eigval_degen           = cp2k_transport_params.eps_eigval_degen;
   transport_params->energy_interval            = cp2k_transport_params.energy_interval;
   transport_params->min_interval               = cp2k_transport_params.min_interval;
   transport_params->temperature                = cp2k_transport_params.temperature;
   transport_params->extra_scf                  = cp2k_transport_params.extra_scf;
   transport_params->cutout                     = cut_l+cut_r;

   std::vector<contact_type> contactvec(transport_params->num_contacts);
   int size_muvec=0;
   for (int i_c=0;i_c<transport_params->num_contacts;i_c++) if (cp2k_transport_params.contacts_data[4+5*i_c]) size_muvec++;
   std::vector<double> muvec(size_muvec);
   int i_m=0;
   int i_n=0;
   for (int i_c=0;i_c<transport_params->num_contacts;i_c++) {
      int i_t;
      if (cp2k_transport_params.contacts_data[4+5*i_c]) {
         i_t=i_m++;
      } else {
         i_t=muvec.size()+i_n++;
      }
      contactvec[i_t].bandwidth = cp2k_transport_params.contacts_data[0+5*i_c];
      contactvec[i_t].ndof      = cp2k_transport_params.contacts_data[1+5*i_c];
      if (cp2k_transport_params.contacts_data[2+5*i_c]>=0) {
          contactvec[i_t].start = cp2k_transport_params.contacts_data[2+5*i_c];
      } else if (!cp2k_transport_params.contacts_data[4+5*i_c]) {
          contactvec[i_t].start = OverlapCut->size_tot/2;
      } else if (cp2k_transport_params.contacts_data[3+5*i_c]==-1) {
          contactvec[i_t].start = OverlapCut->size_tot-contactvec[i_t].ndof*contactvec[i_t].bandwidth;
      } else {
          contactvec[i_t].start = 0;
      }
      contactvec[i_t].inj_sign  = cp2k_transport_params.contacts_data[3+5*i_c];
      contactvec[i_t].n_ele     = cp2k_transport_params.contacts_nelec[i_c];
   }

   switch (transport_params->method) {
      case 0:
         if (!rank) cout << "Writing Matrices" << endl;
         write_matrix(OverlapCut,KohnShamCut,transport_params->n_abscissae,transport_params->n_kpoint,transport_params->num_interval);
         break;
      case 1:
         if (!rank) cout << "Starting ScaLaPackDiag" << endl;
         if (diagscalapack(OverlapCut,KohnShamCut,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
     case 2:
         if (!rank) cout << "Starting CP2K core/valence Hamiltonian + OMEN Poisson local self consistent code" << endl;
         if (semiselfconsistent(OverlapCut,KohnShamCut,muvec,contactvec,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
     case 3:
     case 4:
     default:
         if (!rank) cout << "Starting Transport " << transport_params->method << endl;
         Energyvector energyvector;
         if (energyvector.Execute(OverlapCut,KohnShamCut,muvec,contactvec,transport_params)) throw SCF_Exception(__LINE__,__FILE__);
   }

   if (cut_l+cut_r) {
       cp2kCSR_to_CSR(*P, Overlap);
       c_dscal(Overlap->n_nonzeros,0.5,Overlap->nnz,1);
       Overlap->copy_shifted(OverlapCut,cut_l,Overlap->size_tot-cut_l-cut_r,cut_l,Overlap->size_tot-cut_l-cut_r);
       delete OverlapCut;
       delete KohnShamCut;
   }

   if (!transport_params->extra_scf) CSR_to_cp2kCSR(Overlap, *P);

   delete Overlap;
   delete KohnSham;
   delete transport_params;
}

