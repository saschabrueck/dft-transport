#include "c_scf.H"

int kpointintegration(TCSR<double>*,TCSR<double>*,TCSR<double>*,c_transport_type);
int greensolver(TCSR<double>*,TCSR<double>*,TCSR<double>*);
int diagscalapack(TCSR<double>*,TCSR<double>*,TCSR<double>*,c_transport_type);
int energyvector(TCSR<double>*,TCSR<double>*,TCSR<double>*,c_transport_type);

class SCF_Exception{
public:
    SCF_Exception(const int line,const char* file) {std::cerr<<"Error in line "<<line<<" of file "<<file<<std::endl;}
};

void c_scf_method(c_transport_type transport_env_params, c_DBCSR S, c_DBCSR KS, c_DBCSR * P)
{
   Vector1D<int> row_block_size, col_block_size, local_rows, row_dist, col_dist, 
                 nblkrows_local_all, nrows_local_all, first_rows;
   int sendbuf, nrows_local, num_nodes, rank;
   int* recvbuf;

   row_block_size.assign(S.row_blk_size, S.row_blk_size + S.nblkrows_total);
   col_block_size.assign(S.col_blk_size, S.col_blk_size + S.nblkcols_total);
   local_rows.assign(S.local_rows, S.local_rows + S.nblkrows_local);
   row_dist.assign(S.row_dist, S.row_dist + S.nblkrows_total);
   col_dist.assign(S.col_dist, S.col_dist + S.nblkcols_total);
   nrows_local = std::accumulate(row_block_size.begin()+local_rows.front()-1,row_block_size.begin()+local_rows.back(),0);

   MPI::Intercomm Comm;
   Comm = MPI::COMM_WORLD;
   num_nodes = Comm.Get_size();
   rank = Comm.Get_rank();

   recvbuf = new int[num_nodes];
   std::fill_n(recvbuf, num_nodes, 0);
   sendbuf = S.nblkrows_local;
   Comm.Allgather(&sendbuf, 1, MPI::INT, recvbuf, 1, MPI::INT);
   nblkrows_local_all.assign(recvbuf, recvbuf+num_nodes);

   TCSR<double> *Overlap, *KohnSham, *PMat;
   Overlap = new TCSR<double>(nrows_local, S.n_nze, 0);
   KohnSham = new TCSR<double>(nrows_local, KS.n_nze, 0);
   cDBCSR_to_CSR(S, Overlap);
   cDBCSR_to_CSR(KS, KohnSham);

   Overlap->size_tot=S.fullmatrix_nrows;
   KohnSham->size_tot=KS.fullmatrix_nrows;

   recvbuf = new int[num_nodes];
   std::fill_n(recvbuf, num_nodes, 0);
   sendbuf = nrows_local;
   Comm.Allgather(&sendbuf, 1, MPI::INT, recvbuf, 1, MPI::INT);
   nrows_local_all.reserve(num_nodes);
   nrows_local_all.assign(recvbuf, recvbuf+num_nodes);
   first_rows.reserve(num_nodes);
   partial_sum(nrows_local_all.begin(), nrows_local_all.end(), first_rows.begin());
   if (!rank) Overlap->first_row=0; else Overlap->first_row=first_rows[rank-1]; 
   if (!rank) KohnSham->first_row=0; else KohnSham->first_row=first_rows[rank-1]; 

   PMat = new TCSR<double>(nrows_local, S.n_nze, 0);
   PMat->copy_index(Overlap);

   switch (transport_env_params.method) {
      case 1:
         if (!rank) cout << "Starting ScaLaPackDiag" << endl;
         if (diagscalapack(Overlap,KohnSham,PMat,transport_env_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
      case 2:
         if (!rank) cout << "Starting KPointIntegration" << endl;
         if (kpointintegration(Overlap,KohnSham,PMat,transport_env_params)) throw SCF_Exception(__LINE__,__FILE__);
         break;
     case 3:
         if (!rank) cout << "Starting experimental code" << endl;
         if (greensolver(Overlap,KohnSham,PMat)) throw SCF_Exception(__LINE__,__FILE__);
         break;
     default:
         if (!rank) cout << "Starting Transport" << endl;
         if (energyvector(Overlap,KohnSham,PMat,transport_env_params)) throw SCF_Exception(__LINE__,__FILE__);
   }

   CSR_to_cDBCSR(PMat, *P, row_block_size, col_block_size, row_dist, col_dist, local_rows, nblkrows_local_all);
}
