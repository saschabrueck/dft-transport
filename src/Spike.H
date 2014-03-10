#ifndef SPIKE_H_
#define SPIKE_H_

#include <vector>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <algorithm>

#include <mpi.h>
//#include <mkl.h>
//#include <mkl_lapack.h>
#include "Blas.H"
#include "CSR.H"
#include "Types.H"
#include "Utilities.H"
#include "Pardiso.H"
#include "LinearSolver.H"
#include "mpi.h"

// Timing / benchmarking
#include <ctime>
#include <ratio>
#include <chrono>

// Debugging
#ifdef SPIKE_DEBUG
#include <cstdio>
#include <iostream>
using std::cout;
#define SPIKE_DEBUG_PATH "/users/camauro/transport/c++/src/dbg/"
#endif /* SPIKE_DEBUG */


using namespace std::chrono; 

template <typename T>
class Spike {

 public: 
  Spike(TCSR<T>* matrix, int bandwidth, T* RHS, int RHS_col,
        MPI_Comm communicator, std::vector<int> partition_lines,
        bool diagonal_dense);
  void solve_full(T* sol);
#ifdef SPIKE_DEBUG
  void debug();
#endif /* SPIKE_DEBUG */

 private:        
    
  // Main Part of the alogrithm
  void preprocess();
  void postprocess();

  // Preprocess methods
  void prepare_D();
  void calculate_spikes();
  void calculate_G();
  void delete_D();
  void prepare_send_buffer_spikes();
  void distribute_reduced_spikes();
  void prepare_send_buffer_G();
  void distribute_reduced_G();
  void prepare_sparsity_pattern();
  
  // Postprocess methods
  void solve_reduced_system();
  void distribute_Xr();
  void expand_reduced_system();
  
  // Utils (TODO: most of these belong into LinAlg)
  void get_dense_from_sparse_c(TCSR<T>* matrix_in, 
                               int start_row, int start_col, 
                               int end_row, int end_col, 
                               T* matrix_out);
  void get_dense_from_sparse_f(TCSR<T>* matrix_in,
                               int start_row, int start_col, 
                               int end_row, int end_col,
                               T* matrix_out);
  TCSR<T>* get_sparse_from_sparse(TCSR<T>* matrix_in,
                                  int start_row, int start_col,
                                  int end_row, int end_col);
  
  int  get_max_column(TCSR<T>* matrix_in);
  bool has_right_hand_side();
  void set_mpi_dataype();
  void calculate_lu_decomposition(T* matrix_in, int row, int col);
  void solve_linear_system_dense(T* LU, int row, int col, T* B, int B_cols); 
  void solve_linear_system_sparse(TCSR<T>* M, T* B, int B_cols, T* X);
  void spy(T* matrix, int row, int column);
  void spy(TCSR<T>* matrix, int row, int column);
  void full(T* matrix, int row, int column);
  void full(TCSR<T>* matrix, int row, int col);
  void print_array(T* matrix, int size);
  void merge_matrix_f(T* inout, int ldo, int start_row, 
                      int start_column, T* in, int size, int loi,
                      int length_column);
  void merge_matrix_c(T* inout, int ldo, int start_row, 
                      int start_column, T* in, int size, int loi,
                      int length_column);
  T get_sparse_matrix_value(TCSR<T>* Mat, int i, int j);
  void transpose_dense_matrix(T* mat, int row, int col);
  void trans(T* src, T* dst, const int N, const int M); 
  void MMM_dense_f(int A_rows, int A_cols, int B_cols, 
                   T prefactor_AB, T* A, 
                   int leading_dimension_A, T* B, 
                   int leading_dimension_B, T prefactor_C,
                   T* C, int leading_dimension_C);
  void MMM_dense_c(int A_rows, int A_cols, int B_cols, 
                   T prefactor_AB, T* A,
                   int leading_dimension_A, T* B, 
                   int leading_dimension_B, T prefactor_C, T* C,
                   int leading_dimension_C);
  void xLACPY(char UPLO, int M, int N, T* A, int LDA, T* B, int LDB);
  void xGEMM(char TRANSA, char TRANSB, int M, int N, int K, T ALPHA, T* A,
             int LDA, T* B, int LDB, T BETA, T* C, int LDC);

  
  // Members
  TCSR<T>* M;                     // Partion Matrix
  int M_size, M_rows, M_cols;
  
  T* F;                           // Local part of right hand side (RHS).
  int F_size, F_rows, F_cols;             
  
  T* D_dense;                     // Diagonal Block
  TCSR<T>* D_sparse;
  int D_size, D_rows, D_cols;
  bool D_is_dense;
  
  T* G;                           // Local part of intermediate RHS.
  int G_size, G_rows, G_cols;        
  
  T* V;                           // V spike.
  int V_size, V_rows, V_cols; 
  
  T* W;                           // W spike.
  int W_size, W_rows, W_cols; 
  
  TCSR<T>* Sr;                    // Reduced system
  int Sr_size, Sr_rows, Sr_cols;
    
  T* Gr;                          // Reduced system
  int Gr_size, Gr_rows, Gr_cols;
  
  T* Xr;                          // Reduced system
  int Xr_size, Xr_rows, Xr_cols;
  
  T* Xr_local;                    // Reduced system, Bottom and Top of Xr 
                                  // on local rank  
  int Xr_local_size, Xr_local_rows, Xr_local_cols;
            
  T* Sr_sendbuffer;               // Send buffer for reduced system
  int Sr_sendbuffer_size, Sr_sendbuffer_rows, Sr_sendbuffer_cols;
  
  T* Gr_sendbuffer;               // Send buffer for reduced system
  int Gr_sendbuffer_size, Gr_sendbuffer_rows, Gr_sendbuffer_cols;
  
  T* X;                           // final solution per rank; 
  int X_size, X_rows, X_cols;
  
  int spike_width;
  int partition_height;
  int* ipiv;                      // Pivot vector for factorization.
  std::vector<int> partition_ranges ;
  
  // MPI
  MPI_Comm        comm;           // MPI Communicator
  MPI_Datatype    MPI_data_type;  // Data type 
  MPI_Win         win_Sr, win_Gr, win_Xr, win_X;
  int rank; 
  int num_ranks;
  int master_rank; 
  int last_rank;

  // Debugging (inactive unless SPIKE_DEBUG preprocessor is defined)
  string debug_path;
  void dump(TCSR<T>* matrix, const char* prefix);
  void dump_c(T* matrix, int rows, int cols, const char* prefix);
  void dump_f(T* matrix, int rows, int cols, const char* prefix);
};



template<typename T>
Spike<T>:: Spike(TCSR<T>* matrix, int bandwidth, T* RHS, int RHS_col, 
                 MPI_Comm communicator, std::vector<int> partition_lines,
                 bool diagonal_dense) {

  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &num_ranks);
  set_mpi_dataype();
  partition_ranges  = partition_lines;
  master_rank       = 0;
  last_rank         = num_ranks - 1;
  comm              = communicator;
  partition_height  = partition_lines[rank + 1]-partition_lines[rank];
  int max_column    =  get_max_column(matrix); 
  spike_width       = bandwidth;
  ipiv              = new int[partition_height];
  
  M      = matrix; 
  M_size = matrix->n_nonzeros;  
  M_rows = partition_height; 
  M_cols = max_column;
  
  F      = RHS; 
  F_rows = partition_height; 
  F_cols = RHS_col;
  F_size = F_rows * F_cols;
  
  D_rows = partition_height; 
  D_cols = partition_height;
  D_size = D_rows * D_cols;
  D_is_dense = diagonal_dense;
  if (D_is_dense) {
    D_dense = new T[D_size];
    memset (D_dense, 0, D_size * sizeof(T));
  }
  
  if (rank < last_rank) {
    V_rows = partition_height; 
    V_cols = spike_width;
    V_size = V_rows * V_cols;
    V = new T[V_size];
    memset (V, 0, V_size * sizeof(T));
  }
  
  if (rank > 0) {
    W_rows = partition_height; 
    W_cols = spike_width;
    W_size = W_rows * W_cols;
    W = new T[W_size];
    memset (W, 0, W_size * sizeof(T));
  }
  
  
  G_rows = partition_height; 
  G_cols = RHS_col;
  G_size = G_rows * G_cols;
  if (has_right_hand_side()) {
      G = F;
  }
  
  if (rank == master_rank) {
      Sr_rows = 2 * spike_width * num_ranks;
      Sr_cols = 2 * spike_width * num_ranks;
      Sr_size = (num_ranks - 1) * 4 * spike_width * spike_width + num_ranks
                * 2 * spike_width;
      
      Gr_rows = num_ranks * 2 * spike_width;
      Gr_cols = G_cols;
      Gr_size = Gr_cols * Gr_rows;
      Gr = new T[Gr_size]; 
      memset(Gr, 0, Gr_size * sizeof(T));
  }

#ifdef SPIKE_DEBUG
  debug_path = SPIKE_DEBUG_PATH;
#endif
  
}


template <typename T>
void Spike<T>::solve_full(T* sol) {
  X = sol;
  preprocess();
  postprocess();
}



template <typename T>
void Spike<T>::preprocess() {
#ifdef SPIKE_DEBUG
  debug();
#endif /* SPIKE_DEBUG */
  //printf("rank %d: preparing D\n", rank);
  prepare_D();
  //printf("rank %d: preparing sparsity pattern\n", rank);
  prepare_sparsity_pattern();
  //printf("rank %d: calculating spikes\n", rank);
  calculate_spikes();
  //printf("rank %d: calculating G\n", rank);
  calculate_G();
  delete_D();
  //printf("rank %d: preparing spike send buffer\n", rank);
  prepare_send_buffer_spikes();
  //printf("rank %d: distribute reduced spikes\n", rank);
  distribute_reduced_spikes();
  //printf("rank %d: preparing G send buffer\n", rank);
  prepare_send_buffer_G();
  //printf("rank %d: distribute G\n", rank);
  distribute_reduced_G();
}


template <typename T>
void Spike<T>::postprocess() {

  if (rank == 0) {
    //printf("rank %d: solving reduced system\n", rank);
  }
  solve_reduced_system();

  if (rank == 0) {
    //printf("rank %d: distributing Xr\n", rank);
  }
  distribute_Xr();

  //printf("rank %d: expanding Xr\n", rank);
  expand_reduced_system();
}


#ifdef SPIKE_DEBUG
template <typename T>
void Spike<T>::debug() {

  dump_f(F, F_rows, F_cols, "F");

};
#endif /* SPIKE_DEBUG */


/** \brief Calculation of the final result by using Xr 
 */
template <typename T>
void Spike<T>::expand_reduced_system() {

  X_rows = partition_height; 
  X_cols = F_cols;
  X_size = X_rows * X_cols;

  if (num_ranks == 2) {

    if (rank == master_rank) {

      // X0u = G0u - V0u * X1t    % u -> [1:end-bw,:]
      xLACPY('A', G_rows - spike_width, G_cols, G, G_rows, X, X_rows);

      dump_f(X, X_rows, X_cols, "G0u");

      int X1t_offset = 2 * spike_width;
      xGEMM('N', 'N', V_rows - spike_width, Xr_local_cols, V_cols, T(-1.0), V,
            V_rows, &Xr_local[X1t_offset], Xr_local_rows, T(1.0), X, X_rows);

      dump_f(X, X_rows, X_cols, "X0u");

      // Xb = X0b                % b -> [end-bw+1:end,:]
      int X0b_offset = spike_width;
      int Xb_offset = X_rows - spike_width;
      xLACPY('A', spike_width, Xr_local_cols, &Xr_local[X0b_offset],
             Xr_local_rows, &X[Xb_offset], X_rows);

    } else {

      int Gl_offset = spike_width;
      int Xl_offset = spike_width;

      // X1l = G1l - W1l * X0b    % l -> [bw+1:end,:]
      xLACPY('A', G_rows - spike_width, G_cols, &G[Gl_offset], G_rows,
             &X[Xl_offset], X_rows);

      dump_f(X, X_rows, X_cols, "G1l");

      int Wl_offset = spike_width;
      xGEMM('N', 'N', W_rows - spike_width, Xr_local_cols, W_cols, T(-1.0),
            &W[Wl_offset], W_rows, Xr_local, Xr_local_rows, T(1.0),
            &X[Xl_offset], X_rows);

      dump_f(X, X_rows, X_cols, "X1l");

      // Xt = X1t                % t -> [1:bw,:]
      int X1t_offset = spike_width;
      xLACPY('A', spike_width, Xr_local_cols, &Xr_local[X1t_offset],
             Xr_local_rows, X, X_rows);

    }



  } else {

    // TODO: this is probably not fully working yet

    // Prepare GNprime
    int middle_height = partition_height - 2 * spike_width;
    int middle_size = middle_height * spike_width;
    
    int GNprime_prefactor;
    T* C;
   
    if (has_right_hand_side()) {
      GNprime_prefactor = 1;
      C =  &G[spike_width];
    } else {
      C = new T[middle_size];
      memset(C, 0, middle_size * sizeof(T)); 
      GNprime_prefactor = 0;   // so we don't need to initialize
    }
    
    if (rank == 0) {

        int m = partition_height - 2 * spike_width;
        int n = G_cols;
        int k = spike_width;
         
        T* A = &V[spike_width];
        int lda = partition_height;
           
        T* B = &Xr_local[2 * spike_width];
        int ldb = Xr_local_rows;  
       
        int ldc = G_rows;
          
        MMM_dense_f(m, k, n, T(-1.0), A, lda, B,ldb, T(GNprime_prefactor), C, ldc);
         
        merge_matrix_c(X, partition_height, 0, 0, Xr_local,
                       F_cols * spike_width, Xr_local_rows, spike_width);
        merge_matrix_c(X, partition_height, 0, spike_width, C,
                       F_cols * middle_height, partition_height, middle_height);
        merge_matrix_c(X, partition_height, 0, spike_width + m, &Xr_local[spike_width],
                       F_cols * spike_width, Xr_local_rows, spike_width);
        delete[] V; 

    } else if (rank == last_rank) {
       
        int m = partition_height - 2 * spike_width;
        int n = G_cols;
        int k = spike_width;
        
      
        T* A = &W[spike_width];
        int lda = partition_height;
           
        T* B = &Xr_local[0];
        
        int ldb = Xr_local_rows;
        
        int ldc = G_rows;
          
        MMM_dense_f(m, k, n, T(-1.0), A, lda, B,ldb, T(GNprime_prefactor), C, ldc);
         
        merge_matrix_c(X, partition_height, 0, 0, &Xr_local[spike_width],
                       F_cols * spike_width, Xr_local_rows, spike_width);
        merge_matrix_c(X, partition_height, 0, spike_width, C,
                       F_cols * middle_height, partition_height, middle_height);
        merge_matrix_c(X, partition_height, 0, spike_width + m, &Xr_local[2 * spike_width],
                       F_cols * spike_width, Xr_local_rows, spike_width);
        delete[] W;
        
    } else {

      int m = partition_height - 2 * spike_width;
      int n = G_cols;
      int k = spike_width;
    
      T* A = &W[spike_width];
      int lda = partition_height;
         
      T* B = &Xr_local[0];

      int ldb = Xr_local_rows;
      
      int ldc = G_rows;
        
      MMM_dense_f(m, k, n, T(-1.0), A, lda, B,ldb, T(GNprime_prefactor), C, ldc); 
      
      A = &V[spike_width];
      B = &Xr_local[3 * spike_width];
      
      MMM_dense_f(m, k, n, T(-1.0), A, lda, B,ldb, T(1), C, ldc); 
      
      merge_matrix_c(X, partition_height, 0, 0, &Xr_local[spike_width],
                     F_cols * spike_width, Xr_local_rows, spike_width);
      merge_matrix_c(X, partition_height, 0, spike_width, C,
                     F_cols * middle_height, partition_height, middle_height);
      merge_matrix_c(X, partition_height, 0, spike_width + m, &Xr_local[2 * spike_width],
                     F_cols * spike_width, Xr_local_rows, spike_width);
      delete[] V;
      delete[] W;

    }

  }

  delete[] Xr_local;

  dump_f(X, X_rows, X_cols, "X");

}

template <typename T>
void Spike<T>::distribute_Xr() {

  // NOTE: Xr is sent in row major format since construction of Xr on the 
  // other ranks happens by vertical concatenation (i.e. appending of rows);

  // Let 0 denote the index of the previous rank, 1 the index of the current
  // rank and 2 the index of the next rank. Then, after distribute_Xr(), each
  // rank other than the first and last have
  //
  //          Xr = [X0b;
  //                X1t;
  //                X1b;
  //                X2t;]
  //
  // while the first and last rank have
  //
  //          Xr = [X1t;
  //                X1b;
  //                X2t;
  //                 * ;]
  //             and
  //          Xr = [X0b;
  //                X1t;
  //                X1b;
  //                 * ;]
  //
  // respectively. All matrices are in column major format. Concatenation is
  // according to matlab, i.e. vertical in the above notation.

  transpose_dense_matrix(Xr, Xr_rows, Xr_cols); 
 
  int Xr_offset;
  
  if (rank == master_rank) {
    MPI_Win_create(Xr, Xr_size * sizeof(T), sizeof(T), 
                   MPI_INFO_NULL, comm, &win_Xr);  
  } else {
    MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, comm,
                   &win_Xr);
  }
 
  if (rank == 0 || rank == last_rank) {
    Xr_local_cols = F_cols; 
    Xr_local_rows = 3 * spike_width; 
    Xr_local_size = Xr_local_cols * Xr_local_rows;
    Xr_offset  = 0;
  } else {
    Xr_local_cols = F_cols; 
    Xr_local_rows = 4 * spike_width; 
    Xr_local_size = Xr_local_cols * Xr_local_rows;
  }
  
  if (rank > 0) {
    Xr_offset = spike_width * F_cols + (rank - 1) * (2 * spike_width * F_cols); 
  }
  
  Xr_local = new CPX[Xr_local_size];


  // Start fence
  MPI_Win_fence(0, win_Xr);

  if (rank != master_rank) {
    MPI_Get(Xr_local, Xr_local_size, MPI_data_type, master_rank, Xr_offset,
            Xr_local_size, MPI_data_type, win_Xr);
  }
  if (rank == master_rank) {
    memcpy(Xr_local, Xr + Xr_offset, Xr_local_size * sizeof(T));
  }

  // Stop fence
  MPI_Win_fence(0, win_Xr);
  MPI_Win_free(&win_Xr);
  
  transpose_dense_matrix(Xr_local, Xr_local_rows, Xr_local_cols);
  
  if (rank == master_rank) {  
    delete[] Xr;  
  }

  dump_f(Xr_local, Xr_local_rows, Xr_local_cols, "Xr_local");

}




/** \brief Solve the reduced sparse system
 */
template<typename T>
void Spike<T>::solve_reduced_system() {

  // 2 PARTITION NOTE:
  //  
  //  Faster way:
  //    E   = unity(bandwidth) - Wt*Vb;
  //    X1t = E\(G1t - Wt * G0b);
  //    X0b = G0b - Vb * X1t;
  //
  //  Then send only X0b and X1t to rank1

  if (rank == master_rank) {

    Xr_rows = Gr_cols;
    Xr_cols = Gr_rows;
    Xr_size = Xr_cols * Xr_rows;
    Xr = new T[Xr_size];
    memset (Xr, 0, Xr_size * sizeof(T));
    
    solve_linear_system_sparse(Sr, Gr, Gr_cols, Xr);
    
    delete[] Gr;
    delete Sr;

    dump_f(Xr, Xr_rows, Xr_cols, "Xr");

  }
}

template <typename T>
void Spike<T>::prepare_send_buffer_G() {

  // NOTE: Gr_sendbuffer is constructed in row major format since
  // construction of Gr on the master rank happens by vertical concatenation
  // (i.e. appending of rows);
  

  // 2 PARTITION NOTE:
  //      rank1 only needs to send G1t

  if (has_right_hand_side()) {   

    Gr_sendbuffer_rows = G_cols; 
    Gr_sendbuffer_cols = 2 * spike_width;
    Gr_sendbuffer_size = Gr_sendbuffer_rows * Gr_sendbuffer_cols; 

    // Buffer for the reduced system
    Gr_sendbuffer = new T[Gr_sendbuffer_size];
    
    merge_matrix_c(Gr_sendbuffer, Gr_sendbuffer_cols, 0, 0, G,
                   G_cols * spike_width, partition_height, spike_width);
    merge_matrix_c(Gr_sendbuffer, Gr_sendbuffer_cols, 0, spike_width, 
                   &G[partition_height - spike_width], G_cols * spike_width,
                   partition_height, spike_width);

    transpose_dense_matrix(Gr_sendbuffer, Gr_sendbuffer_rows, Gr_sendbuffer_cols);

    std::swap(Gr_sendbuffer_rows, Gr_sendbuffer_cols);

  }

}


template <typename T>
void Spike<T>::distribute_reduced_G() {

  // 2 PARTITION NOTE:
  //      rank1 only needs to send G1t

  int gr_offset = rank * (Gr_sendbuffer_size);

  if (rank == master_rank) {

    MPI_Win_create(Gr, Gr_size * sizeof(T), sizeof(T), MPI_INFO_NULL,
                   comm, &win_Gr);

  } else {

    MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, comm,
                   &win_Gr);

  }

  // Start fence
  MPI_Win_fence(0, win_Gr);

  if (has_right_hand_side() && rank != master_rank) {

        MPI_Put(Gr_sendbuffer, Gr_sendbuffer_size, MPI_data_type, master_rank,
                 gr_offset, Gr_sendbuffer_size, MPI_data_type, win_Gr);

  } else if (rank == master_rank) {

        memcpy(Gr, Gr_sendbuffer, Gr_sendbuffer_size * sizeof(T));   

  }

  // Stop fence
  MPI_Win_fence(0, win_Gr);
  MPI_Win_free(&win_Gr);

  delete[] Gr_sendbuffer;

  transpose_dense_matrix(Gr, Gr_rows, Gr_cols); 

  if (rank == master_rank) {

    dump_f(Gr, Gr_rows, Gr_cols, "Gr");

  }
}

template <typename T>
void Spike<T>::prepare_send_buffer_spikes() {

  // 2 PARTITION NOTE:
  //      rank1 only needs to send Wt

  CPX* ones = new T[2 * spike_width];
  for (int i = 0; i < 2 * spike_width; ++i) {
     ones[i]=T(1);
  }
  
  if (rank == 0 || rank == last_rank) {
    Sr_sendbuffer_size = 2 * (spike_width * spike_width + spike_width); 
    Sr_sendbuffer_rows = spike_width + 1; 
    Sr_sendbuffer_cols = 2 * spike_width;
  } else {
    Sr_sendbuffer_size = 2 * (2 * spike_width * spike_width + spike_width); 
    Sr_sendbuffer_rows = 2 * spike_width + 1; 
    Sr_sendbuffer_cols = 2 * spike_width;
  }
  
  Sr_sendbuffer = new T[Sr_sendbuffer_size];
  memset(Sr_sendbuffer, 0,Sr_sendbuffer_size * sizeof(T));
   
  if (rank == 0) {

    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, 0, 1, V, 
                   spike_width * spike_width, partition_height, spike_width);
    
    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, spike_width, 1,
                   &V[partition_height  - spike_width],
                   spike_width * spike_width, partition_height, spike_width);
    
    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, 0, 0, &ones[0],
                   Sr_sendbuffer_cols, Sr_sendbuffer_cols, Sr_sendbuffer_cols);
       
  } else if (rank == last_rank) {

    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, 0, 0, W,
                   spike_width * spike_width, partition_height, spike_width);

    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, spike_width, 0,
                   &W[(partition_height - spike_width)],
                   spike_width * spike_width, partition_height, spike_width);

    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, 0, spike_width, &ones[0],
                   Sr_sendbuffer_cols, Sr_sendbuffer_cols, Sr_sendbuffer_cols);
     
  } else {

    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, 0, 0, W,
                   spike_width * spike_width, partition_height, spike_width);

    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, spike_width, 0, 
                   &W[(partition_height - spike_width)],
                   spike_width * spike_width, partition_height, spike_width);     

    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, 0, spike_width + 1, V,
                   spike_width * spike_width, partition_height, spike_width);    

    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, spike_width, spike_width + 1,
                   &V[ (partition_height - spike_width)],
                   spike_width * spike_width, partition_height, spike_width);

    merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_cols, 0, spike_width, &ones[0],
                   Sr_sendbuffer_cols, Sr_sendbuffer_cols, Sr_sendbuffer_cols);   
 
  }
 
}

template <typename T>
void Spike<T>::prepare_D() {

  if (D_is_dense) {

    get_dense_from_sparse_f(M, 0, M_rows * rank, M_rows, M_rows * (rank + 1),
                            D_dense);

    dump_f(D_dense, D_rows, D_cols, "D_dense");

    calculate_lu_decomposition(D_dense, D_rows, D_cols);

  } else {

    D_sparse = get_sparse_from_sparse(M, 0, M_rows * rank, M_rows,
                                      M_rows * (rank + 1));

    dump(D_sparse, "D_sparse");

  }
}

template <typename T>
void Spike<T>::calculate_spikes() {
     
  if (rank < last_rank) {

    if (D_is_dense) {

      get_dense_from_sparse_f(M, 0, M_rows * (rank + 1), M_rows,
                              M_rows * (rank + 1) + V_cols, V);  

      solve_linear_system_dense(D_dense, D_rows, D_cols, V, V_cols);

    } else {

      T* B = new T[V_size];

      get_dense_from_sparse_f(M, 0, M_rows * (rank + 1), M_rows,
                              M_rows * (rank + 1) + V_cols, B);

      //printf("rank %d: start solving for V\n", rank);
      solve_linear_system_sparse(D_sparse, B, V_cols, V);
      //printf("rank %d: done solving for V\n", rank);

      delete[] B;

    }

    dump_f(V, V_rows, V_cols, "V");

  }

  if (rank > 0) {

    if (D_is_dense) {

      get_dense_from_sparse_f(M, 0, M_rows * rank - W_cols, M_rows, M_rows * rank, W);

      solve_linear_system_dense(D_dense, D_rows, D_cols, W, W_cols);        

    } else {

      T* C = new T[W_size];

      get_dense_from_sparse_f(M, 0, M_rows * rank - W_cols, M_rows, M_rows * rank, C);

      //printf("rank %d: start solving for W\n", rank);
      solve_linear_system_sparse(D_sparse, C, W_cols, W);
      //printf("rank %d: stop solving for W\n", rank);

      delete[] C;
    }

    dump_f(W, W_rows, W_cols, "W");

  }
}

template <typename T>
void Spike<T>::calculate_G() {
  if (has_right_hand_side()) {
    if (D_is_dense) {

      solve_linear_system_dense(D_dense, D_rows, D_cols, G, G_cols);   

    } else {

      T* temp = new T[G_size];

      //printf("rank %d: start solving for G\n", rank);
      solve_linear_system_sparse(D_sparse, G, G_cols, temp);
      //printf("rank %d: stop solving for G\n", rank);

      // NOTE:
      // What we use for G is actually 'sol', which is preallocated memory.
      // Thus we don't delete it but swap it with temp instead.  However, if
      // the caller assumes that we don't change this, we're in trouble ...
      std::swap(G, temp);

    }

    dump_f(G, G_rows, G_cols, "G");

  }
}

template <typename T>
void Spike<T>::delete_D() {
  if (D_is_dense) {
    delete[] D_dense;
  } else {
    delete D_sparse;
  }
}

template <typename T>
void Spike<T>::distribute_reduced_spikes() {

  // 2 PARTITION NOTE:
  //      rank1 only needs to send Wt

  // NOTE: since Sr is constructed by adding rows, we change the memory layout
  // to row major

  transpose_dense_matrix(Sr_sendbuffer, Sr_sendbuffer_rows, Sr_sendbuffer_cols);

  int Sr_offset;

  if (rank == master_rank) {

    MPI_Win_create(Sr->nnz, Sr_size * sizeof(T), sizeof(T), MPI_INFO_NULL,
                   comm, &win_Sr);
    Sr_offset = 0;

  } else {

    MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, comm, &win_Sr);
    Sr_offset = 2 * (spike_width * spike_width + spike_width) + 2 * (rank - 1) *
                (2 * spike_width * spike_width + spike_width);

  }

  // Start fence
  MPI_Win_fence(0, win_Sr);

  if (rank != master_rank) {   

    MPI_Put(Sr_sendbuffer, Sr_sendbuffer_size, MPI_data_type, master_rank,
            Sr_offset, Sr_sendbuffer_size, MPI_data_type, win_Sr);

  } else if (rank == master_rank) {

    memcpy(Sr->nnz, Sr_sendbuffer, Sr_sendbuffer_size * sizeof(T));  

  }   

  // Stop fence
  MPI_Win_fence(0, win_Sr);
  MPI_Win_free(&win_Sr);

  delete[] Sr_sendbuffer;

  if (rank == master_rank) {

    dump(Sr, "Sr");

  }

}

template<typename T>
void Spike<T>::prepare_sparsity_pattern() {
  if (rank == master_rank) {
    
    int Sr_num_rows = 2 * spike_width * num_ranks;
    int Sr_num_nnz  = (4 * spike_width * spike_width + 2 * spike_width) *
                      (num_ranks - 2) + 4 * (spike_width * spike_width + 
                      spike_width);
    
    Sr = new TCSR<T>(Sr_num_rows, Sr_num_nnz, 0);
    
    int findx = Sr->findx;
    int num_partitions = num_ranks;


    // Construct sparsity pattern

    // edge_i and index_i
    Sr->edge_i[0] = 0;
    for (int i = 1; i < Sr_rows + 1; ++i) {
      if (i <= 2 * spike_width || i > (num_partitions - 1) * 2 * spike_width) {
        Sr->edge_i[i] = Sr->edge_i[i-1] + 1 + spike_width;
      } else {
        Sr->edge_i[i] = Sr->edge_i[i-1] + 1 + 2 * spike_width;
      }
      Sr->index_i[i - 1] = Sr->edge_i[i] - Sr->edge_i[i - 1];
    }
    Sr->get_row_edge();
    Sr->get_diag_pos();
    
    int index = 0;
    
    // First block row only contains |Vr|
    int Vr_start = 2 * spike_width;
    for (int matrix_row = 0; matrix_row < 2 * spike_width; ++matrix_row) {
      Sr->index_j[index++] = matrix_row + findx;    // diagonal elements
      for (int Vr_column = 0; Vr_column < spike_width; ++Vr_column) {
        Sr->index_j[index++] = Vr_start + Vr_column + findx; // V0t and V0b
      }
    }

    // All other except the last contain |Wr|1|Vr|
    for (int block_row = 1; block_row < num_partitions-1; ++block_row) {
      int Wr_start = (2 * (block_row - 1) + 1) * spike_width;
      int Vr_start = (2 * (block_row + 1)) * spike_width;
      for (int matrix_row = 2 * spike_width * block_row;
           matrix_row < 2 * spike_width * (block_row + 1); ++matrix_row) {
        // Wr
        for (int Wr_column = 0; Wr_column < spike_width; ++Wr_column) {
          Sr->index_j[index++] = Wr_start + Wr_column + findx;
        }
        // 1
        Sr->index_j[index++] = matrix_row + findx;
        // Vr
        for (int Vr_column = 0; Vr_column < spike_width; ++Vr_column) {
          Sr->index_j[index++] = Vr_start + Vr_column + findx;
        }
      }
    }

    // Last block only contains |Wr|
    int Wr_start = (2 * (num_partitions - 2) + 1) * spike_width;
    for (int matrix_row = 2 * spike_width * (num_partitions - 1);
         matrix_row < 2 * spike_width * num_partitions; ++matrix_row) {
      for (int Wr_column = 0; Wr_column < spike_width; ++Wr_column) {
        Sr->index_j[index++]=Wr_start + Wr_column + findx;
      }
      // 1
      Sr->index_j[index++] = matrix_row + findx;
    }
  }
}


/****************************************************************************
 * General routines
 ***************************************************************************/

/** \brief            Extracts a subblock of a sparse matrix by specifying
 *                    the corner indices storing it as dense matrix in C/C++
 *                    convention. It assumes that matrix_out is initialized as
 *                    all zero.
 *
 *  \param[in]        matrix_in
 *                    The CSR matrix to extract the subblock from.
 *
 *  \param[in]        start_row
 *                    Row to start extraction (included, C-numbering).
 *
 *  \param[in]        start_col
 *                    Colum to start extraction (included, C-numbering).
 *
 *  \param[in]        end_row
 *                    Row to stop extraction (excluded, C-numbering).
 *
 *  \param[in]        end_col
 *                    Column to stop extraction (excluded, C-numbering).
 *
 *  \param[out]       matrix_out
 *                    Zero initialized array to contain the subblock.
 */
template<typename T>
void Spike<T>::get_dense_from_sparse_c(TCSR<T>* matrix_in, int start_row, 
                                       int start_col, int end_row, int end_col,
                                       T* matrix_out) {
    int row = end_row - start_row ;
    int column = end_col - start_col;

    int findx = matrix_in->findx; 

    for (int i=start_row; i < end_row ; ++i) {
        for (int j = matrix_in->edge_i[i]; j < matrix_in->edge_i[i + 1]; ++j) {
            int indx = matrix_in->index_j[j];
            if (indx < start_col) continue;
            if (indx >= end_col)  break;
            matrix_out[(i-start_row) * (column) + (indx-start_col)] = 
                matrix_in->nnz[j];
        }
    }
}

/** \brief Extracts a subblock of a sparse matrix by specifying the corner
 *         indices storing it as dense matrix in Fortran convention. It assumes 
 *         that matrix_out is initialized as all zero.
 *
 * \param[in] matrix_in   The CSR matrix to extract the subblock from.
 *
 * \param[in] start_row   Row to start extraction (included, C-numbering).
 *
 * \param[in] start_col   Colum to start extraction (included, C-numbering).
 *
 * \param[in] end_row     Row to stop extraction (excluded, C-numbering).
 *
 * \param[in] end_col     Column to stop extraction (excluded, C-numbering).
 *
 * \param[out] matrix_out Zero initialized array to contain the subblock.
 */
template <typename T>
void Spike<T>::get_dense_from_sparse_f(TCSR<T>* matrix_in, int start_row,
                                       int start_col, int end_row, int end_col,
                                       T* matrix_out) {
  int findx = matrix_in->findx; 

  int rows_out = end_row - start_row;

  for (int row_in = start_row; row_in < end_row; ++row_in) {
    for (int index = matrix_in->edge_i[row_in] - findx;
         index < matrix_in->edge_i[row_in + 1] - findx; ++index) {

      int col_in = matrix_in->index_j[index] - findx;

      if (col_in < start_col) {

        continue;

      } else if (col_in < end_col) {

        int col_out   = col_in - start_col;
        int row_out   = row_in - start_row;
        int array_pos = col_out * rows_out + row_out;

        matrix_out[array_pos] = matrix_in->nnz[index];

      } else {

        break;

      }
    }
  }
}

/** \brief Extracts a subblock of a sparse matrix by specifying the corner
 *         indices and returns the result in CSR format.
 *
 * \param[in] matrix_in   The CSR matrix to extract the subblock from.
 *
 * \param[in] start_row   Row to start extraction (included, C-numbering).
 *
 * \param[in] start_col   Colum to start extraction (included, C-numbering).
 *
 * \param[in] end_row     Row to stop extraction (excluded, C-numbering).
 *
 * \param[in] end_col     Column to stop extraction (excluded, C-numbering).
 */
template <typename T>
TCSR<T>* Spike<T>::get_sparse_from_sparse(TCSR<T>* matrix_in, int start_row,
                                          int start_col, int end_row,
                                          int end_col) {
  T*   nnz_in     = matrix_in->nnz;
  int* edge_i_in  = matrix_in->edge_i;
  int* index_j_in = matrix_in->index_j;
  int  findx      = matrix_in->findx;

  // Count non-zeros
  int nnz_count = 0;
  int start_index = edge_i_in[start_row] - findx;
  int stop_index  = edge_i_in[end_row] - findx;

  for (int index = start_index; index < stop_index; ++index) {

    int col = index_j_in[index] - findx;

    if ((col > start_col - 1) && (col < end_col)) {
      ++nnz_count;
    }

  }

  // Create new CSR matrix
  int rows = end_row - start_row;
  TCSR<T>* matrix_out = new TCSR<T>(rows, nnz_count, findx);

  T*   nnz_out     = matrix_out->nnz;
  int* edge_i_out  = matrix_out->edge_i;
  int* index_j_out = matrix_out->index_j;

  // Copy data
  int out_index = 0;
  for (int in_row = start_row; in_row < end_row + 1; ++in_row) {

    int out_row = in_row - start_row;
    edge_i_out[out_row] = out_index + findx;

    for (int in_index = edge_i_in[in_row] - findx;
         in_index < edge_i_in[in_row + 1] - findx; 
         ++in_index) {

      int in_col = index_j_in[in_index] - findx;

      if (in_col < start_col) {
        continue;
      } else if (in_col < end_col) {

        T   out_val = nnz_in[in_index];
        int out_col = in_col - start_col;

        nnz_out[out_index]     = out_val;
        index_j_out[out_index] = out_col + findx;

        ++out_index;

      } else {
        break;
      }
    }
  }

  // Construct index_i (# elements per line)
  int* index_i_out = matrix_out->index_i;
  for (int out_row = 0; out_row < end_row - start_row + 1; ++out_row) {
    index_i_out[out_row] = edge_i_out[out_row + 1] - edge_i_out[out_row];
  }

  // Put matrix into consistent state
  matrix_out->get_row_edge();
  matrix_out->get_diag_pos();

  return matrix_out;
}

template <typename T>
int Spike<T>::get_max_column(TCSR<T>* matrix_in) {
   return  *std::max_element(
           &matrix_in->index_j[0],
           &matrix_in->index_j[matrix_in->n_nonzeros]);
}

template <typename T>
bool Spike<T>::has_right_hand_side() {
   return (F != NULL);
}

/** \brief Copy one Matrix into an other
 */
template<typename T>
void Spike<T>::merge_matrix_f(T* inout, int ldo, int start_row, 
                            int start_column, T* in, int size, int loi,
                            int length_column) {
    merge_matrix_c(inout, ldo, start_column, start_row, 
                              in, size, loi,
                             length_column);
}

/** \brief Copy one Matrix into an other
 */
template<typename T>
void Spike<T>::merge_matrix_c(T* inout, int ldo, int start_row, 
                            int start_column, T* in, int size, int loi,
                            int length_column) {
  inout += start_row * ldo + start_column; // Go to position
  for (int i = 0, s = 0; s < size ; i += loi, s += length_column) {
    // Copy consecutive column
    memcpy (&inout[0], &in[i], length_column * sizeof(T));
    // Jump to next row in inout matrix
    inout += ldo;    
  }   
}


/** \brief            Solves a linear system using sparse linear algebra
 *
 * \param[in]         M
 *                    Coefficient matrix in CSR format
 *
 * \param[in]         B
 *                    Right hand side of the equation.
 *
 * \param[in]         B_cols
 *                    Number of columns of B.
 *
 * \param[in|out]     X
 *                    The solution.
 */
template <typename T>
void Spike<T>::solve_linear_system_sparse(TCSR<T>* Mat, T* B, int B_cols, T* X) {
  Pardiso::sparse_solve(Mat, B, B_cols, X);
}


/** \brief Get the value of the entry i,j of the sparse matrix
 *
 * \param[in]  Mat   The sparse matrix.
 *
 * \param[in]    i   row number.
 *
 * \param[in]    j   column number.
 *
 */
template <typename T>
T Spike<T>::get_sparse_matrix_value(TCSR<T>* Mat, int i, int j) {
   
  // loop up interval for col_indices
  int start = Mat->edge_i[i];
  int end   = Mat->edge_i[i + 1];

  // serach j in col_indices
  for (int g = start; g < end; ++g) {
    if (Mat->index_j[g] == j) { // if found return data form that pos
      return Mat->nnz[g];
    }
  }
  return T(0, 0);
}

/** \brief            Transpose dense matrix (out of place)
 *
 *  \param[in|out]    matrix
 *                    The matrix to transpose / the transposed matrix.
 *
 *  \param[in]        rows
 *                    Number of rows of matrix before the call.
 *
 *  \param[in]        cols
 *                    Number of columns of matrix before the call.
 *
 */
template <typename T>
void Spike<T>::transpose_dense_matrix(T* matrix, int rows, int cols) {

  int size = rows * cols; 
  T* temp = new T[size];

  // MARKER: why is this needed?
  memcpy(temp, matrix, size * sizeof(T)); 
  
  trans(temp, matrix, rows, cols);
  
  delete[] temp;
}

/** \brief            Transposes a dense matrix src and writes the result into
 *                    dst
 *
 *  \param[in]        src
 *                    The matrix to be transposed.
 *
 *  \param[in]        dst
 *                    The matrix containing the transposed of src after the
 *                    function call.
 *
 *  \param[in]        rows
 *                    Number of row of src.
 *
 *  \param[in]        cols
 *                    Number of columns of src.
 */
template <typename T>
void Spike<T>::trans(T* src, T* dst, const int rows, const int cols) {

  #pragma omp parallel for
  for (int n = 0; n < rows * cols; ++n) {

    int i  = n / rows;
    int j  = n % rows;
    dst[n] = src[cols * j + i];

  }

}


#ifdef SPIKE_DEBUG
template <typename T>
void Spike<T>::dump(TCSR<T>* matrix, const char* prefix) {
  std::cout << "rank " << rank << ": dumping " << prefix << "\n";
  stringstream filename;
  filename << debug_path << prefix << "_rank" << rank << ".csv";
  matrix->write(filename.str().c_str());
};
template <typename T>
void Spike<T>::dump_c(T* matrix, int rows, int cols, const char* prefix) {
  std::cout << "rank " << rank << ": dumping " << prefix << "\n";
  stringstream filename;
  filename << debug_path << prefix << "_rank" << rank << ".csv";
  write_mat_c(matrix, rows, cols, filename.str().c_str());
};
template <typename T>
void Spike<T>::dump_f(T* matrix, int rows, int cols, const char* prefix) {
  std::cout << "rank " << rank << ": dumping " << prefix << "\n";
  stringstream filename;
  filename << debug_path << prefix << "_rank" << rank << ".csv";
  write_mat_f(matrix, rows, cols, filename.str().c_str());
};
#else
template <typename T>
void Spike<T>::dump(TCSR<T>* matrix, const char* prefix) {};
template <typename T>
void Spike<T>::dump_c(T* matrix, int rows, int cols, const char* prefix) {};
template <typename T>
void Spike<T>::dump_f(T* matrix, int rows, int cols, const char* prefix) {};
#endif /* SPIKE_DEBUG */

#endif
