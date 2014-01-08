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

// Timing / benchmarking
#include <ctime>
#include <ratio>
#include <chrono>


using namespace std::chrono; 

template <typename T>
class Spike {

    public: 
        Spike(TCSR<T> *matrix,int bandwidth, T *RHS, int RHS_col, 
                MPI_Comm communicator,  std::vector<int> partition_lines);
        
        void solve_full();
        void get_X(T* sol_X);
        void print_solution_full();
    private:        
        
        // Main Part of the alogrithm
        void preprocess();
        void postprocess();

        //Preprocess methods
        void prepare_D();
        void calculate_spikes();
        void calculate_G();
        void prepare_send_buffer_spikes();
        void distribute_reduced_spikes();
        void prepare_send_buffer_G();
        void distribute_reduced_G();
        void prepare_sparsity_pattern();
        
        //Postprocess methods
        void solve_reduced_system();
        void distribute_Xr();
        void expand_reduced_system();
        
        //Utils
        void get_dense_from_sparse_c(TCSR<T>* matrix_in, 
                                        int start_row,  int start_col, 
                                        int end_row,    int end_col, 
                                        T* matrix_out);
        void get_dense_from_sparse_f(TCSR<T>* matrix_in,
                                        int start_row,  int start_col, 
                                        int end_row,    int end_col,
                                        T* matrix_out);
        
        int  get_max_column(TCSR<T>* matrix_in);
        bool has_right_hand_side();
        void set_mpi_dataype();
        void calculate_lu_decomposition(T* matrix_in,int row, int col);
        void solve_linear_system_dense(T *LU, int row, int col ,T *B, int B_cols); 
        void solve_linear_system_sparse(TCSR<T> *M, T *B, int B_cols, T *X);
        void spy(T * matrix, int row, int column);
        void spy(TCSR<T> * matrix, int row, int column);
        void full(T * matrix, int row, int column);
        void full(TCSR<T> * matrix,int row, int col);
        void print_array(T * matrix, int size);
        void merge_matrix_f(T * inout, int ldo, int start_row, 
                            int start_column, T* in, int size,int loi,
                            int length_column);
        void merge_matrix_c(T * inout, int ldo, int start_row, 
                            int start_column, T* in, int size,int loi,
                            int length_column);
        T get_sparse_matrix_value(TCSR<T> *Mat,int i, int j);
        void transpose_dense_matrix(T* mat,int row, int col);
        void trans(T *src, T *dst, const int N, const int M); 
        void MMM_dense_f(int A_rows, int A_cols, int B_cols, 
                                T prefactor_AB, T *A, 
                                int leading_dimension_A, T *B, 
                                int leading_dimension_B, T prefactor_C,
                                T *C, int leading_dimension_C);
        void MMM_dense_c(int A_rows, int A_cols, int B_cols, 
                             T prefactor_AB, T *A,
                             int leading_dimension_A, T *B, 
                             int leading_dimension_B, T prefactor_C, T *C,
                             int leading_dimension_C);
      
        
        
        //Memebers
        TCSR<T> *M;                      //Partion Matrix
        int M_size,  M_row , M_col;
        
        T *F;                           // Local part of right hand side (RHS).
        int F_size,  F_row , F_col;             
        
        T* D;                           // Diagonal Block
        int D_size,  D_row,  D_col;
        
        T *G;                           // Local part of intermediate RHS.
        int G_size,  G_row , G_col;        
        
        T *V;                           // V spike.
        int V_size,  V_row , V_col; 
        
        T *W;                           // W spike.
        int W_size,  W_row , W_col; 
        
        TCSR<T> *Sr;                    //Reduced system
        int Sr_size,  Sr_row , Sr_col;
          
        T *Gr;                          //Reduced system
        int Gr_size,  Gr_row , Gr_col;
        
        T *Xr;                           //Reduced system
        int Xr_size,  Xr_row , Xr_col;
        
        T *Xr_local;                    //Reduced system, Bottom and Top of Xr 
                                        //on local rank  
        int Xr_local_size,  Xr_local_row , Xr_local_col;
                  
        T *Sr_sendbuffer;                // send buffer for reduced system
        int Sr_sendbuffer_size,  Sr_sendbuffer_row , Sr_sendbuffer_col;
        
        T *Gr_sendbuffer;                // send buffer for reduced system
        int Gr_sendbuffer_size,  Gr_sendbuffer_row , Gr_sendbuffer_col;
        
        T*  X;  //end solution per rank; 
        int X_size,  X_row , X_col;
        
        int spike_width;
        int partition_height;
        int* ipiv;                  // Pivot vector for factorization.
        std::vector<int> partition_ranges ;
        
        //MPI
        MPI_Comm        comm;                  //MPI Communicator
        MPI_Datatype    MPI_data_type;      //The data type 
        MPI_Win         win_Sr,win_Gr,win_Xr,win_X;
        int rank; 
        int world_size;
        int master_rank; 
        int last_rank;
};



template<typename T>
 Spike<T>:: Spike(TCSR<T> *matrix,int bandwidth, T *RHS, int RHS_col, 
                MPI_Comm communicator,  std::vector<int> partition_lines)
{
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator,&world_size);
    set_mpi_dataype();
    partition_ranges = partition_lines;
    master_rank  = 0;
    last_rank    = world_size -1;
    comm = communicator;
    partition_height = partition_lines[rank + 1]-partition_lines[rank];
    int max_column       =  get_max_column(matrix); 
    spike_width = bandwidth;
    ipiv = new int[partition_height];
    
    M      = matrix; 
    M_size = matrix->n_nonzeros;  
    M_row  = partition_height; 
    M_col  = max_column;
    
    F      = RHS; 
    F_size = partition_height*RHS_col;  
    F_row  = partition_height; 
    F_col  = RHS_col;
    
    D_size = partition_height*partition_height;  
    D_row  = partition_height; 
    D_col  = partition_height;
    D = new T[D_size];
    memset (D,0,D_size*sizeof(T));
    
    if(rank < last_rank){
        V_size = partition_height*spike_width;  
        V_row  = partition_height; 
        V_col  = spike_width;
        V = new T[V_size];
        memset (V,0,V_size*sizeof(T));
    }
    
    if(rank > 0){
        W_size = partition_height*spike_width;  
        W_row  = partition_height; 
        W_col  = spike_width;
        W = new T[W_size];
        memset (W,0,W_size*sizeof(T));
    }
    
    
    G_size = partition_height*RHS_col;  
    G_row  = partition_height; 
    G_col  = RHS_col;
    if(has_right_hand_side()){
        G = F;
    }
    
    if(rank == master_rank){
        Sr_row = 2 * spike_width * world_size;
        Sr_col = 2 * spike_width * world_size;
        Sr_size = (world_size -1 )* 4 * spike_width * spike_width + world_size * 2 * spike_width;
        
      
        Gr_row  = world_size * 2 * spike_width;
        Gr_col  = G_col;
        Gr_size = Gr_col * Gr_row;
        Gr = new T[Gr_size]; 
        memset (Gr,0,Gr_size*sizeof(T));
    }
    
}


template <typename T>
void Spike<T>::solve_full(){
  preprocess();
  postprocess();
}



template <typename T>
void Spike<T>::preprocess(){
   prepare_D();
   prepare_sparsity_pattern();
   calculate_spikes();
   calculate_G();
   prepare_send_buffer_spikes();
   distribute_reduced_spikes();
   prepare_send_buffer_G();
   distribute_reduced_G();
}


template <typename T>
void Spike<T>::postprocess(){
   solve_reduced_system();
   distribute_Xr();
   expand_reduced_system();
}

template <typename T>
void Spike<T>::get_X(T* sol_X){
    sol_X = X;
}

template <typename T>
void Spike<T>::print_solution_full(){
  T* full_X; 
  int full_X_size;
  if (rank == master_rank) {
   full_X_size = F_col*partition_ranges[world_size];
   full_X = new T[full_X_size];
    
    MPI_Win_create(full_X, full_X_size*sizeof(T), sizeof(T), 
                   MPI_INFO_NULL, comm, &win_X);  
  } else {
    MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, comm,
                   &win_X);
  }
  
  transpose_dense_matrix(X, X_col, X_row);
  
  int Xr_offset = partition_ranges[rank]*F_col; 

  MPI_Win_fence(0, win_X); // start fence
  if (rank != master_rank) {
        MPI_Put(X, X_size, MPI_data_type, master_rank, Xr_offset,
                X_size, MPI_data_type, win_X);
  }
  if (rank == master_rank) {
        memcpy(full_X+Xr_offset,X, X_size*sizeof(T));
  }
  MPI_Win_fence(0, win_X); // stop fence
  MPI_Win_free(&win_X);
  
  if (rank == master_rank) {
    full(full_X,partition_ranges[world_size],F_col);
    delete   full_X;  
  }

}


/** \brief Calculation of the final result by using Xr 
 */
template <typename T>
void Spike<T>::expand_reduced_system(){
  X_row = partition_height; 
  X_col = F_col;
  X_size = X_row*X_col;

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
    memset ( C, 0, middle_size*sizeof(T) ); 
    GNprime_prefactor = 0;   // so we don't need to initialize
  }
  
  X = new T[X_size];  
  
  if (rank == 0) {

      int m = partition_height - 2*spike_width;
      int n = G_col;
      int k = spike_width;
       
      T * A = &V[spike_width];
      int lda = partition_height;
         
      T * B = &Xr_local[2*spike_width];
      int ldb = Xr_local_row;  
     
      int ldc = G_row;
        
      MMM_dense_f(m, k, n,T(-1.0),A, lda,B,ldb,T(GNprime_prefactor), C, ldc);
       
      merge_matrix_c(X, partition_height, 0, 0, Xr_local,
                     F_col * spike_width, Xr_local_row, spike_width);
      merge_matrix_c(X, partition_height, 0,spike_width , C,
                     F_col * middle_height, partition_height, middle_height);
      merge_matrix_c(X, partition_height, 0, spike_width+m, &Xr_local[spike_width],
                     F_col * spike_width, Xr_local_row, spike_width);
      delete[] V; 
  } 
  else if (rank == last_rank) {
     
      int m = partition_height - 2*spike_width;
      int n = G_col;
      int k = spike_width;
      
    
      T * A = &W[spike_width];
      int lda = partition_height;
         
      T * B = &Xr_local[0];
      
      int ldb = Xr_local_row;
      
      int ldc = G_row;
        
      MMM_dense_f(m, k, n,T(-1.0),A, lda,B,ldb,T(GNprime_prefactor), C, ldc);
       
      merge_matrix_c(X, partition_height, 0, 0, &Xr_local[spike_width],
                     F_col * spike_width, Xr_local_row, spike_width);
      merge_matrix_c(X, partition_height, 0,spike_width , C,
                     F_col * middle_height, partition_height, middle_height);
      merge_matrix_c(X, partition_height, 0, spike_width+m, &Xr_local[2*spike_width],
                     F_col * spike_width, Xr_local_row, spike_width);
       delete[] W;
  } 
  else {
      int m = partition_height - 2*spike_width;
      int n = G_col;
      int k = spike_width;
    
      T * A = &W[spike_width];
      int lda = partition_height;
         
      T * B = &Xr_local[0];

      int ldb = Xr_local_row;
      
      int ldc = G_row;
        
      MMM_dense_f(m, k, n,T(-1.0),A, lda,B,ldb,T(GNprime_prefactor), C, ldc); 
      
      A = &V[spike_width];
      B = &Xr_local[3*spike_width];
      
      MMM_dense_f(m, k, n,T(-1.0),A, lda,B,ldb,T(1), C, ldc); 
      
      merge_matrix_c(X, partition_height, 0, 0, &Xr_local[spike_width],
                     F_col * spike_width, Xr_local_row, spike_width);
      merge_matrix_c(X, partition_height, 0,spike_width , C,
                     F_col * middle_height, partition_height, middle_height);
      merge_matrix_c(X, partition_height, 0, spike_width+m, &Xr_local[2*spike_width],
                     F_col * spike_width, Xr_local_row, spike_width);
     delete[] V;
     delete[] W;
  }

    delete[] Xr_local;
}

template <typename T>
void Spike<T>::distribute_Xr() {
 
  int Xr_offset;
  
  if (rank == master_rank) {
    MPI_Win_create(Xr, Xr_size*sizeof(T), sizeof(T), 
                   MPI_INFO_NULL, comm, &win_Xr);  
  } else {
    MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, comm,
                   &win_Xr);
  }
 
  if (rank == 0 || rank == last_rank) {
      Xr_local_col = F_col; 
      Xr_local_row = 3 * spike_width; 
      Xr_local_size = Xr_local_col*Xr_local_row;
      Xr_offset  = 0;
  } else{
      Xr_local_col = F_col; 
      Xr_local_row = 4 * spike_width; 
      Xr_local_size = Xr_local_col*Xr_local_row;
  }
  
  if(rank > 0){
      Xr_offset = spike_width*F_col+(rank - 1) * (2*spike_width*F_col); 
  }
  
  Xr_local = new CPX[Xr_local_size];

  MPI_Win_fence(0, win_Xr); // start fence
  if (rank != master_rank) {
        MPI_Get(Xr_local, Xr_local_size, MPI_data_type, master_rank, Xr_offset,
                Xr_local_size, MPI_data_type, win_Xr);
  }
  if (rank == master_rank) {
      memcpy(Xr_local, Xr+Xr_offset, Xr_local_size*sizeof(T));
  }
  MPI_Win_fence(0, win_Xr); // stop fence
  MPI_Win_free(&win_Xr);
  
  // transpose Xr_local to put it back in F format
  transpose_dense_matrix(Xr_local, Xr_local_row, Xr_local_col);
  
  if (rank == master_rank) {  
    delete[]   Xr;  
  }
}




/** \brief Solve the reduced sparse system
 */
template<typename T>
void Spike<T>::solve_reduced_system(){
  // We dont need D anymore
  delete[] D;
  if(rank == master_rank){
        Xr_row  = Gr_col;
        Xr_col  = Gr_row;
        Xr_size = Xr_col * Xr_row;
        Xr = new T[Xr_size];
        memset ( Xr, 0, Xr_size*sizeof(T) );
    
        transpose_dense_matrix(Gr, Gr_row, Gr_col); 
       
        solve_linear_system_sparse(Sr, Gr, Gr_col, Xr);
       
        transpose_dense_matrix(Xr, Xr_row, Xr_col); 
        
       // full(Xr, Xr_col, Xr_row);
        delete[] Gr;
        delete Sr;
  }
}

template <typename T>
void Spike<T>::prepare_send_buffer_G(){
    if(has_right_hand_side()){   
        Gr_sendbuffer_row  =  G_col; 
        Gr_sendbuffer_col  =  2*spike_width;
        Gr_sendbuffer_size = Gr_sendbuffer_row * Gr_sendbuffer_col; 
        Gr_sendbuffer = new T[Gr_sendbuffer_size]; // Buffer of the reduced system
        
        merge_matrix_c(Gr_sendbuffer, Gr_sendbuffer_col, 0, 0, G,
                     G_col * spike_width, partition_height, spike_width);
        merge_matrix_c(Gr_sendbuffer, Gr_sendbuffer_col, 0, spike_width, 
                   &G[partition_height - spike_width], 
                   G_col * spike_width, partition_height, spike_width);
//        
       transpose_dense_matrix(Gr_sendbuffer, Gr_sendbuffer_row, Gr_sendbuffer_col);
       // swap 
       int tmp = Gr_sendbuffer_col; 
       Gr_sendbuffer_col = Gr_sendbuffer_row;
       Gr_sendbuffer_row = tmp;
     
    }
}


template <typename T>
void Spike<T>::distribute_reduced_G(){
  int gr_offset = rank * (Gr_sendbuffer_size);
  if (rank == master_rank) {
    MPI_Win_create(Gr, Gr_size*sizeof(T), sizeof(T), MPI_INFO_NULL,
                   comm, &win_Gr);
  } else {
    MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, comm,
                   &win_Gr);
  }

  MPI_Win_fence(0, win_Gr); // start fence
  if (has_right_hand_side() && rank != master_rank) {
        MPI_Put(Gr_sendbuffer,Gr_sendbuffer_size, MPI_data_type, master_rank,
                 gr_offset,Gr_sendbuffer_size ,MPI_data_type, win_Gr);
  }else if(rank == master_rank){
        memcpy(Gr, Gr_sendbuffer, Gr_sendbuffer_size*sizeof(T));   
  }
  MPI_Win_fence(0, win_Gr); // stop fence
  MPI_Win_free(&win_Gr);
  delete[]   Gr_sendbuffer;
}

template <typename T>
void Spike<T>::prepare_send_buffer_spikes(){
     CPX *ones = new T[2 * spike_width];
     for (int i = 0; i < 2 * spike_width; ++i) {
        ones[i]=T(1);
     }
    
    if(rank == 0 || rank == last_rank){
        Sr_sendbuffer_size = 2 * (spike_width * spike_width + spike_width); 
        Sr_sendbuffer_row  = spike_width + 1; 
        Sr_sendbuffer_col  = 2 * spike_width;
    }else{
        Sr_sendbuffer_size = 2 * (2 * spike_width * spike_width + spike_width); 
        Sr_sendbuffer_row  = 2 * spike_width + 1; 
        Sr_sendbuffer_col  = 2 * spike_width;
    }
    
   Sr_sendbuffer = new T[Sr_sendbuffer_size];
   memset(Sr_sendbuffer,0,Sr_sendbuffer_size*sizeof(T));
    
   if(rank == 0 ){
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col , 0, 1, V, 
                   spike_width * spike_width, partition_height ,spike_width);
      
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col, spike_width, 1,
                   &V[partition_height  - spike_width],
                   spike_width * spike_width,partition_height, spike_width);
      
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col, 0, 0, &ones[0],
                   Sr_sendbuffer_col, Sr_sendbuffer_col,Sr_sendbuffer_col);
        
   }else if(rank == last_rank){
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col, 0, 0, W,
                   spike_width * spike_width, partition_height, spike_width);
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col, spike_width, 0,
                   &W[(partition_height - spike_width)],
                   spike_width * spike_width, partition_height, spike_width);
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col, 0, spike_width, &ones[0],
                   Sr_sendbuffer_col, Sr_sendbuffer_col, Sr_sendbuffer_col);
      
   }else{
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col, 0, 0, W, spike_width * spike_width,partition_height,
                   spike_width);
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col, spike_width, 0, 
                   &W[ (partition_height - spike_width)],
                   spike_width * spike_width,partition_height, spike_width);     
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col, 0,spike_width + 1, V,
                   spike_width * spike_width,partition_height, spike_width);    
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col, spike_width ,spike_width + 1,
                   &V[ (partition_height - spike_width)],
                   spike_width * spike_width,partition_height, spike_width);
      merge_matrix_f(Sr_sendbuffer, Sr_sendbuffer_col,  0 ,spike_width, &ones[0],
                   Sr_sendbuffer_col, Sr_sendbuffer_col, Sr_sendbuffer_col);   
     
   }
   // we need to transpose other wise we have wrong regeions when using just one put
    transpose_dense_matrix(Sr_sendbuffer, Sr_sendbuffer_row, Sr_sendbuffer_col);
    
}

template <typename T>
void Spike<T>::prepare_D(){

    get_dense_from_sparse_f(M, 0, M_row * rank,
                            M_row, M_row * (rank + 1),
                            D);
     
    calculate_lu_decomposition(D,D_row,D_col);
}

template <typename T>
void Spike<T>::calculate_spikes(){
     

    if(rank < last_rank){
      // V
      get_dense_from_sparse_f(M, 0, M_row * (rank+1),
                            M_row,  M_row * (rank+1)+V_col,
                            V);  
          
      solve_linear_system_dense(D, D_row, D_col ,V, V_col);     
    }

    if(rank > 0){
      // W
      get_dense_from_sparse_f(M, 0, M_row * rank - W_col,
                               M_row, M_row * rank,
                              W);    
       
      solve_linear_system_dense(D, D_row, D_col ,W, W_col);        
    }
    
}

template <typename T>
void Spike<T>::calculate_G(){
    if(has_right_hand_side()){
   
      solve_linear_system_dense(D, D_row, D_col ,G, G_col);   
        // if(rank == 0)  full(G,G_col,G_row);
    }
}

template <typename T>
void Spike<T>::distribute_reduced_spikes(){
  int Sr_offset;
  if (rank == master_rank) {
    MPI_Win_create(Sr->nnz, Sr_size*sizeof(T), sizeof(T), MPI_INFO_NULL,
                   comm, &win_Sr);
    Sr_offset = 0;
  } else {
    MPI_Win_create(MPI_BOTTOM, 0, 1, MPI_INFO_NULL, comm, &win_Sr);
    Sr_offset = 2 * (spike_width * spike_width + spike_width) + 2 * (rank - 1) *
                (2 * spike_width * spike_width + spike_width);
  }

  MPI_Win_fence(0, win_Sr); // start fence
     if(rank != master_rank){   
        MPI_Put(Sr_sendbuffer, Sr_sendbuffer_size, MPI_data_type, master_rank,
          Sr_offset, Sr_sendbuffer_size, MPI_data_type, win_Sr);
     }else if(rank == master_rank){
          memcpy(Sr->nnz, Sr_sendbuffer, Sr_sendbuffer_size*sizeof(T));  
     }   
  MPI_Win_fence(0, win_Sr); // stop fence
  MPI_Win_free(&win_Sr);
  delete[]   Sr_sendbuffer;
}

template<typename T>
void Spike<T>::prepare_sparsity_pattern() {
    if(rank == master_rank){
        
        int   Sr_num_rows = 2 * spike_width * world_size;
        int Sr_num_nnz = (4 * spike_width * spike_width + 2 * spike_width) * 
               (world_size - 2) + 
               4 * (spike_width * spike_width + spike_width);
        
        Sr = new TCSR<T>(Sr_num_rows, Sr_num_nnz,0);
        
        int findx=Sr->findx;
        int num_partitions = world_size;
        // Construct sparsity pattern

        // The rows
        Sr->edge_i[0] = 0;
        for(int i = 1; i <= Sr_row ; i++ ){
          if(i <= 2*spike_width || i > (num_partitions-1)*2*spike_width)
               Sr->edge_i[i] = Sr->edge_i[i-1]+1+spike_width;
          else
               Sr->edge_i[i] =  Sr->edge_i[i-1]+1+2*spike_width;
        }
        
        int index = 0;
        

//        // First block row only contains Vr
        int Vr_start = 2 * spike_width;
        for (int matrix_row = 0; matrix_row < 2 * spike_width; ++matrix_row) {
          Sr->index_j[index++] = matrix_row + findx;    // diagonal elements
          for (int Vr_column = 0; Vr_column < spike_width; ++Vr_column) {
            Sr->index_j[index++] = Vr_start + Vr_column + findx; // V0t and V0b
          }
        }

        // All other except the last contain Wr 1 Vr
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

        // Last block only contains Wr
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

/** \brief Extracts a subblock of a sparse matrix by specifying the corner
 *         indices storing it as dense matrix in C/C++ convention. It assumes 
 *         that matrix_out is initialized as all zero.
 *
 *  \param[in] matrix_in   The CSR matrix to extract the subblock from.
 *
 *  \param[in] start_row   Row to start extraction (included, C-numbering).
 *
 *  \param[in] start_col   Colum to start extraction (included, C-numbering).
 *
 *  \param[in] end_row     Row to stop extraction (excluded, C-numbering).
 *
 *  \param[in] end_col     Column to stop extraction (excluded, C-numbering).
 *
 *  \param[out] matrix_out Zero initialized array to contain the subblock.
 */
template<typename T>
void Spike<T>::get_dense_from_sparse_c(TCSR<T> *matrix_in, int start_row, 
                                       int start_col, int end_row, int end_col,
                                       T* matrix_out)
{
    int row = end_row - start_row ;
    int column = end_col - start_col;

    int findx = matrix_in->findx; 

    int i,j;

    for(int i=start_row; i < end_row ; ++i){
        for(int j = matrix_in->edge_i[i]; j < matrix_in->edge_i[i+1]; ++j){
            int indx = matrix_in->index_j[j];
            if( indx < start_col)  continue;
            if( indx >= end_col)    break;
            // printf("i = %d ;  edge[i=%d]=%d ; j=%d  ; index=%d  a = %d   
            //        length_r = %d\n  ",i,i,edge_i[i],j,indx,a,r_length);
            matrix_out[(i-start_row) * (column) + (indx-start_col)] = 
                matrix_in->nnz[j];
        }
    }
}

/** \brief Extracts a subblock of a sparse matrix by specifying the corner
 *         indices storing it as dense matrix in Fortran convention. It assumes 
 *         that matrix_out is initialized as all zero.
 *
 *  \param[in] matrix_in   The CSR matrix to extract the subblock from.
 *
 *  \param[in] start_row   Row to start extraction (included, C-numbering).
 *
 *  \param[in] start_col   Colum to start extraction (included, C-numbering).
 *
 *  \param[in] end_row     Row to stop extraction (excluded, C-numbering).
 *
 *  \param[in] end_col     Column to stop extraction (excluded, C-numbering).
 *
 *  \param[out] matrix_out Zero initialized array to contain the subblock.
 */
template <typename T>
void Spike<T>::get_dense_from_sparse_f(TCSR<T> *matrix_in, int start_row,
                                       int start_col, int end_row, int end_col,
                                       T* matrix_out)
{
    int row = end_row - start_row;

    int findx = matrix_in->findx; 

    for(int i = start_row; i < end_row; ++i){
        for(int j = matrix_in->edge_i[i]-findx; j < matrix_in->edge_i[i+1]-findx; ++j){
            int indx = matrix_in->index_j[j]-findx;
            if( indx < start_col)  continue;
            if( indx >= end_col)    break;
            // printf("i = %d ;  edge[i=%d]=%d ; j=%d  ; index=%d  \n",i,i,matrix_in->edge_i[i],j,indx);
            matrix_out[(i-start_row) + (indx - start_col) * row] = 
                matrix_in->nnz[j];
            
            //printf("nnz = %f   %d %d \n ", real(matrix_in->nnz[j]),i,indx);
        }
    }
}

template <typename T>
int Spike<T>::get_max_column(TCSR<T>* matrix_in){
   return  *std::max_element(
           &matrix_in->index_j[0],
           &matrix_in->index_j[matrix_in->n_nonzeros]);
}

template <typename T>
bool Spike<T>::has_right_hand_side(){
   return (F != NULL);
}

/** \brief Copy one Matrix into an other
 */
template<typename T>
void Spike<T>::merge_matrix_f(T * inout, int ldo, int start_row, 
                            int start_column, T* in, int size, int loi,
                            int length_column){
    merge_matrix_c(inout,ldo, start_column, start_row, 
                              in,  size, loi,
                             length_column);
}

/** \brief Copy one Matrix into an other
 */
template<typename T>
void Spike<T>::merge_matrix_c(T * inout, int ldo, int start_row, 
                            int start_column, T* in, int size,int loi,
                            int length_column){
  inout += start_row*ldo + start_column; // Go to position
  for(int i = 0 , s = 0; s < size ; i += loi, s+=length_column){
    // Copy consecutive column
    memcpy (&inout[0],&in[i],length_column*sizeof(T) );
    // Jump to next row in inout matrix
    inout += ldo;    
  }   
}


/** \brief Solves a linear system using sparse linear algebra
 *
 *  \param[in]  M   Coefficient matrix in CSR format
 *
 *  \param[in]  B   Right hand side of the equation.
 *
 *  \param[in]  B_cols  Number of columns of B.
 *
 *  \param[in|out] X  The solution.
 */
template <typename T>
void Spike<T>::solve_linear_system_sparse(TCSR<T> *Mat, T *B, int B_cols, T *X) {
  Pardiso::sparse_solve(Mat, B, B_cols, X);
}


/** \brief Get the value of the entry i,j of the sparse matrix
 *
 *  \param[in]  Mat   The sparse matrix.
 *
 *  \param[in]    i   row number.
 *
 *  \param[in]    j   column number.
 *
 */
template <typename T>
T Spike<T>::get_sparse_matrix_value(TCSR<T> *Mat,int i, int j){
   
    // loop up interval for col_indices
    int start = Mat->edge_i[i];
    int end   = Mat->edge_i[i+1];

    // serach j in col_indices
    for(int g = start; g < end; g++){
              if(Mat->index_j[g] == j) // if found return data form that pos
                      return Mat->nnz[g];
    }
    return T(0,0);
}

/** \brief Transpose dense matrix 
 *
 *  \param[in|out] mat  The trsposed matrix.
 *
 *  \param[in]    row   Number of row of mat.
 *
 *  \param[in]    col   Number of columns of mat.
 *
 */
template <typename T>
void Spike<T>::transpose_dense_matrix(T* mat,int row, int col){
    int size = row*col; 
    T* temp = new T[size];
    memcpy(temp,mat,size*sizeof(T)); 
    
    trans(temp, mat,  row, col);
    
    delete[] temp; 
}

/** \brief Transposes a dense matrix src and writes the result in dst
 *
 *  \param[in]  src   Dense source  Matrix
 *
 *  \param[in]  dst   Dense destination Matrix
 *
 *  \param[in]  N      Number of row of src.
 *
 *  \param[in]  M     Number of columns of src.
 */
template <typename T>
void Spike<T>::trans(T *src, T *dst, const int N, const int M) {
    #pragma omp parallel for
    for(int n = 0; n<N*M; n++) {
        int i = n/N;
        int j = n%N;
        dst[n] = src[M*j + i];
    }
}

#endif