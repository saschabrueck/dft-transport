#include "Spike.H"
#include <iomanip>
#include <sstream>
#include <string>

template <>
void Spike<CPX>::set_mpi_dataype() {
  MPI_data_type = MPI_COMPLEX16;
}

template <>
void Spike<double>::set_mpi_dataype() {
  MPI_data_type = MPI_DOUBLE;   
}

template <>
void Spike<CPX>::calculate_lu_decomposition(CPX* m, int rows, int cols) {
  int info; 
  c_zgetrf(rows, cols, m, rows, ipiv, &info);
}

template <>
void Spike<double>::calculate_lu_decomposition(double* m, int rows, int cols) {
  int info;
  c_dgetrf(rows, cols, m, rows, ipiv, &info); 
}

/** \brief Solves a linear system using dense linear algebra
 *
 *  Solves the system D_dense for a given right hand side B.
 *
 *  \param[in|out]  B   Right hand side of the equation.
 *
 *  \param[in]      B_cols  Number of columns of B
 *
 * M_TODO: since our diagonal blocks are actually symmetric,
 * this could be sped up by using a corresponding LAPACK
 * routine. But it would need a change in our extraction
 * routine.
 */
template <>
void Spike<CPX>::solve_linear_system_dense(CPX *LU, int rows, int cols ,CPX *B, int B_cols) {
  int info;
  c_zgetrs('N', rows, B_cols, LU, rows, ipiv, B,  rows, &info);
}

template <>
void Spike<double>::solve_linear_system_dense(double *LU, int rows, int cols,double *B, int B_cols) {
  int info;
  c_dgetrs('N', rows, B_cols, LU, rows, ipiv, B,  rows, &info);
}



template <>
void Spike<CPX>::spy(CPX* matrix, int rows, int columns) {
  std::cout << "\n";
  for(int i = 0; i < rows; ++i) {
    for(int j = 0; j < columns; ++j) {
      int position = (i * columns) + j;
      if (real(matrix[position]) != 0 || imag(matrix[position]) != 0) {
        std::cout << "\033[0;31m" << "*" << "\033[0;30m";
      } else {
         std::cout << "*";     
      }
    }
    std::cout << std::endl; 
  }
}

template <>
void Spike<double>::spy(double* matrix, int rows, int columns) {
  std::cout << "\n";
  for(int i = 0; i < rows; ++i) {
    for(int j = 0; j < columns; ++j) {
      int position = (i * columns) + j;
      if (matrix[position] != 0) {
        std::cout << "\033[0;31m" << "*" << "\033[0;30m";
      } else {
        std::cout << "*";     
      }
    }
    std::cout << std::endl; 
  }
}

template <>
void Spike<CPX>::spy(TCSR<CPX>* matrix, int rows, int cols){
  for(int i = 0; i < rows; ++i){
    for(int j = 0; j < cols; ++j){
      CPX f = get_sparse_matrix_value(matrix,i,j);
      if (real(f) != 0 || imag(f) != 0) {
        std::cout << "\033[0;31m" << "*" << "\033[0;30m";
      } else {
        std::cout << "*";
      }
    }
    std::cout << std::endl; 
  }
}

template <>
void Spike<CPX>::full(TCSR<CPX>* matrix, int rows, int cols){
  for(int i = 0; i < rows; ++i){
    for(int j = 0; j < cols; ++j){
      CPX f = get_sparse_matrix_value(matrix,i,j);
      if (real(f) != 0 || imag(f) != 0) {
        std::cout << "\033[0;31m" << f << "\033[0;30m";
      } else {
        std::cout << "0";     
      }
    }
    std::cout << std::endl; 
  }
}


template <>
void Spike<CPX>::full(CPX* matrix, int rows, int columns) {
  std::cout << "\n[";
  int field_width = 25;
  for(int i = 0; i < rows; ++i) {
    for(int j = 0; j < columns; ++j) {
      int position = (i * columns) + j;
      if (real(matrix[position]) != 0 || imag(matrix[position]) != 0) {
        std::stringstream number_to_print;
        number_to_print << real(matrix[position]) << "+" 
                        << imag(matrix[position]) << "i";
        std::cout << "\033[0;31m" << std::setw(field_width)
                  << number_to_print.str() << "\033[0;30m";
      } else {
        std::cout <<std::setw(field_width)<< "0";     
      }
    }
    std::cout <<";" <<std::endl; 
  }
  std::cout << "]\n";
}

template <>
void Spike<double>::full(double* matrix, int rows, int columns) {
  std::cout << "\n";
  int field_width = 15;
  for(int i = 0; i < rows; ++i) {
    for(int j = 0; j < columns; ++j) {
      int position = (i * columns) + j;
      if (matrix[position] != 0) {
        std::cout << "\033[0;31m" << std::setw(field_width) 
                  << matrix[position] << "\033[0;30m";
      } else {
        std::cout <<std::setw(field_width) << "0";     
      }
    }
    std::cout << std::endl; 
  }
}


template <>
void Spike<CPX>::print_array(CPX* matrix, int size) {
  std::cout << "\n";
  int field_width = 10;
  for(int i = 0; i < size; ++i) {
    std::stringstream number_to_print;
    number_to_print << real(matrix[i]) << "+" 
                    << imag(matrix[i]) << "i";
    std::cout << "\033[0;31m" << std::setw(field_width)
              << number_to_print.str() << "\033[0;30m";
  }
}

/** \brief Templated routines for dense MMM
  *
  * The leading dimension specifies the number of storage locations between
  * elements in the same row (if the matrix is stored in Fortran format),
  * unless the matrix is the second argument B in which case it is the
  * number of columns.
  *
  *
  * For matrices without padding
  *
  *   leading_dimension_A = A_rows / C_rows
  *   leading_dimension_B = B_cols / C_cols (? but seems to be true)
  *   leading_dimension_C = C_rows / A_rows
  *
  * Note: C_rows = A_rows,
  *       A_cols = B_rows,
  *       C_cols = B_cols
  */
template <>
void Spike<CPX>::MMM_dense_f(int A_rows, int A_cols, int B_cols, 
                             CPX prefactor_AB, CPX *A,
                             int leading_dimension_A, CPX *B, 
                             int leading_dimension_B, CPX prefactor_C, CPX *C,
                             int leading_dimension_C) {
   
  c_zgemm('N', 'N', A_rows, B_cols, A_cols, prefactor_AB, A,
          leading_dimension_A, B, leading_dimension_B, prefactor_C, C,
          leading_dimension_C);
  
}
template <>
void Spike<CPX>::MMM_dense_c(int A_rows, int A_cols, int B_cols, 
                             CPX prefactor_AB, CPX *A,
                             int leading_dimension_A, CPX *B, 
                             int leading_dimension_B, CPX prefactor_C, CPX *C,
                             int leading_dimension_C) {
  // A_rows = m, A_cols = k, B_cols = n
  // c-convention: all matrices need transposition so the arguments
  //               have to be reversed.
  // TODO: have to think about that again, probably
  int B_rows = A_cols;
  c_zgemm('T', 'T', B_rows, A_cols, B_cols, prefactor_AB, B,
          leading_dimension_B, A, leading_dimension_A, prefactor_C, C,
          leading_dimension_C);
}
template <>
void Spike<double>::MMM_dense_f(int A_rows, int A_cols, int B_cols, 
                                double prefactor_AB, double *A, 
                                int leading_dimension_A, double *B, 
                                int leading_dimension_B, double prefactor_C,
                                double *C, int leading_dimension_C) {
  c_dgemm('N', 'N', A_rows, B_cols, A_cols, prefactor_AB, A,
          leading_dimension_A, B, leading_dimension_B, prefactor_C, C,
          leading_dimension_C);
}
template <>
void Spike<double>::MMM_dense_c(int A_rows, int A_cols, int B_cols, 
                                double prefactor_AB, double *A, 
                                int leading_dimension_A, double *B, 
                                int leading_dimension_B, double prefactor_C,
                                double *C, int leading_dimension_C) {
  // TODO: have to think about that again, probably
  int B_rows = A_cols;
  c_dgemm('T', 'T', B_rows, A_cols, B_cols, prefactor_AB, B,
          leading_dimension_B, A, leading_dimension_A, prefactor_C, C,
          leading_dimension_C);
}

template <>
void Spike<CPX>::xLACPY(char UPLO, int M, int N, CPX* A, int LDA, CPX* B,
                        int LDB) {
  c_zlacpy(UPLO, M, N, A, LDA, B, LDB);
};

template <>
void Spike<double>::xLACPY(char UPLO, int M, int N, double* A, int LDA, double* B,
                           int LDB) {
  c_dlacpy(UPLO, M, N, A, LDA, B, LDB);
};

template <>
void Spike<CPX>::xGEMM(char TRANSA, char TRANSB, int M, int N, int K,
                       CPX ALPHA, CPX* A, int LDA, CPX* B, int LDB, CPX BETA,
                       CPX* C, int LDC) {
  c_zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
};

template <>
void Spike<double>::xGEMM(char TRANSA, char TRANSB, int M, int N, int K,
                          double ALPHA, double* A, int LDA, double* B, int LDB,
                          double BETA, double* C, int LDC) {
  c_dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
};
