/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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

/* templated interfaces to BLAS/LAPACK routines */
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
