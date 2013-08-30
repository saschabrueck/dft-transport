#include "tmprGF.H"
#include <cstddef>
#include <vector>
#include <algorithm>    // swap in C++98
#include <utility>      // swap in C++11


namespace tmprGF {

/** \brief Function to invert a CSR matrix using rGF
 * 
 *  This function inverts sparse block tridiagonal matrices using the rGF
 *  algorithm.
 *  The sparsity pattern in the input is preserved in the output.
 *
 *  \param *matrix The matrix to be inverted (in-place)
 *  \pre The given matrix is a non-singular matrix (complex or real) 
 *       in CSR format
 *  \post The input matrix has been replaced by its inverse
 */
void sparse_invert(TCSR<CPX> *matrix, std::vector<int> Bmax) {

  std::vector<int> Bmin(Bmax.size(), 0);

  for (unsigned int i = 0; i < Bmin.size() - 1; ++i) {
    Bmin[i + 1] = Bmax[i] + 1;
  }

  // calculate memory requirements and positions:
  int num_blocks = Bmin.size();
  std::vector<int> GR_start_index(num_blocks + 1, Bmin[0]);
  std::vector<int> GRNNp1_start_index(num_blocks, Bmin[0]);
  int GR_size = 0;
  int GRNNp1_size = 0;
  int current_diagonal = 0;
  int next_diagonal = 0;
  int largest_diagonal = 0;
  for (int i = 0; i < num_blocks - 1; ++i) {
    current_diagonal = Bmax[i] - Bmin[i] + 1;
    next_diagonal = Bmax[i+1] - Bmin[i+1] + 1;
    GR_size += current_diagonal * current_diagonal;
    GRNNp1_size += current_diagonal * next_diagonal;
    GR_start_index[i+1] = GR_start_index[i] +
                          (current_diagonal * current_diagonal);
    GRNNp1_start_index[i+1] = GRNNp1_start_index[i] +
                              (current_diagonal * next_diagonal);
    if (current_diagonal > largest_diagonal) {
      largest_diagonal = current_diagonal;
    }
  }
  GR_size += next_diagonal * next_diagonal;
  GR_start_index[num_blocks] = GR_start_index[num_blocks-1] +
                               (next_diagonal * next_diagonal);
  if (current_diagonal > largest_diagonal) {
    largest_diagonal = current_diagonal;
  }

  // solve
  CPX *GR = new CPX[GR_size];
  CPX *GRNNp1 = new CPX[GRNNp1_size];
  rGF *myRGF = new rGF(matrix);
  myRGF->solve_blocks(Bmin, Bmax, GR, GRNNp1);

  // this here serves as a cache for the lower diagonal blocks
  CPX *T10_cur = new CPX[largest_diagonal*largest_diagonal];
  CPX *T10_next = new CPX[largest_diagonal*largest_diagonal];

  // append indices to Bmin to avoid special handling of last blocks in
  // the following code.
  Bmin.push_back(Bmax[num_blocks - 1] + 1);
  Bmin.push_back(Bmax[num_blocks - 1] + 1);

  std::stringstream filename("");

  for (int block_row = 0; block_row < num_blocks; ++block_row) {

    std::cout << block_row << ": Going from row " << Bmin[block_row] << " to "
              << Bmax[block_row] << ", up to column " << Bmin[block_row + 2]
              << "\n";

    for (int row = Bmin[block_row]; row <= Bmax[block_row]; ++row) {
      int nnz_start = matrix->edge_i[row] - matrix->findx;
      int nnz_end = matrix->edge_i[row+1] - matrix->findx;
      for (int innz = nnz_start; innz < nnz_end; ++innz) {
        int column = matrix->index_j[innz] - matrix->findx;

        // Fetch the element at position (row, column)
        // and possibly store it for later use

        // M_TODO: possibly test if we're before the left block
        // and raise an exception.

        if (column < Bmin[block_row]) {
          // left to diagonal
          int j_T10_cur = column - Bmin[block_row - 1];
          int i_T10_cur = row - Bmin[block_row];
          int T10_cur_row_length = Bmin[block_row] - Bmin[block_row - 1];
          int T10_cur_pos = (i_T10_cur * T10_cur_row_length) + j_T10_cur;
          matrix->nnz[innz] = T10_cur[T10_cur_pos];
        }
        else if (column < Bmin[block_row + 1]) {
          // diagonal
          int j_GR = column - Bmin[block_row];
          int i_GR = row - Bmin[block_row];
          int GR_row_length = Bmin[block_row+1] - Bmin[block_row];
          int GR_pos = GR_start_index[block_row] + (i_GR * GR_row_length) +
                                        j_GR;
          matrix->nnz[innz] = GR[GR_pos];
        }
        else if (column < Bmin[block_row + 2]) {
          // right to diagonal
          int j_GRNNp1 = column - Bmin[block_row + 1];
          int i_GRNNp1 = row - Bmin[block_row];
          int GRNNp1_row_length = Bmin[block_row + 2] - Bmin[block_row + 1];
          int GRNNp1_pos = GRNNp1_start_index[block_row] +
                           (i_GRNNp1 * GRNNp1_row_length) + j_GRNNp1;
          matrix->nnz[innz] = GRNNp1[GRNNp1_pos];
          // store this element in T10_next for the following
          // block_row (transposed)
          int T10_next_row_length = Bmin[block_row + 1] - Bmin[block_row];
          int T10_next_pos = (T10_next_row_length * j_GRNNp1) + i_GRNNp1;
          T10_next[T10_next_pos] = GRNNp1[GRNNp1_pos];
        }
        else {
          // M_TODO: outside right block, raise exception
        }
      }
    }
    std::swap(T10_cur, T10_next);
  }




  delete[] GR;
  delete[] GRNNp1;
  delete[] T10_cur;
  delete[] T10_next;

}

} // namespace
