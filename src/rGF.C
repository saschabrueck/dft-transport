#include <vector>
#include <assert.h>
#include <string.h>
#include <sstream>
#include <algorithm>

#include "rGF.H"

/** \brief Constructor
 *  
 *  Constructor for the rGF solver.
 *
 *  \param e_minus_h E (energy) times unity matrix minus the Hamiltonian.
 */
rGF::rGF(TCSR<CPX>* e_minus_h) {
    matrix = e_minus_h;
    fortran_index = 0; // fortran index
}

/// \brief Destructor
rGF::~rGF() {}

/** \brief rGF solving algorithm
 *
 *  Member function to solve for the Green's function in a block by block
 *  format (not in one big sparse matrix).
 *
 *  \param Bmin Vector containing the index of the first element of each diagonal
 *              block in the Hamiltonian (input).
 *
 *  \param Bmax Vector containing the index of the corresponding last element.
 *              NOTE: This are the lines corresponding to the first and last
 *              orbital of the first and last atoms respectively of the 
 *              indexed cell (input).
 *
 *  \param GR Array containing the diagonal blocks of GR (output).
 *
 *
 *  M_TODO:
 *    - consider writing a pardiso class like this rGF class
 *    - consider creating a parent class (like Mathieu)
 *
 *  Potential optimizations:
 *    - having Bmin, Bmax, X_start_index and as gR members
 *    - getting rid of GR and GRNNp1 by writing to the sparse matrix directly
 *    - multithreaded BLAS
 *    - separate thread for GRNNp1 calculation
 *    - multithreaded mat_vec_mult and trans_vec_mat_mult
 *
 */

void rGF::solve_blocks(std::vector<int> Bmin, std::vector<int> Bmax, CPX *GR,
                       CPX *GRNNp1)
{
  int num_blocks = Bmin.size();
  std::vector<int> GR_start_index(num_blocks + 1, Bmin[0]);
  std::vector<int> GRNNp1_start_index(num_blocks, Bmin[0]);

  // Determining block sizes, corresponding start and end indices and
  // biggest block's size
  int current_diagonal=0;
  int next_diagonal=0;
  int largest_diagonal=0;
  for (int i = 0; i < num_blocks-1; ++i){
    current_diagonal = Bmax[i] - Bmin[i] + 1;
    next_diagonal = Bmax[i+1] - Bmin[i+1] + 1;
    GR_start_index[i+1] = GR_start_index[i] +
                          (current_diagonal * current_diagonal);
    GRNNp1_start_index[i+1] = GRNNp1_start_index[i] +
                              (current_diagonal * next_diagonal);
    if (current_diagonal > largest_diagonal) {
      largest_diagonal = current_diagonal;
    }
  }
  GR_start_index[num_blocks] = GR_start_index[num_blocks-1] +
                               (next_diagonal * next_diagonal);
  if (next_diagonal > largest_diagonal) {
    largest_diagonal = next_diagonal;
  }

  //double time_start = get_time(0.0);

  // Allocate working memory
  sparse_CSR = new TCSR<CPX>(largest_diagonal, largest_diagonal *
                              largest_diagonal, fortran_index);
  sparse_CSC = new TCSC<CPX,int>(largest_diagonal, largest_diagonal *
                                  largest_diagonal, fortran_index);
  tmp0 = new CPX[largest_diagonal * largest_diagonal];
  tmp1 = new CPX[largest_diagonal * largest_diagonal];
  tmp2 = new CPX[largest_diagonal * largest_diagonal];
  CPX *gR = new CPX[GR_start_index[num_blocks]];

  // rGF, STAGE 1
  // last element of the diagonal
  calculate_gR_init(num_blocks - 1, Bmin, Bmax, GR_start_index, gR);

  // Recursion upwards along the diagonal:
  // gR_{i-1} = (E - H_{i-1} - T_{i-1,i} * gR_{i} * T_{i,i-1})^{-1}
  for (int block = num_blocks - 2; block > 0; --block) {
    //double sub_time_start = get_time(0.0);
    calculate_sigmaR(block, Bmin, Bmax, GR_start_index, gR);
    calculate_gR_rec(block, Bmin, Bmax, GR_start_index, gR);
    //std::cout << "gR_total_" << block << ": " << get_time(sub_time_start) << "\n";
  }

  // gR_0 goes to tmp2 goes to GR, tmp0<->tmp2 for the later recursion
  int diagonal_length = Bmax[0] - Bmin[0] + 1;
  int diagonal_block_size = diagonal_length * diagonal_length;
  calculate_sigmaR(0, Bmin, Bmax, GR_start_index, gR);
  calculate_gR_rec(0, Bmin, Bmax, GR_start_index, tmp2);
  c_zcopy(diagonal_block_size, tmp2, 1, GR, 1);
  std::swap(tmp0, tmp2);

  //write_mat_c(GR, diagonal_length, diagonal_length, "tests/rGF/GR_0.csv");

  //std::cout << "rGF, stage1 total: " << get_time(time_start) << "\n";
  //time_start = get_time(0.0);



  // rGF STAGE 2
  //
  // calculated GR_{0} = (E-H_{0} - T_{0,1} * gR_{1} * T_{1,0})^{-1} = gR_{0}

  // M_TODO: technically speaking gR is overallocated by the first block as
  // gR[0] gR[GR_start_index[1]] is never actually used.


  /*
  // GRN1_current = 1*GRr_current + GRN1_current
  // M_TODO: since we're not reusing GRN1_current we could allocate tightly
  CPX GRN1_current = new CPX[largest_diagonal * largest_diagonal]; // overallocation
  set_to_zero(largest_diagonal * largest_diagonal, GRN1_current);
  c_zcopy(diagonal_length * diagonal_length, GRr_current, 1, GRN1_current, 1);
  
  
  // tmp = GRN1_current * Gammal
  CPX tmp1 = new CPX[largest_diagonal * largest_diagonal];
  c_zgemm('N', 'N', diagonal_length, diagonal_length, diagonal_length,
          CPX(1.0, 0.0), GRN1_current, diagonal_length, Gammal, 
          diagonal_length, CPX(0.0, 0.0), tmp1, diagonal_length);
  
  // Zl = tmp1 * complexconjugate(GRN1_current);
  CPX Zl = new CPX[largest_diagonal * largest_diagonal];
  c_zgemm('N', 'C', diagonal_length, diagonal_length, diagonal_length,
          CPX(1.0, 0.0), M, diagonal_length, GRN1_current, diagonal_length,
          CPX(0.0, 0.0), Zl, diagonal_length);
  */
   

  // recursion along the diagonal (downwards)
  //
  // calculate
  // GR_{i} = gR_{i} + gR_{i} * T_{i,i-1} * GR_{i-1} * T_{i-1,i} * gR_{i}
  // GR_{i-1,i} = -1 * GR_{i-1} * T_{i-1,i} * gR_{i}
  //
  // 1 <= i <= num_blocks
  
  for (int block = 1; block < num_blocks; ++block) {
    diagonal_length = Bmax[block] - Bmin[block] + 1;
    diagonal_block_size = GR_start_index[block] - GR_start_index[block - 1];
    calculate_GR(block, Bmin, Bmax, &gR[GR_start_index[block]],
                 &GR[GR_start_index[block]],
                 &GRNNp1[GRNNp1_start_index[block -1]]);

    /*std::stringstream filename;
    filename.str("");
    filename << "tests/rGF/GR_" << block << ".csv";
    write_mat_c(&GR[GR_start_index[block]], diagonal_length, diagonal_length, 
              filename.str());
    filename.str("");
    filename << "tests/rGF/GRNNp1_" << block << ".csv";
    write_mat_c(&GRNNp1[GRNNp1_start_index[block - 1]], Bmax[block - 1] - 
              Bmin[block - 1] + 1, diagonal_length, filename.str());*/
  }

  delete sparse_CSR;
  delete sparse_CSC;
  delete[] tmp0;
  delete[] tmp1;
  delete[] tmp2;
  delete[] gR;

  //std::cout << "rGF, stage2 total: " << get_time(time_start) << "\n";

}



/** \brief Extraction of a diagonal block from the main sparse matrix
 *
 * This routine extracts a diagonal block from the instance's matrix 
 * (this->matrix) and stores it in diagonal_block as full matrix.
 *
 * \param block Index of the block to be extracted (numbering starts with 0)
 *              (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param diagonal_block The full matrix in which the diagonal block is to be
 *                       stored (output).
 */
void rGF::get_diagonal_block(int block, std::vector<int> Bmin,
                                 std::vector<int> Bmax, CPX *diagonal_block)
{
  int diagonal_length = Bmax[block] - Bmin[block] + 1;
  assert(diagonal_length <= Bmax[block] + 1);

  set_to_zero(diagonal_length * diagonal_length, diagonal_block);

  for (int i_absolute = Bmin[block]; i_absolute <= Bmax[block]; 
       ++i_absolute) {
    int i_relative = i_absolute - Bmin[block];
    for (int j_absolute = matrix->edge_i[i_absolute] - matrix->findx;
         j_absolute < matrix->edge_i[i_absolute+1] - matrix->findx;
         ++j_absolute) { 
      int j_relative = matrix->index_j[j_absolute] - matrix->findx - 
                                Bmin[block];
      if ((j_relative >=  0) && (j_relative < diagonal_length)) {
        diagonal_block[(i_relative * diagonal_length) + j_relative] =
                                                         matrix->nnz[j_absolute];
      }
    }
  }
}


/** \brief Extraction of two offdiagonal blocks from the main sparse matrix
 *
 * This routine extracts the two offdiagonal blocks T_{i,i+1} and T_{i+1,i}
 * from the instance's matrix (this->matrix) and stores them in the sparse
 * matrices T_iip1 (CSC) and T_ip1i (CSR) respectively. T_{i,i+1} is the 
 * conjugate of T_{i+1,i}, T_{i+1,i} will be taken as source for construction 
 * of both blocks.
 *
 * \param block Index i for extraction of T_{i,i+1} and T_{i+1,i} (numbering
 *              starts with 0) (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param T_iip1 Sparse (CSC) matrix to contain T_{i,i+1} (output).
 *
 * \param T_iip1 Sparse (CSR) matrix to contain T_{i+1,i} (output).
 */
void rGF::get_Tip1i_csc_Tiip1_csr(int block, std::vector<int> Bmin,
                                  std::vector<int> Bmax, TCSC<CPX,int> *T_ip1i, 
                                  TCSR<CPX> *T_iip1)
{
  assert(block < (int)(Bmin.size() - 1));
  assert(block >= 0);

  int current_diagonal = Bmax[block] - Bmin[block] + 1;

  T_iip1->size = current_diagonal; // CSR->size = number of rows
  T_ip1i->size = current_diagonal; // CSC->size = number of columns

  int new_values_index = 0;

  for (int row_absolute = Bmin[block]; row_absolute <= Bmax[block] ; 
      ++row_absolute) {

    int row = row_absolute - Bmin[block];
    T_iip1->index_i[row] = 0;
    T_ip1i->index_j[row] = 0;

    int row_start = matrix->edge_i[row_absolute] - matrix->findx;
    int next_row_start = matrix->edge_i[row_absolute + 1] - matrix->findx;

    for (int value_index = row_start; value_index < next_row_start;
         ++value_index) {

      int values_column = matrix->index_j[value_index] - matrix->findx;

      if ((Bmin[block+1] <= values_column) && (values_column <= Bmax[block+1])) {

        int new_values_column = values_column - Bmin[block+1];

        T_iip1->index_j[new_values_index] = new_values_column + fortran_index;
        T_iip1->nnz[new_values_index] = matrix->nnz[value_index];
        T_iip1->index_i[row]++;

        T_ip1i->index_i[new_values_index] = new_values_column + fortran_index;
        T_ip1i->nnz[new_values_index] = matrix->nnz[value_index];
        T_ip1i->index_j[row]++;

        new_values_index++;
      } 
      else if (values_column > Bmax[block+1]) {
        break;
      }
    }
  }

  T_iip1->n_nonzeros = new_values_index;
  T_ip1i->n_nonzeros = new_values_index;
  T_iip1->get_row_edge();
  T_ip1i->get_column_edge();
}


/** \brief Extraction of two offdiagonal blocks from the main sparse matrix
 *
 * This routine extracts the two offdiagonal blocks T_{i,i+1} and T_{i+1,i}
 * from the instance's matrix (this->matrix) and stores them in the sparse
 * matrices T_ip1i (CSR) and T_iip1 (CSC) respectively. T_{i+1,i} is the 
 * conjugate of T_{i,i+1}. T_{i,i+1} will be taken as source for construction
 * of both blocks.
 *
 * \param block Index i for extraction of T_{i,i+1} and T_{i+1,i} (numbering
 *              starts with 0) (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param T_ip1i Sparse (CSR) matrix to contain T_{i,i+1} (output).
 *
 * \param T_iip1 Sparse (CSC) matrix to contain T_{i+1,i} (output).
 */
void rGF::get_Tip1i_csr_Tiip1_csc(int block, std::vector<int> Bmin,
                                  std::vector<int> Bmax, TCSR<CPX> *T_ip1i, 
                                  TCSC<CPX,int> *T_iip1)
{
  assert(block < (int)(Bmin.size() - 1));
  assert(block >= 0);

  int next_diagonal = Bmax[block+1] - Bmin[block+1] + 1;

  T_ip1i->size = next_diagonal; // CSR->size = number of rows     
  T_iip1->size = next_diagonal; // CSC->size = number of columns

  int new_values_index = 0;

  for (int row_absolute = Bmin[block+1]; row_absolute <= Bmax[block+1];
       ++row_absolute) {

    int row = row_absolute - Bmin[block+1];
    T_iip1->index_j[row] = 0;
    T_ip1i->index_i[row] = 0;

    int row_start = matrix->edge_i[row_absolute] - matrix->findx;
    int next_row_start = matrix->edge_i[row_absolute + 1] - matrix->findx;

    for (int value_index = row_start; value_index < next_row_start;
         value_index++) {

      int values_column = matrix->index_j[value_index] - matrix->findx;

      if ((Bmin[block] <= values_column) && (values_column <= Bmax[block])) {

        int new_values_column = values_column - Bmin[block];

        T_ip1i->index_j[new_values_index] = new_values_column + fortran_index;
        T_ip1i->nnz[new_values_index] = matrix->nnz[value_index];
        T_ip1i->index_i[row]++;

        T_iip1->index_i[new_values_index] = new_values_column + fortran_index;
        T_iip1->nnz[new_values_index] = matrix->nnz[value_index];
        T_iip1->index_j[row]++;

        ++new_values_index;
      }
      else if (values_column > Bmax[block+1]) {
        break;
      }
    }
  }

  T_iip1->n_nonzeros = new_values_index;
  T_ip1i->n_nonzeros = new_values_index;
  T_iip1->get_column_edge();
  T_ip1i->get_row_edge();
}


/** \brief Function to calculate a block of gR
 *
 * Function to calcuate gR = (matrix_{i,i}-SigR)^{-1} for a given subblock i,i
 * of the main matrix.
 *
 * \pre sigmaR is stored in tmp0.
 *
 * \param block The index of the block to calculate (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param GR_start_index Vector of block start indices within the matrix
 *                       array (input).
 *
 * \param gR Array containing gR (output).
 *
 */
void rGF::calculate_gR_rec(int block, std::vector<int> Bmin, 
                           std::vector<int> Bmax, 
                           std::vector<int> GR_start_index, CPX *gR)
{

  //double start_time = get_time(0.0);

  int diagonal_length = Bmax[block] - Bmin[block] + 1;
  int diagonal_block_size = diagonal_length * diagonal_length;

  CPX *sigmaR = tmp0;
  CPX *inverted_gR = tmp1;
  get_diagonal_block(block, Bmin, Bmax, inverted_gR);

  // inverted_gR += -1 * sigmaR:
  c_zaxpy(diagonal_block_size, CPX(-1.0, 0.0), sigmaR, 1, inverted_gR, 1);

  // gR = inverted_gR^{-1}:
  int *pivot_vector = new int[diagonal_length];
  int status;
  c_zgetrf(diagonal_length, diagonal_length, inverted_gR, diagonal_length,
           pivot_vector, &status);
  set_to_unity(diagonal_length, &gR[GR_start_index[block]]);
  c_zgetrs('N', diagonal_length, diagonal_length, inverted_gR, diagonal_length,
           pivot_vector, &gR[GR_start_index[block]], diagonal_length, &status);

  /*std::stringstream gR_filename;
  gR_filename << "tests/rGF/gR_" << block << ".csv";
  write_mat_c(&gR[GR_start_index[block]], diagonal_length, diagonal_length,
            gR_filename.str());*/

  delete[] pivot_vector;

  //std::cout << "gR_inv_" << block << ": " << get_time(start_time) << "\n";
}

/** \brief Function to calculate a block of gR
 *
 * Function to calcuate gR = (matrix_{i,i})^{-1} for a given subblock i,i
 * of the main matrix. To be used for the lowest block of gR if the self
 * energies are already applied to the matrix.
 *
 * \param block The index of the block to calculate (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param GR_start_index Vector of block start indices within the matrix
 *                       array (input).
 *
 * \param gR Array containing gR (output).
 *
 */
void rGF::calculate_gR_init(int block, std::vector<int> Bmin, 
                            std::vector<int> Bmax, 
                            std::vector<int> GR_start_index, CPX *gR)
{

  //double start_time = get_time(0.0);

  int diagonal_length = Bmax[block] - Bmin[block] + 1;

  CPX *inverted_gR = tmp0;
  get_diagonal_block(block, Bmin, Bmax, inverted_gR);

  // gR = inverted_gR^{-1}:
  int *pivot_vector = new int[diagonal_length];
  int status;
  c_zgetrf(diagonal_length, diagonal_length, inverted_gR, diagonal_length,
           pivot_vector, &status);
  set_to_unity(diagonal_length, &gR[GR_start_index[block]]);
  c_zgetrs('N', diagonal_length, diagonal_length, inverted_gR, diagonal_length,
           pivot_vector, &gR[GR_start_index[block]], diagonal_length, &status);

  /*std::stringstream gR_filename;
  gR_filename << "tests/rGF/gR_" << block << ".csv";
  write_mat_c(&gR[GR_start_index[block]], diagonal_length, diagonal_length,
            gR_filename.str());*/

  delete[] pivot_vector;

  //std::cout << "gR_inv_" << block << ": " << get_time(start_time) << "\n";
}


/** \brief Function to calculate sigmaR
 *
 * This function calculates sigmaR as
 *
 *  sigmaR_{i} = T_{i,i+1} * gR_{i+1} * T_{i+1,i}
 *
 * \post sigmaR will be stored in tmp0 after return.
 *
 * \param block Index of the block of gR to calculate (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param GR_start_index Vector of block start indices within the matrix
 *                       array (input).
 *
 * \param gR The matrix containing the already calculated blocks of gR
 *           (input).
 *
 */
void rGF::calculate_sigmaR(int block, std::vector<int> Bmin,
                           std::vector<int> Bmax,
                           std::vector<int> GR_start_index, CPX *gR)
{

  int result_size;
  int precursor_size;

  //double start_time = get_time(0.0);
  
  assert(block < (int)(Bmin.size() - 1));
  assert(block >= 0);

  result_size = Bmax[block] - Bmin[block] + 1;
  precursor_size = Bmax[block+1] - Bmin[block+1] + 1;
  int oversize = result_size;
  if ( oversize < precursor_size) oversize = precursor_size;

  TCSR<CPX> *T_iip1 = sparse_CSR;
  TCSC<CPX,int> *T_ip1i = sparse_CSC;
  CPX *sigmaR = tmp0;
  CPX *T_iip1_gR_ip1 = tmp1;

  get_Tip1i_csc_Tiip1_csr(block, Bmin, Bmax, T_ip1i, T_iip1);
  // T_iip1_gR_ip1 =  T_{i,i+1} * gR_{i+1,i+1}
  transpose_full_matrix(&gR[GR_start_index[block+1]], precursor_size, 
                        precursor_size, sigmaR);
  T_iip1->trans_mat_vec_mult(sigmaR, T_iip1_gR_ip1,
                            precursor_size, precursor_size);
  // sigmaR = T_iip1_gR_ip1 * T_ip1i
  T_ip1i->vec_mat_mult(T_iip1_gR_ip1, sigmaR, result_size);

  //std::cout << "sigmaR_" << block << ": " << get_time(start_time) << "\n";

}


/** \brief Function to calculate GR_{i}
 *
 * This function calculates GR_{i} as 
 *
 *   GR_{i} = gR_{i} + gR_{i} * T_{i,i-1} * GR_{i-1} * T{i-1,i} * gR_{i}
 *
 * and GR_{i-1,i} as
 *
 *   GR_{i-1,i} = -1 * GR_{i-1} * T_{i-1,i} * gR_{i}
 *
 * \pre tmp0 contains GR_{i-1}
 *
 * \post tmp0 contains GR_{i}
 *
 * \post tmp1 contains GR_{i-1, i}
 *
 * \param block Index i of the block GR_{i} to calculate (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param gR The matrix containing gR_{i} (input).
 *
 * \note: There's commented out code to calculate blocks of the first column of
 * GR.
 *
 */
void rGF::calculate_GR(int block, std::vector<int> Bmin, std::vector<int> Bmax, 
                       CPX *gR, CPX *GR, CPX *GRNNp1)
{

  //double start_time = get_time(0.0);

  assert(block > 0 && block < (int)Bmin.size());
  int result_size = Bmax[block] - Bmin[block] + 1;
  int precursor_size = Bmax[block-1] - Bmin[block-1] + 1;
  int oversize = result_size;
  if (oversize < precursor_size) {oversize = precursor_size;}

  TCSC<CPX,int> *T_im1i = sparse_CSC;
  TCSR<CPX> *T_iim1 = sparse_CSR;
  CPX *GR_block = tmp0;
  CPX *tmp = tmp1;

  get_Tip1i_csr_Tiip1_csc(block - 1, Bmin, Bmax, T_iim1, T_im1i);

  /* Calculation of GRN1_block
  *if (we_want_to_calc_GRN1) {
  *  int first_size = Bmax[0] - Bmin[0] + 1;
  *  transpose_full_matrix(GRN1_block, precursor_size, first_size, 
  *                        GRN1_transposed);
  *  T_iim1->trans_mat_vec_mult(GRN1_block_transposed, M_tmp, first_size, 
  *                            precursor_size);
  *  c_zgemm('N', 'N', result_size, first_size, result_size, CPX(-1.0, 0.0),
  *          gR, result_size, M_tmp, result_size, CPX(0.0, 0.0), GRN1_block,
  *          result_size);
  *}
  */

  // GR_{i} = gR_{i} + gR_{i} * T_{i,i-1} * GR_{i-1} * T_{i-1,i} * gR_{i}
  // GR_{i-1,i} =                      -1 * GR_{i-1} * T_{i-1,i} * gR_{i}

  // change to fortran memory layout
  transpose_full_matrix(GR_block, precursor_size, precursor_size, tmp);

  // Check preconditions
  /*std::stringstream filename;
  filename.str(""); filename << "tests/rGF/GR_" << block-1 << "_it" << block << ".csv";
  write_mat_c(tmp, precursor_size, precursor_size, filename.str());
  filename.str(""); filename << "tests/rGF/gR_" << block << "_it" << block << ".csv";
  write_mat_c(gR, result_size, result_size, filename.str());
  filename.str(""); filename << "tests/rGF/T_" << block-1 << block << ".csv";
  const std::string &filename_str = filename.str();
  const char *cfilename = filename_str.c_str();
  T_im1i->write(cfilename);*/

  // GR_block = tmp * T_{i-1,i}, dim: precursor_size*result_size (but we over-
  // allocate). The routine produces Fortran-style matrices, i.e. the dimensions
  // are result_size*precursor_size
  // NOTE: vec_mat_mult(in, out, number_of_rows for 'in')
  T_im1i->vec_mat_mult(tmp, GR_block, precursor_size);
  /*filename.str(""); filename << "tests/rGF/A_" << block << ".csv";
  write_mat_c(GR_block, result_size, precursor_size, filename.str());*/

  // GRNNp1[..] = GR_block * gR_{i} = GR_{i-1} * T_{i-1,i} * gR_{i}, 
  // dim: precursor_size*result_size
  // NOTE: LDA is effectively the stride of the matrix A in BLAS (result_size
  // in our case)
  c_zgemm('N', 'T', result_size, precursor_size, result_size, CPX(1.0, 0.0),
          gR, result_size, GR_block, precursor_size, CPX(0.0, 0.0), GRNNp1,
          result_size);
  /*filename.str(""); filename << "tests/rGF/B_" << block << ".csv";
  write_mat_c(GRNNp1, precursor_size, result_size, filename.str());*/

  // tmp = T_{i,i-1} * GRNNp1 = T_{i,i-1} * GR_block * T_{i-1,i} * gR_{i},
  // dim: result_size*result_size
  // NOTE: trans_mat_vec_mult(in, out, cols_of_non_transposed_in, <ignored>)
  T_iim1->trans_mat_vec_mult(GRNNp1, tmp, result_size, result_size);
  /*filename.str(""); filename << "tests/rGF/C_" << block << ".csv";
  write_mat_c(tmp, result_size, result_size, filename.str());*/

  // GRNNp1 *= -1, dim: precursor_size*result_size
  // M_TODO: we could also account for the factor at result construction time
  c_zscal(precursor_size * result_size, CPX(-1.0, 0.0), GRNNp1, 1);
  /*filename.str(""); filename << "tests/rGF/D_" << block << ".csv";
  write_mat_c(GRNNp1, precursor_size, result_size, filename.str());*/

  // GR_block = gR * tmp = gR_{i} * T_{i,i-1} * GR_block * T_{i-1,i} * gR_{i},
  // dim: result_size*result_size
  c_zgemm('T', 'N', result_size, result_size, result_size, CPX(1.0, 0.0),
          tmp, result_size, gR, result_size, CPX(0.0, 0.0), GR_block, 
          result_size);
  /*filename.str(""); filename << "tests/rGF/E_" << block << ".csv";
  write_mat_c(GR_block, result_size, result_size, filename.str());*/

  // GR_block += gR = gR_{i} + gR_{i} * T_{i,i-1} * GR_block * T_{i-1,i} * gR_{i},
  // dim: result_size*result_size
  c_zaxpy(result_size * result_size, CPX(1.0, 0.0), gR, 1, GR_block, 1);
  c_zcopy(result_size * result_size, GR_block, 1, GR, 1);
  /*filename.str(""); filename << "tests/rGF/F_" << block << ".csv";
  write_mat_c(GR_block, result_size, result_size, filename.str());*/

  //std::cout << "GR_" << block << ": " << get_time(start_time) << "\n"; */

  return;

}


/** \brief Function to transpose (dense) matrices
 *
 * \param normal The rows*cols matrix to be transposed (input).
 *
 * \param rows The number of rows of the matrix (input).
 *
 * \param cols The number of columns of the matrix (input).
 *
 * \param transposed The matrix in transposed form (a cols*rows matrix) (output).
 */
void rGF::transpose_full_matrix(CPX *normal, int rows, int cols, CPX *transposed)
{
  for (int j=0; j < cols; ++j) {
    c_zcopy(rows, &normal[j*rows], 1, &transposed[j], cols);
  }
}


/** \brief Set a full matrix to unity
 *
 * \param diagonal_length The diagonal length of the square matrix (input)
 *
 * \param unity The array containing the matrix to set to unity (input/output)
 */
void rGF::set_to_unity(int diagonal_length, CPX *unity)
{
  set_to_zero(diagonal_length * diagonal_length, unity);
  CPX *ones = new CPX[diagonal_length];
  for (int i = 0; i < diagonal_length; ++i) {
    ones[i] = CPX(1.0, 0.0);
  }
  int N = diagonal_length;
  c_zcopy(N, ones, 1, unity, N+1);
}
