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
 *  CHECKPOINT:
 *    - consider writing a pardiso class like this rGF class
 *    - consider creating a parent class (like Mathieu)
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
  int largest_block=0;
  for (int i = 0; i < num_blocks-1; ++i){
    current_diagonal = Bmax[i] - Bmin[i] + 1;
    next_diagonal = Bmax[i+1] - Bmin[i+1] + 1;
    GR_start_index[i+1] = GR_start_index[i] +
                          (current_diagonal * current_diagonal);
    GRNNp1_start_index[i+1] = GRNNp1_start_index[i] +
                              (current_diagonal * next_diagonal);
    if (current_diagonal > largest_block) {
      largest_block = current_diagonal;
    }
  }
  GR_start_index[num_blocks] = GR_start_index[num_blocks-1] +
                               (next_diagonal * next_diagonal);
  if (next_diagonal > largest_block) {
    largest_block = next_diagonal;
  }

  //double time_start = get_time(0.0);

  // rGF, STAGE 1
  CPX *gR = new CPX[GR_start_index[num_blocks]];

  // last element of the diagonal
  calculate_gR(num_blocks - 1, Bmin, Bmax, GR_start_index, gR);

  // Recursion upwards along the diagonal:
  // gR_{i-1} = (E - H_{i-1} - T_{i-1,i} * gR_{i} * T_{i,i-1})^{-1}
  CPX *sigmaR = new CPX[largest_block * largest_block];
  for (int block = num_blocks - 2; block > 0; --block) {
    //double sub_time_start = get_time(0.0);
    calculate_sigmaR(block, (char *)"r", Bmin, Bmax, GR_start_index, gR, sigmaR);
    calculate_gR(block, Bmin, Bmax, GR_start_index, sigmaR, gR);
    //std::cout << "gR_total_" << block << ": " << get_time(sub_time_start) << "\n";
  }


  // rGF STAGE 2
  //
  // calculate GR_{0} = (E-H_{0} - T_{0,1} * gR_{1} * T_{1,0})^{-1}

  // 1: calculate sigmaR[0]
  int diagonal_length = Bmax[0] - Bmin[0] + 1;
  int diagonal_block_size = diagonal_length * diagonal_length;
  calculate_sigmaR(0, (char *)"r", Bmin, Bmax, GR_start_index, gR, sigmaR);

  // 2: GR_current = gR[0]
  CPX *GR_current = new CPX[largest_block * largest_block]; // overallocation
  set_to_unity(largest_block, GR_current);
  calculate_gR(0, Bmin, Bmax, GR_start_index, sigmaR, GR_current);

  //std::cout << "rGF, stage1 total: " << get_time(time_start) << "\n";
  //time_start = get_time(0.0);

  // M_TODO: technically speaking gR is overallocated by the first block as
  // gR[0] gR[GR_start_index[1]] is never actually used.

  // 6: GR_{0} = GR_current
  c_zcopy(diagonal_block_size, GR_current, 1, GR, 1);

  // recursion along the diagonal (downwards)
  //
  // calculate
  // GR_{i} = gR_{i} + gR_{i} * T_{i,i-1} * GR_{i-1} * T_{i-1,i} * gR_{i}
  // GR_{i-1,i} = -1 * GR_{i-1} * T_{i-1,i} * gR_{i}
  //
  // 1 <= i <= num_blocks
  
  CPX *GRNNp1_current = new CPX[largest_block * largest_block]; // overallocation
  for (int block = 1; block < num_blocks; ++block) {
  //for (int block = 1; block < num_blocks-1; ++block) {
    diagonal_length = Bmax[block] - Bmin[block] + 1;
    diagonal_block_size = GR_start_index[block] - GR_start_index[block - 1];
    int upper_off_diagonal_size = GRNNp1_start_index[block] - 
                                  GRNNp1_start_index[block - 1];
    calculate_GR(block, Bmin, Bmax, &gR[GR_start_index[block]], GR_current, 
                 GRNNp1_current);

    c_zcopy(diagonal_block_size, GR_current, 1,
            &GR[GR_start_index[block]], 1);
    c_zcopy(upper_off_diagonal_size, GRNNp1_current, 1,
            &GRNNp1[GRNNp1_start_index[block - 1]], 1);

  }

  delete[] gR;
  delete[] sigmaR;
  delete[] GR_current;
  delete[] GRNNp1_current;

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
                                 std::vector<int> Bmax,
                                 CPX *diagonal_block)
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
                                     std::vector<int> Bmax,
                                     TCSC<CPX,int> *T_ip1i, TCSR<CPX> *T_iip1)
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
        //this works for SmallNW:
        //T_ip1i->nnz[new_values_index] = conj(matrix->nnz[value_index]);
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
                                     std::vector<int> Bmax,
                                     TCSR<CPX> *T_ip1i, TCSC<CPX,int> *T_iip1)
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
        //MARKER: this works for SmallNW:
        //T_iip1->nnz[new_values_index] = conj(matrix->nnz[value_index]);
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
 * \param block The index of the block to calculate (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param GR_start_index Vector of block start indices within the matrix
 *                       array (input).
 *
 * \param sigmaR Boundary self-energy (optional: if not given, it is assumed
 *               that the main matrix (this->matrix) contains the sigmaR
 *               already) (input).
 *
 * \param gR Array containing gR (output).
 *
 */
void rGF::calculate_gR(int block, std::vector<int> Bmin, std::vector<int> Bmax, 
                       std::vector<int> GR_start_index, CPX *sigmaR, CPX *gR)
{

  //double start_time = get_time(0.0);

  int diagonal_length = Bmax[block] - Bmin[block] + 1;
  int diagonal_block_size = diagonal_length * diagonal_length;

  CPX *inverted_gR = new CPX[diagonal_block_size];
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

  delete[] inverted_gR;
  delete[] pivot_vector;

  //std::cout << "gR_inv_" << block << ": " << get_time(start_time) << "\n";
}

void rGF::calculate_gR(int block, std::vector<int> Bmin, std::vector<int> Bmax, 
                       std::vector<int> GR_start_index, CPX *gR)
{

  //double start_time = get_time(0.0);

  int diagonal_length = Bmax[block] - Bmin[block] + 1;
  int diagonal_block_size = diagonal_length * diagonal_length;

  CPX *inverted_gR = new CPX[diagonal_block_size];
  get_diagonal_block(block, Bmin, Bmax, inverted_gR);

  // gR = inverted_gR^{-1}:
  int *pivot_vector = new int[diagonal_length];
  int status;
  c_zgetrf(diagonal_length, diagonal_length, inverted_gR, diagonal_length,
           pivot_vector, &status);
  set_to_unity(diagonal_length, &gR[GR_start_index[block]]);
  c_zgetrs('N', diagonal_length, diagonal_length, inverted_gR, diagonal_length,
           pivot_vector, &gR[GR_start_index[block]], diagonal_length, &status);

  delete[] inverted_gR;
  delete[] pivot_vector;

  //std::cout << "gR_inv_" << block << ": " << get_time(start_time) << "\n";
}


/** \brief Function to calculate sigmaR
 *
 * This function calculates sigmaR either from the right or from the left
 * side.
 *
 * For the calculation from the right it calculates
 *
 *  sigmaR_{i} = T_{i,i+1} * gR_{i+1} * T_{i+1,i}
 *
 * whereas from the left side
 *
 *  sigmaR_{i} = T_{i,i-1} * gR_{i-1} * T_{i-1,i}
 *
 *
 * \param block Index of the block of gR to calculate (input).
 *
 * \param side Indicator of the side to calculate gR from, valid values
 *             are "r" for right or "l" for left (input).
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
 * \param sigmaR The matrix containing sigmaR (output).
 *
 */
void rGF::calculate_sigmaR(int block, char *side,
                           std::vector<int> Bmin,
                           std::vector<int> Bmax,
                           std::vector<int> GR_start_index, CPX *gR, CPX *sigmaR)
{
  // It seems as if T_iip1 and T_ip1i would not be used afterwards
  // and therefore can be allocated in this function.

  int result_size;
  int precursor_size; // size of the predecessor in the recusion (size of
                      // block i+1 for side='r' and size of block i-1 for 
                      // side='l')

  //double start_time = get_time(0.0);
  
  if (!strcmp(side, "r")) {

    assert(block < (int)(Bmin.size() - 1));
    assert(block >= 0);

    result_size = Bmax[block] - Bmin[block] + 1;
    precursor_size = Bmax[block+1] - Bmin[block+1] + 1;
    int oversize = result_size;
    if ( oversize < precursor_size) oversize = precursor_size;

    TCSR<CPX> *T_iip1 = new TCSR<CPX>(oversize, oversize * oversize,
                                      fortran_index);
    TCSC<CPX,int> *T_ip1i = new TCSC<CPX,int>(oversize , oversize * oversize, 
                                              fortran_index);
    CPX *gR_ip1_transposed = new CPX[oversize * oversize];
    CPX *T_iip1_gR_ip1 = new CPX[oversize * oversize];

    get_Tip1i_csc_Tiip1_csr(block, Bmin, Bmax, T_ip1i, T_iip1);
    // T_iip1_gR_ip1 =  T_{i,i+1} * gR_{i+1,i+1}
    transpose_full_matrix(&gR[GR_start_index[block+1]], precursor_size, 
                          precursor_size, gR_ip1_transposed);
    T_iip1->trans_mat_vec_mult(gR_ip1_transposed, T_iip1_gR_ip1,
                              precursor_size, precursor_size);
    // sigmaR = T_iip1_gR_ip1 * T_ip1i
    T_ip1i->vec_mat_mult(T_iip1_gR_ip1, sigmaR, result_size);

    delete T_iip1;
    delete T_ip1i;
    delete[] gR_ip1_transposed;
    delete[] T_iip1_gR_ip1;

  }
  else if (!strcmp(side, "l")) {

    assert(block < (int)(Bmin.size()));
    assert(block > 0);

    result_size = Bmax[block] - Bmin[block] + 1;
    precursor_size = Bmax[block-1] - Bmin[block-1] + 1;
    int oversize = result_size;
    if ( oversize < precursor_size) oversize = precursor_size;

    TCSR<CPX> *T_iim1 = new TCSR<CPX>(oversize, oversize * oversize,
                                      fortran_index);
    TCSC<CPX,int> *T_im1i = new TCSC<CPX,int> (oversize, oversize * oversize, 
                                               fortran_index);
    CPX *gR_im1_transposed = new CPX[oversize * oversize];
    CPX *T_iim1_gR_im1 = new CPX[oversize * oversize];

    get_Tip1i_csr_Tiip1_csc(block - 1, Bmin, Bmax, T_iim1, T_im1i);
    // T_iim1_gR_im1 = T_{i,i-1} * gR_{i-1,i-1}
    transpose_full_matrix(&gR[GR_start_index[block-1]],precursor_size,
                          precursor_size, gR_im1_transposed);
    T_iim1->trans_mat_vec_mult(gR_im1_transposed, T_iim1_gR_im1,
                               precursor_size, precursor_size);
    // sigmaR = T_iim1_gR_im1 * T_im1i
    T_im1i->vec_mat_mult(T_iim1_gR_im1, sigmaR, result_size);

    delete T_iim1;
    delete T_im1i;
    delete[] gR_im1_transposed;
    delete[] T_iim1_gR_im1;

  }
  else {
    // M_TODO: throw exception
  }

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
 * \param block Index i of the block GR_{i} to calculate (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param gR The matrix containing gR_{i} (input).
 *
 * \param GR_block On start: the matrix containing GR_{i-1} (input).
 *                 On exit: the matrix containing GR_{i} (output).
 *
 * \param GRNNp1_block The matrix containing GR_{i-1,i} (output).
 *
 * \note: There's commented out code to calculate blocks of the first column of
 * GR.
 *
 */
void rGF::calculate_GR(int block, std::vector<int> Bmin, std::vector<int> Bmax, 
                       CPX *gR, CPX *GR_block, CPX *GRNNp1_block)
{

  //double start_time = get_time(0.0);

  assert(block > 0 && block < (int)Bmin.size());
  int result_size = Bmax[block] - Bmin[block] + 1;
  int precursor_size = Bmax[block-1] - Bmin[block-1] + 1;
  int oversize = result_size;
  if (oversize < precursor_size) {oversize = precursor_size;}

  // CSR->size = #rows
  // CSC->size = #cols
  TCSC<CPX,int> *T_im1i = new TCSC<CPX,int>(oversize, oversize * oversize, 
                                            fortran_index);
  TCSR<CPX> *T_iim1 = new TCSR<CPX>(oversize, oversize * oversize, fortran_index);
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


  std::stringstream filename;

  // change to fortran memory layout
  CPX *tmp = new CPX[oversize * oversize];
  transpose_full_matrix(GR_block, precursor_size, precursor_size, tmp);

  // tmp = GR_block * T_{i-1,i}, dim: precursor_size*result_size (but we over-
  // allocate). The routine produces Fortran-style matrices, i.e. the dimensions
  // are result_size*precursor_size
  // NOTE: vec_mat_mult(in, out, number_of_rows for 'in')
  T_im1i->vec_mat_mult(tmp, GR_block, precursor_size);

  // GRNNp1_block = GR_block * tmp1 = GR_{i-1} * T_{i-1,i} * gR_{i}, 
  // dim: precursor_size*result_size
  // NOTE: LDA is effectively the stride of the matrix A in BLAS (result_size
  // in our case)
  //c_zgemm('T', 'T', result_size, precursor_size, result_size, CPX(1.0, 0.0),
  c_zgemm('N', 'T', result_size, precursor_size, result_size, CPX(1.0, 0.0),
          gR, result_size, GR_block, result_size, CPX(0.0, 0.0), GRNNp1_block, 
          result_size);

  // tmp = T_{i,i-1} * GRNNp1_block = T_{i,i-1} * GR_block * T_{i-1,i} * gR_{i},
  // dim: result_size*result_size
  // NOTE: trans_mat_vec_mult(in, out, cols_of_non_transposed_in, <ignored>)
  T_iim1->trans_mat_vec_mult(GRNNp1_block, tmp, result_size, result_size);

  // GRNNp1_block *= -1, dim: precursor_size*result_size
  // M_TODO: we could also account for the factor at result construction time
  c_zscal(precursor_size * result_size, CPX(-1.0, 0.0), GRNNp1_block, 1);

  // GR_block = gR * tmp = gR_{i} * T_{i,i-1} * GR_block * T_{i-1,i} * gR_{i},
  // dim: result_size*result_size
  c_zgemm('T', 'T', result_size, result_size, result_size, CPX(1.0, 0.0),
          tmp, result_size, gR, result_size, CPX(0.0, 0.0), GR_block, 
          result_size);

  // GR_block += gR = gR_{i} + gR_{i} * T_{i,i-1} * GR_block * T_{i-1,i} * gR_{i},
  // dim: result_size*result_size
  c_zaxpy(result_size*result_size, CPX(1.0, 0.0), gR, 1, GR_block, 1);


  delete[] tmp;
  delete T_im1i;
  delete T_iim1;

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
