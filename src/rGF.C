#include <vector>
#include <assert.h>
#include <string.h>

#include "rGF.H"

/** \brief Constructor
 *  
 *  Constructor for the rGF solver.
 *
 *  \param e_minus_h E (energy) times unity matrix minus the Hamiltonian.
 */
rGF::rGF(TCSR<CPX>* e_minus_h) {
    matrix = e_minus_h;
    findx  = 0; // fortran index
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
 *  \param sigmaRl Array containing the retarded boundary self-energy of the left
 *                 contact (input).
 *
 *  \param gammal Array containing the broadening function of the left contact
 *                (input).
 *
 *  \param sigmaRr Array containing the retarded boundary self-energy of the 
 *                 right contact (input).
 *
 *  \param gammar Array containing the broadening function of the right contact
 *                (input).
 *
 *  \param GR Array containing the diagonal blocks of GR (output).
 *
 *  \param GRNNp1 Array containing the first off-diagonal blocks of GR (output).
 *
 *
 *  CHECKPOINT:
 *    - consider writing a pardiso class like this rGF class
 *    - consider creating a parent class (like Mathieu)
 *
 */
void rGF::solve_blocks(std::vector<int> Bmin, 
                       std::vector<int> Bmax,
                       CPX *sigmaRl, CPX *gammal, CPX *sigmaRr, 
                       CPX *gammar, CPX *GR, CPX *GRNNp1)
{
  int num_blocks = Bmin.size();
  std::vector<int> GR_start_index(num_blocks + 1, Bmin[0]);
  std::vector<int> GRNNp1_start_index(num_blocks, Bmin[0]);

  // Determining block sizes, corresponding start and end indices and
  // biggest block's size
  int current_diagonal, next_diagonal;
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

  // rGF, STAGE 1
  CPX *gR = new CPX[GR_start_index[num_blocks]];

  // last element of the diagonal
  calculate_gR(num_blocks - 1, Bmin, Bmax, GR_start_index, sigmaRr, gR);

  // Recursion upwards along the diagonal:
  // gR_{i-1} = (E - H_{i-1} - T_{i-1,i} * gR_{i} * T_{i,i-1})^{-1}
  CPX *sigmaR = new CPX[largest_block * largest_block];
  for (int block = num_blocks - 2; block > 0; --block) {
    calculate_sigmaR(block, (char *)"r", Bmin, Bmax, GR_start_index, gR, sigmaR);
    calculate_gR(block, Bmin, Bmax, GR_start_index, sigmaR, gR);
  }

  
  // rGF STAGE 2
  //
  // calculate GR_{0} = (E-H_{0} - T_{0,1} * gR_{1} * T_{1,0} - sigmaRl)^{-1}

  // 1: calculate sigmaR[0]
  int block_diagonal = Bmax[0] - Bmin[0] + 1;
  int block_size = block_diagonal * block_diagonal;
  calculate_sigmaR(0, (char *)"r", Bmin, Bmax, GR_start_index, gR, sigmaR);

  // 2: sigmaR = 1*sigmaR_{0} + sigmaRl
  c_zaxpy(block_size, CPX(1.0, 0.0), sigmaRl, 1, sigmaR, 1);

  // 3: GR_current = gR[0]
  CPX *GR_current = new CPX[largest_block * largest_block]; // overallocation
  set_to_unity(largest_block, GR_current);
  calculate_gR(0, Bmin, Bmax, GR_start_index, sigmaR, GR_current);


  /* M_TODO: discuss if this is really not needed
  // 4: GRN1_current = 1*GRr_current + GRN1_current
  // M_TODO: since we're not reusing GRN1_current we could allocate tightly
  CPX GRN1_current = new CPX[largest_block * largest_block]; // overallocation
  set_to_zero(largest_block * largest_block, GRN1_current);
  c_zcopy(block_size, GRr_current, 1, GRN1_current, 1);
  
  
  // 5: M_tmp = GRN1_current * Gammal
  CPX M_tmp = new CPX[largest_block * largest_block];
  c_zgemm('N', 'N', block_diagonal, block_diagonal, block_diagonal,
          CPX(1.0, 0.0), GRN1_current, block_diagonal, Gammal, 
          block_diagonal, CPX(0.0, 0.0), M_tmp, block_diagonal);
  
  // 6: Zl = M_tmp * complexconjugate(GRN1_current);
  CPX Zl = new CPX[largest_block * largest_block];
  c_zgemm('N', 'C', block_diagonal, block_diagonal, block_diagonal,
          CPX(1.0, 0.0), M, block_diagonal, GRN1_current, block_diagonal,
          CPX(0.0, 0.0), Zl, block_diagonal);
  */
   
  // 7: GR_{0} = GRr_current
  c_zcopy(block_diagonal, GR_current, block_diagonal + 1, GR, 1);

  //M_TODO: calculate GR_{12}

  
  // recursion along the diagonal (downwards)
  //
  // calculate
  // GR_{i} = gR_{i} + gR_{i} * T_{i,i-1} * GR_{i-1} * T_{i-1,i} * gR_{i}
  //
  // 1 <= i <= num_blocks
  
  for (int block = 1; block < num_blocks; ++block) {
    block_diagonal = Bmax[block] - Bmin[block] + 1;
    calculate_GR(block, (char *)"r", Bmin, Bmax, gR, GR_current);
    c_zcopy(block_diagonal, GR_current, block_diagonal+1, 
            &GR[GR_start_index[block]], 1);
  }

  //MTODO: calculate GR_{ii+1}

  delete[] gR;
  delete[] sigmaR;
  delete[] GR_current;
  /*
  delete[] GRN1_current;
  delete[] M_tmp;
  delete[] Zl;*/
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
 *
 * \note Status: reviewed in detail (not yet tested).
 */
void rGF::get_diagonal_block(int block, std::vector<int> Bmin,
                                 std::vector<int> Bmax,
                                 CPX *diagonal_block)
{
  int diagonal_length = Bmax[block] - Bmin[block] + 1;
  assert(diagonal_length < Bmax[block]);

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
        diagonal_block[i_relative + (j_relative * diagonal_length)] =
                                                         matrix->nnz[j_absolute];
      }
    }
  }
}


/** \brief Extraction of two offdiagonal blocks from the main sparse matrix
 *
 * This routine extracts the two offdiagonal blocks T_{i,i+1} and T_{i+1,i}
 * from the instance's matrix (this->matrix) and stores them in the sparse
 * matrices T_iip1 (CSR) and T_ip1i (CSC) respectively. T_{i+1,i} is the 
 * conjugate of T_{i,i+1}.
 *
 * \param block Index i for extraction of T_{i,i+1} and T_{i+1,i} (numbering
 *              starts with 0) (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param T_iip1 Sparse (CSR) matrix to contain T_{i,i+1} (output).
 *
 * \param T_iip1 Sparse (CSC) matrix to contain T_{i+1,i} (output).
 *
 * \note Status: reviewed in detail (not yet tested).
 */
void rGF::get_Tiip1_csr_Tip1i_csc(int block, std::vector<int> Bmin,
                                     std::vector<int> Bmax,
                                     TCSR<CPX> *T_iip1, TCSC<CPX,int> *T_ip1i)
{
  assert(block < (int)(Bmin.size() - 1));
  assert(block >= 0);

  int current_diagonal = Bmax[block] - Bmin[block] + 1;

  T_iip1->size = current_diagonal; // CSR->size = number of rows
  T_ip1i->size = current_diagonal; // CSC->size = number of columns

  int num_nonzeros = 0;
  for (int i_absolute = Bmin[block]; i_absolute <= Bmax[block] ; 
      ++i_absolute) {

    int i_relative = i_absolute - Bmin[block];
    T_iip1->index_i[i_relative] = 0;
    T_ip1i->index_j[i_relative] = 0;

    for (int j_absolute = matrix->edge_i[i_absolute] - matrix->findx;
         j_absolute < matrix->edge_i[i_absolute+1] - matrix->findx; 
         ++j_absolute) {

      int j_relative = matrix->index_j[j_absolute] - matrix->findx -
                       Bmin[block+1];
      if (j_relative >= 0) {
        T_iip1->index_j[num_nonzeros] = j_relative + findx;
        T_ip1i->index_i[num_nonzeros] = j_relative + findx;
        T_iip1->nnz[num_nonzeros] = matrix->nnz[j_absolute];
        T_ip1i->nnz[num_nonzeros] = conj(matrix->nnz[j_absolute]);
        T_iip1->index_i[i_relative]++;
        T_ip1i->index_j[i_relative]++;
        num_nonzeros++;
      } 
    }
  }

  T_iip1->n_nonzeros = num_nonzeros;
  T_ip1i->n_nonzeros = num_nonzeros;
  T_iip1->get_row_edge();
  T_ip1i->get_column_edge();
}


/** \brief Extraction of two offdiagonal blocks from the main sparse matrix
 *
 * This routine extracts the two offdiagonal blocks T_{i,i+1} and T_{i+1,i}
 * from the instance's matrix (this->matrix) and stores them in the sparse
 * matrices T_iip1 (CSC) and T_ip1i (CSR) respectively. T_{i+1,i} is the 
 * conjugate of T_{i,i+1}.
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
 *
 * \note Status: reviewed in detail (not yet tested).
 */
void rGF::get_Tiip1_csc_Tip1i_csr(int block, std::vector<int> Bmin,
                                     std::vector<int> Bmax,
                                     TCSC<CPX,int> *T_iip1, TCSR<CPX> *T_ip1i)
{
  assert(block < (int)(Bmin.size() - 1));
  assert(block >= 0);

  // rows of T_iim1, cols of T_im1i
  int next_diagonal = Bmax[block+1] - Bmin[block+1] + 1;

  T_iip1->size = next_diagonal; // CSR->size = number of rows     
  T_ip1i->size = next_diagonal; // CSC->size = number of columns

  int num_nonzeros = 0;

  for (int i_absolute = Bmin[block+1]; i_absolute <= Bmax[block+1];
       ++i_absolute) {

    int i_relative = i_absolute - Bmin[block+1];
    T_iip1->index_i[i_relative] = 0;
    T_ip1i->index_j[i_relative] = 0;

    for (int j_absolute = matrix->edge_i[i_absolute] - matrix->findx;
         j_absolute < matrix->edge_i[i_absolute+1] - matrix->findx; 
         ++j_absolute) {

      int j_relative = matrix->index_j[j_absolute] - matrix->findx -
                       Bmin[block];
      if (j_relative < (Bmax[block] - Bmin[block] + 1)) {
        T_iip1->index_i[num_nonzeros] = j_relative+findx;
        T_ip1i->index_j[num_nonzeros] = j_relative+findx;
        T_iip1->nnz[num_nonzeros] = conj(matrix->nnz[j_absolute]);
        T_ip1i->nnz[num_nonzeros] = matrix->nnz[j_absolute];
        T_iip1->index_j[i_relative]++;
        T_ip1i->index_i[i_relative]++;
        num_nonzeros++;
      }
    }
  }

  T_iip1->n_nonzeros = num_nonzeros;
  T_ip1i->n_nonzeros = num_nonzeros;
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
 * \param sigmaR Boundary self-energy (input).
 *
 * \param gR Array containing gR (output).
 *
 */
void rGF::calculate_gR(int block, std::vector<int> Bmin,
                       std::vector<int> Bmax, 
                       std::vector<int> GR_start_index,
                       CPX *sigmaR, CPX *gR)
{

  int diagonal_length = Bmax[block] - Bmin[block] + 1;
  int block_size = diagonal_length * diagonal_length;

  CPX *inverted_gR = new CPX[block_size];
  get_diagonal_block(block, Bmin, Bmax, inverted_gR);

  // inverted_gR += -1 * sigmaR:
  c_zaxpy(block_size, CPX(-1.0, 0.0), sigmaR, 1, inverted_gR, 1);

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
  
  if (!strcmp(side, "r")) {
    // M_TODO: check range of block argument to be valid
    result_size = Bmax[block] - Bmin[block] + 1;
    precursor_size = Bmax[block+1] - Bmin[block+1] + 1;
    int oversize = result_size;
    if ( oversize < precursor_size) oversize = precursor_size;

    // M_TODO: these here are over-allocated
    TCSR<CPX> *T_iip1 = new TCSR<CPX>(oversize, oversize * oversize, findx);
    TCSC<CPX,int> *T_ip1i = new TCSC<CPX,int>(oversize , oversize * oversize, 
                                              findx);
    CPX *gR_ip1_transposed = new CPX[oversize * oversize];
    CPX *T_iip1_gR_ip1 = new CPX[oversize * oversize];

    get_Tiip1_csr_Tip1i_csc(block, Bmin, Bmax, T_iip1, T_ip1i);
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
    // M_TODO: check range of block argument to be valid
    result_size = Bmax[block] - Bmin[block] + 1;
    precursor_size = Bmax[block-1] - Bmin[block-1] + 1;
    int oversize = result_size;
    if ( oversize < precursor_size) oversize = precursor_size;

    // M_TODO: these here are over-allocated
    TCSR<CPX> *T_iim1 = new TCSR<CPX>(oversize, oversize * oversize, findx);
    TCSC<CPX,int> *T_im1i = new TCSC<CPX,int> (oversize, oversize * oversize, 
                                               findx);
    CPX *gR_im1_transposed = new CPX[oversize * oversize];
    CPX *T_iim1_gR_im1 = new CPX[oversize * oversize];

    get_Tiip1_csc_Tip1i_csr(block - 1, Bmin, Bmax, T_im1i, T_iim1);
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

}


/** \brief Function to calculate GR_{i}
 *
 * This function calculates GR_{i} as 
 *
 *   GR_{i} = gR_{i} + gR_{i} * T_{i,i+1} * GR_{i+1} * T{i,i+1} * gR_{i}
 *
 * in case of side = 'r' and
 *
 *   GR_{i} = gR_{i} + gR_{i} * T_{i,i-1} * GR_{i-1} * T{i,i-1} * gR_{i}
 *
 * \param block Index i of the block GR_{i} to calculate (input).
 *
 * \param side Indicator of the side to calculate gR from, valid values
 *             are "r" for right or "l" for left (input).
 *
 * \param Bmin Vector of block start indices along the diagonal (input).
 *
 * \param Bmax Vector of block end indices along the diagonal (input).
 *
 * \param gR The matrix containing the already calculated blocks of gR
 *           (input).
 *
 * \param GR_block The matrix containing the i'th block of GR (output).
 *
 * \note: There's commented out code to calculate blocks of the first column of
 * GR.
 *
 */
void rGF::calculate_GR(int block, char *side, std::vector<int> Bmin, 
                       std::vector<int> Bmax, CPX *gR, CPX *GR_block)
{
  int result_size;
  int precursor_size; // size of the predecessor in the recusion (size of
                      // block i+1 for side='r' and size of block i-1 for 
                      // side='l')
  int block_to_get;
  int oversize;


  if (!strcmp(side, "r")){
    assert(block > 0 && block < (int)Bmin.size());
    result_size = Bmax[block] - Bmin[block] + 1;
    precursor_size = Bmax[block+1] - Bmin[block+1] + 1;
    oversize = result_size;
    if ( oversize < precursor_size) oversize = precursor_size;
    block_to_get = block;
  }
  else if (!strcmp(side, "l")) {
    assert(block >= 0 && block < (int)(Bmin.size() - 1));
    result_size = Bmax[block] - Bmin[block] + 1;
    precursor_size = Bmax[block-1] - Bmin[block-1] + 1;
    oversize = result_size;
    if ( oversize < precursor_size) oversize = precursor_size;
    block_to_get = block - 1;
  }
  else {
    // M_TODO: throw exception
  }

  // M_TODO: these here are over-allocated
  TCSC<CPX,int> *T_iip1 = new TCSC<CPX,int>(oversize, oversize * oversize, findx);
  TCSR<CPX> *T_ip1i = new TCSR<CPX>(oversize, oversize * oversize, findx);
  get_Tiip1_csc_Tip1i_csr(block_to_get, Bmin, Bmax, T_iip1, T_ip1i);

  /* Calculation of GRN1_block
  *if (!strcomp(side, "r") && we_want_to_calc_GRN1) {
  *  int first_size = Bmax[0] - Bmin[0] + 1;
  *  transpose_full_matrix(GRN1_block, precursor_size, first_size, 
  *                        GRN1_transposed);
  *  T_ip1i->trans_mat_vec_mult(GRN1_block_transposed, M_tmp, first_size, 
  *                            precursor_size);
  *  c_zgemm('N', 'N', result_size, first_size, result_size, CPX(-1.0, 0.0),
  *          gR, result_size, M_tmp, result_size, CPX(0.0, 0.0), GRN1_block,
  *          result_size);
  *}
  */

  // GR_block = T_{i+1,i} * GR_block * T_{i,i+1}
  CPX *GR_block_transposed = new CPX[result_size * result_size];
  transpose_full_matrix(GR_block, precursor_size, precursor_size, 
                        GR_block_transposed); 
  CPX *M_tmp = new CPX[result_size * result_size];
  T_ip1i->trans_mat_vec_mult(GR_block_transposed, M_tmp, precursor_size, 
                            precursor_size);
  T_iip1->vec_mat_mult(M_tmp, GR_block, result_size);

  // M_tmp = GR_block * gR
  c_zgemm('N', 'N', result_size, result_size, result_size, CPX(1.0, 0.0),
          GR_block, result_size, gR, result_size, CPX(0.0, 0.0), M_tmp, 
          result_size);
  // GR_block = gR * M_tmp = gR * GR_block * gR
  c_zgemm('N', 'N', result_size, result_size, result_size, CPX(1.0, 0.0),
          gR, result_size, M_tmp, result_size, CPX(0.0, 0.0), GR_block,
          result_size);

  // GR_block += gR
  c_zaxpy(result_size*result_size, CPX(1.0, 0.0), gR, 1, GR_block, 1);

  delete T_iip1;
  delete T_ip1i;
  delete[] GR_block_transposed;
  delete[] M_tmp;
}


/** \brief Function to transpose (dense) matrices
 *
 * \param normal The matrix to be transposed (input).
 *
 * \param rows The number of rows of the matrix (input).
 *
 * \param cols The number of columns of the matrix (input).
 *
 * \param transposed The matrix in transposed form (output).
 */
void rGF::transpose_full_matrix(CPX *normal, int rows, int cols, CPX *transposed)
{
  for (int i=0; i < cols; ++i) {
    c_zcopy(rows, &normal[i*rows], 1, &transposed[i], cols);
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


/** \brief Function to write a double complex matrix to a file
 *
 *  \param NR Number of rows.
 *
 *  \param matrix Array containing the matrix to write to the file.
 *  
 *  \param NC Number of columns.
 *
 *  \param filename Name of the file to write to.
 */
void rGF::write_matrix(const char *filename,CPX *matrix,int NR,int NC)
{
    int IC,IR;
    ofstream myfile;
    
    myfile.open(filename);
    myfile.precision(8);
    for(IR=0;IR<NR;IR++){
        for(IC=0;IC<NC;IC++){
            myfile<<real(matrix[IR+IC*NR])<<" "<<imag(matrix[IR+IC*NR])<<" ";
        }
        myfile<<"\n";
    }
    myfile.close();
}

