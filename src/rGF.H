#ifndef CP2K_OMENDFT_TRANSPORT_RGF_H_
#define CP2K_OMENDFT_TRANSPORT_RGF_H_

#include "CSR.H"
#include "CSC.H"
#include "Types.H"
#include "Utilities.H"

/**
 * M_TODO: write documentation
 */
class rGF {
	
public:

  rGF(TCSR<CPX>* mat);
  ~rGF();
  void solve_blocks(std::vector<int>, std::vector<int>, CPX*, CPX*);

private:

  /// E * 1 - H
  TCSR<CPX>* matrix;

  int fortran_index;

  void get_diagonal_block(int, std::vector<int>, std::vector<int>, CPX*);
  void get_Tip1i_csr_Tiip1_csc(int, std::vector<int>, std::vector<int>, 
                               TCSR<CPX>*, TCSC<CPX,int>*);
  void get_Tip1i_csc_Tiip1_csr(int, std::vector<int>, std::vector<int>,
                               TCSC<CPX,int>*, TCSR<CPX>*);
  void calculate_gR(int, std::vector<int>, std::vector<int>,
                    std::vector<int>, CPX*, CPX*);
  void calculate_gR(int, std::vector<int>, std::vector<int>,
                    std::vector<int>, CPX*);
  void calculate_sigmaR(int, char*, std::vector<int>,
                        std::vector<int>, std::vector<int>, CPX*, CPX*);
  void calculate_GR(int, std::vector<int>, std::vector<int>, CPX*, CPX*,
                    CPX*);
  void transpose_full_matrix(CPX*, int, int, CPX*);
  void set_to_unity(int, CPX*);
  void test_x(TCSC<CPX,int>*, TCSR<CPX>*, int);
};



#endif
