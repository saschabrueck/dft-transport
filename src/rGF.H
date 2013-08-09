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
  void solve_blocks(std::vector<int>, std::vector<int>, CPX*, CPX*, CPX*,
                    CPX*, CPX*, CPX*);
				
private:

  /// E * 1 - H
  TCSR<CPX>* matrix;

  int findx;

  void get_diagonal_block(int, std::vector<int>, std::vector<int>, CPX*);
  void get_Tiip1_csr_Tip1i_csc(int, std::vector<int>, std::vector<int>, 
                               TCSR<CPX>*, TCSC<CPX,int>*);
  void get_Tiip1_csc_Tip1i_csr(int, std::vector<int>, std::vector<int>,
                               TCSC<CPX,int>*, TCSR<CPX>*);
  void calculate_gR(int, std::vector<int>, std::vector<int>,
                    std::vector<int>, CPX*, CPX*);
  void calculate_sigmaR(int, char*, std::vector<int>,
                        std::vector<int>, std::vector<int>, CPX*, CPX*);
  void calculate_GR(int, char*, std::vector<int>, std::vector<int>,
                    CPX*, CPX*);
  void transpose_full_matrix(CPX*, int, int, CPX*);
  void set_to_unity(int, CPX*);
  void write_matrix(const char*, CPX*, int, int);

  template <typename T>
  void set_to_zero(int, T*);
};


/** 
 * M_TODO: write documentation
 */
template <typename T>
void rGF::set_to_zero(int length, T *array) {
  for (int i=0; i < length; ++i) {
    array[i] = (T)0;
  }
}


#endif
