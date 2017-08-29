/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

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
  TCSR<CPX> *matrix;

  int fortran_index;

  void get_diagonal_block(int, std::vector<int>, std::vector<int>, CPX*);
  void get_Tip1i_csr_Tiip1_csc(int, std::vector<int>, std::vector<int>, 
                               TCSR<CPX>*, TCSC<CPX,int>*);
  void get_Tip1i_csc_Tiip1_csr(int, std::vector<int>, std::vector<int>,
                               TCSC<CPX,int>*, TCSR<CPX>*);
  void calculate_gR_rec(int, std::vector<int>, std::vector<int>,
                    std::vector<int>, CPX*);
  void calculate_gR_init(int, std::vector<int>, std::vector<int>,
                    std::vector<int>, CPX*);
  void calculate_sigmaR(int,  std::vector<int>, std::vector<int>, 
                        std::vector<int>, CPX*);
  void calculate_GR(int, std::vector<int>, std::vector<int>, CPX*, CPX*, 
                    CPX*);
  void transpose_full_matrix(CPX*, int, int, CPX*);
  void set_to_unity(int, CPX*);
  void test_x(TCSC<CPX,int>*, TCSR<CPX>*, int);

  // working memory
  TCSR<CPX> *sparse_CSR;
  TCSC<CPX,int> *sparse_CSC;
  CPX* tmp0;
  CPX* tmp1;
  CPX* tmp2;
};



#endif
