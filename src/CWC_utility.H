/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <cuda_runtime.h>
#include <stdio.h>
#include "Types.H"

extern "C" {

  void set_gpu(int,char*);
  void cublas_init(void**);
  void cublas_finalize(void*);
  void cusparse_init(void**);
  void cusparse_finalize(void*);
  size_t allocate_data_on_device(void**,size_t);
  void deallocate_data_on_device(void*);
  size_t deallocate_data_on_dev(void*,size_t);
  void copy_data_to_device(void*,void*,int,int,size_t);
  void memcpy_to_device(void*,void*,size_t);
  void copy_data_to_host(void*,void*,int,int,size_t);
  void memcpy_to_host(void*,void*,size_t);
  void dgemm_on_dev(void*,char,char,int,int,int,double,double*,\
		    int,double*,int,double,double*,int);
  void zgemm_on_dev(void*,char,char,int,int,int,CPX,CPX*,\
		    int,CPX*,int,CPX,CPX*,int);
  void zaxpy_on_dev(void*,int,CPX,CPX*,int,CPX*,int);
  void d_copy_csr_to_device(int,int,int*,int*,double*,int*,int*,double*);
  void z_copy_csr_to_device(int,int,int*,int*,CPX*,int*,int*,CPX*);
  void d_csr_mult_f(void*,int,int,int,int,int*,int*,double*,double,\
		    double*,double,double*);
  void z_csr_mult_f(void*,int,int,int,int,int*,int*,CPX*,CPX,\
		    CPX*,CPX,CPX*);
  void z_csr_mult_fo(void*,int,int,int,int,int*,int*,CPX*,CPX,\
		     CPX*,CPX,CPX*);
  void d_transpose_matrix(double*,double*,int,int);
  void z_transpose_matrix(CPX*,CPX*,int,int);
  void d_init_var_on_dev(double*,int,cudaStream_t);
  void d_init_eye_on_dev(double*,int,cudaStream_t);
  void z_init_var_on_dev(CPX*,int,cudaStream_t);
  void z_init_eye_on_dev(CPX*,int,cudaStream_t);
  void correct_diag_on_dev(CPX*,int,cudaStream_t);
  void change_var_type_on_dev(double*,CPX*,int,cudaStream_t);
  void change_sign_imag_on_dev(CPX*,int);
  void d_extract_diag_on_dev(double*,int*,int*,CPX*,int,int,int,int,int,\
			     cudaStream_t);
  void d_extract_not_diag_on_dev(double*,int*,int*,CPX*,int,int,int,int,\
				 int,int,int,cudaStream_t); 
  void d_symmetrize_matrix(double*,int,cudaStream_t);
  void z_extract_diag_on_dev(CPX*,int*,int*,CPX*,int,int,int,int,int,\
			     cudaStream_t);
  void z_extract_not_diag_on_dev(CPX*,int*,int*,CPX*,int,int,int,int,\
				 int,int,int,cudaStream_t); 
  void z_symmetrize_matrix(CPX*,int,cudaStream_t);
  void z_symmetrize_matrix_2(CPX*,int,cudaStream_t);
}
