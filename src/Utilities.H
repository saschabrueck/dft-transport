/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __UTILITIES
#define __UTILITIES

#include <mpi.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <math.h>
#include <cmath>
#include "LinkedList.H"

#ifdef CUDA
//#include <cuda.h>
//#include <cuda_runtime_api.h>
//#include <cublas.h>
//#include <cblas.h>
#endif

bool sortx(const XYZPOS&, const XYZPOS&);

bool sorty(const XYZPOS&, const XYZPOS&);

bool sortz(const XYZPOS&, const XYZPOS&);

bool sortn(const CONNEC&, const CONNEC&);

bool sorti(const IJPOS&, const IJPOS&);

bool sortj(const IJPOS&, const IJPOS&);

bool my_isnan(double*,int);

void dreshape(int, int, double*, int*);

void icopy(int,int*,int*);

int Round(double);

void init_CSR(CSR**,int,int,int);

void del_CSR(CSR**,int);

double min_vec(double*,int,int);

double max_sign_abs_vec(double*,int,int);

void change_sign(double*,int,int,int);

void sort_vec(double*,int);

void sort_vec(double*,int*,int);

void sort_abs_imag(CPX*,int*,int);

void sort_abs_vec(CPX*,int*,int);

void check_mpi(int,int,int,int,int,int,int,int);

double get_time(double);

double randn();

void domain_decomposition(int,int,int,int,int*,int*,int,int,int,int);

int get_cpu_color(int,int,int,int);

int get_number_of_blocks(int,int,int,int*);

int get_msize(int,int,int,int*);

int get_max_orb(int*,int);

int sum_vec(int,int*,int*);

int min_active_vec(int,int*,int*);

int max_active_vec(int,int*,int*);

void cuda_initialize();

void cuda_finalize();

/************************************************************************************************/

template <typename T>
void init_var(T *var,int N);

template <typename T>
void init_var(T *var,int N)
{
    for(int i=0;i<N;i++){
        var[i] = (T)0;
    }
}

/************************************************************************************************/

template <typename T>
T sum_vec_comp(T *vec,int N,int inc);

template <typename T>
T sum_vec_comp(T *vec,int N,int inc)
{
  T sum_comp = 0;
  int i;

  for(i=0;i<N;i++){
      sum_comp = sum_comp+vec[inc*i];
  }
  return sum_comp;
}

/************************************************************************************************/

template <typename T>
T max_vec(T *vec,int N);

template <typename T>
T max_vec(T *vec,int N)
{
    T max_val = -INF;
    int i;
    
    for(i=0;i<N;i++){
        if(vec[i]>max_val){
            max_val = vec[i];
        }
    }
    return max_val;
}

/************************************************************************************************/

template <typename T>
T max_vec(T *vec,int N,int inc);

template <typename T>
T max_vec(T *vec,int N,int inc)
{
  T max_val = -INF;
  int i;

  for(i=0;i<N;i++){
    if(vec[inc*i]>max_val){
      max_val = vec[inc*i];
    }
  }
  return max_val;
}

/************************************************************************************************/

template <typename T>
T max_vec(T *vec,int *ind,int N);

template <typename T>
T max_vec(T *vec,int *ind,int N)
{
    T max_val = -INF;
    int i;
    
    for(i=0;i<N;i++){
        if(vec[ind[i]]>max_val){
            max_val = vec[ind[i]];
        }
    }
    return max_val;
}

/************************************************************************************************/

template <typename T,typename W>
T convert(W val);

template<>
inline double convert(CPX val)
{
    return real(val);
}

template<>
inline CPX convert(CPX val)
{
    return val;
}

template<>
inline double convert(double val)
{
    return val;
}

template<>
inline CPX convert(double val)
{
    return CPX(val,0.0);
}

/************************************************************************************************/

#endif

