#ifndef __UTILITIES
#define __UTILITIES

using namespace std;

#include <mpi.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctime>
#include <math.h>
#include <cmath>
#include "LinkedList.H"

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

void domain_decomposition(int,int,int,int,int*,int*,int);

int get_number_of_blocks(int,int,int,int*);

int get_msize(int,int,int,int*);

int get_max_orb(int*,int);

int sum_vec(int,int*,int*);

int min_active_vec(int,int*,int*);

int max_active_vec(int,int*,int*);

template <typename T>
void set_to_zero(int length, T *array) {
  for (int i=0; i < length; ++i) {
    array[i] = (T)0;
  }
}

void set_random(int, int, CPX*);

void write_matrix_to_file(const char*, CPX*, int, int);

int get_csr_rows(const char*);

void read_csr(const char*, int, int, CPX*);

void print_mat(TCSR<CPX>*, int, int);

void print_mat(CPX*, int, int, int);

void print_mat(CPX*, int, int, int, int, int);
 
void print_mat_spy(CPX *, int , int , int );

void print_mat_spy(TCSR<CPX> *, int , int );

void write_mat_c(CPX*, int, int, const char*);

void write_mat_c(CPX*, int, int, std::basic_string<char>);

void write_mat_f(CPX*, int, int, const char*);

void write_mat_f(CPX*, int, int, std::basic_string<char>);

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

#ifdef _AIX

template<>
inline double convert<double,CPX>(CPX val)
{
    return real(val);
}

template<>
inline CPX convert<CPX,CPX>(CPX val)
{
    return val;
}

template<>
inline double convert<double,double>(double val)
{
    return val;
}

template<>
inline CPX convert<CPX,double>(double val)
{
    return CPX(val,0.0);
}

#else

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

#endif
/************************************************************************************************/

template <typename T>
T fermi(T E,double mu,double temp,int mode)
{
    double k_b=K_BOLTZMANN;
    if (mode==0) {
        if (temp<=0.0) {
            if (real(E)<=mu) {
                return 1.0;
            } else {
                return 0.0;
            }
        } else {
            return 1.0/(1.0+exp((E-mu)/(k_b*temp)));
        }
    } else if (mode==1) {
        if (temp<=0.0) {
            return 0.0;
        } else {
            return log(1.0+exp(-(E-mu)/(k_b*temp)));
        }
    } else if (mode==2) {
        if (temp<=0.0 || abs(E-mu)>10.0) { // CHECK WHY -nan FOR >10.0
            return 0.0;
        } else {
            return exp((E-mu)/(k_b*temp))/(k_b*temp)/pow(1.0+exp((E-mu)/(k_b*temp)),2);
        }
    } else if (mode==3) {
        if (temp<=0.0) {
            return 0.0;
        } else {
            return exp((E-mu)/(k_b*temp))/(k_b*temp)/(1.0+exp(-(E-mu)/(k_b*temp)));
        }
    } else {
        return 0.0;
    }
}

/************************************************************************************************/

#endif

