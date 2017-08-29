/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __POISSON
#define __POISSON

#include <mpi.h>
#include "WireGenerator.H"
#include "Utilities.H"
#include "Blas.H"
#include "Types.H"
#include "FEMGrid.H"
#include "CSR.H"
#include "Hypre.H"
#include "Fermi.H"

class Poisson{
public:

    TCSR<double> *P,*P1,*P2,*RP2;
    
    Poisson();
    ~Poisson();
    
    void init(WireGenerator*,WireStructure*,FEMGrid*,int,int,MPI_Comm);
    void solve(double*,double*,double*,double*,double,FermiStructure*,FEMGrid*,\
	       WireGenerator*,double,double*,double,double,double,double*,double,\
               int,int,int,MPI_Comm,MPI_Comm,int,MPI_Comm,int);
    
private:
    
    int N3D,mpi_rank,mpi_size,findx;
    int del_PMatrix,del_P1DMatrix;
    double Eps0,e,kB,Vol;
      
    void solve_FDM(double*,double*,double*,double*,double,FermiStructure*,FEMGrid*,\
		   WireGenerator*,double,double*,double,double*,double,int,int,\
		   MPI_Comm,MPI_Comm);
    void solve_FEM(double*,double*,double*,double*,double,FEMGrid*,WireGenerator*,\
		   double,double*,double,double,double,double*,double,int,int,\
                   MPI_Comm,MPI_Comm,int,MPI_Comm,int);
    void geometry(double*,double*,double*,double*,double*,double*);
    void create_FEM_matrix(int*,int,int*,double*,int,double*,int,double*);
    void create_FDM_matrix(FEMGrid*);
    void get_points(double*,int*,double*,int);
    double get_permitivity(double*,int*);
    void sort_ij(IJPOS*,int*,int*,int*,int);
    void assemble_matrix(int,double*,double*,IJPOS*,int*,int,int);
    void reduce_P(int*,int);
    void init_j_shift(int*,int*,int);
    void init_V(double*,double*,int,int**,int*,int,double*,int*,int,double,int*,int,\
		double,int*,int,double);
    void iterate_V(double*,double*,int*,int);
    void find_fermi(double*,double*,double*,double*,double,int*,int,double);
    void find_1D_fermi(double*,double*,double*,double*,double*,Fermi*,FermiStructure*,\
		       double*,double*,int,double,double*,int,int,double);
    void init_doping(double*,double*,int*,int,int,int);
    void calc_cd_cdd(double*,double*,double*,double*,double*,double,int*,int*,\
                     int*,int*,int*,int,double,int,int);
    void copy_cd_cdd(double*,double*,double*,double*,int*,int*,int*,int*,int,\
		     int,int);
    void calc_1D_cd_cdd(double*,double*,double*,double*,double*,Fermi*,\
			FermiStructure*,double,double,int,int);
    void generate_output(const char*,double*,int,int);
    void generate_output(const char*,int*,int,int);
};

#endif
