/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __FEMGRID
#define __FEMGRID

#include <mpi.h>
#include <string.h>
#include <algorithm>
#include "WireGenerator.H"
#include "Utilities.H"
#include "Blas.H"
#include "Types.H"

class FEMGrid{
public:

    double *grid,*tetraVol,*doping,*Eps;
    int *tetrahedron,*atom_index,*real_at_index,*ratom_index;
    int *channel_index,*poisson_index,**gate_index,*ground_index;
    int *source_index,*drain_index;    
    int NAtom,NGrid,NTetra,NChannel,NPoisson,*NGate,NGround;
    int no_diff_gate;
    int NSource,NDrain;    
    
    FEMGrid();
    ~FEMGrid();
    
    void execute_task(WireGenerator*,WireStructure*,int,int,MPI_Comm,MPI_Comm);
    void execute_reduced_task(WireGenerator*,WireStructure*,int,MPI_Comm,MPI_Comm);
    void write_grid(const char*); 
    void write_tetrahedron(const char*);
    void write_gate(const char*);
    void delete_tetrahedron();
    
private:

    int NB,N3D,mpi_size,mpi_rank;
    int del_all_var;

    void make_grid(WireGenerator*,int);
    void make_1D_grid(WireGenerator*);
    void sort_grid_point(WireGenerator*,double*,int,int,int,double);
    void sort_xyz(int,XYZPOS*,int*,int*,int*,int*);
    
    void write_permittivity(const char*);
    void delaunay(int,char*,int,int,int,MPI_Comm,MPI_Comm);
    double calc_tri_area(double*,double*,double*,int);
    void calc_vec(double*,double*,double*,double*);
    double calc_tri_volume(double*,double*,double*,double*);
    double calc_tri_volume(double,double*,double*,double*);
    double calc_max_distance(double*,double*,double*,double*);
    void check_dimension(WireGenerator*,WireStructure*);
    void find_contact(WireStructure*);
    void find_doping(WireStructure*,double,int);
    void find_1D_doping(WireStructure*);
    void add_doping(WireStructure*,double,int);
    void find_permittivity(WireStructure*);
    void find_1D_permittivity(WireStructure*);
    void reduce_atom_index();
    int is_in_gate(GATE*,double*,double);
    void calc_doping_info(WireStructure*);
    double calc_quad_area(double*,double*,double*,double*,int);
    double calc_quad_volume(DOPING*);
    int is_in_quad3D(DOPING*,double*);
    int is_in_quad3D(MAT*,double*);
    int is_in_quad2D(GATE*,double*);
    double check_volume(MAT*,double*);
    
};

#endif
