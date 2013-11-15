#ifndef __FEMGRID
#define __FEMGRID

#include <mpi.h>
#include <algorithm>
#include "WireGenerator.H"
#include "Utilities.H"
#include "Blas.H"
#include "Types.H"

class FEMGrid{
public:

    double *grid,*tetraVol,*doping,*Eps;
    int *tetrahedron,*atom_index,*real_at_index,*ratom_index;
    int *channel_index,*poisson_index,*gate_index,*ground_index;
    int *source_index,*drain_index;    
    int NAtom,NGrid,NTetra,NChannel,NPoisson,NGate,NGround;
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
