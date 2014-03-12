#ifndef __WIREGENERATOR
#define __WIREGENERATOR

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <string.h>

#include "Utilities.H"
#include "Material.H"
#include "Blas.H"
//#include "AtomStrain.H"

class WireGenerator{

public:
    
    WireGenerator(int,int);
    ~WireGenerator();
    
    void execute_task(WireStructure*,Material*,int,MPI_Comm comm);
    void write_LM(const char*);
    void write_TLM(const char*);
    void write_AM(const char*);
    void write_positions(const char *);
    void get_positions(double*);
    void get_atomtype(double*);
    void get_neighbors(double*);
    void get_uc_positions(double*);
    void get_uc_atomtype(double*);
    void get_uc_neighbors(double*);
    int get_no_bond();
    int get_no_atom();	   
 
    Unit_Cell *unit_cell;
    double *Layer_Matrix,*Around_Matrix,cell_width,cell_area;
    double y_width,Ly,Lz,volume,volume_tot;
    double lenx, leny, lenz;
    double Vol_atom;
    double *EMidGap,*EGap,EMGL,EMGR;
    double *atomic_mass;
    double EVmaxL,EVmaxR,ECminL,ECminR;
    int NLayer,*Lmin,*Lmax,NSlab,*Smin,*Smax,*Neigh_Matrix;
    int *Smin_tot,*Smax_tot;
    int Channel_tot,Around_tot;
    int No_Atom,Around_Atom,Rough_Atom;
    int *neighbor_layer;
    int layer_per_slab;
    int l_per_slab_ref;
    int LB_size,RB_size;
    int **index_boundary,*index_channel,**bound_pos,*ch_pos,**bound_conv,*ch_conv;
    int *index_arch,*arch_pos,*arch_conv;
    int *b_shift,*orb_per_at;
    int QMfirst,QMlast;
    int freedim,poisson_dim;
    int NA,NB,SP,SLM;
    int NBlock;
    int lattice_type;
    int max_neighbors;
    WireStructure* nwire;
    BOUNDARY_ENN** Boundary;
    BOUNDARY_ENN** Enn;
    
protected:

    void delete_variables();
    void check_dimension(WireStructure*);
    void calc_info(WireStructure*);
    void make_unit_cell(const char*,double,double,double);
    void change_orientation(double*,double*,double*,double);
    void make_strained_unit_cell(Strain*);
    void add_strain(Strain*,int,StrainDomain**,int);
    void strain_coordinate(double*,int,int,double*,double,double,double,\
			   double,double,double,int);
    void get_strain_matrix(double*,double*,double*,double*,double*,double*,\
			   double*,StrainDomain*);
    void make_wire(WireStructure*);
    void make_bulk_layer_matrix(WireStructure*,int*);    
    void update_atom_pos(WireStructure*,Material*);
    void get_orbital_per_atom(int*);
    void get_mid_gap_energy(Material*);
    void get_bulk_orbital_per_atom(int*);
    void generate_roughness(WireStructure*,int);
    void generate_alloy_disorder(WireStructure*);
    void replicate_unit_cell(int,WireStructure*);
    void read_atom_position(int,WireStructure*);
    void get_face_info(int*,int*,int*,int,double*,int,int,WireStructure*);
    int is_in_face(double*,MAT*,int,char*,double);
    void generate_random_surface(double*,int*,int,int,int*,WireStructure*);
    void get_face_center(double*,int,int,int,int*);
    void get_acv_value(double*,double*,double*,double*,Roughness*);
    void eig(double*,double*,int);
    void separate_dimension(WireStructure*);
    void convert_position(WireStructure*);
    void bcast_wire_info(int,int,int);
    void bcast_oxide_info();
    void exchange_wire_info(int,double*,double*,int);
    void exchange_oxide_info(int,int,double*,int);
    int make_slab(WireStructure*);
    void make_connections(WireStructure*);
    void find_3D_boundary(int*,int*,int*,int*,int*,int*,WireStructure*);
    void find_2D_boundary(int*,int*,int*,int*,WireStructure*);
    void check_uc_axis(int*);
    double get_length(int);
    double calc_quad_area(double*,double*,double*,double*,int);
    double calc_tri_area(double*,double*,double*,int);
    void calc_vec(double*,double*,double*,double*);
    double calc_tri_volume(double,double*,double*,double*);
    double calc_quad_volume(MAT*);
    int is_in_quad3D(MAT*,double*);
    int is_in_quad2D(SURFACE*,double*);
    void sort_xyz(XYZPOS*,int*,int*,int*,int*);
    void init_variables(WireStructure*);
    void get_cell_width();
    void cut_boundary_slab(WireStructure*);
    void update_boundary_slab(WireStructure*);
    void cut_layer(WireStructure*);
    void boundary_dimension(int*,int*,int,int);
    void get_boundary_size(WireStructure*);
    void get_qm_region(WireStructure*);
    void roll_cnt();
    void check_periodicity(WireStructure*);
    void update_periodicity(WireStructure*);
    void hydrogen_passivation(WireStructure*);
    void adapt_parameters(Material*);
    void get_d_matrix(double*,double,double,double);

    int tb,N3D,N2D,N1D,no_element,NA_per_slice;
    int hxmin,hxmax,hymin,hymax,hzmin,hzmax; //set in make_wire
    int mpi_size,mpi_rank;
    int **inv_bound_pos,*inv_ch_pos;
    int NyFold,NzFold;
    int transport_type;
    double *XYZH, *HXYZ;
    double cnt_radius;
    MPI_Comm wg_comm;
};

#endif

