#ifndef __TYPES
#define __TYPES

using namespace std;

#include <math.h>
#include <complex>
#include <algorithm>

#ifndef PI
#define PI 3.14159265358979
#endif

// k in [eV/K]
#ifndef K_BOLTZMANN
#define K_BOLTZMANN 8.617332478e-5
#endif

#ifndef E_ELECTRON
#define E_ELECTRON 1.6022e-19
#endif

#ifndef H_BAR
#define H_BAR 1.0546e-34
#endif

#ifndef ROWBLOCK
#define ROWBLOCK 100
#endif

#ifndef COLBLOCK
#define COLBLOCK 100
#endif

#ifndef MAX_COMM
#define MAX_COMM 20
#endif

#ifndef INF
#define INF 1E8
#endif

#ifndef SEP_TASK
#define SEP_TASK 50
#endif

#ifndef tollim
#define tollim 1E-10
#endif

#ifndef scale_norm
#define scale_norm 1.0
#endif

#ifndef max_diff_gate
#define max_diff_gate 4
#endif

#ifdef Add_
#define fortran_name(x,y) (x ## _)
#endif

#ifdef NoChange
#define fortran_name(x,y) (x)
#endif

#ifdef UpCase
#define fortran_name(x,y) (y)
#endif

typedef complex<double> CPLXType;
typedef CPLXType CPX;
typedef CPX *CPXp;
typedef CPX (*CPXpfn)(CPX);

struct POINT3D
{
    double coord[3];
};

struct POINT2D
{
    double coord[2];
};

struct Axis
{
    double vec[3];
};

struct Bond
{
    double vec[3];
    double SP3[16];
    double d0;
    int neigh;
};

struct Atom
{
    int type;
    double coord[3];
    double int_disp[3];
};

struct at_type
{
    Bond** bond;
};

struct Unit_Cell
{
    at_type** type;
    at_type** type_original;
    Atom** atom;
    Axis** axis;
    POINT3D** BZPoint;    
    double strain_matrix[9];
    double bond_pass;
    double passivation[4];
    int BZNP;    
};

struct Strain
{
    int on;
    double Eps_xx,Eps_yy,Eps_zz;
    double Eps_xy,Eps_xz,Eps_yz;
    double zeta;
};

struct StrainDomain
{
    double xmin,xmax;
    double ymin,ymax;
    double zmin,zmax;
    double Eps_vec[6];
};

struct Contact
{
    int active;
    int type[2];    
    double barrier[2];
    double virtual_CB[2];    //energy of metal conduction band edge below the Fermi level    
};

struct Roughness
{
    int on,seed;
    char *type;
    double rms,Lms;
};

struct Alloy
{ 
    int on,seed;
    double composition;
};

struct MAT
{
    POINT3D** p;
    double radius[2];
    double face_area[6];
    double p_ref[6][3];
    double vec_dir[6][3];
    double volume;
    int id_number;
    char *type;
    char *cross_section;
};

struct SURFACE
{
    POINT2D** p;
    double radius[2];
    double face_area;
    int id_number;
    char *type;
    char *cross_section;
};

struct DOPING
{
    POINT3D** p;
    double radius[2];
    double face_area[6];
    double p_ref[6][3];
    double vec_dir[6][3];
    double volume;
    char *type;
    double NA,ND;
};

struct DopingDomain
{
    double xmin,xmax;
    double ymin,ymax;
    double zmin,zmax;
    double slope[3];
    double conc;
};

struct GATE
{
    POINT3D** p;
    double radius[2];
    double angle[2];
    double face_area;
    char *type;
};

struct WireStructure
{
    char *first_atom;
    char *at_file;
    char *dop_file;
    char *perm_file;
    char *vtot_file;
    char *vact_file;
    char *phiy_file;
    char *phiz_file;
    char *tetra_file;
    char *fit_file;
    char *grid_file;
    char *energy_file;
    char *ph_energy_file;
    char *ph_mode_file;
    int restart,IG_start,IS_start,ID_start;
    int NDim;
    int NxFold,NyFold,NzFold;
    int no_element;
    int no_ch_element;
    int no_ox_element;
    int no_ro_element;
    int no_strain_domain;
    int no_doping_domain;
    double a0,c0,u0;
    double x[3],y[3],z[3];
    double max_bond_deformation;
    double bs_target;
    MAT **mat;
    SURFACE **surf;
    Strain *strain;
    StrainDomain **strain_domain;
    Roughness *rough;
    Alloy *alloy;
    Contact *Schottky;    
    int QMregion;
    int update_bs_target;
    int update_fermi;
    int update_at;
    int update_dop;
    int update_perm;
    int update_tetra;
    int update_fitness;
    int update_grid;
    int update_energy;
    int no_ph_energy;
    int no_ph_k;
    int ph_file_type;
    int sc_max_iteration;
    int sc_restart;
    int sc_step;
    double sc_id_crit;
    double sc_rh_crit;
    int sc_dist_dep;
    int sc_diag_def;
    int sc_k_coupling;
    int sc_e_ph_coupling;
    double sc_scale_fact;
    double sc_vbound;
    double sc_memory_fact;
    int incoherent_injection;
    int full_current;
    int robust_numerics;
    int convergence_booster;
    double Ef0;
    double QMstart,QMstop;
    double dsp3;
    int hydrogen_passivation;
    int grid_accuracy;
    double Eps_wire[4];
    double Eps_ox[4];
    double Xi_wire;
    double phi_m;
    int no_gate;
    int no_doping;
    int no_ground;
    int no_diff_gate;
    GATE **gate;
    DOPING **doping;
    DopingDomain **doping_domain;
    GATE **ground;
    double Ls,Lc,Ld;
    double ND_S,NA_S;
    double ND_D,NA_D;
    double VFF_alpha;
    double VFF_beta;
    int relax_atoms;
    int strain_bc;
    int charge_average;
    int charge_transfer;
    int replicate_unit_cell;
    int read_atom_pos;
    int open_system;
    int periodic_system;
    int read_hamiltonian;
    int bulk_mat_id;
};

struct XYZPOS
{
    double x;
    double y;
    double z;
    int index;
};

struct IJPOS
{
    int i;
    int j;
    int index;
};

struct BOUNDARY_ENN
{
    double *UNN, *UNNp, *UNNm;
    double cell_width, cell_area;
    int *UNNneigh, NA, NC, NA_tot;
};

struct CONNEC
{
    int neigh;
    int index;
};

struct ENERGY
{
    double Elimit;        //interval after a mode to go from dE_in to dE_f
    double Emin_tail;     //energy offset below first mode
    double EOffset;       //energy offset above first mode
    double dE_in;         //energy interval after a mode
    double dE_in_tot[2];  //[0] for electrons [1] for phonons
    double dE_f;          //energy interval after mode energy + Elimit
    double dE_f_tot[2];   //[0] for electrons [1] for phonons
    double dE_sep;        //minimum distance between mode and discretization point
    double dE_sep_tot[2]; //[0] for electrons [1] for phonons
    double EExt;          //energy added after Emax (only for phonons)
    int NEmax;            //maximum number of energy points
    int regular_mesh;
};

struct PARAM
{
    char *mat_name;
    char *table_file;
    int strain_model;
    double mat_binary_x[3];
    int tb;
    int last_first;
    int n_of_modes;
    int Nk,Nky,Nkz;
    int rot_sym;
    int bs_solver;
    double eta_res;
    double eta;
    int CPU_per_sample;
    int CPU_per_temp_point;
    int CPU_per_vd_point;
    int CPU_per_vg_point;
    int CPU_per_kz_point;
    int CPU_ppoint;
    int CPU_per_wire;
    int CPU_per_bc;
    int NPROW;
    int NPCOL;
    int NPCS;
    int spec_decomp;
    char **command;
    int no_comm;
    double Temp;
    int poisson_solver;
    double poisson_criterion;
    int poisson_iteration;
    double poisson_inner_criterion;
    int poisson_inner_iteration;
    int max_proc_poisson;
    char *directory;
    int plot_all_k;
    int lattype;
    int transport_type;
    int injection_type;
    int NM_for_feast;
    double electron_weight;
};

struct FermiStructure
{
    double *Ekl,*Ekr;
    double dkx,*dky,*dkz;
    int Nkx,Nky,Nkz;
    int n_of_modes;
    double Temp;
    double spin_factor;
    double cell_area;
    int derivative;
};

struct VOLTAGE
{
    int NVG;
    double Vgmin[max_diff_gate];
    double Vgmax[max_diff_gate];
    int NVS;
    double Vsmin;
    double Vsmax;
    int NVD;
    double Vdmin;
    double Vdmax;
    int NTEMP;
    double Tmin;
    double Tmax; 	
};

struct transport_parameters 
{
    int    n_occ;
    int    n_atoms;
    int    method;
    int    injection_method;
    int    linear_solver;
    int    n_abscissae;
    int    n_kpoint;
    int    num_interval;
    int    num_contacts;
    int    tasks_per_point;
    int    extra_scf;
    double colzero_threshold;
    double eps_limit;
    double eps_decay;
    double eps_singularity_curvatures;
    double eps_mu;
    double eps_eigval_degen;
    double energy_interval;
    double min_interval;
    double temperature;
};

struct contact_type
{
    int    bandwidth, ndof, start, inj_sign;
    double n_ele;
};

namespace transport_methods {
enum transport_method {
  EQ=0,
  NEGF=1,
  GF=2,
  WF=3,
};
} // namespace

#endif





