/*********************************************************************************************************************************
Definition of a wire structure:

1) the order of the parameters must not be changed
2) the points in the structure mat_coord(i,j) and ox_coord(i,j) must follow the order described below
3) comments must remain on one line or start again with //

Commands:

CB_Bandstructure	= conduction bandstructure of contacts
VB_Bandstructure	= valence bandstructure of contacts
CB_Transmission_UWF	= conduction band transmission with UMFPACK
VB_Transmission_UWF	= valence band transmission with UMFPACK
CB_Transmission_SWF	= conduction band transmission with SuperLU_DIST
VB_Transmission_SWF	= valence band transmission with SuperLU_DIST
CB_Transmission_MWF	= conduction band transmission with MUMPS
VB_Transmission_MWF	= valence band transmission with MUMPS
CB_Transmission_RGF	= conduction band transmission with recursive GF
VB_Transmission_RGF	= valence band transmission with recursive GF
EL_SC_UWF		= self-consistent electron simulation with UMFPACK
HO_SC_UWF		= self-consistent hole simulation with UMFPACK
EL_SC_SWF		= self-consistent electron simulation with SuperLU_DIST
HO_SC_SWF		= self-consistent hole simulation SuperLU_DIST
EL_SC_MWF		= self-consistent electron simulation with MUMPS
HO_SC_MWF		= self-consistent hole simulation MUMPS
EL_SC_RGF		= self-consistent electron simulation with recursive GF
HO_SC_RGF		= self-consistent hole simulation with recursive GF
Write_Layer_Matrix      = write file with atom positions + connections
Write_Grid_Matrix       = write file with grid points + file with index of atom position

*********************************************************************************************************************************/
/*Parameters*/

mat_name            	= mat_par;
lattice_type		= cnt;
a0                      = 0.142;                	//lattice constant
first_atom              = cation;               	//atom situated at [0 0 0]

poisson_iteration	= 10;
poisson_criterion	= 1e-3;
poisson_inner_criterion	= 5e-3;
max_proc_poisson	= 32;

read_hamiltonian	= 1;

//injection_type		= [2 90];

tb			= 10;                   	//tight-binding order
dsp3			= 30;                   	//passivation energy [eV]

Temp           		= 300;				//operation temperature

n_of_modes		= 16;              		//number of modes for bandstructure calculation
Nk			= 101;                           //number of k points in bandstructure calculation. High Nk important for good energy grid in pMOS
bs_solver		= full;

max_bond_def		= 0.1;         			//maximum relative bond deformation (should only be changed if very large strain)

last_first		= 1;                    	//last super cell is equal to first super cell if last_first=1

x                       = [0 1 0];			//transport direction
y                       = [1 0 0];			//direction of confinement
z     			= [0 0 1];              	//direction of confinement

eta_res			= 0;				//imaginary part of the energy in the reservoir (if 0 less time to get BC), eta_res<=i*1e-6
eta			= 0;			//imaginary part of the energy in the device

Elimit			= 50e-3;                	//energy interval after a mode where the energy grid is finer (should not be changed)
Emin_tail		= 4.0e-3;			//energy interval below the lowest band (should not be changed)
EOffset                 = 25*UT;			//Emax = Emin + EOffset or Emax = max(Efl,Efr) + EOffset
dE_in			= 5.0e-4;			//smallest energy interval (should not be changed)
dE_f			= 1.0e-3;			//largest energy interval (should not be changed)
dE_sep			= 1.0e-4;			//distance between a channel turn-on and the following energy point (should not be changed)
NEmax			= 1000;				//maximum number of calculated energy point

//CPU_ppoint		= 8;				//number of CPUs per energy point, must be a divider of mpi_size
//CPU_per_bc		= 8;
//spec_decomp		= 1;

strain.on		= 0;                   		//strain.on=1 => strain present
strain.Eps_xx		= 0.014;                	//strain in x direction (transport)
strain.Eps_yy		= -0.025;             		//strain in y direction
strain.Eps_zz		= -0.025;             		//strain in z direction

Eps_wire		= 3.9;                 	//Wire permitivity
Eps_ox			= [1.0 3.9];                  	//Oxide permitivity
Xi_wire			= 4.50;                 	//Wire affinity
phi_m			= 4.50;                 	//Gate work function

NVG			= 1;				//number of gate voltages Vg=Vgmin:(Vgmax-Vgmin)/(NVG-1):Vgmax
Vgmin                   = 0.0;				//absolute minimum gate potential
Vgmax                   = 0.0;				//absolute maximum gate potential

NVS			= 1;				//number of source voltages Vs=Vsmin:(Vsmax-Vsmin)/(NVS-1):Vsmax
Vsmin                   = 0.0;				//absolute minimum source potential
Vsmax                   = 0.0;				//absolute maximum source potential

NVD			= 1;				//number of drain voltages Vd=Vdmin:(Vdmax-Vdmin)/(NVD-1):Vdmax
Vdmin                   = 0.0;				//absolute minimum drain potential
Vdmax                   = 0.0;				//absolute maximum drain potential

update_fitness		= 0;
fitness_file		= FitMat_dat;

update_energy		= 0;
energy_file		= E_dat;

//restart			= [3 8 0 0];
vtot_file		= vtot_dat;
vact_file		= vact_dat;

//directory		= 

/*********************************************************************************************************************************/
/*Structure*/

grid_accuracy		= 2;                    	//number of grid points added between two neighbor atoms (2 is a good value)

no_mat			= 3;				//number of pieces that form the nanowire (channel + oxide)
no_channel_mat          = 1;                    	//number of pieces that form the nanowire channel
no_oxide_mat            = 2;                    	//number of pieces that form the oxide around the wire  

Lc			= 3.54;     			//channel length
Ls			= 2.42;  			//source length
Ld			= 2.42;				//drain length
yc			= 0.0;				//y coord. of circle center
zc			= 0.0;				//z coord. of circle center

mat_type(1)		= circle;
mat_cs(1)		= yes;
mat_radius(1)		= 0.3;
mat_coord(1,1)	        = [0.0 yc zc];
mat_coord(1,2)		= [Ls+Lc+Ld yc zc];

ox_type(1)		= circle;
ox_id(1)		= 1;   
ox_cs(1)		= no;
ox_radius(1)		= 0.35;
ox_coord(1,1)		= [0.0 yc zc];
ox_coord(1,2)		= [Ls+Lc+Ld yc zc];

ox_type(2)		= circle;
ox_id(2)		= 2;   
ox_cs(2)		= no;
ox_radius(2)		= 1.35;
ox_coord(2,1)		= [0.0 yc zc];
ox_coord(2,2)		= [Ls+Lc+Ld yc zc];

no_gate			= 1;				//number of gates that control the wire

gate_type(1)		= circle;
gate_radius(1)		= 1.35;
gate_coord(1,1)		= [Ls yc zc];
gate_coord(1,2)		= [Ls+Lc yc zc];

no_doping		= 2;				//number of doping regions in the wire

ND_S			= 1.2e27;				//donator concentration in source  [m^{-3}]
NA_S			= 1e12;				//acceptor concentration in source [m^{-3}]
ND_D			= 1.2e27;				//donator concentration in drain   [m^{-3}]
NA_D			= 1e12;				//acceptor concentration in drain  [m^{-3}]

doping_type(1)		= circle;
doping_radius(1)	= 0.4;
doping_ND(1)		= ND_S;
doping_NA(1)		= NA_S;
doping_coord(1,1)	= [0.0 yc zc];
doping_coord(1,2)	= [Ls yc zc];

doping_type(2)		= circle;
doping_radius(2)	= 0.4;
doping_ND(2)		= ND_D;
doping_NA(2)		= NA_D;
doping_coord(2,1)	= [Ls+Lc yc zc];
doping_coord(2,2)	= [Ls+Lc+Ld yc zc];

/*********************************************************************************************************************************/
/*Commands*/

command(1)            	= Write_Grid_Matrix;
command(2)		= EL_SC_MWF;
