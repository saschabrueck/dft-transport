%{
#include <stdio.h>
#include <stdlib.h>    
#include "InputParameter.H"

extern int yylineno;
extern char *yytext;

int yyerror (const char *s) 
{
  fprintf(stderr, "GreenSolver: parse error at line %d  (token = '%s')\n",
          yylineno, yytext );

  return 1;
}
double kB = 1.38e-23;
double e  = 1.6022e-19;
double Temp;
double Ls,Lc,Ld,yc,zc;
double ts,tc,td,hs,hc,hd,tox,hox,tground;
double Lqw,Lb,Lb1,Lb2,Lsp,Lsp1,Lsp2; 
double x0,y0,z0;
double ND_S,ND_D,NA_S,NA_D;
int ncm,nom;
%}

%union {
  char* c;
  int i;
  double d;
}

%token <i>INTEGER
%token <d>DOUBLE
%token <c>CHAR

%token SEMICOLON
%token NUMBER VAR
%token X Y Z
%token PLOTALLK
%token TRANSTYPE
%token PSOLVER PCRITERION PITERATION PMAXPROC PINCRITERION PINITERATION CHAVG CHTRS
%token NDIM NKY NKZ NXFOLD NYFOLD NZFOLD PHIYFILE PHIZFILE ROTSYM BSSOLVER
%token EQUAL PLUS MINUS MULT LBRACE RBRACE LPAR RPAR COMMA POINT
%token DSP3 HPASS TB NM NK LF NAME BINARYX ETAR ETA MBD LATTICE FIRSTATOM
%token NPCS SPDEC NPR NPC CPU CPUKZ CPUSAMP CPUVG CPUVD CPUTEMP CPUWI CPUBC ALL
%token ALPHA BETA RELAX_ATOMS STRAIN_BC
%token SCHOTTKY ACTIVE TYPE BARRIER VIRTUALCB
%token ELIMIT EMIN EOFFSET DEIN DEF DESEP EEXT NE REMESH
%token UPDATEEN ENERGYFILE
%token STRAINMOD STRAIN ON EXX EYY EZZ EXY EXZ EYZ ZETA
%token NOSTRAIND SDOM XMIN XMAX YMIN YMAX ZMIN ZMAX EPSVEC
%token NODOPD DDOM SLOPE DCONC
%token ROUGHNESS SEED FORM RMS LMS
%token REPUC READATPOS UPDATEAT ATOMFILE UPDATEDOP DOPFILE UPDATEPERM PERMFILE
%token UPDATEFERM FERMIL
%token PHENFILE NOPHEN PHMOFILE PHFILETY NOPHK 
%token SCMAXIT SCRESTART SCIDCRIT SCRHCRIT SCDISTDEP SCDIAGDEF SCKCOUP SCSCALE SCVBOUND INCINJ
%token QMSTART QMSTOP
%token GRID UPDATETETRA TETRAFILE UPDATEFIT FITFILE UPDATEGRID GRIDFILE
%token NOMAT NOCMAT NOOMAT LS LC LD YC ZC MCOORD OCOORD MTYPE OTYPE MID OID MCS OCS MRADIUS ORADIUS
%token NORMAT RCOORD RTYPE RCS RRADIUS
%token TS TC TD HS HC HD TOX TGR HOX X0 Y0 Z0
%token LQW LB LB1 LB2 LSP LSP1 LSP2
%token EPSWIRE EPSOX XIWIRE PHIM
%token NDS NAS NDD NAD
%token NOGATE NODOPING GCOORD DCOORD GTYPE DTYPE DND DNA DRADIUS GRADIUS GANGLE
%token NOGROUND GRCOORD GRTYPE GRRADIUS
%token UT TEMP
%token NVG VGMIN VGMAX NVS VSMIN VSMAX NVD VDMIN VDMAX NTEMP TMIN TMAX
%token DIRECTORY COMMAND
%token RESTART VTFILE VAFILE
%token FULLCURRENT 
%token ROBUSTNUM CBOOST
%token LATTYPE

%type <d>double
%type <d>expr
%type <d>exprdef

%%

input		: command separator
		| input command separator
		;

command		: inittbpar
                | initplotopt
                | initbspar
                | initcpu
                | initpoisson
                | initxyz
                | initenergy
                | initstrain
                | initstdomain
                | initroughness
                | initqmregion
                | initupdateat
                | initupdatedop
                | initupdateperm
                | initupdatefermi
                | initscdata
                | initnomat
                | initmtype
                | initmcs
                | initmat
                | initradius
                | initdimvar
                | initgate
                | initground
                | initdoping
                | initdopdomain
                | initgrid
                | initvolt
                | initrestart
                | initdirectory
                | initcommand
                | initrelax
                | initphonon
                | initcurrent
                | initschottky
                | initnumerics
		;

inittbpar       : DSP3 EQUAL double {init_dsp3($3);}
                | HPASS EQUAL INTEGER {init_hpass($3);}
                | TRANSTYPE EQUAL INTEGER {init_transport_type($3);}
                | TB EQUAL INTEGER {init_tb($3);}
                | MBD EQUAL double {init_mbd($3);}
                | FIRSTATOM EQUAL CHAR {init_first_atom($3);}
                | LATTICE EQUAL double {init_a0_c0_u0($3,$3,0.0);}
                | LATTICE EQUAL LBRACE double double RBRACE {init_a0_c0_u0($4,$5,0.0);}
                | LATTICE EQUAL LBRACE double double double RBRACE {init_a0_c0_u0($4,$5,$6);}
                | LF EQUAL INTEGER {init_last_first($3);}
                | NDIM EQUAL INTEGER {init_ndim($3);}
                | NAME EQUAL CHAR {init_mat_name($3);} 
                | BINARYX EQUAL DOUBLE {init_mat_binary_x($3,0.0,0.0);}
                | BINARYX EQUAL LBRACE DOUBLE RBRACE {init_mat_binary_x($4,0.0,0.0);}
                | BINARYX EQUAL LBRACE DOUBLE DOUBLE RBRACE {init_mat_binary_x($4,$5,0.0);}
                | BINARYX EQUAL LBRACE DOUBLE DOUBLE DOUBLE RBRACE {init_mat_binary_x($4,$5,$6);}
                | LATTYPE EQUAL CHAR {init_lattype($3);}
		;

initplotopt     : PLOTALLK EQUAL INTEGER {init_plot_all_k($3);}
                ;

initbspar       : NM EQUAL INTEGER {init_n_of_modes($3);}
	        | NK EQUAL INTEGER {init_nk($3);}
		| NKY EQUAL INTEGER {init_nky($3);}
		| NKZ EQUAL INTEGER {init_nkz($3);}
                | PHIYFILE EQUAL CHAR {init_phiy_file($3);}
                | PHIZFILE EQUAL CHAR {init_phiz_file($3);}
                | ROTSYM EQUAL INTEGER {init_rot_sym($3);}
		| NXFOLD EQUAL INTEGER {init_nxfold($3);}
		| NYFOLD EQUAL INTEGER {init_nyfold($3);}
		| NZFOLD EQUAL INTEGER {init_nzfold($3);}
                | BSSOLVER EQUAL CHAR {init_bs_solver($3);}
		;

initcpu         : NPR EQUAL INTEGER {init_NPROW($3);}
                | NPC EQUAL INTEGER {init_NPCOL($3);}
                | NPCS EQUAL INTEGER {init_NPCS($3);}
                | SPDEC EQUAL INTEGER {init_spec_decomp($3);}
                | CPU EQUAL INTEGER {init_CPU_ppoint($3);}
                | CPUKZ EQUAL INTEGER {init_CPU_per_kz_point($3);}
                | CPUKZ EQUAL ALL {init_CPU_per_kz_point(-1);}
                | CPUSAMP EQUAL INTEGER {init_CPU_per_sample($3);}
                | CPUSAMP EQUAL ALL {init_CPU_per_sample(-1);}
                | CPUVG EQUAL INTEGER {init_CPU_per_vg_point($3);}
                | CPUVG EQUAL ALL {init_CPU_per_vg_point(-1);}
                | CPUVD EQUAL INTEGER {init_CPU_per_vd_point($3);}
                | CPUVD EQUAL ALL {init_CPU_per_vd_point(-1);}
                | CPUTEMP EQUAL INTEGER {init_CPU_per_temp_point($3);}
                | CPUTEMP EQUAL ALL {init_CPU_per_temp_point(-1);}
                | CPUWI EQUAL INTEGER {init_CPU_per_wire($3);}
                | CPUWI EQUAL ALL {init_CPU_per_wire(-1);}
                | CPUBC EQUAL INTEGER {init_CPU_per_bc($3);}
		;

initpoisson     : PSOLVER EQUAL INTEGER {init_poisson_solver($3);}
                | PCRITERION EQUAL double {init_poisson_criterion($3);}
                | PITERATION EQUAL INTEGER {init_poisson_iteration($3);}
                | PMAXPROC EQUAL INTEGER {init_max_proc_poisson($3);}
                | PINCRITERION EQUAL double {init_poisson_inner_criterion($3);}
                | PINITERATION EQUAL INTEGER {init_poisson_inner_iteration($3);}
                | CHAVG EQUAL INTEGER {init_charge_average($3);}
                | CHTRS EQUAL INTEGER {init_charge_transfer($3);}
                ;

initxyz         : X EQUAL LBRACE double double double RBRACE {init_x($4,$5,$6);}
                | Y EQUAL LBRACE double double double RBRACE {init_y($4,$5,$6);}
                | Z EQUAL LBRACE double double double RBRACE {init_z($4,$5,$6);}
		;

initenergy      : ELIMIT EQUAL double {init_Elimit($3);}
                | EMIN EQUAL double {init_Emin_tail($3);}
                | EOFFSET EQUAL expr {init_EOffset($3);}
                | DEIN EQUAL double {init_dE_in($3);}
                | DEF EQUAL double {init_dE_f($3);}
                | DESEP EQUAL double {init_dE_sep($3);}
                | EEXT EQUAL expr {init_EExt($3);}
                | NE EQUAL INTEGER {init_NE($3);}
                | REMESH EQUAL INTEGER {init_reg_mesh($3);}
                | ETAR EQUAL double {init_eta_res($3);}
                | ETA EQUAL double {init_eta($3);}
                | UPDATEEN EQUAL INTEGER {init_update_energy($3);}
                | ENERGYFILE EQUAL CHAR {init_energy_file($3);}
		;

initphonon      : PHENFILE EQUAL CHAR {init_ph_energy_file($3);}
                | NOPHEN EQUAL INTEGER {init_no_ph_energy($3);}
                | PHMOFILE EQUAL CHAR {init_ph_mode_file($3);}
                | PHFILETY EQUAL INTEGER {init_ph_file_type($3);}
                | NOPHK EQUAL INTEGER {init_no_ph_k($3);}
                | SCMAXIT EQUAL INTEGER {init_sc_max_iteration($3);}
                | SCRESTART EQUAL INTEGER {init_sc_restart($3,-1);}
                | SCRESTART EQUAL LBRACE INTEGER INTEGER RBRACE {init_sc_restart($4,$5);}
                | SCIDCRIT EQUAL double {init_sc_id_crit($3);}
                | SCRHCRIT EQUAL double {init_sc_rh_crit($3);}
                | SCDISTDEP EQUAL INTEGER {init_sc_dist_dep($3);}
                | SCDIAGDEF EQUAL INTEGER {init_sc_diag_def($3);}
                | SCKCOUP EQUAL INTEGER {init_sc_k_coupling($3);}
                | SCSCALE EQUAL double {init_sc_scale_fact($3);}
                | SCVBOUND EQUAL double {init_sc_vbound($3);}
                | INCINJ EQUAL INTEGER {init_incoherent_injection($3);}
                ;

initstrain      : STRAIN POINT ON EQUAL INTEGER {init_strain($5);}
                | STRAIN POINT EXX EQUAL double {init_Eps_xx($5);}
                | STRAIN POINT EYY EQUAL double {init_Eps_yy($5);}
                | STRAIN POINT EZZ EQUAL double {init_Eps_zz($5);}
                | STRAIN POINT EXY EQUAL double {init_Eps_xy($5);}
                | STRAIN POINT EXZ EQUAL double {init_Eps_xz($5);}
                | STRAIN POINT EYZ EQUAL double {init_Eps_yz($5);}
                | STRAIN POINT ZETA EQUAL double {init_zeta($5);}
                | STRAINMOD EQUAL INTEGER {init_strain_model($3);}
                | ALPHA EQUAL double {init_alpha($3);}
                | BETA EQUAL double {init_beta($3);}
                | STRAIN_BC EQUAL INTEGER {init_strain_bc($3);}
		;

initschottky    : SCHOTTKY POINT ACTIVE EQUAL INTEGER {init_schottky($5);}
                | SCHOTTKY POINT TYPE EQUAL INTEGER {init_schottky_type($5,$5);}
                | SCHOTTKY POINT TYPE EQUAL LBRACE INTEGER INTEGER RBRACE {init_schottky_type($6,$7);}
                | SCHOTTKY POINT BARRIER EQUAL double {init_schottky_barrier($5,$5);}
                | SCHOTTKY POINT BARRIER EQUAL LBRACE double double RBRACE {init_schottky_barrier($6,$7);}
                | SCHOTTKY POINT VIRTUALCB EQUAL double {init_schottky_virtual_cb($5,$5);}
                | SCHOTTKY POINT VIRTUALCB EQUAL LBRACE double double RBRACE {init_schottky_virtual_cb($6,$7);}
		;

initstdomain    : NOSTRAIND EQUAL INTEGER {init_no_strain_domain($3);}
                | SDOM LPAR INTEGER RPAR POINT XMIN EQUAL expr {init_stdomain_coord($3-1,$8,0);}
                | SDOM LPAR INTEGER RPAR POINT XMAX EQUAL expr {init_stdomain_coord($3-1,$8,1);}
                | SDOM LPAR INTEGER RPAR POINT YMIN EQUAL expr {init_stdomain_coord($3-1,$8,2);}
                | SDOM LPAR INTEGER RPAR POINT YMAX EQUAL expr {init_stdomain_coord($3-1,$8,3);}
                | SDOM LPAR INTEGER RPAR POINT ZMIN EQUAL expr {init_stdomain_coord($3-1,$8,4);}
                | SDOM LPAR INTEGER RPAR POINT ZMAX EQUAL expr {init_stdomain_coord($3-1,$8,5);}
                | SDOM LPAR INTEGER RPAR POINT EPSVEC EQUAL LBRACE double double double RBRACE {init_stdomain_eps_vec($3-1,$9,$10,$11,0,0,0);}
                | SDOM LPAR INTEGER RPAR POINT EPSVEC EQUAL LBRACE double double double double double double RBRACE {init_stdomain_eps_vec($3-1,$9,$10,$11,$12,$13,$14);}
                ;

initroughness   : ROUGHNESS POINT ON EQUAL INTEGER {init_roughness($5);}
                | ROUGHNESS POINT SEED EQUAL INTEGER {init_roughness_seed($5);}
                | ROUGHNESS POINT FORM EQUAL CHAR {init_roughness_type($5);}
                | ROUGHNESS POINT RMS EQUAL double {init_roughness_rms($5);}
                | ROUGHNESS POINT LMS EQUAL double {init_roughness_lms($5);}
		;

initqmregion    : QMSTART EQUAL expr {init_qmstart($3);}
                | QMSTOP EQUAL expr {init_qmstop($3);}
		;

initcurrent     : FULLCURRENT EQUAL INTEGER {init_full_current($3);}
                ;

initnumerics    : ROBUSTNUM EQUAL INTEGER {init_robust_numerics($3);}
                | CBOOST EQUAL INTEGER {init_convergence_booster($3);}
                ;

initupdatefermi : UPDATEFERM EQUAL INTEGER {init_update_fermi($3);}
                | FERMIL EQUAL double {init_fermi_level($3);}
		;

initupdateat    : UPDATEAT EQUAL INTEGER {init_update_at($3);}
                | ATOMFILE EQUAL CHAR {init_at_file($3);}
                | REPUC EQUAL INTEGER {init_replicate_unit_cell($3);} 
                | READATPOS EQUAL INTEGER {init_read_atom_pos($3);} 
		;

initupdatedop   : UPDATEDOP EQUAL INTEGER {init_update_dop($3);}
                | DOPFILE EQUAL CHAR {init_dop_file($3);}
		;

initupdateperm  : UPDATEPERM EQUAL INTEGER {init_update_perm($3);}
                | PERMFILE EQUAL CHAR {init_perm_file($3);}
		;

initscdata      : EPSWIRE EQUAL double {init_Eps_wire($3,0,0,0);}
                | EPSOX EQUAL double {init_Eps_ox($3,0,0,0);}
                | EPSWIRE EQUAL LBRACE double double RBRACE {init_Eps_wire($4,$5,0,0);}
                | EPSOX EQUAL LBRACE double double RBRACE {init_Eps_ox($4,$5,0,0);}
                | EPSWIRE EQUAL LBRACE double double double RBRACE {init_Eps_wire($4,$5,$6,0);}
                | EPSOX EQUAL LBRACE double double double RBRACE {init_Eps_ox($4,$5,$6,0);}
                | EPSWIRE EQUAL LBRACE double double double double RBRACE {init_Eps_wire($4,$5,$6,$7);}
                | EPSOX EQUAL LBRACE double double double double RBRACE {init_Eps_ox($4,$5,$6,$7);}
                | XIWIRE EQUAL double {init_Xi_wire($3);}
                | PHIM EQUAL double {init_phi_m($3);}
                | TEMP EQUAL double {Temp=$3;init_temp($3);}
		;

initnomat       : NOMAT EQUAL INTEGER {init_no_mat($3);}
                | NOCMAT EQUAL INTEGER {ncm = $3;init_no_ch_mat($3);}
                | NOOMAT EQUAL INTEGER {nom = $3;init_no_ox_mat($3)}
                | NORMAT EQUAL INTEGER {init_no_ro_mat($3)}
                ;

initmtype       : MTYPE LPAR INTEGER RPAR EQUAL CHAR {init_mat_type($3-1,$6);}
                | OTYPE LPAR INTEGER RPAR EQUAL CHAR {init_mat_type(ncm+$3-1,$6);}
                | RTYPE LPAR INTEGER RPAR EQUAL CHAR {init_mat_type(ncm+nom+$3-1,$6);}
                | MID LPAR INTEGER RPAR EQUAL INTEGER {init_mat_id($3-1,$6);}
                | OID LPAR INTEGER RPAR EQUAL INTEGER {init_mat_id(ncm+$3-1,$6);}
		;

initmcs         : MCS LPAR INTEGER RPAR EQUAL CHAR {init_mat_cs($3-1,$6);}
                | OCS LPAR INTEGER RPAR EQUAL CHAR {init_mat_cs(ncm+$3-1,$6);}
                | RCS LPAR INTEGER RPAR EQUAL CHAR {init_mat_cs(ncm+nom+$3-1,$6);}
		;

initmat         : MCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr expr RBRACE {init_mat($3-1,$5-1,$9,$10,$11);}
                | OCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr expr RBRACE {init_mat(ncm+$3-1,$5-1,$9,$10,$11);}
                | RCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr expr RBRACE {init_mat(ncm+nom+$3-1,$5-1,$9,$10,$11);}
                | MCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_mat($3-1,$5-1,$9,$10,0.0);}
                | OCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_mat(ncm+$3-1,$5-1,$9,$10,0.0);}
                | RCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_mat(ncm+nom+$3-1,$5-1,$9,$10,0.0);}
                | MCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL expr {init_mat($3-1,$5-1,$8,0.0,0.0);}
                | OCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL expr {init_mat(ncm+$3-1,$5-1,$8,0.0,0.0);}
                | RCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL expr {init_mat(ncm+nom+$3-1,$5-1,$8,0.0,0.0);}
		;

initradius      : MRADIUS LPAR INTEGER RPAR EQUAL expr {init_mradius($3-1,$6,$6);}
                | ORADIUS LPAR INTEGER RPAR EQUAL expr {init_mradius(ncm+$3-1,$6,$6);}
                | RRADIUS LPAR INTEGER RPAR EQUAL expr {init_mradius(ncm+nom+$3-1,$6,$6);}
                | DRADIUS LPAR INTEGER RPAR EQUAL expr {init_dradius($3-1,$6,$6);}
                | GRADIUS LPAR INTEGER RPAR EQUAL expr {init_gradius($3-1,$6,$6);}
                | GRRADIUS LPAR INTEGER RPAR EQUAL expr {init_grradius($3-1,$6);}
                | MRADIUS LPAR INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_mradius($3-1,$7,$8);}
                | ORADIUS LPAR INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_mradius(ncm+$3-1,$7,$8);}
                | RRADIUS LPAR INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_mradius(ncm+nom+$3-1,$7,$8);}
                | DRADIUS LPAR INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_dradius($3-1,$7,$8);}
                | GRADIUS LPAR INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_gradius($3-1,$7,$8);}
                | GANGLE LPAR INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_gangle($3-1,$7,$8);}
		;

initdimvar      : LS EQUAL expr {Ls=$3;init_Ls($3);}
                | LC EQUAL expr {Lc=$3;init_Lc($3);}
                | LD EQUAL expr {Ld=$3;init_Ld($3);} 
                | TS EQUAL double {ts=$3;}
                | TC EQUAL double {tc=$3;}
                | TD EQUAL double {td=$3;}
                | HS EQUAL double {hs=$3;}
                | HC EQUAL double {hc=$3;}
                | HD EQUAL double {hd=$3;}
                | TOX EQUAL double {tox=$3;}
                | TGR EQUAL double {tground=$3;}
                | HOX EQUAL double {hox=$3;}
                | X0 EQUAL double {x0=$3;}
                | Y0 EQUAL double {y0=$3;}
                | Z0 EQUAL double {z0=$3;}
                | YC EQUAL double {yc=$3;}
                | ZC EQUAL double {zc=$3;}
                | LQW EQUAL double {Lqw=$3;}
                | LB EQUAL double {Lb=$3;} 
                | LB1 EQUAL double {Lb1=$3;}
                | LB2 EQUAL double {Lb2=$3;}
                | LSP EQUAL double {Lsp=$3;}
                | LSP1 EQUAL double {Lsp1=$3;}
                | LSP2 EQUAL double {Lsp2=$3;}
		;

initgate        : NOGATE EQUAL INTEGER {init_no_gate($3);}
                | GCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr expr RBRACE {init_gate($3-1,$5-1,$9,$10,$11);}
                | GCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_gate($3-1,$5-1,$9,$10,0.0);}
                | GTYPE LPAR INTEGER RPAR EQUAL CHAR {init_gate_type($3-1,$6);}
		;

initground      : NOGROUND EQUAL INTEGER {init_no_ground($3);}
                | GRCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr expr RBRACE {init_ground($3-1,$5-1,$9,$10,$11);}
                | GRCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_ground($3-1,$5-1,$9,$10,0.0);}
                | GRTYPE LPAR INTEGER RPAR EQUAL CHAR {init_ground_type($3-1,$6);}
		;

initdoping      : NODOPING EQUAL INTEGER {init_no_doping($3);}
                | DCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr expr RBRACE {init_doping($3-1,$5-1,$9,$10,$11);}
                | DCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL LBRACE expr expr RBRACE {init_doping($3-1,$5-1,$9,$10,0.0);}
                | DCOORD LPAR INTEGER COMMA INTEGER RPAR EQUAL expr {init_doping($3-1,$5-1,$8,0.0,0.0);}
                | DTYPE LPAR INTEGER RPAR EQUAL CHAR {init_doping_type($3-1,$6);}
                | DND LPAR INTEGER RPAR EQUAL expr {init_doping_nd($3-1,$6);}
                | DNA LPAR INTEGER RPAR EQUAL expr {init_doping_na($3-1,$6);}
                | NDS EQUAL double {ND_S=$3;init_NDS($3);}
                | NAS EQUAL double {NA_S=$3;init_NAS($3);}
                | NDD EQUAL double {ND_D=$3;init_NDD($3);}
                | NAD EQUAL double {NA_D=$3;init_NAD($3);}
		;

initdopdomain   : NODOPD EQUAL INTEGER {init_no_doping_domain($3);}
                | DDOM LPAR INTEGER RPAR POINT XMIN EQUAL expr {init_dop_domain_coord($3-1,$8,0);}
                | DDOM LPAR INTEGER RPAR POINT XMAX EQUAL expr {init_dop_domain_coord($3-1,$8,1);}
                | DDOM LPAR INTEGER RPAR POINT YMIN EQUAL expr {init_dop_domain_coord($3-1,$8,2);}
                | DDOM LPAR INTEGER RPAR POINT YMAX EQUAL expr {init_dop_domain_coord($3-1,$8,3);}
                | DDOM LPAR INTEGER RPAR POINT ZMIN EQUAL expr {init_dop_domain_coord($3-1,$8,4);}
                | DDOM LPAR INTEGER RPAR POINT ZMAX EQUAL expr {init_dop_domain_coord($3-1,$8,5);}
                | DDOM LPAR INTEGER RPAR POINT SLOPE EQUAL double {init_dop_domain_slope($3-1,$8,1e10,1e10);}
                | DDOM LPAR INTEGER RPAR POINT SLOPE EQUAL LBRACE double double RBRACE {init_dop_domain_slope($3-1,$9,$10,1e10);}
                | DDOM LPAR INTEGER RPAR POINT SLOPE EQUAL LBRACE double double double RBRACE {init_dop_domain_slope($3-1,$9,$10,$11);}
                | DDOM LPAR INTEGER RPAR POINT DCONC EQUAL expr {init_dop_domain_conc($3-1,$8);}
                ;

initgrid        : GRID EQUAL INTEGER {init_grid($3);}
                | UPDATETETRA EQUAL INTEGER {init_update_tetra($3);}
                | TETRAFILE EQUAL CHAR {init_tetra_file($3);}       
                | UPDATEFIT EQUAL INTEGER {init_update_fitness($3);}
                | FITFILE EQUAL CHAR {init_fitness_file($3);}       
                | UPDATEGRID EQUAL INTEGER {init_update_grid($3);}
       		| GRIDFILE EQUAL CHAR {init_grid_file($3);}       
		;

initvolt        : NVG EQUAL INTEGER {init_nvg($3);}
                | NVS EQUAL INTEGER {init_nvs($3);}
                | NVD EQUAL INTEGER {init_nvd($3);}
                | NTEMP EQUAL INTEGER {init_ntemp($3);}
                | VGMIN EQUAL double {init_vgmin($3);}
                | VGMAX EQUAL double {init_vgmax($3);}
                | VSMIN EQUAL double {init_vsmin($3);}
                | VSMAX EQUAL double {init_vsmax($3);}
                | VDMIN EQUAL double {init_vdmin($3);}
                | VDMAX EQUAL double {init_vdmax($3);}
                | TMIN EQUAL double {init_tmin($3);}
                | TMAX EQUAL double {init_tmax($3);}
		;

initrestart     : RESTART EQUAL LBRACE INTEGER INTEGER INTEGER INTEGER RBRACE {init_restart($4,$5,$6,$7);}
                | RESTART EQUAL INTEGER {init_restart($3,0,0,0);}
                | VTFILE EQUAL CHAR {init_vtot_file($3);}
                | VAFILE EQUAL CHAR {init_vact_file($3);}
                ;

initdirectory   : DIRECTORY EQUAL CHAR {init_directory($3);}
		;

initcommand     : COMMAND LPAR INTEGER RPAR EQUAL CHAR {init_command($6);}
		;

initrelax       : RELAX_ATOMS EQUAL INTEGER {init_relax($3);}
		;

double		: INTEGER { $$ = (double) $1; }
		| DOUBLE
		;

expr            : exprdef { $$ = $1; }
                | expr MINUS exprdef {$$ = $1 - $3;}
                | expr PLUS exprdef { $$ = $1 + $3; }
                | expr MULT exprdef { $$ = $1 * $3; }
                ;

exprdef         : double { $$ = $1; }
                | LS { $$ = Ls; }
                | LC { $$ = Lc; }
                | LD { $$ = Ld; }
                | TS { $$ = ts; }
                | TC { $$ = tc; }
                | TD { $$ = td; }
                | HS { $$ = hs; }
                | HC { $$ = hc; }
                | HD { $$ = hd; }
                | X0 { $$ = x0; }
                | Y0 { $$ = y0; }
                | Z0 { $$ = z0; }
                | TOX { $$ = tox; }
                | TGR { $$ = tground; }
                | HOX { $$ = hox; }
                | YC { $$ = yc; }
                | ZC { $$ = zc; }
                | LQW { $$ = Lqw; }
                | LB { $$ = Lb; }
                | LB1 { $$ = Lb1; }
                | LB2 { $$ = Lb2; }
                | LSP { $$ = Lsp; }
                | LSP1 { $$ = Lsp1; }
                | LSP2 { $$ = Lsp2; }
                | NDS { $$ = ND_S; }
                | NAS { $$ = NA_S; }
                | NDD { $$ = ND_D; }
                | NAD { $$ = NA_D; }
                | UT  { $$ = kB*Temp/e; }
                ;

separator	: SEMICOLON
		;

%%


