%{

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
    
#include "y.tab.h"

/* avoid isatty() function */
#define YY_NEVER_INTERACTIVE 1

int comment_caller;

%}

sign      [-]
int       [0-9]+
signed    ({sign})?{int}
exp       [DdEe]{signed}
float     ((({signed})?\.{signed}|{signed}\.({signed})?)({exp})?)|{signed}{exp}
variable  [a-zA-Z_/][a-zA-Z0-9_/]*

%x comment
%option noyywrap
%option yylineno

%%

"/*"                { comment_caller = YY_START; BEGIN(comment); } /* BEGIN(comment); */
<comment>[^*\n]*    /* eat anything that's not a '*' */
<comment>"*"+[^*/\n]* /* eat up '*'s not followed by '/'s */
<comment>\n         /* eat anything that's not a '*' */
<comment>"*"+"/"    { BEGIN(comment_caller); }


"//".*\n	    /* C++ comments */

transport_type                 { return TRANSTYPE; }
injection_type                 { return INJTYPE; }
electron_weight                { return EWEIGHT; }
open_system                    { return OPENSYS; }
periodic_system                { return PERIODICSYS; }
read_hamiltonian               { return READHAM; }
dsp3                           { return DSP3; }
h_passivation                  { return HPASS; }
tb                             { return TB; }
n_of_modes                     { return NM; }
Nk                             { return NK; }
bs_solver                      { return BSSOLVER; }
update_bs_target               { return UPDATEBST; }
bs_target                      { return BSTARGET; }
last_first                     { return LF; }
mat_name                       { return NAME; }
bg_data                        { return BGDATA; }
mat_binary_x                   { return BINARYX; }
eta_res                        { return ETAR; }
eta                            { return ETA; }
alpha                          { return ALPHA; }
beta                           { return BETA; }
relax_atoms                    { return RELAX_ATOMS; }
strain_bc                      { return STRAIN_BC; }
max_bond_def                   { return MBD; }
a0                             { return LATTICE; }
first_atom                     { return FIRSTATOM; }
lattice_type                   { return LATTYPE; }

plot_all_k                     { return PLOTALLK;}

NPROW                          { return NPR; }
NPCOL                          { return NPC; }
NPCS                           { return NPCS; }
spec_decomp                    { return SPDEC; }
CPU_ppoint                     { return CPU; }
CPU_per_kz_point               { return CPUKZ; }
CPU_per_sample                 { return CPUSAMP; }
CPU_per_vg_point               { return CPUVG; }
CPU_per_vd_point               { return CPUVD; }
CPU_per_temp_point             { return CPUTEMP; }
CPU_per_wire                   { return CPUWI; }
CPU_per_bc                     { return CPUBC; }
all                            { return ALL; } 

poisson_solver                 { return PSOLVER; }
poisson_criterion              { return PCRITERION; }
poisson_iteration              { return PITERATION; }
max_proc_poisson               { return PMAXPROC; }
poisson_inner_criterion        { return PINCRITERION; }
poisson_inner_iteration        { return PINITERATION; }
charge_average                 { return CHAVG; }
charge_transfer                { return CHTRS; }

NDim                           { return NDIM; }
Nky                            { return NKY; }
Nkz                            { return NKZ; }
phiy_file                      { return PHIYFILE; }
phiz_file                      { return PHIZFILE; }
rot_sym                        { return ROTSYM; }

NxFold                         { return NXFOLD; } 
NyFold                         { return NYFOLD; } 
NzFold                         { return NZFOLD; } 

x                              { return X; }
y                              { return Y; }
z                              { return Z; }

Elimit                         { return ELIMIT; }
Emin_tail                      { return EMIN; }
EOffset                        { return EOFFSET; } 
dE_in                          { return DEIN; }
dE_f                           { return DEF; }
dE_sep                         { return DESEP; }
EExt                           { return EEXT; }
NEmax                          { return NE; }
update_energy                  { return UPDATEEN; }
energy_file                    { return ENERGYFILE; }
regular_mesh                   { return REMESH; }

strain_model                   { return STRAINMOD; }

strain                         { return STRAIN; }
on                             { return ON; }
Eps_xx                         { return EXX; }
Eps_yy                         { return EYY; }
Eps_zz                         { return EZZ; }
Eps_xy                         { return EXY; }
Eps_xz                         { return EXZ; }
Eps_yz                         { return EYZ; }
zeta                           { return ZETA; }

no_strain_domain               { return NOSTRAIND; }
strain_domain                  { return SDOM; }
xmin                           { return XMIN; }
xmax                           { return XMAX; }
ymin                           { return YMIN; }
ymax                           { return YMAX; }
zmin                           { return ZMIN; }
zmax                           { return ZMAX; }
Eps_vec                        { return EPSVEC; }

no_doping_domain               { return NODOPD; }
doping_domain                  { return DDOM; }
slope                          { return SLOPE; }
ND_NA                          { return DCONC; }

roughness                      { return ROUGHNESS; }
seed                           { return SEED; }
form                           { return FORM; }
rms                            { return RMS; }
Lms                            { return LMS; }

alloy                          { return ALDIS; }
composition                    { return COMP; }

update_fermi                   { return UPDATEFERM; }
fermi_level                    { return FERMIL; }

Schottky                       { return SCHOTTKY; }
active                         { return ACTIVE; }
type                           { return TYPE; }
barrier                        { return BARRIER; }       
virtual_CB                     { return VIRTUALCB; }

ph_energy_file                 { return PHENFILE; }
no_ph_energy                   { return NOPHEN; }
ph_mode_file                   { return PHMOFILE; }
ph_file_type                   { return PHFILETY; }
no_ph_k                        { return NOPHK; }
sc_max_iteration               { return SCMAXIT; }
sc_restart                     { return SCRESTART; }
sc_id_crit                     { return SCIDCRIT; }
sc_rh_crit                     { return SCRHCRIT; }
sc_dist_dep                    { return SCDISTDEP; }
sc_diag_def                    { return SCDIAGDEF; }
sc_k_coupling                  { return SCKCOUP; }
sc_e_ph_coupling               { return SCEPHCOUP; }
sc_scale_fact                  { return SCSCALE; }
sc_vbound                      { return SCVBOUND; }
sc_memory_fact                 { return SCMEM; }
incoherent_injection           { return INCINJ; }

full_current                   { return FULLCURRENT; }
robust_numerics                { return ROBUSTNUM; }
convergence_booster            { return CBOOST; }

QMstart                        { return QMSTART; }
QMstop                         { return QMSTOP; }

replicate_unit_cell            { return REPUC; }
read_atom_pos                  { return READATPOS; }
update_atoms                   { return UPDATEAT; }
atom_file                      { return ATOMFILE; }
update_doping                  { return UPDATEDOP; }
doping_file                    { return DOPFILE; }
update_perm                    { return UPDATEPERM; }
perm_file                      { return PERMFILE; }

Eps_wire                       { return EPSWIRE; }
Eps_ox                         { return EPSOX; }
Xi_wire                        { return XIWIRE; }
phi_m                          { return PHIM; }

grid_accuracy                  { return GRID; }
update_tetra                   { return UPDATETETRA; }
tetra_file                     { return TETRAFILE; }
update_fitness                 { return UPDATEFIT; }
fitness_file                   { return FITFILE; }
update_grid                    { return UPDATEGRID; }
grid_file                      { return GRIDFILE; }

no_mat                         { return NOMAT; }
no_channel_mat                 { return NOCMAT; }
no_oxide_mat                   { return NOOMAT; }
no_rough_mat                   { return NORMAT; }
mat_coord                      { return MCOORD; }
ox_coord                       { return OCOORD; }
ro_coord                       { return RCOORD; }
mat_type                       { return MTYPE; }
mat_id                         { return MID; }
ox_type                        { return OTYPE; }
ox_id                          { return OID; }
ro_type                        { return RTYPE; }
mat_cs                         { return MCS; }
ox_cs                          { return OCS; }
ro_cs                          { return RCS; }
mat_radius                     { return MRADIUS; }
ox_radius                      { return ORADIUS; }
ro_radius                      { return RRADIUS; }

Ls                             { return LS; }
Lc                             { return LC; }
Ld                             { return LD; }

ts                             { return TS; }
tc                             { return TC; }
td                             { return TD; }
hs                             { return HS; }
hc                             { return HC; }
hd                             { return HD; }
tox                            { return TOX; }
tground                        { return TGR; }
hox                            { return HOX; }
x0                             { return X0; }
y0                             { return Y0; }
z0                             { return Z0; }

yc                             { return YC; }
zc                             { return ZC; }

Lqw                            { return LQW; }
Lb                             { return LB; }
Lb1                            { return LB1; }
Lb2                            { return LB2; }
Lsp                            { return LSP; }
Lsp1                           { return LSP1; }
Lsp2                           { return LSP2; }

ND_S                           { return NDS; }
NA_S                           { return NAS; }
ND_D                           { return NDD; }
NA_D                           { return NAD; }

no_gate                        { return NOGATE; }
no_ground                      { return NOGROUND; }
no_doping                      { return NODOPING; }
gate_coord                     { return GCOORD; }
ground_coord                   { return GRCOORD; }
doping_coord                   { return DCOORD; }
gate_type                      { return GTYPE; }
ground_type                    { return GRTYPE; }
doping_type                    { return DTYPE; }
doping_ND                      { return DND; }
doping_NA                      { return DNA; }
doping_radius                  { return DRADIUS; }
gate_radius                    { return GRADIUS; }
ground_radius                  { return GRRADIUS; }
gate_angle                     { return GANGLE; }

NVG                            { return NVG; }
Vgmin                          { return VGMIN; }
Vgmax                          { return VGMAX; }
NVS                            { return NVS; }
Vsmin                          { return VSMIN; }
Vsmax                          { return VSMAX; }
NVD                            { return NVD; }
Vdmin                          { return VDMIN; }
Vdmax                          { return VDMAX; }
NTemp                          { return NTEMP; }
Tmin                           { return TMIN; }
Tmax                           { return TMAX; }

restart                        { return RESTART; }
vtot_file                      { return VTFILE; }
vact_file                      { return VAFILE; }

directory                      { return DIRECTORY; }
command                        { return COMMAND; }

Temp                           { return TEMP; }
UT                             { return UT; }

"="		               { return EQUAL; }
"+"		               { return PLUS; }
"-"                            { return MINUS; }
"*"		               { return MULT; }
"["		               { return LBRACE; }
"]"		               { return RBRACE; }
"("		               { return LPAR; }
")"		               { return RPAR; }
","		               { return COMMA; }
"."		               { return POINT; }
";"		               { BEGIN(INITIAL); return SEMICOLON; }

[ \t\n]+            /* ignore whitespace */

{signed}            { yylval.i = atoi(yytext); 
                      return INTEGER;
		    }

{float}		    { /*  accept FORTRAN DOUBLE PRECISION exponent 
		          notation with 'd' or 'D'                 */
 
		      char* cp;
		      for (cp = yytext; *cp != 0; cp++) 
		          if (*cp == 'D' || *cp == 'd') *cp = 'e';

	              yylval.d = atof(yytext); 
		      return DOUBLE; 
                     }

{variable}           { yylval.c = yytext; return CHAR; }


<<EOF>>             {
                      yy_delete_buffer( YY_CURRENT_BUFFER );
                      yyterminate();
                    }
%%

