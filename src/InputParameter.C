#include <iostream>
#include "InputParameter.H"
#include "Types.H"
#include "AtomStrain.H"

extern PARAM *parameter;
extern WireStructure *nanowire;
extern ENERGY *En;
extern VOLTAGE *voltage;
extern Strain *strain;

void init_parameters()
{

    parameter                = new PARAM();
    nanowire                 = new WireStructure();
    En                       = new ENERGY();
    voltage                  = new VOLTAGE();
	
    nanowire->strain         = new Strain();
    nanowire->rough          = new Roughness();
    nanowire->rough->type    = new char[255];
    nanowire->first_atom     = new char[255];
    nanowire->at_file        = new char[255];
    nanowire->dop_file       = new char[255];
    nanowire->perm_file      = new char[255];
    nanowire->vtot_file      = new char[255];
    nanowire->vact_file      = new char[255];
    nanowire->phiy_file      = new char[255];
    nanowire->phiz_file      = new char[255];
    nanowire->tetra_file     = new char[255];
    nanowire->fit_file       = new char[255];
    nanowire->grid_file      = new char[255];
    nanowire->energy_file    = new char[255];
    nanowire->ph_energy_file = new char[255];
    nanowire->ph_mode_file   = new char[255];
    parameter->command       = new char*[MAX_COMM];
    parameter->mat_name      = new char[255];
    parameter->directory     = new char[255];

    default_parameters();
}

/************************************************************************************************/

void delete_parameters()
{
    int i,j;

    for(i=0;i<nanowire->no_element;i++){
        for(j=0;j<8;j++){
            delete nanowire->mat[i]->p[j];
        }
        delete nanowire->mat[i]->p;
        delete[] nanowire->mat[i]->type;
        delete[] nanowire->mat[i]->cross_section;
        delete nanowire->mat[i];

        for(j=0;j<4;j++){
            delete nanowire->surf[i]->p[j];
        }
        delete nanowire->surf[i]->p;
        delete[] nanowire->surf[i]->type;
        delete[] nanowire->surf[i]->cross_section;
        delete nanowire->surf[i];
    }   // end for

    if (nanowire->mat!=NULL){
        delete nanowire->mat;
    }

    if (nanowire->surf!=NULL){
        delete nanowire->surf;
    }

    for(i=0;i<nanowire->no_gate;i++){
        for(j=0;j<4;j++){
	    delete nanowire->gate[i]->p[j];
        }
        delete nanowire->gate[i]->p;
        delete[] nanowire->gate[i]->type;
        delete nanowire->gate[i];
    }  // end for

    if (nanowire->gate!=NULL){
        delete nanowire->gate;
    }

    for(i=0;i<nanowire->no_ground;i++){
        for(j=0;j<4;j++){
            delete nanowire->ground[i]->p[j];
        }
        delete nanowire->ground[i]->p;
        delete[] nanowire->ground[i]->type;
        delete nanowire->ground[i];
    }  // end for

    if (nanowire->ground!=NULL){
        delete nanowire->ground;
    }

    for(i=0;i<nanowire->no_doping;i++){
        for(j=0;j<8;j++){
            delete nanowire->doping[i]->p[j];
        }
        delete nanowire->doping[i]->p;
        delete[] nanowire->doping[i]->type;
        delete nanowire->doping[i];
    }

    if (nanowire->doping!=NULL){
        delete nanowire->doping;
    }

    for(i=0;i<nanowire->no_strain_domain;i++){
        delete nanowire->strain_domain[i];
    }

    if(nanowire->strain_domain!=NULL){
        delete nanowire->strain_domain;
    }

    for(i=0;i<nanowire->no_doping_domain;i++){
        delete nanowire->doping_domain[i];
    }

    if(nanowire->doping_domain!=NULL){
        delete nanowire->doping_domain;
    }

    if (nanowire->first_atom!=NULL){
        delete[] nanowire->first_atom;
    }

    if (nanowire->at_file!=NULL){
        delete[] nanowire->at_file;
    }

    if (nanowire->dop_file!=NULL){
        delete[] nanowire->dop_file;
    }

    if (nanowire->perm_file!=NULL){
        delete[] nanowire->perm_file;
    }

    if (nanowire->vtot_file!=NULL){
        delete[] nanowire->vtot_file;
    }

    if (nanowire->vact_file!=NULL){
        delete[] nanowire->vact_file;
    }

    if (nanowire->phiy_file!=NULL){
	delete[] nanowire->phiy_file;
    }

    if (nanowire->phiz_file!=NULL){
	delete[] nanowire->phiz_file;
    }

    if (nanowire->tetra_file!=NULL){
	delete[] nanowire->tetra_file;
    }

    if (nanowire->fit_file!=NULL){   
        delete[] nanowire->fit_file; 
    }                                 

    if (nanowire->grid_file!=NULL){    
        delete[] nanowire->grid_file;  
    }

    if (nanowire->energy_file!=NULL){    
        delete[] nanowire->energy_file;  
    }  

    if (nanowire->ph_energy_file!=NULL){    
        delete[] nanowire->ph_energy_file;  
    } 

    if (nanowire->ph_mode_file!=NULL){    
        delete[] nanowire->ph_mode_file;  
    } 

    if (nanowire->strain!=NULL){
        delete nanowire->strain;
    }

    if (nanowire->rough->type!=NULL){
        delete nanowire->rough->type;
    }

    if (nanowire->rough!=NULL){
        delete nanowire->rough;
    }

    if (nanowire->Schottky!=NULL){
        delete nanowire->Schottky;
    }
    
    if(nanowire!=NULL){delete nanowire;}

    for(i=0;i<parameter->no_comm;i++){
        delete[] parameter->command[i];
    }

    if(parameter->command!=NULL){
    delete[] parameter->command;
    }

    if(parameter->mat_name!=NULL){
        delete[] parameter->mat_name;
    }

    if(parameter->directory!=NULL){
        delete[] parameter->directory;
    }

    if(parameter!=NULL){delete parameter;}
    if(En!=NULL){delete En;}
    if(voltage!=NULL){delete voltage;}

}

/************************************************************************************************/

void default_parameters()
{
    strcpy(nanowire->first_atom,"cation");
    strcpy(nanowire->rough->type,"square");
    nanowire->restart                  = 0;
    nanowire->IG_start                 = 0;
    nanowire->IS_start                 = 0;
    nanowire->ID_start                 = 0;
    nanowire->no_element               = 0;
    nanowire->no_ch_element            = 0;
    nanowire->no_ox_element            = 0;
    nanowire->no_ro_element            = 0;
    nanowire->no_strain_domain         = 0;
    nanowire->no_doping_domain         = 0;
    nanowire->a0                       = 0.543;
    nanowire->c0                       = 0.543;
    nanowire->u0                       = 0.0;
    nanowire->x[0]                     = 1;
    nanowire->x[1]                     = 0;
    nanowire->x[2]                     = 0;
    nanowire->y[0]                     = 0;
    nanowire->y[1]                     = 1;
    nanowire->y[2]                     = 0;
    nanowire->z[0]                     = 0;
    nanowire->z[1]                     = 0;
    nanowire->z[2]                     = 1;
    nanowire->max_bond_deformation     = 0.1;
    nanowire->NDim                     = 3;
    nanowire->NxFold                   = 1;
    nanowire->NyFold                   = 1;
    nanowire->NzFold                   = 1;
    nanowire->replicate_unit_cell      = 0;
    nanowire->read_atom_pos            = 0;
    nanowire->strain->on               = 0;
    nanowire->strain->Eps_xx           = 0.0;
    nanowire->strain->Eps_yy           = 0.0;
    nanowire->strain->Eps_zz           = 0.0;
    nanowire->strain->Eps_xy           = 0.0;
    nanowire->strain->Eps_xz           = 0.0;
    nanowire->strain->Eps_yz           = 0.0;
    nanowire->strain->zeta             = 0.0;   
    nanowire->rough->on                = 0;
    nanowire->rough->seed              = -1;
    nanowire->rough->rms               = 0.0;
    nanowire->rough->Lms               = 0.0;
    nanowire->QMregion                 = 0;
    nanowire->update_energy            = 0;
    nanowire->update_fermi             = 0;
    nanowire->update_at                = 0;
    nanowire->update_dop               = 0;
    nanowire->update_perm              = 0;
    nanowire->dsp3                     = 30;
    nanowire->hydrogen_passivation     = 0;
    nanowire->grid_accuracy            = 2;
    nanowire->update_tetra             = 0;
    nanowire->update_fitness           = 0;
    nanowire->update_grid              = 0;
    nanowire->Eps_wire[0]              = 11.9;
    nanowire->Eps_ox[0]                = 3.9;
    nanowire->Xi_wire                  = 4.05;
    nanowire->phi_m                    = 4.25;
    nanowire->no_gate                  = 0;
    nanowire->no_doping                = 0;
    nanowire->no_ground                = 0;
    nanowire->Ls                       = 10.0;
    nanowire->Lc                       = 10.0;
    nanowire->Ld                       = 10.0;
    nanowire->ND_S                     = 1e26;
    nanowire->NA_S                     = 1e12;
    nanowire->ND_D                     = 1e26;
    nanowire->NA_D                     = 1e12;
    nanowire->VFF_alpha                = 0.0;
    nanowire->VFF_beta                 = 0.0;
    nanowire->relax_atoms              = 0;
    nanowire->strain_bc                = 0;
    nanowire->charge_average           = 0;
    nanowire->charge_transfer          = 1;
    nanowire->no_ph_energy             = 30;
    nanowire->no_ph_k                  = -1;
    nanowire->ph_file_type             = 1;
    nanowire->sc_max_iteration         = 100;
    nanowire->sc_restart               = 0;
    nanowire->sc_step                  = 20;
    nanowire->sc_id_crit               = 1e-2;
    nanowire->sc_rh_crit               = 1e-3;
    nanowire->sc_dist_dep              = 0;
    nanowire->sc_diag_def              = 0;
    nanowire->sc_k_coupling            = 1;
    nanowire->sc_scale_fact            = 1.0;
    nanowire->sc_vbound                = 0.0;
    nanowire->incoherent_injection     = 0;
    nanowire->full_current             = 0;
    nanowire->robust_numerics          = 0;
    nanowire->convergence_booster      = 0;
    nanowire->mat                      = new MAT*[max(nanowire->no_element,1)];
    nanowire->surf                     = new SURFACE*[max(nanowire->no_element,1)];
    nanowire->gate                     = new GATE*[max(nanowire->no_gate,1)];
    nanowire->doping                   = new DOPING*[max(nanowire->no_doping,1)];
    nanowire->ground                   = new GATE*[max(nanowire->no_ground,1)];
    nanowire->strain_domain            = new StrainDomain*[max(nanowire->no_strain_domain,1)];
    nanowire->doping_domain            = new DopingDomain*[max(nanowire->no_doping_domain,1)];
    nanowire->Schottky                 = new Contact();
    nanowire->Schottky->active         = 0;
    nanowire->Schottky->type[0]        = 0;
    nanowire->Schottky->barrier[0]     = 0.2;
    nanowire->Schottky->virtual_CB[0]  = -0.2;
    nanowire->Schottky->type[1]        = 1;
    nanowire->Schottky->barrier[1]     = 0.2;
    nanowire->Schottky->virtual_CB[1]  = -0.2;

    En->Elimit                         = 50e-3;
    En->Emin_tail                      = 4e-3;
    En->EOffset                        = 0.5;
    En->dE_in                          = 1e-3;
    En->dE_f                           = 5e-3;
    En->dE_sep                         = 1e-4;
    En->EExt                           = 0.0;
    En->NEmax                          = 100;
    En->regular_mesh                   = 0;
    
    strcpy(parameter->mat_name,"Si");
    parameter->strain_model            = 1;
    parameter->mat_binary_x[0]         = 0.0;
    parameter->mat_binary_x[1]         = 0.0;
    parameter->mat_binary_x[2]         = 0.0;
    parameter->tb                      = 10;
    parameter->last_first              = 1;
    parameter->n_of_modes              = 16;
    parameter->Nk                      = 101;
    parameter->Nky                     = 1;
    parameter->Nkz                     = 1;
    parameter->bs_solver               = 0;
    parameter->rot_sym                 = 0;
    parameter->eta_res                 = 0.0;
    parameter->eta                     = 0.0;
    parameter->plot_all_k              = 0;
    parameter->CPU_per_sample          = -1;
    parameter->CPU_per_temp_point      = -1;
    parameter->CPU_per_vd_point        = -1;
    parameter->CPU_per_vg_point        = -1;
    parameter->CPU_per_kz_point        = -1;
    parameter->CPU_ppoint              = 1;
    parameter->CPU_per_wire            = 512;
    parameter->CPU_per_bc              = 2;
    parameter->NPROW                   = 1;
    parameter->NPCOL                   = 1;
    parameter->NPCS                    = 1;
    parameter->spec_decomp             = 0;
    parameter->no_comm                 = 0;
    parameter->Temp                    = 300;
    parameter->poisson_solver          = 1;
    parameter->poisson_criterion       = 1e-3;
    parameter->poisson_iteration       = 15;
    parameter->max_proc_poisson        = 64;
    parameter->poisson_inner_criterion = 1e-3;
    parameter->poisson_inner_iteration = 15;
    parameter->lattype                 = 1;//ZB by default
    parameter->transport_type          = 0;

    strcpy(parameter->directory,".");	

    voltage->NVG                       = 1;
    voltage->Vgmin                     = 0.0;
    voltage->Vgmax                     = 0.0;
    voltage->NVS                       = 1;
    voltage->Vsmin                     = 0.0;
    voltage->Vsmax                     = 0.0;
    voltage->NVD                       = 1;
    voltage->Vdmin                     = 0.0;
    voltage->Vdmax                     = 0.0;
    voltage->NTEMP                     = 1;
    voltage->Tmin                      = 300.0;
    voltage->Tmax                      = 300.0;

}

/************************************************************************************************/

extern "C"{
    
    void init_restart(int restart,int IG_start,int IS_start,int ID_start)
    {
        nanowire->restart  = restart;
        nanowire->IG_start = IG_start;
        nanowire->IS_start = IS_start;
        nanowire->ID_start = ID_start;
    }

    void init_vtot_file(char *vtot_file)
    {
        strcpy(nanowire->vtot_file,vtot_file);
    }

    void init_vact_file(char *vact_file)
    {
        strcpy(nanowire->vact_file,vact_file);
    }
    void init_tb(int tb)
    {
        parameter->tb = tb;
    }
    
    void init_dsp3(double dsp3)
    {
        nanowire->dsp3 = dsp3;
    }

    void init_hpass(int hpass)
    {
        nanowire->hydrogen_passivation = hpass;
    }
    
    void init_a0_c0_u0(double a0,double c0,double u0)
    {
        nanowire->a0 = a0;
	nanowire->c0 = c0;
	nanowire->u0 = u0;
    }
    
    void init_first_atom(char *first_atom)
    {
        strcpy(nanowire->first_atom,first_atom);
    }
    
    void init_mbd(double mbd)
    {
        nanowire->max_bond_deformation = mbd;
    }
    
    void init_n_of_modes(int n_of_modes)
    {
        parameter->n_of_modes = n_of_modes;
    }

    void init_nk(int Nk)
    {
        parameter->Nk = Nk;
    }
    
    void init_last_first(int last_first)
    {
        parameter->last_first = last_first;
    }

    void init_ndim(int NDim)
    {
        nanowire->NDim = NDim;
    }

    void init_nxfold(int NxFold)
    {
        nanowire->NxFold = NxFold;
    }

    void init_nyfold(int NyFold)
    {
        nanowire->NyFold = NyFold;
    }

    void init_nzfold(int NzFold)
    {
        nanowire->NzFold = NzFold;
    }

    void init_nky(int Nky)
    {
        parameter->Nky = Nky;
    }

    void init_nkz(int Nkz)
    {
        parameter->Nkz = Nkz;
    } 

    void init_phiy_file(char *phiy_file)
    {
        strcpy(nanowire->phiy_file,phiy_file);
    }

    void init_phiz_file(char *phiz_file)
    {
        strcpy(nanowire->phiz_file,phiz_file);
    }

    void init_replicate_unit_cell(int replicate_unit_cell)
    {
        nanowire->replicate_unit_cell = replicate_unit_cell;
    }

    void init_read_atom_pos(int read_atom_pos)
    {
        nanowire->read_atom_pos = read_atom_pos;
    }

    void init_rot_sym(int rot_sym)
    {
        parameter->rot_sym = rot_sym;
    }
    
    void init_bs_solver(char *bs_solver)
    {
        if(!strcmp(bs_solver,"sparse")){
	    parameter->bs_solver = 0;
	}
	if(!strcmp(bs_solver,"full")){
	    parameter->bs_solver = 1;
	}
    }

    void init_mat_name(char *mat_name)
    {
        strcpy(parameter->mat_name,mat_name);
    }

    void init_strain_model(int strain_model)
    {
        parameter->strain_model = strain_model;
    }

    void init_mat_binary_x(double x1,double x2,double x3)
    {
        parameter->mat_binary_x[0] = x1;
	parameter->mat_binary_x[1] = x2;
	parameter->mat_binary_x[2] = x3;
    }
    
    void init_eta_res(double eta_res)
    {
        parameter->eta_res = eta_res;
    }
    
    void init_eta(double eta)
    {
        parameter->eta = eta;
    }

    void init_plot_all_k(int plot_all_k)
    {
        parameter->plot_all_k = plot_all_k;
    }
    
    void init_NPROW(int NPROW)
    {
        parameter->NPROW = NPROW;
    }
    
    void init_NPCOL(int NPCOL)
    {
        parameter->NPCOL = NPCOL;
    }

    void init_NPCS(int NPCS)
    {
        parameter->NPCS = NPCS;
    }

    void init_spec_decomp(int spec_decomp)
    {
        parameter->spec_decomp = spec_decomp;
    }
    
    void init_CPU_ppoint(int CPU_ppoint)
    {
        parameter->CPU_ppoint = CPU_ppoint;
    }

    void init_CPU_per_kz_point(int CPU_per_kz_point)
    {
        parameter->CPU_per_kz_point = CPU_per_kz_point;
    }

    void init_CPU_per_sample(int CPU_per_sample)
    {
        parameter->CPU_per_sample = CPU_per_sample;
    }

    void init_CPU_per_vg_point(int CPU_per_vg_point)
    {
        parameter->CPU_per_vg_point = CPU_per_vg_point;
    }
    
    void init_CPU_per_vd_point(int CPU_per_vd_point)
    {
        parameter->CPU_per_vd_point = CPU_per_vd_point;
    }

    void init_CPU_per_temp_point(int CPU_per_temp_point)
    {
        parameter->CPU_per_temp_point = CPU_per_temp_point;
    }

    void init_CPU_per_wire(int CPU_per_wire)
    {
        parameter->CPU_per_wire = CPU_per_wire;
    }

    void init_CPU_per_bc(int CPU_per_bc)
    {
        parameter->CPU_per_bc = CPU_per_bc;
    }

    void init_poisson_solver(int poisson_solver)
    {
        parameter->poisson_solver = poisson_solver;
    }

    void init_poisson_criterion(double poisson_criterion)
    {
        parameter->poisson_criterion = poisson_criterion;
    }

    void init_poisson_iteration(int poisson_iteration)
    {
        parameter->poisson_iteration = poisson_iteration;   
    }
    
    void init_max_proc_poisson(int max_proc_poisson)
    {
        parameter->max_proc_poisson = max_proc_poisson;   
    }

    void init_poisson_inner_criterion(double poisson_inner_criterion)
    {
        parameter->poisson_inner_criterion = poisson_inner_criterion;
    }

    void init_poisson_inner_iteration(int poisson_inner_iteration)
    {
        parameter->poisson_inner_iteration = poisson_inner_iteration;   
    }

    void init_charge_average(int charge_average)
    {
        nanowire->charge_average = charge_average;   
    }

    void init_charge_transfer(int charge_transfer)
    {
        nanowire->charge_transfer = charge_transfer;   
    }
    
    void init_x(double x1,double x2,double x3)
    {
        nanowire->x[0] = x1;
        nanowire->x[1] = x2;
        nanowire->x[2] = x3;
    }
    
    void init_y(double y1,double y2,double y3)
    {
        nanowire->y[0] = y1;
        nanowire->y[1] = y2;
        nanowire->y[2] = y3;
    }
    
    void init_z(double z1,double z2,double z3)
    {
        nanowire->z[0] = z1;
        nanowire->z[1] = z2;
        nanowire->z[2] = z3;
    }
    
    void init_Elimit(double Elimit)
    {
        En->Elimit = Elimit;
    }
    
    void init_Emin_tail(double Emin_tail)
    {
        En->Emin_tail = Emin_tail;
    }
    
    void init_EOffset(double EOffset)
    {
        En->EOffset = EOffset;
    }
    
    void init_dE_in(double dE_in)
    {
        En->dE_in = dE_in;
    }
    
    void init_dE_f(double dE_f)
    {
        En->dE_f = dE_f;
    }
    
    void init_dE_sep(double dE_sep)
    {
        En->dE_sep = dE_sep;
    }

    void init_EExt(double EExt)
    {
        En->EExt = EExt;
    }
    
    void init_NE(int NE)
    {
        En->NEmax = NE;
    }

    void init_reg_mesh(int regular_mesh)
    {
        En->regular_mesh = regular_mesh;
    }
    
    void init_strain(int on)
    {
        nanowire->strain->on = on;
    }
    
    void init_Eps_xx(double Eps_xx)
    {
        nanowire->strain->Eps_xx = Eps_xx;
    }
    
    void init_Eps_yy(double Eps_yy)
    {
        nanowire->strain->Eps_yy = Eps_yy;
    }
    
    void init_Eps_zz(double Eps_zz)
    {
        nanowire->strain->Eps_zz = Eps_zz;
    }

    void init_Eps_xy(double Eps_xy)
    {
        nanowire->strain->Eps_xy = Eps_xy;
    }
    
    void init_Eps_xz(double Eps_xz)
    {
        nanowire->strain->Eps_xz = Eps_xz;
    }
    
    void init_Eps_yz(double Eps_yz)
    {
        nanowire->strain->Eps_yz = Eps_yz;
    }

    void init_zeta(double zeta)
    {
        nanowire->strain->zeta = zeta;
    }

    void init_no_strain_domain(int no_strain_domain)
    {
        int i,j;

        if (nanowire->strain_domain!=NULL){
	    delete nanowire->strain_domain;
        }
        
        nanowire->no_strain_domain = no_strain_domain;
        
        nanowire->strain_domain    = new StrainDomain*[max(no_strain_domain,1)];

        for(i=0;i<no_strain_domain;i++){
            nanowire->strain_domain[i]       = new StrainDomain();
	    nanowire->strain_domain[i]->xmin = -INF;
	    nanowire->strain_domain[i]->xmax = INF;
	    nanowire->strain_domain[i]->ymin = -INF;
	    nanowire->strain_domain[i]->ymax = INF;
	    nanowire->strain_domain[i]->zmin = -INF;
	    nanowire->strain_domain[i]->zmax = INF;
	    for(j=0;j<6;j++){
	        nanowire->strain_domain[i]->Eps_vec[j] = 0.0;
	    }
        }
    }

    void init_stdomain_coord(int index,double coord,int type)
    {
        switch(type){
	case 0:
	    nanowire->strain_domain[index]->xmin = coord;
	    break;
	case 1:
	    nanowire->strain_domain[index]->xmax = coord;
	    break;
	case 2:
	    nanowire->strain_domain[index]->ymin = coord;
	    break;
        case 3:
	    nanowire->strain_domain[index]->ymax = coord;
	    break;
	case 4:
	    nanowire->strain_domain[index]->zmin = coord;
	    break;
	case 5:
	    nanowire->strain_domain[index]->zmax = coord;
	    break;  
	}
    }

    void init_stdomain_eps_vec(int index,double Eps_xx,double Eps_yy,double Eps_zz,\
			       double Eps_xy,double Eps_xz,double Eps_yz)
    {
        nanowire->strain_domain[index]->Eps_vec[0] = Eps_xx;
	nanowire->strain_domain[index]->Eps_vec[1] = Eps_yy;
	nanowire->strain_domain[index]->Eps_vec[2] = Eps_zz;
	nanowire->strain_domain[index]->Eps_vec[3] = Eps_xy;
	nanowire->strain_domain[index]->Eps_vec[4] = Eps_xz;
	nanowire->strain_domain[index]->Eps_vec[5] = Eps_yz;
    }

    void init_schottky(int active)    
    {
        nanowire->Schottky->active = active;    
    }

    void init_schottky_type(int type1,int type2)    
    {
        nanowire->Schottky->type[0] = type1;
        nanowire->Schottky->type[1] = type2;
    }

    void init_schottky_barrier(double barrier1,double barrier2)    
    {
        nanowire->Schottky->barrier[0] = barrier1;
        nanowire->Schottky->barrier[1] = barrier2;
    }

    void init_schottky_virtual_cb(double vcb1,double vcb2)    
    {
        nanowire->Schottky->virtual_CB[0] = vcb1;
        nanowire->Schottky->virtual_CB[1] = vcb2;
    }    

    void init_roughness(int on)
    {
        nanowire->rough->on = on;
    }

    void init_roughness_seed(int seed)
    {
        nanowire->rough->seed = seed;
    }
    
    void init_roughness_type(char *type)
    {
        strcpy(nanowire->rough->type,type);
    }
    
    void init_roughness_rms(double rms)
    {
        nanowire->rough->rms = rms;
    }
    
    void init_roughness_lms(double Lms)
    {
        nanowire->rough->Lms = Lms;
    }

    void init_qmstart(double start)
    {
        nanowire->QMstart = start;
	nanowire->QMregion++;
    }

    void init_qmstop(double stop)
    {
        nanowire->QMstop = stop;
	nanowire->QMregion++;
    }

    void init_update_fermi(int update_fermi)
    {
        nanowire->update_fermi = update_fermi;
    }

    void init_fermi_level(double Ef0)
    {
        nanowire->Ef0 = Ef0;
    }

    void init_update_at(int update_at)
    {
        nanowire->update_at = update_at;
    }

    void init_update_dop(int update_dop)
    {
        nanowire->update_dop = update_dop;
    }

    void init_update_perm(int update_perm)
    {
        nanowire->update_perm = update_perm;
    }

    void init_Eps_wire(double epsw_r1,double epsw_r2,double epsw_r3,double epsw_r4)
    {
        nanowire->Eps_wire[0] = epsw_r1;
	nanowire->Eps_wire[1] = epsw_r2;
	nanowire->Eps_wire[2] = epsw_r3;
	nanowire->Eps_wire[3] = epsw_r4;
    }
    
    void init_Eps_ox(double epso_r1,double epso_r2,double epso_r3,double epso_r4)
    {
        nanowire->Eps_ox[0] = epso_r1;
	nanowire->Eps_ox[1] = epso_r2;
	nanowire->Eps_ox[2] = epso_r3;
	nanowire->Eps_ox[3] = epso_r4;
    }
    
    void init_Xi_wire(double Xi_wire)
    {
        nanowire->Xi_wire = Xi_wire;
    }
    
    void init_phi_m(double phi_m)
    {
        nanowire->phi_m = phi_m;
    }

    void init_grid(int grid_accuracy)
    {
        nanowire->grid_accuracy = grid_accuracy;
    }

    void init_update_tetra(int update_tetra)
    {
        nanowire->update_tetra = update_tetra;
    }

    void init_tetra_file(char *tetra_file)
    {
        strcpy(nanowire->tetra_file,tetra_file);
    }                                            
    
    void init_update_fitness(int update_fitness)
    {
        nanowire->update_fitness = update_fitness;
    }
                                                 
    void init_fitness_file(char *fit_file)       
    {                                            
        strcpy(nanowire->fit_file,fit_file); 
    }                                                

    void init_update_grid(int update_grid)     
    {                                                
        nanowire->update_grid = update_grid;   
    }                                                
                                                 
    void init_grid_file(char *grid_file)           
    {                                                
        strcpy(nanowire->grid_file,grid_file);         
    }

    void init_update_energy(int update_energy)     
    {                                                
        nanowire->update_energy = update_energy;   
    }                                                
                                                 
    void init_energy_file(char *energy_file)           
    {                                                
        strcpy(nanowire->energy_file,energy_file);         
    }

    void init_ph_energy_file(char *ph_energy_file)           
    {                                                
        strcpy(nanowire->ph_energy_file,ph_energy_file);         
    }

    void init_no_ph_energy(int no_ph_energy)           
    {                                                
        nanowire->no_ph_energy = no_ph_energy;         
    }

    void init_ph_mode_file(char *ph_mode_file)           
    {                                                
        strcpy(nanowire->ph_mode_file,ph_mode_file);         
    }

    void init_ph_file_type(int ph_file_type)           
    {                                                
        nanowire->ph_file_type = ph_file_type;         
    }

    void init_no_ph_k(int no_ph_k)           
    {                                                
        nanowire->no_ph_k = no_ph_k;         
    }

    void init_sc_max_iteration(int sc_max_iteration)           
    {                                                
        nanowire->sc_max_iteration = sc_max_iteration;         
    }

    void init_sc_restart(int sc_restart,int sc_step)           
    {                                                
        nanowire->sc_restart = sc_restart;
	
	if(sc_step>0){
	    nanowire->sc_step = sc_step;
	}
    }

    void init_sc_id_crit(double sc_id_crit)           
    {                                                
        nanowire->sc_id_crit = sc_id_crit;         
    }

    void init_sc_rh_crit(double sc_rh_crit)           
    {                                                
        nanowire->sc_rh_crit = sc_rh_crit;         
    }

    void init_sc_dist_dep(int sc_dist_dep)           
    {                                                
        nanowire->sc_dist_dep = sc_dist_dep;         
    }

    void init_sc_diag_def(int sc_diag_def)           
    {                                                
        nanowire->sc_diag_def = sc_diag_def;         
    }

    void init_sc_k_coupling(int sc_k_coupling)           
    {                                                
        nanowire->sc_k_coupling = sc_k_coupling;         
    }

    void init_sc_scale_fact(double sc_scale_fact)           
    {                                                
        nanowire->sc_scale_fact = sc_scale_fact;         
    }

    void init_sc_vbound(double sc_vbound)           
    {                                                
        nanowire->sc_vbound = sc_vbound;         
    }

    void init_incoherent_injection(int incoherent_injection)           
    {                                                
        nanowire->incoherent_injection = incoherent_injection;         
    }

    void init_full_current(int full_current)           
    {                                                
        nanowire->full_current = full_current;         
    }

    void init_robust_numerics(int robust_numerics)           
    {                                                
        nanowire->robust_numerics = robust_numerics;         
    }

    void init_convergence_booster(int convergence_booster)           
    {                                                
        nanowire->convergence_booster = convergence_booster;         
    }

    void init_no_mat(int no_mat)
    {
        int i,j;

        if(nanowire->mat!=NULL){
            delete nanowire->mat;
        }

        if(nanowire->surf!=NULL){
            delete nanowire->surf;
        }

        nanowire->no_element = no_mat;
        
        nanowire->mat        = new MAT*[max(no_mat,1)];
        nanowire->surf       = new SURFACE*[max(no_mat,1)];

        for(i=0;i<no_mat;i++){
            nanowire->mat[i]     = new MAT();
            nanowire->surf[i]    = new SURFACE();
            nanowire->mat[i]->p  = new POINT3D*[8];
            nanowire->surf[i]->p = new POINT2D*[4];

            for(j=0;j<8;j++){
                nanowire->mat[i]->p[j]  = new POINT3D();
            }
            for(j=0;j<4;j++){
                nanowire->surf[i]->p[j] = new POINT2D();
            }
        }
    }
    
    void init_no_ch_mat(int no_ch_mat)
    {
        nanowire->no_ch_element = no_ch_mat;
    }
    
    void init_no_ox_mat(int no_ox_mat)
    {
        nanowire->no_ox_element = no_ox_mat;
    }

    void init_no_ro_mat(int no_ro_mat)
    {
        nanowire->no_ro_element = no_ro_mat;
    }
    
    void init_mat_type(int mat_index,char* type)
    {
        if(mat_index<nanowire->no_element){

            nanowire->mat[mat_index]->type  = new char[255];
            nanowire->surf[mat_index]->type = new char[255];

            strcpy(nanowire->mat[mat_index]->type,type);
            strcpy(nanowire->surf[mat_index]->type,type);

	    nanowire->mat[mat_index]->id_number  = 1;
	    nanowire->surf[mat_index]->id_number = 1;
        }
    }

    void init_mat_id(int mat_index,int id_number)
    {
        nanowire->mat[mat_index]->id_number  = id_number;
        nanowire->surf[mat_index]->id_number = id_number;
    }

    void init_mat_cs(int mat_index,char* cs)
    {
        if(mat_index<nanowire->no_element){

            nanowire->mat[mat_index]->cross_section  = new char[255];
            nanowire->surf[mat_index]->cross_section = new char[255];

            strcpy(nanowire->mat[mat_index]->cross_section,cs);
            strcpy(nanowire->surf[mat_index]->cross_section,cs);
        }
    }
    
    void init_mat(int mat_index,int p_index,double x,double y,double z)
    {
        int surf_index[8] = {0,-1,-1,1,3,-1,-1,2};

        if(mat_index<nanowire->no_element){
            nanowire->mat[mat_index]->p[p_index]->coord[0] = x;
            nanowire->mat[mat_index]->p[p_index]->coord[1] = y;
            nanowire->mat[mat_index]->p[p_index]->coord[2] = z;

            if((p_index == 0)||(p_index == 3)||(p_index == 4)||(p_index == 7)){
                nanowire->surf[mat_index]->p[surf_index[p_index]]->coord[0] = y;
                nanowire->surf[mat_index]->p[surf_index[p_index]]->coord[1] = z;
            }
        }else{
            printf("Warning: you defined more material regions than no_mat\n");
        }
    }

    void init_mradius(int mat_index,double radius_1,double radius_2)
    {

        if(mat_index<nanowire->no_element){
            
            nanowire->mat[mat_index]->radius[0]  = radius_1;
	    nanowire->mat[mat_index]->radius[1]  = radius_2;
            nanowire->surf[mat_index]->radius[0] = radius_1;
	    nanowire->surf[mat_index]->radius[1] = radius_1;
            
        }else{
            printf("Warning: you defined more material regions than no_mat\n");
        }
    }

    void init_no_gate(int no_gate)
    {
        int i,j;

        if (nanowire->gate!=NULL){
	    delete nanowire->gate;
        }
        
        nanowire->no_gate = no_gate;
        
        nanowire->gate    = new GATE*[max(no_gate,1)];

        for(i=0;i<no_gate;i++){

            nanowire->gate[i]     = new GATE();
            nanowire->gate[i]->p  = new POINT3D*[4];

            for(j=0;j<4;j++){
                nanowire->gate[i]->p[j]  = new POINT3D();
            }
        }
    }

    void init_no_ground(int no_ground)
    {
        int i,j;

        if (nanowire->ground!=NULL){
	    delete nanowire->ground;
        }
        
        nanowire->no_ground = no_ground;
        
        nanowire->ground    = new GATE*[max(no_ground,1)];

        for(i=0;i<no_ground;i++){

            nanowire->ground[i]     = new GATE();
            nanowire->ground[i]->p  = new POINT3D*[4];

            for(j=0;j<4;j++){
                nanowire->ground[i]->p[j]  = new POINT3D();
            }
        }
    }

    void init_no_doping(int no_doping)
    {
        int i,j;

        if (nanowire->doping!=NULL){
	    delete nanowire->doping;
        }
        
        nanowire->no_doping = no_doping;
        
        nanowire->doping    = new DOPING*[max(no_doping,1)];

        for(i=0;i<no_doping;i++){
            nanowire->doping[i]     = new DOPING();
            nanowire->doping[i]->p  = new POINT3D*[8];
            for(j=0;j<8;j++){
                nanowire->doping[i]->p[j]  = new POINT3D();
            }
        }
    }

    void init_gate(int gate_index,int p_index,double x,double y,double z)
    {

        if(gate_index<nanowire->no_gate){
            nanowire->gate[gate_index]->p[p_index]->coord[0] = x;
            nanowire->gate[gate_index]->p[p_index]->coord[1] = y;
            nanowire->gate[gate_index]->p[p_index]->coord[2] = z;
        }else{
            printf("Warning: you defined more gate regions than no_gate\n");
        }
    }

    void init_gradius(int gate_index,double radius_1,double radius_2)
    {

        if(gate_index<nanowire->no_gate){
            nanowire->gate[gate_index]->radius[0] = radius_1;
	    nanowire->gate[gate_index]->radius[1] = radius_2;
        }else{
            printf("Warning: you defined more gate regions than no_gate\n");
        }
    }

    void init_gangle(int gate_index,double angle_1,double angle_2)
    {

        if(gate_index<nanowire->no_gate){
            nanowire->gate[gate_index]->angle[0] = angle_1;
	    nanowire->gate[gate_index]->angle[1] = angle_2;
        }else{
            printf("Warning: you defined more gate regions than no_gate\n");
        }
    }

    void init_ground(int ground_index,int p_index,double x,double y,double z)
    {

        if(ground_index<nanowire->no_ground){
            nanowire->ground[ground_index]->p[p_index]->coord[0] = x;
            nanowire->ground[ground_index]->p[p_index]->coord[1] = y;
            nanowire->ground[ground_index]->p[p_index]->coord[2] = z;
        }else{
            printf("Warning: you defined more ground regions than no_ground\n");
        }
    }

    void init_grradius(int ground_index,double radius)
    {

        if(ground_index<nanowire->no_ground){
            nanowire->ground[ground_index]->radius[0] = radius;
        }else{
            printf("Warning: you defined more ground regions than no_ground\n");
        }
    }

    void init_doping(int doping_index,int p_index,double x,double y,double z)
    {

        if(doping_index<nanowire->no_doping){
            nanowire->doping[doping_index]->p[p_index]->coord[0] = x;
            nanowire->doping[doping_index]->p[p_index]->coord[1] = y;
            nanowire->doping[doping_index]->p[p_index]->coord[2] = z;
        }else{
            printf("Warning: you defined more doping regions than no_doping\n");
        }
    }

    void init_dradius(int doping_index,double radius_1,double radius_2)
    {

        if(doping_index<nanowire->no_doping){
            nanowire->doping[doping_index]->radius[0] = radius_1;
	    nanowire->doping[doping_index]->radius[1] = radius_2;
        }else{
            printf("Warning: you defined more doping regions than no_doping\n");
        }
    }

    void init_gate_type(int gate_index,char* type)
    {
        if(gate_index<nanowire->no_gate){
            nanowire->gate[gate_index]->type  = new char[255];
            strcpy(nanowire->gate[gate_index]->type,type);
        }
    }

    void init_ground_type(int ground_index,char* type)
    {
        if(ground_index<nanowire->no_ground){
            nanowire->ground[ground_index]->type  = new char[255];
            strcpy(nanowire->ground[ground_index]->type,type);
        }
    }

    void init_doping_type(int doping_index,char* type)
    {
        if(doping_index<nanowire->no_doping){
            nanowire->doping[doping_index]->type  = new char[255];
            strcpy(nanowire->doping[doping_index]->type,type);
        }
    }

    void init_doping_nd(int doping_index,double ND)
    {
        if(doping_index<nanowire->no_doping){
            nanowire->doping[doping_index]->ND  = ND;
        }
    }

    void init_doping_na(int doping_index,double NA)
    {
        if(doping_index<nanowire->no_doping){
            nanowire->doping[doping_index]->NA  = NA;
        }
    }

    void init_no_doping_domain(int no_doping_domain)
    {
        int i,j;

        if (nanowire->doping_domain!=NULL){
	    delete nanowire->doping_domain;
        }
        
        nanowire->no_doping_domain = no_doping_domain;
        
        nanowire->doping_domain    = new DopingDomain*[max(no_doping_domain,1)];

        for(i=0;i<no_doping_domain;i++){
            nanowire->doping_domain[i]       = new DopingDomain();
	    nanowire->doping_domain[i]->xmin = -1000.0;
	    nanowire->doping_domain[i]->xmax = 1000.0;
	    nanowire->doping_domain[i]->ymin = -1000.0;
	    nanowire->doping_domain[i]->ymax = 1000.0;
	    nanowire->doping_domain[i]->zmin = -1000.0;
	    nanowire->doping_domain[i]->zmax = 1000.0;
	    for(j=0;j<3;j++){
	        nanowire->doping_domain[i]->slope[j] = INF;
	    }
        }
    }

    void init_dop_domain_coord(int index,double coord,int type)
    {
        switch(type){
	case 0:
	    nanowire->doping_domain[index]->xmin = coord;
	    break;
	case 1:
	    nanowire->doping_domain[index]->xmax = coord;
	    break;
	case 2:
	    nanowire->doping_domain[index]->ymin = coord;
	    break;
        case 3:
	    nanowire->doping_domain[index]->ymax = coord;
	    break;
	case 4:
	    nanowire->doping_domain[index]->zmin = coord;
	    break;
	case 5:
	    nanowire->doping_domain[index]->zmax = coord;
	    break;  
	}
    }

    void init_dop_domain_slope(int index,double Lx,double Ly,double Lz)
    {
        nanowire->doping_domain[index]->slope[0] = Lx;
	nanowire->doping_domain[index]->slope[1] = Ly;
	nanowire->doping_domain[index]->slope[2] = Lz;
    }

    void init_dop_domain_conc(int index,double conc)
    {
        nanowire->doping_domain[index]->conc = conc;
    }

    void init_temp(double Temp)
    {
        parameter->Temp = Temp;
    }

    void init_Ls(double Ls)
    {
        nanowire->Ls = Ls;
    }

    void init_Lc(double Lc)
    {
        nanowire->Lc = Lc;
    }

    void init_Ld(double Ld)
    {
        nanowire->Ld = Ld;
    }

    void init_alpha(double alpha)
    {
        nanowire->VFF_alpha = alpha;
    }

    void init_beta(double beta)
    {
        nanowire->VFF_beta = beta;
    }

    void init_relax(int relax)
    {
        nanowire->relax_atoms = relax;
    }

    void init_strain_bc(int strain_bc)
    {
        nanowire->strain_bc = strain_bc;
    }
        
    void init_NDS(double ND_S)
    {
        nanowire->ND_S = ND_S;
    }

    void init_NAS(double NA_S)
    {
        nanowire->NA_S = NA_S;
    }

    void init_NDD(double ND_D)
    {
        nanowire->ND_D = ND_D;
    }

    void init_NAD(double NA_D)
    {
        nanowire->NA_D = NA_D;
    }

    void init_nvg(int NVG)
    {
        voltage->NVG = NVG;
    }

    void init_vgmin(double Vgmin)
    {
        voltage->Vgmin = Vgmin;
    }

    void init_vgmax(double Vgmax)
    {
        voltage->Vgmax = Vgmax;
    }

    void init_nvs(int NVS)
    {
        voltage->NVS = NVS;
    }

    void init_vsmin(double Vsmin)
    {
        voltage->Vsmin = Vsmin;
    }

    void init_vsmax(double Vsmax)
    {
        voltage->Vsmax = Vsmax;
    }

    void init_nvd(int NVD)
    {
        voltage->NVD = NVD;
    }

    void init_vdmin(double Vdmin)
    {
        voltage->Vdmin = Vdmin;
    }

    void init_vdmax(double Vdmax)
    {
        voltage->Vdmax = Vdmax;
    }

    void init_ntemp(int NTEMP)
    {
        voltage->NTEMP = NTEMP;
    }

    void init_tmin(double Tmin)
    {
        voltage->Tmin = Tmin;
    }

    void init_tmax(double Tmax)
    {
        voltage->Tmax = Tmax;
    }

    void init_at_file(char *at_file)
    {
        strcpy(nanowire->at_file,at_file);
    }

    void init_dop_file(char *dop_file)
    {
        strcpy(nanowire->dop_file,dop_file);
    }

    void init_perm_file(char *perm_file)
    {
        strcpy(nanowire->perm_file,perm_file);
    }

    void init_directory(char *directory)
    {
        //check if the directory exits?
        //if yes then files will be written...
        //else create the directory first.
	strcpy(parameter->directory,directory);
    }

    void init_transport_type(int transport_type)
    {
        parameter->transport_type = transport_type;
    }
    
    void init_command(char *command)
    {
        if(parameter->no_comm<MAX_COMM){
            parameter->command[parameter->no_comm] = new char[255];
            strcpy(parameter->command[parameter->no_comm],command);
            parameter->no_comm++;
        }
    }

    void init_lattype(const char *latname)
    {
      
        if(strcmp(latname,"ZB")==0 || strcmp(latname,"zb")==0 || strcmp(latname,"ZincBlende")==0 || strcmp(latname,"zincblende")==0 ){
	    parameter->lattype = 1;
	}
	if(strcmp(latname,"CU")==0 ||	strcmp(latname,"cu")==0 || strcmp(latname,"Cubic")==0 || strcmp(latname,"cubic")==0 ){
	    parameter->lattype = 2;
	}
	if(strcmp(latname,"WZ")==0 ||	strcmp(latname,"wz")==0 || strcmp(latname,"Wurtzite")==0 || strcmp(latname,"wurtzite")==0 ){
	    parameter->lattype = 3;
	}
	if(strcmp(latname,"hex")==0 ||	strcmp(latname,"HEX")==0 || strcmp(latname,"Graphene")==0 || strcmp(latname,"graphene")==0 ){
	    parameter->lattype = 4;
	}
	if(strcmp(latname,"rhombo")==0 || strcmp(latname,"rhombohedral")==0){
	    parameter->lattype = 5;
	}
	if(strcmp(latname,"cnt")==0 ||	strcmp(latname,"CNT")==0 || strcmp(latname,"Nanotube")==0 || strcmp(latname,"nanotube")==0 ){
	    parameter->lattype = 6;
	}
	if(strcmp(latname,"bilayer_graphene")==0 || strcmp(latname,"bilayer")==0){
	    parameter->lattype = 7;
	}
        if(strcmp(latname,"multilayer_graphene")==0 || strcmp(latname,"multilayer")==0){
	    parameter->lattype = 8;
	}

    }
/******************************************************************/

}
