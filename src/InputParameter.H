/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __INPUT
#define __INPUT

#ifdef __cplusplus
extern "C" {
#endif

  void init_parameters();
  void default_parameters();
  void init_tb(int);
  void init_dsp3(double);
  void init_hpass(int);
  void init_a0_c0_u0(double,double,double);
  void init_first_atom(char*);
  void init_mbd(double);
  void init_n_of_modes(int);
  void init_nk(int);
  void init_last_first(int);
  void init_mat_name(char*);
  void init_bg_data(char*);
  void init_strain_model(int);
  void init_mat_binary_x(double,double,double);
  void init_eta_res(double);
  void init_eta(double);
  void init_plot_all_k(int);
  void init_NPROW(int);
  void init_NPCOL(int);
  void init_CPU_ppoint(int);
  void init_x(double,double,double);
  void init_y(double,double,double);
  void init_z(double,double,double);
  void init_Elimit(double);
  void init_Emin_tail(double);
  void init_EOffset(double);
  void init_dE_in(double,double);
  void init_dE_f(double,double);
  void init_dE_sep(double,double);
  void init_EExt(double);
  void init_NE(int);
  void init_reg_mesh(int);
  void init_strain(int);
  void init_Eps_xx(double);
  void init_Eps_yy(double);
  void init_Eps_zz(double);
  void init_Eps_xy(double);
  void init_Eps_xz(double);
  void init_Eps_yz(double);
  void init_shottky(int);      
  void init_schottky_type(int,int);      
  void init_schottky_barrier(double,double);      
  void init_schottky_virtual_cb(double,double);      
  void init_grid(int);
  void init_no_mat(int);
  void init_no_ch_mat(int);
  void init_no_ox_mat(int);
  void init_mat_type(int,char*);
  void init_mat_cs(int,char*);
  void init_mat(int,int,double,double,double);
  void init_mradius(int,double,double);
  void init_Eps_wire(double,double,double,double);
  void init_Eps_ox(double,double,double,double);
  void init_Xi_wire(double);
  void init_phi_m(double);
  void init_no_gate(int);
  void init_gradius(int,double,double);
  void init_gangle(int,double,double);
  void init_no_doping(int);
  void init_dradius(int,double,double);
  void init_gate(int,int,double,double,double);
  void init_doping(int,int,double,double,double);
  void init_gate_type(int,char*);
  void init_doping_type(int,char*);
  void init_doping_nd(int,double);
  void init_doping_na(int,double);
  void init_no_doping_domain(int);
  void init_dop_domain_coord(int,double,int);
  void init_dop_domain_slope(int,double,double,double);
  void init_dop_domain_conc(int,double);
  void init_temp(double);
  void init_command(char*);
  void init_Ls(double);
  void init_Lc(double);
  void init_Ld(double);
  void init_NDS(double);
  void init_NAS(double);
  void init_NDD(double);
  void init_NAD(double);
  void init_nvg(int,int);
  void init_vgmin(double,double,double,double);
  void init_vgmax(double,double,double,double);
  void init_nvs(int);
  void init_vsmin(double);
  void init_vsmax(double);
  void init_nvd(int);
  void init_vdmin(double);
  void init_vdmax(double);
  void init_ntemp(int);
  void init_tmin(double);
  void init_tmax(double);
  void init_directory(char*);
  void delete_parameters();
  void init_ndim(int);
  void init_nky(int);
  void init_nkz(int);
  void init_phiy_file(char*);
  void init_phiz_file(char*);
  void replicate_unit_cell(int);
  void read_atom_pos(int);
  void init_rot_sym(int);
  void init_nxfold(int);
  void init_nyfold(int);
  void init_nzfold(int);
  void init_bs_solver(char *);
  void init_update_bs_target(int);
  void init_bs_target(double);
  void init_NPCS(int);
  void init_spec_decomp(int);
  void init_CPU_per_kz_point(int);
  void init_CPU_per_sample(int);
  void init_CPU_per_vg_point(int);
  void init_CPU_per_vd_point(int);
  void init_CPU_per_temp_point(int);
  void init_CPU_per_wire(int);
  void init_CPU_per_bc(int);
  void init_poisson_solver(int);
  void init_poisson_criterion(double);
  void init_poisson_iteration(int);
  void init_max_proc_poisson(int);
  void init_poisson_inner_criterion(double);
  void init_poisson_inner_iteration(int);
  void init_charge_average(int);
  void init_charge_transfer(int);
  void init_zeta(double);
  void init_no_strain_domain(int);
  void init_stdomain_coord(int,double,int);
  void init_stdomain_eps_vec(int,double,double,double,\
			     double,double,double);
  void init_roughness(int);
  void init_roughness_seed(int);
  void init_roughness_type(char*);
  void init_roughness_rms(double);
  void init_roughness_lms(double);
  void init_alloy_disorder(int);
  void init_alloy_seed(int);
  void init_alloy_composition(double);
  void init_qmstart(double);
  void init_qmstop(double);
  void init_update_fermi(int);
  void init_fermi_level(double);
  void init_update_at(int);
  void init_at_file(char*);
  void init_update_dop(int);
  void init_dop_file(char*);
  void init_update_perm(int);
  void init_perm_file(char*);
  void init_no_ro_mat(int);
  void init_mat_id(int,int);
  void init_bulk_mat_id(int);
  void init_grradius(int,double);
  void init_no_ground(int);
  void init_ground(int,int,double,double,double);
  void init_ground_type(int,char*);
  void init_update_tetra(int);
  void init_tetra_file(char*);
  void init_update_fitness(int);
  void init_fitness_file(char*);
  void init_update_grid(int);
  void init_grid_file(char*);
  void init_update_energy(int);
  void init_energy_file(char*);
  void init_ph_energy_file(char*);
  void init_no_ph_energy(int);
  void init_ph_mode_file(char*);
  void init_no_ph_k(int);
  void init_ph_file_type(int);
  void init_sc_max_iteration(int);
  void init_sc_restart(int,int);
  void init_sc_id_crit(double);
  void init_sc_rh_crit(double);
  void init_sc_dist_dep(int);
  void init_sc_diag_def(int);
  void init_sc_k_coupling(int);
  void init_sc_e_ph_coupling(int);
  void init_sc_scale_fact(double);
  void init_sc_vbound(double);
  void init_sc_memory_fact(double);
  void init_incoherent_injection(int);
  void init_full_current(int);
  void init_robust_numerics(int);
  void init_convergence_booster(int);
  void init_restart(int,int,int,int);
  void init_vtot_file(char*);
  void init_vact_file(char*);
  void init_transport_type(int);
  void init_injection_type(int,int);
  void init_electron_weight(double);
  void init_open_system(int);
  void init_periodic_system(int);
  void init_read_hamiltonian(int);
  void init_lattype(const char*);
#ifdef __cplusplus
}
#endif

#endif

