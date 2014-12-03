/**  \brief CP2K library routines 
 *   \author Mohammad Hossein Bani-Hashemian
 */

#ifndef __libcp2k_H
#define __libcp2k_H

#include "c_scf.H"

extern "C" {

   void cp_c_init_cp2k(int *init_mpi, int *ierr);
   void cp_c_finalize_cp2k(int *finalize_mpi, int *ierr);
   void cp_c_create_fenv(int *new_env_id, char *input_file_path, char *output_file_path, int *ierr);
   void cp_c_create_fenv_comm(int *new_env_id, char *input_file_path, char *output_file_path, int *mpi_comm, int *ierr);
   void cp_c_destroy_fenv(int *env_id, int *ierr);

   void cp_c_get_natom(int *env_id, int *natom, int *ierr);
   void cp_c_get_nparticle(int *env_id, int *nparticle, int *ierr);
   void cp_c_get_energy(int *env_id, double *e_pot, int *ierr);
   void cp_c_get_force(int *env_id, double force[], int *n_el, int *ierr);
   void cp_c_get_pos(int *env_id, double pos[], int *n_el, int *ierr);

   void cp_c_calc_energy_force(int *env_id, int *calc_force, int *ierr);

   void cp_c_run_input(char *input_file_path, char *output_file_path, int *ierr);
   void cp_c_run_input_comm(char *input_file_path, char *output_file_path, int *mpi_comm, int *ierr);

   void cp_c_ext_method_set_ptr(int *f_env_id, ptr2function, int *ierr);

}

#endif
