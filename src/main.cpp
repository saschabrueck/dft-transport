#include <string.h>
#include <mpi.h>
#include <assert.h>
#include "c_scf.H"

#define OUTPUT_FILE "./output.out"

extern "C"
{
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

   void cp_c_ext_scf_set_ptr(int *f_env_id, ptr2function, int *ierr);
}

int main (int argc, char **argv)
{
   int error, f_env_id, natom, calc_force, n_el, n_el_pos, n_el_force;
   int fcomm, init_mpi, finalize_mpi;
   double e_pot;
   char *input_file;
   char *output_file;
   double *pos, *force;
//   int myid, nproc;
   MPI::Intercomm mpi_comm;

   std::string inpfile_path = argv[1];
//   std::string outfile_path = argv[2];
   std::string outfile_path = OUTPUT_FILE;
   input_file=new char[inpfile_path.size()+1];
   input_file[inpfile_path.size()]=0;
   memcpy(input_file,inpfile_path.c_str(),inpfile_path.size());
   output_file=new char[outfile_path.size()+1];
   output_file[outfile_path.size()]=0;
   memcpy(output_file,outfile_path.c_str(),outfile_path.size());

   init_mpi = 0;
   finalize_mpi = 0;
   calc_force = 1;

   MPI::Init(argc, argv);
   MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
//   myid = MPI::COMM_WORLD.Get_rank();
//   nproc = MPI::COMM_WORLD.Get_size();
   mpi_comm = MPI::COMM_WORLD;
   fcomm = MPI_Comm_c2f(mpi_comm);
 
   cp_c_init_cp2k(&init_mpi, &error);
   cp_c_create_fenv_comm(&f_env_id, input_file, output_file, &fcomm, &error);
   cp_c_ext_scf_set_ptr(&f_env_id, &c_scf_method, &error);
   cp_c_get_natom(&f_env_id, &natom, &error);

   n_el = natom*3;
   n_el_pos = natom*3;
   n_el_force = natom*3;
   pos = new double[n_el_pos];
   force = new double[n_el_force];
   std::fill_n(pos, n_el_pos, 0);
   std::fill_n(force, n_el_force, 0);
   cp_c_calc_energy_force(&f_env_id, &calc_force, &error);
   cp_c_get_energy(&f_env_id, &e_pot, &error);
   cp_c_get_force(&f_env_id, force, &n_el_force, &error);
   cp_c_get_pos(&f_env_id, pos, &n_el, &error);
   cp_c_destroy_fenv(&f_env_id, &error);
   cp_c_finalize_cp2k(&finalize_mpi, &error);
   
   delete [] pos;
   delete [] force;
   delete [] input_file;
   delete [] output_file;

   MPI::Finalize();
   return 0;
}
