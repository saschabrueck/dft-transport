#include <string.h>
#include <mpi.h>
#include <assert.h>
#include "c_scf.H"

// OMEN INPUT
#include "Types.H"
#include "InputParameter.H"
#include "Material.H"
#include "WireGenerator.H"
#include "FEMGrid.H"
#include "Poisson.H"
// WRITE IN MATRIX
#include "CSR.H"

int semiselfconsistent(TCSR<double>*,TCSR<double>*,c_transport_type);

extern "C" {
    void yyrestart(FILE *);
    void yyparse();
    extern FILE *yyin;
}

PARAM *parameter;
WireStructure *nanowire;
ENERGY *En;
VOLTAGE *voltage;

WireGenerator* Wire;
FEMGrid *FEM;
Poisson *OMEN_Poisson_Solver;

// CP2K LIBRARY ROUTINES
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
   char *command_file;
   double *pos, *force;

   std::string inpfile_path = argv[1]+std::string(".inp");
   input_file=new char[inpfile_path.size()+1];
   input_file[inpfile_path.size()]=0;
   memcpy(input_file,inpfile_path.c_str(),inpfile_path.size());
   std::string outfile_path = argv[1]+std::string(".out");
   output_file=new char[outfile_path.size()+1];
   output_file[outfile_path.size()]=0;
   memcpy(output_file,outfile_path.c_str(),outfile_path.size());
   std::string commandfile_path = argv[1]+std::string(".cmd");
   command_file=new char[commandfile_path.size()+1];
   command_file[commandfile_path.size()]=0;
   memcpy(command_file,commandfile_path.c_str(),commandfile_path.size());

   MPI::Intercomm mpi_comm;

   init_mpi = 0;
   finalize_mpi = 0;
   calc_force = 1;

   MPI::Init(argc, argv);
   MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
   mpi_comm = MPI::COMM_WORLD;
   fcomm = MPI_Comm_c2f(mpi_comm);

   int run_cp2k=1;
   if (FILE *trpafile = fopen("TransportParams","r")) {
      run_cp2k=0;
      fclose(trpafile);
   }
   if (run_cp2k) {
      cp_c_init_cp2k(&init_mpi, &error);
      cp_c_create_fenv_comm(&f_env_id, input_file, output_file, &fcomm, &error);
      cp_c_ext_scf_set_ptr(&f_env_id, &c_scf_method, &error);
      cp_c_get_natom(&f_env_id, &natom, &error);
   }

   n_el = natom*3;
   n_el_pos = natom*3;
   n_el_force = natom*3;
   pos = new double[n_el_pos];
   force = new double[n_el_force];
   std::fill_n(pos, n_el_pos, 0);
   std::fill_n(force, n_el_force, 0);

   int do_omen_poisson=0;
   if ((yyin = fopen(command_file,"r")) != NULL) do_omen_poisson=1;
   if (do_omen_poisson) {
      init_parameters();
      yyrestart(yyin);
      yyparse();
      fclose(yyin);
      Material* material = new Material(parameter->mat_name,parameter->table_file,nanowire->read_hamiltonian,parameter->mat_binary_x,parameter->strain_model,parameter->Temp);
      material->initialize(nanowire->sc_dist_dep);
      Wire = new WireGenerator(parameter->lattype,0);
      Wire->execute_task(nanowire,material,0,MPI::COMM_WORLD);
//      WireGenerator* Wire = new WireGenerator();
//      Wire->execute_task(nanowire);
      FEM = new FEMGrid();
      MPI_Comm newcomm;
      int iam;
      MPI_Comm_rank(MPI_COMM_WORLD,&iam);
      MPI_Comm_split(MPI_COMM_WORLD,iam,iam,&newcomm);
      FEM->execute_task(Wire,nanowire,1,1,newcomm,MPI::COMM_WORLD);
//      FEM->execute_task(Wire,nanowire);
      OMEN_Poisson_Solver = new Poisson();
      OMEN_Poisson_Solver->init(Wire,nanowire,FEM,1,1,MPI::COMM_WORLD);
//      OMEN_Poisson_Solver->init(nanowire,FEM);
   }

   int w_size, w_rank;
   MPI_Comm_size(MPI_COMM_WORLD,&w_size);
   MPI_Comm_rank(MPI_COMM_WORLD,&w_rank);
   if (run_cp2k) {
      if (!w_rank) cout << "Starting CP2K" << endl;
      cp_c_calc_energy_force(&f_env_id, &calc_force, &error);
      cp_c_get_energy(&f_env_id, &e_pot, &error);
      cp_c_get_force(&f_env_id, force, &n_el_force, &error);
      cp_c_get_pos(&f_env_id, pos, &n_el, &error);
      cp_c_destroy_fenv(&f_env_id, &error);
      cp_c_finalize_cp2k(&finalize_mpi, &error);
   } else {
      c_transport_type transport_env_params;
      ifstream paraminfile("TransportParams");
      paraminfile >> transport_env_params.method;
      paraminfile >> transport_env_params.bandwidth;
      paraminfile >> transport_env_params.n_cells;
      paraminfile >> transport_env_params.n_occ;
      paraminfile >> transport_env_params.n_abscissae;
      paraminfile >> transport_env_params.n_kpoint;
      paraminfile >> transport_env_params.extra_int_param1;
      paraminfile >> transport_env_params.extra_int_param2;
      paraminfile >> transport_env_params.extra_int_param3;
      paraminfile >> transport_env_params.evoltfactor;
      paraminfile >> transport_env_params.colzero_threshold;
      paraminfile >> transport_env_params.eps_limit;
      paraminfile >> transport_env_params.eps_decay;
      paraminfile >> transport_env_params.eps_singularities;
      paraminfile >> transport_env_params.extra_param1;
      paraminfile >> transport_env_params.extra_param2;
      paraminfile >> transport_env_params.extra_param3;
      paraminfile.close();
      TCSR<double>* KohnSham = new TCSR<double>("KohnSham",w_size,w_rank);
      c_dscal(KohnSham->n_nonzeros,1.0/transport_env_params.evoltfactor,KohnSham->nnz,1);
      TCSR<double>* Overlap;
      if (FILE *ovlfile = fopen("Overlap","r")) {
         fclose(ovlfile);
         Overlap = new TCSR<double>("Overlap",w_size,w_rank);
      } else {
         TCSR<double>* OverlapRoot = NULL;
         if (!w_rank) {
            OverlapRoot = new TCSR<double>(KohnSham->size_tot,KohnSham->size_tot,KohnSham->findx);
            OverlapRoot->set_to_id();
         }
         Overlap = new TCSR<double>(KohnSham,OverlapRoot,0,MPI_COMM_WORLD);
         if (!w_rank) delete OverlapRoot;
      }
      semiselfconsistent(Overlap,KohnSham,transport_env_params);
      delete Overlap;
      delete KohnSham;
   }

   if (do_omen_poisson) {
      delete OMEN_Poisson_Solver;
      delete FEM;
      delete Wire;
      delete_parameters();
   }
   
   delete [] pos;
   delete [] force;
   delete [] input_file;
   delete [] output_file;
   delete [] command_file;

   MPI::Finalize();
   return 0;
}
