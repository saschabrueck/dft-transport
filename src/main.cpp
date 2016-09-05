#include <string.h>
#include <mpi.h>
#include <assert.h>
#include <fstream>
#include <iostream>
#include "libcp2k.h"

#ifdef HAVE_OMEN_POISSON
#include "Types.H"
#include "InputParameter.H"
#include "Material.H"
#include "WireGenerator.H"
#include "FEMGrid.H"
#include "Poisson.H"

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
#endif

int main (int argc, char **argv)
{
   char *input_file;
   char *output_file;
   char *command_file;
#ifdef libcp2k
   force_env_t force_env;
   int natom, n_el, n_el_pos, n_el_force;
   double e_pot;
   double *pos, *force;
#endif

   std::cout.precision(12);

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

   MPI::Init(argc, argv);
   MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);

#ifdef libcp2k
   cp2k_init_without_mpi();
   cp2k_create_force_env_comm(&force_env, input_file, output_file, MPI_Comm_c2f(MPI::COMM_WORLD));
   cp2k_transport_set_callback(force_env, &c_scf_method);
   cp2k_get_natom(force_env, &natom);

   n_el = natom*3;
   n_el_pos = natom*3;
   n_el_force = natom*3;
   pos = new double[n_el_pos];
   force = new double[n_el_force];
   std::fill_n(pos, n_el_pos, 0);
   std::fill_n(force, n_el_force, 0);
#endif

#ifdef HAVE_OMEN_POISSON
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
      FEM->execute_task(Wire,nanowire,1,1,MPI::COMM_SELF,MPI::COMM_WORLD);
//      FEM->execute_task(Wire,nanowire);
      OMEN_Poisson_Solver = new Poisson();
      OMEN_Poisson_Solver->init(Wire,nanowire,FEM,1,1,MPI::COMM_WORLD);
//      OMEN_Poisson_Solver->init(nanowire,FEM);
   }
#endif

#ifdef libcp2k
   cp2k_calc_energy_force(force_env);
   cp2k_get_potential_energy(force_env, &e_pot);
   cp2k_get_forces(force_env, force, n_el_force);
   cp2k_get_positions(force_env, pos, n_el);
   cp2k_destroy_force_env(force_env);
   cp2k_finalize_without_mpi();
#endif

#ifdef HAVE_OMEN_POISSON
   if (do_omen_poisson) {
      delete OMEN_Poisson_Solver;
      delete FEM;
      delete Wire;
      delete_parameters();
   }
#endif

#ifdef libcp2k   
   delete [] pos;
   delete [] force;
#endif
   delete [] input_file;
   delete [] output_file;
   delete [] command_file;

   MPI::Finalize();
   return 0;
}
