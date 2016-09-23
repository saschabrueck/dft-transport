#include <mpi.h>
#include <vector>
#include <assert.h>
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
   MPI_Init(&argc,&argv);

#ifdef HAVE_OMEN_POISSON
   if (argc>2) {
      yyin = fopen(argv[2],"r");
      init_parameters();
      yyrestart(yyin);
      yyparse();
      fclose(yyin);
      Material* material = new Material(parameter->mat_name,parameter->table_file,nanowire->read_hamiltonian,parameter->mat_binary_x,parameter->strain_model,parameter->Temp);
      material->initialize(nanowire->sc_dist_dep);
      Wire = new WireGenerator(parameter->lattype,0);
      Wire->execute_task(nanowire,material,0,MPI_COMM_WORLD);
//      WireGenerator* Wire = new WireGenerator();
//      Wire->execute_task(nanowire);
      FEM = new FEMGrid();
      FEM->execute_task(Wire,nanowire,1,1,MPI_COMM_SELF,MPI_COMM_WORLD);
//      FEM->execute_task(Wire,nanowire);
      OMEN_Poisson_Solver = new Poisson();
      OMEN_Poisson_Solver->init(Wire,nanowire,FEM,1,1,MPI_COMM_WORLD);
//      OMEN_Poisson_Solver->init(nanowire,FEM);
   }
#endif

#ifdef libcp2k
   cp2k_init_without_mpi();
   force_env_t force_env;
   cp2k_create_force_env_comm(&force_env, argv[1], "__STD_OUT__", MPI_Comm_c2f(MPI_COMM_WORLD));
   cp2k_transport_set_callback(force_env, &c_scf_method);
   int natom;
   cp2k_get_natom(force_env, &natom);
   cp2k_calc_energy_force(force_env);
   double e_pot;
   cp2k_get_potential_energy(force_env, &e_pot);
   std::vector<double> force(natom*3,0.0);
   cp2k_get_forces(force_env, &force[0], force.size());
   std::vector<double> pos(natom*3,0.0);
   cp2k_get_positions(force_env, &pos[0], pos.size());
   cp2k_destroy_force_env(force_env);
   cp2k_finalize_without_mpi();
#endif

#ifdef HAVE_OMEN_POISSON
   if (argc>2) {
      delete OMEN_Poisson_Solver;
      delete FEM;
      delete Wire;
      delete_parameters();
   }
#endif

   MPI_Finalize();
   return 0;
}
