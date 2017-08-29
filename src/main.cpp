/*
Copyright (c) 2017 ETH Zurich
Sascha Brueck, Mauro Calderara, Mohammad Hossein Bani-Hashemian, and Mathieu Luisier

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <mpi.h>
#include <vector>
#include <assert.h>
#include "libcp2k.h"

#ifdef HAVE_OMEN_POISSON
#include "Types.H"
#include "InputParameter.H"
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
      Wire = new WireGenerator(parameter->lattype,0);
      Wire->execute_simple(nanowire,MPI_COMM_WORLD);
      FEM = new FEMGrid();
      int worldsize;
      MPI_Comm_size(MPI_COMM_WORLD,&worldsize);
      FEM->execute_task(Wire,nanowire,worldsize,1,MPI_COMM_WORLD,MPI_COMM_WORLD);
      OMEN_Poisson_Solver = new Poisson();
      OMEN_Poisson_Solver->init(Wire,nanowire,FEM,worldsize,1,MPI_COMM_WORLD);
   }
#endif

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
