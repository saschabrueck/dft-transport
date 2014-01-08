
#include <iostream>
#include <sstream>
#include "CSR.H"
#include "Utilities.H"
#include "Blas.H"

// MPI
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "Spike.H"
#include <omp.h>


int main(int argc, char** argv) {

  //MPI_Status mpistatus;
  int rank,world_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
 
  std::stringstream matrixfile; 
  matrixfile << "/Users/ottobibartiu/Documents/Repository/SPIKE/input/random/rank"<< rank << ".csr";
 //matrixfile << "/Users/ottobibartiu/Documents/Repository/SPIKE/prototype/m3000p2b300/rank"<< rank << ".csr";
  
  TCSR<CPX> *matrix = new TCSR<CPX>(matrixfile.str().c_str());
  
  CPX *RHS = NULL;
  // generate RHS
  if( rank == 0 || rank == world_size-1 ){
        RHS = new CPX[180];
        for (int i = 0; i < 60 ; ++i) {
          RHS[i] = {1,2};
        }
         for (int i = 60; i < 120 ; ++i) {
          RHS[i] = {3,4};
        }
        for (int i = 120; i < 180 ; ++i) {
          RHS[i] = {5,6};
        }
  }
  int rhs_col = 3;
  std::vector<int> partitions = {0,60,120,180,240,300};
//  std::vector<int> partitions = {0,1500, 3000};
  
  Spike<CPX> *spike = new Spike<CPX>(matrix,11 ,RHS, rhs_col, MPI_COMM_WORLD,
                                     partitions); 
  
  spike->solve_full();
 
  spike->print_solution_full();
  
  
  MPI_Finalize();
  return 0;
}
