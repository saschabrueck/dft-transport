#include <iostream>
#include <sstream>
#include "CSR.H"
#include "Utilities.H"
//#include "Blas.H"

// MPI
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mkl.h>
#include "intel_Spike.H"
#include <omp.h>


using namespace std;

int main(int argc, char** argv) {

  //MPI_Status mpistatus;
  int rank,world_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);


  //printf("Hello World from rank = %d\n", rank);
  /**
   * LOAD MATRIX BLOCKS
   * 
   */
  
  std::stringstream matrixfile; 
  //matrixfile << "/Users/ottobibartiu/Documents/Repository/SPIKE/input/random/rank"<< rank << ".csr";
  matrixfile << "../input/random/rank" << rank << ".csr";
  
  TCSR<double> *matrix = new TCSR<double>(matrixfile.str().c_str());
  
  double *RHS = NULL;
  // generate RHS
  if( rank == 0 || rank == world_size-1 ){
        RHS = new double[60*2];
        for (int i = 0; i < 120 ; ++i) {
          RHS[i] = 1.0;
        }
  }
  int rhs_col = 2;
  std::vector<int> partitions = {0, 60, 120, 180, 240, 300};
  
  intel_SPIKE *spike = new intel_SPIKE(matrix, RHS, rhs_col,
                                      MPI_COMM_WORLD, partitions, 0); 
  spike->solve_full();
 
  
  MPI_Finalize();

  return 0;
}
