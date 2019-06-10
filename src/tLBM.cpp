#include <iostream>
#include <mpi.h>
#include "Problem.h"
#include "TLBM_Partition.h"

int main(int argc, char* argv[]){
  
  int rank, size;
  MPI_Init(&ragc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);



  MPI_Finalize();
  return 0;
}
