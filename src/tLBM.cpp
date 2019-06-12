#include <iostream>
#include <cstdio> // sometimes C-style printf works better while using MPI
#include <mpi.h>
#include "Problem.h"
#include "TLBM_Partition.h"

int main(int argc, char* argv[]){
  
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  printf("Constructing a partition object for rank %d of %d\n",rank,size);
  TLBM_Partition myPart(rank,size);

  printf("Rank %d complete!\n",rank);

  MPI_Finalize();
  return 0;
}
