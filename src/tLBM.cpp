#include <iostream>
#include <cstdio>
#include <mpi.h>
#include "Problem.h"
#include "TLBM_Partition.h"

int main(int argc, char* argv[]){
  
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (rank == 0)
  {
	  printf("Commencing test with %d processes\n",size);
  }

  TLBM_Partition myPart(rank,size);

//  printf("Rank %d, cut-size: %d \n",rank,myPart.get_cut_size());

  if (rank == 0)
  {
	  printf("Test complete.");
  }

  MPI_Finalize();
  return 0;
}
