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

  TLBM_Partition myPart(rank,size,MPI_COMM_WORLD);

//  printf("Rank %d, cut-size: %d \n",rank,myPart.get_cut_size());

  // note: tLBM needs to ask a partition how many time steps there are and
  // when to write, etc...
  int numTs, tsRepFreq, plotFreq;

  numTs = myPart.get_num_ts();
  tsRepFreq = myPart.get_ts_rep_freq();
  plotFreq = myPart.get_plot_freq();

  double timeStart, timeEnd, execTime;

  timeStart = MPI_Wtime();
  for(int ts = 0; ts < numTs; ++ts)
  {
	  if((rank == 0) & ((ts+1)%tsRepFreq==0))
	  {
		  // say something comforting
		  printf("Executing time step: %d \n",ts+1);
	  }

      myPart.take_LBM_time_step(ts%2);

	  if((ts+1)%plotFreq == 0)
	  {
		  // plot the data
	  }

  }
  timeEnd = MPI_Wtime();
  execTime = timeEnd - timeStart;

  if (rank == 0)
  {
	  printf("Test complete.\n");
	  printf("Elapsed time: %g seconds \n",execTime);
  }

  MPI_Finalize();
  return 0;
}
