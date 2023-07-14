#include <iostream>
#include <cstdio>
#include <mpi.h>
#include "Problem.h"
#include "TLBM_Partition.h"

int main(int argc, char* argv[]){

	// initialize MPI environment
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if (rank == 0)
	{
		printf("Commencing test with %d processes\n",size);
	}

	// initialize the partition.  TLBM_Partition reads input files to obtain problem data.
	TLBM_Partition myPart(rank,size,MPI_COMM_WORLD);


	int numTs, tsRepFreq, warmupTs, plotFreq, timeAvgFlag;


	numTs = myPart.get_num_ts();


	tsRepFreq = myPart.get_ts_rep_freq();


	warmupTs = myPart.get_warmupTs();

	plotFreq = myPart.get_plot_freq();

	timeAvgFlag = myPart.get_time_avg_flag();

	double timeStart, timeEnd, execTime;

	// carry out time-stepping process
	timeStart = MPI_Wtime();

        // commmence time stepping
	for(int ts = 0; ts < numTs; ++ts)
	{
		if((rank == 0) & ((ts+1)%tsRepFreq==0))
		{
			// periodically report progress
			printf("Executing time step: %d \n",ts+1);
		}

		// carry-out a single LBM time step (including data communication)
		myPart.take_LBM_time_step(ts%2);

		if(((ts+1)%plotFreq == 0) & ((ts+1) > warmupTs) )
		{
			// plot the data at specified intervals
			myPart.write_data();
		}

	}
	timeEnd = MPI_Wtime();
	execTime = timeEnd - timeStart;

	if (rank == 0)
	{
		// report basic performance parameters
		printf("Test complete.\n");
		printf("Elapsed time: %g seconds \n",execTime);
		double numNodes = static_cast<double>(myPart.get_num_global_nodes());
		double LPUs = numNodes*(static_cast<double>(numTs))/execTime;
		printf("Estimated LPU/s = %g \n",LPUs);
	}

	if (timeAvgFlag == 1)
	{
		if (rank == 0)
			printf("Writing time average data.\n");

		myPart.write_time_avg_data();
	}


	// clean up MPI environment
	MPI_Finalize();
	return 0;
}
