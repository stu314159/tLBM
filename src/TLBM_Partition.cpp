#include "TLBM_Partition.h"
#include <stdexcept>
#include <string>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <exception>


TLBM_Partition::TLBM_Partition(int r, int s, MPI_Comm c):
rank(r),size(s), comm(c)
{
//  printf("rank %d entering constructor \n",rank);
  tlbm_initialize();

}

TLBM_Partition::~TLBM_Partition(){
  delete myLattice;
  delete [] adjMatrix;
  delete [] fEven;
  delete [] fOdd;
  delete [] ndType;

  delete [] ux;
  delete [] uy;
  delete [] uz;
  delete [] rho;

}

int TLBM_Partition::get_num_ts()
{
	return thisProblem.numTs;
}

int TLBM_Partition::get_ts_rep_freq()
{
	return thisProblem.tsRepFreq;
}

int TLBM_Partition::get_plot_freq()
{
	return thisProblem.plotFreq;
}

int TLBM_Partition::tlbm_initialize(){

  thisProblem.load_input();
  // construct appropriate lattice type
  if (thisProblem.latticeType == std::string("D3Q15"))
  {
	 myLattice = new D3Q15LatticeStructure<real>;
  }else if(thisProblem.latticeType == std::string("D3Q19"))
  {
	  myLattice = new D3Q19LatticeStructure<real>;
  }else if (thisProblem.latticeType == std::string("D3Q27"))
  {
	  myLattice = new D3Q27LatticeStructure<real>;
  } else {
	  throw std::invalid_argument("Invalid Lattice Structure!");
  }

  load_parts();

  create_adj_matrix();

  compute_halo_data();

  make_adj_matrix_local();

  allocate_arrays();

  load_ndType();

  initialize_data_arrays();

  write_node_ordering();

  make_interior_node_list();

  return 0;
}

void TLBM_Partition::load_ndType()
{
	std::ifstream ndtype("ndType.lbm");
	int nt;
	int gNdInd = 0;
	int localNdInd;

	while (ndtype >> nt){
		if (partsG[gNdInd] == rank)// if gNdInd is a local node
		{
			localNdInd = globalToLocal.at(gNdInd); // get the local node index
			ndType[localNdInd] = nt; // assign node type to the ndType array
		}

		++gNdInd;
	}
	ndtype.close();
}

void TLBM_Partition::load_parts(){
 // read parts.lbm and obtain information about lattice points in my partition as well
// as my neighbors.

   std::ifstream parts("parts.lbm");
   int p;
   int gNdInd = 0; // global node index of current lattice point
   int localNdInd = 0; // local node index counter
   partSizes = std::vector<int>(size,0);
   int nNodes = thisProblem.nx*thisProblem.ny*thisProblem.nz;
   partsG = std::vector<int> (nNodes,0);
   std::pair<std::map<int,int>::iterator,bool> ret;

   while (parts >> p){
	   partSizes[p]+=1; // increment the # LPs in partition p
	   if (p == rank)
	   {
	     localNdList.push_back(gNdInd); // add to the local node list
	     // register nodes in the global to local and local to global maps
	     ret = localToGlobal.insert( std::pair<int,int>(localNdInd,gNdInd) );
	     if (ret.second == false){
	    	 throw "local to global key already existed!";
	     }
	     ret = globalToLocal.insert( std::pair<int,int>(gNdInd,localNdInd));
	     if (ret.second == false){
	    	 throw "global to local key already existed!";
	     }
	     localNdInd += 1; // increment the local node index
	   }
	   partsG[gNdInd]=p; // load partition number into partsG vector
	   gNdInd+=1;
   }
   parts.close(); // needed?

   // compute write offsets
    writeOffset = 0;
    for (auto i = partSizes.begin(); i < partSizes.begin()+rank;++i){
    	writeOffset += *i;
    }

    numLnodes = partSizes[rank]; // just to make this easier

}

void TLBM_Partition::create_adj_matrix(){
  unsigned int numSpd = myLattice->get_numSpd();
//  printf("Rank %d, numSpd: %d\n",rank,numSpd);
  adjMatrix = new int[numLnodes*numSpd];
  const int * ex = myLattice->get_ex();
  const int * ey = myLattice->get_ey();
  const int * ez = myLattice->get_ez();
  int tgt, node;
  for (int nd = 0; nd < numLnodes;nd++){
	  node = localNdList[nd];// note this is a global node number of the local node
	  for(unsigned spd = 0; spd < numSpd; spd++){
		  tgt = get_tgt_index(node,ex[spd],ey[spd],ez[spd]);
//		  printf("Rank %d, node %d, spd %d, tgt node: %d \n",
//				  rank,node,spd,tgt);
		  adjMatrix[getIDx(numSpd,nd,spd)] = tgt;
	  }
  }
}

void TLBM_Partition::compute_halo_data()
{
	int numSpd = myLattice->get_numSpd();
	int tgtNd, tgtP;
	for (int nd = 0; nd < numLnodes; nd++)
	{
		// iterate through the adjacency list for all of my local nodes
		for(auto spd = 0; spd < numSpd; spd++)
		{
			tgtNd = adjMatrix[getIDx(numSpd,nd,spd)]; // get global node number of tgt node
			tgtP = partsG[tgtNd];// get partition number of tgt node
			// if tgtP not equal to rank, then nd is in the set of boundary nodes
			if (tgtP != rank)
			{
				boundaryNdList.insert(nd); // local node number of boundary nodes
				// Create a Halo Data Object for in/out comms from tgtP
				// within the nodes Halo Data Organizer
				ngbSet.insert(tgtP);

			}

		}

	}
//	printf("Rank %u has %lu nds on bnl \n",rank,boundaryNdList.size());
//	printf("Rank %u has %lu neighbors \n",rank,ngbSet.size());

	// add set of neighbors to the HaloDataOrganizer for data in and out.
    for(const auto & ngbIt : ngbSet)
    {
    	HDO_out.add_neighbor(ngbIt);
    	HDO_in.add_neighbor(ngbIt);
    }

    // iterate over boundary nodes, traverse the boundary node
    // adjacency list, and add nodes from the halo into the HDOs.

    for (const auto & bnlIt : boundaryNdList)
    {

    	int nd = bnlIt; // bnl is already the local node number of boundary nodes
    	for(int spd = 0; spd < numSpd; ++spd)
    	{
    		const int * bbSpd = myLattice->get_bbSpd();
    		tgtNd = adjMatrix[getIDx(numSpd,nd,spd)]; // get global node number of tgt node
    		tgtP = partsG[tgtNd];
    		if ( tgtP != rank)
    		{
    			HDO_out[tgtP].insert_item(tgtNd,spd);
    			HDO_in[tgtP].insert_item(nd,bbSpd[spd]);
    			// HDO_out now recorded the global node number (in order) and
    			// speed of all data to be transferred to each neighbor.
    			//
    			// HDO_in has the local node number and speeds of all data that
    			// will be received from neighboring partitions.
    			//
    			// above is okay since local node numbers are generated in order
    			// of increasing global node number.
    		}

    	}
    }

    // compute the total number of halo nodes.
    numHaloNodes = HDO_out.get_num_halo_nodes();
    //printf("Rank %d, num halo nodes: %d \n",rank,numHaloNodes);
    totalNodes  = numLnodes + numHaloNodes;

    // make a list of the halo nodes (by global node number)
    haloNodes = HDO_out.get_halo_nodes(); // set of global node numbers for halo nodes

   // generate local nodes for the halo nodes and add to the local2global node map.
    int lNd = numLnodes; // initialize to the next local node number
    for (const auto & hnIt : haloNodes)
    {
    	localNdList.push_back(hnIt);
    	globalToLocal[hnIt] = lNd;
    	localToGlobal[lNd] = hnIt;
    	++lNd;
    }

}

void TLBM_Partition::make_adj_matrix_local()
{
	// now that all entries in the (global) adjacency matrix
	// have a local node identity, the adjacency list can be
	// converted to local node numbers
	int numSpd = myLattice->get_numSpd();
	int tgtNd;
	for (int nd = 0; nd < numLnodes; nd++)
	{
		// iterate through the adjacency list for all of my local nodes
		for(auto spd = 0; spd < numSpd; spd++)
		{
		  tgtNd = adjMatrix[getIDx(numSpd,nd,spd)]; // get global node number of tgt node
		  adjMatrix[getIDx(numSpd,nd,spd)] = globalToLocal.at(tgtNd);
		  // use map::at to generate an exception if tgtNd is not in the map
		}
	}
}

void TLBM_Partition::allocate_arrays()
{
	int numSpd = myLattice->get_numSpd();
	fEven = new real[numSpd*numLnodes];
	fOdd = new real[numSpd*numLnodes];
	ndType = new int[numLnodes];

	ux = new real[numLnodes];
	uy = new real[numLnodes];
	uz = new real[numLnodes];
	rho = new real[numLnodes];

}

void TLBM_Partition::initialize_data_arrays()
{
	int numSpd = myLattice->get_numSpd();
	const real * w = myLattice->get_w();
	real rho = thisProblem.rhoLBM;
	for(auto nd = 0; nd<numLnodes; ++nd)
	{
		for(auto spd = 0; spd < numSpd; ++spd)
		{
			fEven[getIDx(numSpd,nd,spd)] = w[spd]*rho;
			fOdd[getIDx(numSpd,nd,spd)] = w[spd]*rho;
		}
	}

}

void TLBM_Partition::write_node_ordering()
{
	MPI_File fh;
	MPI_Status status;
	int rc;

	rc = MPI_File_open(comm,"ordering.b_dat",
			MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);

	//int offset_s = firstSlice*Nx*Ny*sizeof(int);
	//MPI_File_write_at(fh_snl,offset_s,snl+HALO*Nx*Ny,numEntries,MPI_INT,&mpi_s1);
	int offset = writeOffset*sizeof(int);
	MPI_File_write_at(fh,offset,localNdList.data(),numLnodes,MPI_INT,&status);

	MPI_File_close(&fh);

}

int TLBM_Partition::get_cut_size()
{

	return HDO_out.get_cut_size();
}

int TLBM_Partition::get_tgt_index(int gInd, int ex, int ey, int ez){

	LatticeIndex tgt; //ret.X = 0; ret.Y = 0; ret.Z = 0;

	LatticeIndex src =  get_xyz_index(gInd);

	tgt.X = (src.X + ex)%thisProblem.nx;
	tgt.Y = (src.Y + ey)%thisProblem.ny;
	tgt.Z = (src.Z + ez)%thisProblem.nz;

	// get correct "wrap-around" behavior
	if (tgt.X < 0)
		tgt.X += thisProblem.nx;

	if (tgt.Y < 0)
		tgt.Y += thisProblem.ny;

	if (tgt.Z < 0)
		tgt.Z += thisProblem.nz;

	return get_gInd(tgt);
}



LatticeIndex TLBM_Partition::get_xyz_index(int gInd){
//	z = g_nd/(self.Nx*self.Ny)
//	y = (g_nd - z*self.Nx*self.Ny)/self.Nx
//	x = (g_nd - z*self.Nx*self.Ny - y*self.Nx)
//	return (x,y,z)
	LatticeIndex rVal;
	rVal.Z = gInd/(thisProblem.nx*thisProblem.ny);
	rVal.Y = (gInd - rVal.Z*thisProblem.nx*thisProblem.ny)/thisProblem.nx;
	rVal.X = (gInd - rVal.Z*thisProblem.nx*thisProblem.ny - rVal.Y*thisProblem.nx);

	return rVal;

}

int TLBM_Partition::get_gInd(int x, int y, int z){
//	return x+y*self.Nx + z*self.Nx*self.Ny
	return x + y*thisProblem.nx + z*thisProblem.nx*thisProblem.ny;
}

int TLBM_Partition::get_gInd(LatticeIndex myXYZ){
	return myXYZ.X + myXYZ.Y*thisProblem.nx + myXYZ.Z*thisProblem.nx*thisProblem.ny;
}
void TLBM_Partition::process_node_list(real * fOut, const real * fIn,
		const std::set<int>& nodeList)
{
	// create scratch data arrays
	const int numSpd = myLattice->get_numSpd();
	for(auto const & nd : nodeList)
	{
		myLattice->compute_macroscopic_data(ux,uy,uz,rho,fIn,nd);

	}
}


void TLBM_Partition::make_interior_node_list()
{
	for(int nd = 0; nd < numLnodes; ++nd)
	{
		if (!boundaryNdList.count(nd))
		{ // not in the boundary node list
			interiorNdList.insert(nd);
		}
	}

}
void TLBM_Partition::take_LBM_time_step(bool isEven)
{
	// set fIn and fOut
	if (isEven)
	{
		fIn = fEven;
		fOut = fOdd;
	}else{
		fIn = fOdd;
		fOut = fEven;
	}

	// process boundary nodes
	process_node_list(fOut,fIn,boundaryNdList);

	// extract halo data

	// initiate MPI Isend/Irecv

	// process interior nodes

	// ensure MPI comms are complete

	// distribute incoming halo data

}

