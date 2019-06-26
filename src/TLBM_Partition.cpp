#include "TLBM_Partition.h"
#include <stdexcept>
#include <string>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <exception>


TLBM_Partition::TLBM_Partition(int r, int s, MPI_Comm c):
rank(r),size(s), comm(c),dataWriteNum(0)
{
//  printf("rank %d entering constructor \n",rank);
  tlbm_initialize();

}

TLBM_Partition::~TLBM_Partition(){
  delete myLattice;
  delete [] adjMatrix;
  delete [] fEven;
  delete [] fOdd;
  delete [] fEq;
  delete [] ndType;

  delete [] ux;
  delete [] uy;
  delete [] uz;
  delete [] rho;

  delete [] mpiInRequest;
  delete [] mpiOutRequest;
  delete [] mpiStatus;

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

  finalize_halo_data_arrays();

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

    // construct an array of MPI_Request(s) and Status(s)
    mpiInRequest = new MPI_Request[ngbSet.size()];
    mpiOutRequest = new MPI_Request[ngbSet.size()];
    mpiStatus = new MPI_Status[ngbSet.size()];

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
    			HDO_in[tgtP].insert_item(localToGlobal[nd],bbSpd[spd]);
    			// HDO_out now recorded the global node number (in order) and
    			// speed of all data to be transferred to each neighbor.
    			//
    			// HDO_in has the global node number and speeds of all data that
    			// will be received from neighboring partitions.
    			// (will be converted back to local node number when the HDO allocates and fills low-level arrays)
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

void TLBM_Partition::finalize_halo_data_arrays()
{
	//cycle through Halo Data Objects contained within each Halo Data Organizer (both in and out)
	// a) allocate arrays
	// b) load the nodeNums and spds arrays with the *local* node number and speed corresponding to
	//    each entry to be made in the buffer
	HDO_in.allocate_halo_arrays();
	HDO_out.allocate_halo_arrays();

	HDO_in.fill_arrays(globalToLocal);
	HDO_out.fill_arrays(globalToLocal);



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
	fEq = new real[numSpd*numLnodes];
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

	rc = MPI_File_open(comm,(char *)"ordering.b_dat",
			MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	if(rc)
	{
		throw "Error opening file to write node ordering";
	}

	//int offset_s = firstSlice*Nx*Ny*sizeof(int);
	//MPI_File_write_at(fh_snl,offset_s,snl+HALO*Nx*Ny,numEntries,MPI_INT,&mpi_s1);
	int offset = writeOffset*sizeof(int);
	MPI_File_write_at(fh,offset,localNdList.data(),numLnodes,MPI_INT,&status);

	MPI_File_close(&fh);

}

void TLBM_Partition::write_data()
{
	MPI_File fh;
	MPI_Status status;
	int rc;
	std::string fileName = "ux"+std::to_string(dataWriteNum)+".b_dat";
	int offset = writeOffset*sizeof(real);

	rc = MPI_File_open(comm,fileName.c_str(),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,ux,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write ux";
	}

	fileName = "uy"+std::to_string(dataWriteNum)+".b_dat";
	rc = MPI_File_open(comm,fileName.c_str(),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,uy,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write uy";
	}

	fileName = "uz"+std::to_string(dataWriteNum)+".b_dat";
	rc = MPI_File_open(comm,fileName.c_str(),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,uz,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write uz";
	}

	fileName = "density"+std::to_string(dataWriteNum)+".b_dat";
	rc = MPI_File_open(comm,fileName.c_str(),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,rho,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write rho";
	}
	++dataWriteNum;


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

	for(auto const & nd : nodeList)
	{
		myLattice->compute_macroscopic_data(ux,uy,uz,rho,fIn,nd);
		if (ndType[nd] == 1) // 1 is a solid node
		{
			ux[nd] = 0; uy[nd] = 0; uz[nd] = 0; // set macroscopic speed to zero
			myLattice->bounce_back(fOut,fIn,nd);
			continue; // skip to next iteration of the for-loop
		}

		// other than solid nodes, all nodes will need to compute equilibrium.
		myLattice->compute_equilibrium(fEq,ux,uy,uz,rho,nd);

		// node type 2 = inlet velocity node
		// node type 3 = outlet pressure node
		// node type 5 = specified u_z node
		if (ndType[nd] == 2)
		{
			myLattice->set_inlet_bc_macro(fIn,ux, uy, uz,rho,
					thisProblem.uLBM,nd);
			myLattice->set_inlet_bc_micro(const_cast<real *>(fIn),fEq,nd);
		}
		if (ndType[nd] == 3)
		{
			myLattice->set_outlet_bc_macro(fIn,ux,rho,thisProblem.rhoLBM,nd);
			myLattice->set_outlet_bc_micro(const_cast<real *>(fIn),fEq,nd);
		}
		if (ndType[nd] == 5)
		{
			myLattice->set_uz_bc(const_cast<real *>(fIn),ux,uy,uz,rho,thisProblem.uLBM,nd);
		}

		real omega = thisProblem.omega;
		if (thisProblem.cs > 0)
		{
			// create data structures needed for each individual lattice point
			real S[9] = {0,0,0,0,0,0,0,0,0};
			myLattice->compute_strain_tensor(S,fIn, fEq, nd);
			myLattice->apply_turbulence_model(omega,S,thisProblem.cs);
		}

		// pick between relaxation methodologies
		switch(thisProblem.dynamics)
		{
		case 1:
			myLattice->relax(fOut,fIn,fEq,omega,nd); break;

		case 2: // for now, do not do this.
			real piFlat[9] = {0,0,0,0,0,0,0,0,0};
			myLattice->compute_piflat(piFlat,fIn,fEq,nd);

			break;
//		case 3:

		}

		stream_node_data(fOut, fIn, nd);


	}
}

void TLBM_Partition::stream_node_data(real * fOut, const real * fIn, const int nd)
{
	int numSpd = myLattice->get_numSpd();
	for(int spd = 0; spd<numSpd;++spd)
	{
		int idx = getIDx(numSpd,nd,spd);
		int t_idx = adjMatrix[idx];
		fOut[t_idx] = fIn[idx];
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

int TLBM_Partition::get_num_global_nodes()
{
	int Nx = thisProblem.nx;
	int Ny = thisProblem.ny;
	int Nz = thisProblem.nz;

	return Nx*Ny*Nz;
}

real TLBM_Partition::get_data_member(const real * f, const int nd, const int spd)
{
	int numSpd = myLattice->get_numSpd();
	return f[getIDx(numSpd,nd,spd)];
}

void TLBM_Partition::set_data_member(real * f,const real val, const int nd, const int spd)
{
	int numSpd = myLattice->get_numSpd();
	f[getIDx(numSpd,nd,spd)] = val;
}

void TLBM_Partition::extract_halo_data(real * fOut)
{
	int numSpd = myLattice->get_numSpd();
	HDO_out.extract_halo_data(fOut,numSpd);

}

void TLBM_Partition::insert_halo_data(real * fOut)
{
	//printf("Rank %d inserting halo data\n",rank);
	int numSpd = myLattice->get_numSpd();
	HDO_in.insert_halo_data(fOut,numSpd);

}

void TLBM_Partition::initiate_data_exchange()
{

	int ngbIndex = 0;
	int count;
	real * in_buff;
	real * out_buff;


	for(auto & ngbIt : ngbSet)
	{
		out_buff = HDO_out[ngbIt].get_buffer();
		in_buff = HDO_in[ngbIt].get_buffer();
		count = HDO_out[ngbIt].get_num_items();
		MPI_Isend(static_cast<void *>(out_buff),count,MPI_DTYPE,ngbIt,ngbIndex,comm,mpiOutRequest+ngbIndex);
		MPI_Irecv(static_cast<void*>(in_buff),count,MPI_DTYPE,ngbIt,MPI_ANY_TAG,comm,mpiInRequest+ngbIndex);
		++ngbIndex;

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
	extract_halo_data(fOut);

	// initiate MPI Isend/Irecv
	initiate_data_exchange();

	// process interior nodes
	process_node_list(fOut,fIn,interiorNdList);


		// ensure MPI comms are complete
	MPI_Waitall(ngbSet.size(),mpiInRequest,mpiStatus);


	// distribute incoming halo data
	insert_halo_data(fOut);

}

