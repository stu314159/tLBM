#include "TLBM_Partition.h"
#include <stdexcept>
#include <string>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <exception>

using namespace H5;
TLBM_Partition::TLBM_Partition(int r, int s, MPI_Comm c):
rank(r),size(s), comm(c),dataWriteNum(0)
{
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

  delete [] Fx;
  delete [] Fy;
  delete [] Fz;

  if (timeAvg == 1)
  {
	delete [] uAvg;
	delete [] vAvg;
	delete [] wAvg;
	delete [] rhoAvg;
  }

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

int TLBM_Partition::get_warmupTs()
{
	return thisProblem.warmupTs;
}

int TLBM_Partition::is_restart()
{
    return thisProblem.restartFlag;
}

int TLBM_Partition::get_time_avg_flag()
{
	return thisProblem.timeAvgFlag;
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



  make_force_calc_map();

  // swap HALO information to ensure consistency


  swap_halo_node_data();

  swap_halo_speed_data();

  // the HDO_in objects should now have the comp_ndNums and comp_spds arrays populated
  // with information passed from their neighboring partitions
  check_halo_data();




  return 0;
}

void TLBM_Partition::load_ndType()
{
	std::ifstream ndtype("ndType.lbm");
	int nt;
	int gNdInd = 0;
    int localNdInd;
	int nNodes = thisProblem.nx*thisProblem.ny*thisProblem.nz;
	std::vector<int> allNdType = std::vector<int> (nNodes,0);

	while (ndtype >> nt)
	{
		allNdType[gNdInd] = nt;
		++gNdInd;
	}
	ndtype.close();

	for (const auto & ndit : localNdList)// remember: localNdList has *global* node numbers for local nodes (including halo)
	{
		localNdInd = globalToLocal.at(ndit); // get the local node number
		ndType[localNdInd] = allNdType[ndit]; // assign node type to local ndType array
	}
}

void TLBM_Partition::make_force_calc_map()
{
	// iterate through all local nodes (not including the halo) and identify which nodes are on a fluid/surface boundary.
	// the local node number will be the key for the forceCalcMap and the set of all lattice speeds pointing to fluid
	// nodes is the value.

	int numSpd = myLattice->get_numSpd();
	int idx, ngbNd, ngbNd_type;

	for (int nd = 0; nd < numLnodes; nd++)
	{
		if ( ndType[nd] == 1 ) // if a solid node
		{ // iterate through neighbors.
			for ( int spd=0; spd < numSpd; spd++)
			{
			   	idx =  getIDx(numLnodes,nd,spd);
			   	ngbNd = adjMatrix[idx];
			   	ngbNd_type = ndType[ngbNd];

			   	if (ngbNd_type == 0) // if a neighbor node type is 0 (fluid), then add to the force calc map
			   	{
			   		forceCalcMap[nd].insert(spd);
			   	}
			}
		}
	}

//	printf("Rank %d has %lu surface nodes for the force calculation.\n ",rank,forceCalcMap.size());

}

void TLBM_Partition::calc_force()
{
	// initialize all of the force values to zero
	for(int nd = 0; nd < numLnodes; nd++)
	{
	  Fx[nd] = 0; Fy[nd] = 0; Fz[nd] = 0;
	}

	// get lattice information
//	int numSpd = myLattice->get_numSpd();
	const int * ex = myLattice->get_ex();
    const int * ey = myLattice->get_ey();
	const int * ez = myLattice->get_ez();
	const int * bb_spds = myLattice->get_bbSpd();

	// iterate through all of the keys of the forceCalcMap
	for( auto & fnd_pairs : forceCalcMap )
	{
		int fnd = fnd_pairs.first; // local node number of LP on fluid/surface interface
		for ( auto & spd : fnd_pairs.second ) // iterate through speeds of fnd that point to fluid nodes
		{
			int idx = getIDx(totalNodes,fnd,spd);
			int bb_spd = bb_spds[spd];
//			int ngbNd = adjMatrix[idx];
			real f1 = fOut[idx]; // get fOut in speed towards fluid neighbor
//			int idx_bb = getIDx(numSpd,ngbNd,bb_spd);
			int idx_bb = getIDx(totalNodes,fnd,bb_spd);
			real f2 = fOut[idx_bb]; // get fOut from neighbor node, in speed towards fnd

			Fx[fnd] += ex[bb_spd] * ( f1 + f2 );
			Fy[fnd] += ey[bb_spd] * ( f1 + f2 );
			Fz[fnd] += ez[bb_spd] * ( f1 + f2 );

		}
	}
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
	  node = localNdList[nd];// note this (node) is a global node number of the local node (nd)
	  for(unsigned spd = 0; spd < numSpd; spd++){
		  tgt = get_tgt_index(node,ex[spd],ey[spd],ez[spd]);
//		  printf("Rank %d, node %d, spd %d, tgt node: %d \n",
//				  rank,node,spd,tgt);
		  adjMatrix[getIDx(numLnodes,nd,spd)] = tgt;
	  }
  }
}

void TLBM_Partition::print_adjacency(const int nd){
	unsigned int numSpd = myLattice->get_numSpd();
	std::cout << "adjacency for nd " << nd << ": ";
	for(unsigned int s=0; s<numSpd; s++)
	{
		std::cout << adjMatrix[getIDx(numLnodes,nd,s)] << ", ";
	}
	std::cout << std::endl;

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
			tgtNd = adjMatrix[getIDx(numLnodes,nd,spd)]; // get global node number of tgt node
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

    const int * bbSpd = myLattice->get_bbSpd();
    for (const auto & bnlIt : boundaryNdList)
    {

    	int nd = bnlIt; // bnl is already the local node number of boundary nodes
    	for(int spd = 0; spd < numSpd; ++spd)
    	{

    		tgtNd = adjMatrix[getIDx(numLnodes,nd,spd)]; // get global node number of tgt node
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
    			// of increasing global node number. (both for local and halo nodes)
    		}

    	}
    }

    // compute the total number of halo nodes.
    numHaloNodes = HDO_out.get_num_halo_nodes();
    //printf("Rank %d, num halo nodes: %d \n",rank,numHaloNodes);
    totalNodes  = numLnodes + numHaloNodes;

    // make a list of the halo nodes (by global node number)
    haloNodes = HDO_out.get_halo_nodes(); // set of global node numbers for halo nodes (in global node # order)

   // generate local nodes for the halo nodes and add to the local2global node map.
    int lNd = numLnodes; // initialize to the next local node number

    std::pair<std::map<int,int>::iterator,bool> ret;
    for (const auto & hnIt : haloNodes)
    {
    	localNdList.push_back(hnIt);
    	ret = globalToLocal.insert(std::pair<int,int>(hnIt,lNd));
    	if (ret.second == false)
    	{
    	  	 throw "global to local key already existed!";
    	}

    	ret = localToGlobal.insert(std::pair<int,int>(lNd,hnIt));
    	if (ret.second == false)
    	{
    	  	 throw "local to global key already existed!";
    	}

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

void TLBM_Partition::check_halo_data()
{
	int numHaloErrors = HDO_in.check_ndNums_and_spds(globalToLocal);

	if(numHaloErrors > 0)
	{
		printf("Rank %d reports %d Halo Errors \n", rank, numHaloErrors);
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
		  tgtNd = adjMatrix[getIDx(numLnodes,nd,spd)]; // get global node number of tgt node
		  adjMatrix[getIDx(numLnodes,nd,spd)] = globalToLocal.at(tgtNd);
		  // use map::at to generate an exception if tgtNd is not in the map
		}
	}
}

void TLBM_Partition::allocate_arrays()
{
	int numSpd = myLattice->get_numSpd();
	fEven = new real[numSpd*totalNodes];//larger to store halo data (as stream target) as well.
	fOdd = new real[numSpd*totalNodes];// larger to store halo data (as stream target) as well
	fEq = new real[numSpd*numLnodes];
	ndType = new int[totalNodes];

	ux = new real[numLnodes];
	uy = new real[numLnodes];
	uz = new real[numLnodes];
	rho = new real[numLnodes];

	Fx = new real[numLnodes];
	Fy = new real[numLnodes];
	Fz = new real[numLnodes];

	if (timeAvg == 1)
	{
		uAvg = new real[numLnodes];
		vAvg = new real[numLnodes];
		wAvg = new real[numLnodes];
		rhoAvg = new real[numLnodes];
	}

}

void TLBM_Partition::initialize_data_arrays()
{
	int numSpd = myLattice->get_numSpd();
	const real * w = myLattice->get_w();
	real rho = thisProblem.rhoLBM;
	if (is_restart())
	{
		load_restart_data();
	}else
	{
		for(auto nd = 0; nd<totalNodes; ++nd)
		{
			for(auto spd = 0; spd < numSpd; ++spd)
			{
				fEven[getIDx(totalNodes,nd,spd)] = w[spd]*rho;
				fOdd[getIDx(totalNodes,nd,spd)] = w[spd]*rho;
			}
		}
	}

	if (timeAvg == 1)
	{
		for(auto nd = 0; nd<numLnodes; ++nd)
		{
			uAvg[nd] = 0;
			vAvg[nd] = 0;
			wAvg[nd] = 0;
			rhoAvg[nd] = 0;
		}
	}

}

void TLBM_Partition::load_restart_data()
{
	if (rank == 0)
	{
		printf("Loading Restart Data \n");
	}
	// declare arrays and allocate memory to hold restart data
	real * ux_r;
	real * uy_r;
	real * uz_r;
	real * rho_r;
	int t_nodes = (thisProblem.nx)*(thisProblem.ny)*(thisProblem.nz);
	ux_r = new real[t_nodes];
	uy_r = new real[t_nodes];
	uz_r = new real[t_nodes];
	rho_r = new real[t_nodes];

	// get data from the HDF5 formatted "restart.h5"
	const H5std_string FILE_NAME("restart.h5");
	H5File file(FILE_NAME,H5F_ACC_RDONLY);
	const H5std_string PRESSURE("density/rho");
	const H5std_string X_VELO("velocity/x");
	const H5std_string Y_VELO("velocity/y");
	const H5std_string Z_VELO("velocity/z");

	DataSet p_dataset = file.openDataSet(PRESSURE);
	DataSet x_dataset = file.openDataSet(X_VELO);
	DataSet y_dataset = file.openDataSet(Y_VELO);
	DataSet z_dataset = file.openDataSet(Z_VELO);

	// read data into the arrays
	p_dataset.read(rho_r,PredType::NATIVE_FLOAT);
	x_dataset.read(ux_r,PredType::NATIVE_FLOAT);
	y_dataset.read(uy_r,PredType::NATIVE_FLOAT);
	z_dataset.read(uz_r,PredType::NATIVE_FLOAT);

	// load data from arrays ux_r,uy_r,uz_r,rho_r into the local partition arrays
	for(auto nd = 0; nd<numLnodes; ++nd)
	{
		int gnd = localToGlobal[nd];
		ux[nd] = ux_r[gnd];
		uy[nd] = uy_r[gnd];
		uz[nd] = uz_r[gnd];
		rho[nd] = rho_r[gnd];
	}

	// using the macroscopic data arrays, initialize fEven and fOdd
	for(auto nd = 0; nd<numLnodes; ++nd)
	{
		myLattice->compute_equilibrium(fEven,ux,uy,uz,rho,nd, totalNodes);
		myLattice->compute_equilibrium(fOdd,ux,uy,uz,rho,nd, totalNodes);
	}

	// free memory allocated for restart data
	delete [] ux_r;
	delete [] uy_r;
	delete [] uz_r;
	delete [] rho_r;

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

void TLBM_Partition::write_time_avg_data()
{
	MPI_File fh;
	MPI_Status status;
	int rc;
	std::string filename = "uAvg.b_dat";
	int offset = writeOffset*sizeof(real);

	rc = MPI_File_open(comm,(char *)(filename.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,uAvg,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write uAvg";
	}

	filename = "vAvg.b_dat";
	rc = MPI_File_open(comm,(char *)(filename.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,vAvg,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write vAvg";
	}

	filename = "wAvg.b_dat";
	rc = MPI_File_open(comm,(char *)(filename.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,wAvg,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write wAvg";
	}

	filename = "rhoAvg.b_dat";
	rc = MPI_File_open(comm,(char *)(filename.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,rhoAvg,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write rhoAvg";
	}
}

void TLBM_Partition::write_data()
{
	MPI_File fh;
	MPI_Status status;
	int rc;
	std::string fileName = "ux"+std::to_string(dataWriteNum)+".b_dat";
	int offset = writeOffset*sizeof(real);

	rc = MPI_File_open(comm,(char *)(fileName.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,ux,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write ux";
	}

	fileName = "uy"+std::to_string(dataWriteNum)+".b_dat";
	rc = MPI_File_open(comm,(char *)(fileName.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,uy,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write uy";
	}

	fileName = "uz"+std::to_string(dataWriteNum)+".b_dat";
	rc = MPI_File_open(comm,(char *)(fileName.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,uz,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write uz";
	}

	fileName = "density"+std::to_string(dataWriteNum)+".b_dat";
	rc = MPI_File_open(comm,(char *)(fileName.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	MPI_File_write_at(fh,offset,rho,numLnodes,MPI_DTYPE,&status);
	MPI_File_close(&fh);
	if(rc)
	{
		throw "Error opening file to write rho";
	}

    calc_force();

    fileName = "Fx"+std::to_string(dataWriteNum)+".b_dat";
    rc = MPI_File_open(comm,(char *)(fileName.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
    MPI_File_write_at(fh,offset,Fx,numLnodes,MPI_DTYPE,&status);
    MPI_File_close(&fh);
    if(rc)
    {
    	throw "Error opening file to write Fx";
    }

    fileName = "Fy"+std::to_string(dataWriteNum)+".b_dat";
    rc = MPI_File_open(comm,(char *)(fileName.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
    MPI_File_write_at(fh,offset,Fy,numLnodes,MPI_DTYPE,&status);
    MPI_File_close(&fh);
    if(rc)
    {
    	throw "Error opening file to write Fy";
    }

    fileName = "Fz"+std::to_string(dataWriteNum)+".b_dat";
    rc = MPI_File_open(comm,(char *)(fileName.c_str()),MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
    MPI_File_write_at(fh,offset,Fz,numLnodes,MPI_DTYPE,&status);
    MPI_File_close(&fh);
    if(rc)
    {
    	throw "Error opening file to write Fz";
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
void TLBM_Partition::process_node_list(real * fOut, real * fIn,
		const std::set<int>& nodeList)
{

	for(auto const & nd : nodeList)
	{
		myLattice->compute_macroscopic_data(ux,uy,uz,rho,fIn,nd,totalNodes);
		if (ndType[nd] == 0)
		{
			myLattice->compute_equilibrium(fEq,ux,uy,uz,rho,nd,numLnodes);

		}
		if (ndType[nd] == 1) // 1 is a solid node
		{
			ux[nd] = 0; uy[nd] = 0; uz[nd] = 0; // set macroscopic speed to zero
			myLattice->compute_equilibrium(fEq,ux,uy,uz,rho,nd,numLnodes);
			myLattice->bounce_back(fIn,nd,totalNodes);
			//continue; // skip to next iteration of the for-loop
		}


		if (ndType[nd] != 1) {

			// node type 2 = inlet (west) velocity node
			// node type 3 = outlet (east) pressure node
			// node type 5 = specified u_z node
			if (ndType[nd] == 2) {
				myLattice->set_inletW_bc_macro(fIn, ux, uy, uz, rho,
						thisProblem.uLBM, nd,totalNodes);
				// other than solid nodes, all nodes will need to compute equilibrium.
				myLattice->compute_equilibrium(fEq, ux, uy, uz, rho, nd,numLnodes);
				myLattice->set_inletW_bc_micro(fIn, fEq, nd,totalNodes, numLnodes);
			}
			if (ndType[nd] == 3) {
				myLattice->set_outletE_bc_macro(fIn, uz, rho, thisProblem.rhoLBM,
						nd,totalNodes);
				// other than solid nodes, all nodes will need to compute equilibrium.
				myLattice->compute_equilibrium(fEq, ux, uy, uz, rho, nd, numLnodes);
				myLattice->set_outletE_bc_micro(fIn, fEq, nd,totalNodes, numLnodes);
			}
			if (ndType[nd] == 5) {
				myLattice->set_uz_bc(fIn, ux, uy, uz, rho, thisProblem.uLBM,
						nd,totalNodes);
				// other than solid nodes, all nodes will need to compute equilibrium.
				myLattice->compute_equilibrium(fEq, ux, uy, uz, rho, nd, numLnodes);
			}

			real omega = thisProblem.omega;
			if (thisProblem.cs > 0) {
				// create data structures needed for each individual lattice point
				real S[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
				myLattice->compute_strain_tensor(S, fIn, fEq, nd,totalNodes, numLnodes);
				omega = myLattice->apply_turbulence_model(omega, S,
						thisProblem.cs);
			}

			// pick between relaxation methodologies
			switch (thisProblem.dynamics) {
			case 1:
				myLattice->relax(fIn, fEq, omega, nd,totalNodes, numLnodes);
				break;

			case 2: // for now, do not do this.
				real piFlat[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
				myLattice->compute_piflat(piFlat, fIn, fEq, nd, totalNodes);

				break;
				//		case 3:

			}
		}
		stream_node_data(fOut, fIn, nd);

		
	}
}

void TLBM_Partition::update_time_avg()
{
	for (auto nd=0; nd<numLnodes; ++nd)
	{
		uAvg[nd]+=ux[nd];
		vAvg[nd]+=uy[nd];
		wAvg[nd]+=uz[nd];
		rhoAvg[nd]+=rho[nd];
	}
}

void TLBM_Partition::stream_node_data(real * fOut, const real * fIn, const int nd)
{
	int numSpd = myLattice->get_numSpd();
	for(int spd = 0; spd<numSpd;++spd)
	{
		int idx = getIDx(totalNodes,nd,spd);
		int t_nd = adjMatrix[getIDx(numLnodes,nd,spd)];
		int t_idx = getIDx(totalNodes,t_nd,spd);
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
//	int numSpd = myLattice->get_numSpd();
	return f[getIDx(totalNodes,nd,spd)];
}

void TLBM_Partition::set_data_member(real * f,const real val, const int nd, const int spd)
{
//	int numSpd = myLattice->get_numSpd();
	f[getIDx(totalNodes,nd,spd)] = val;
}

void TLBM_Partition::extract_halo_data(real * fOut)
{
//	int numSpd = myLattice->get_numSpd();
	HDO_out.extract_halo_data(fOut,totalNodes);

}

void TLBM_Partition::insert_halo_data(real * fOut)
{
	//printf("Rank %d inserting halo data\n",rank);
//	int numSpd = myLattice->get_numSpd();
	HDO_in.insert_halo_data(fOut,totalNodes);

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

void TLBM_Partition::swap_halo_node_data()
{
	int ngbIndex = 0;
	int count;

	int * in_buff;
	int * out_buff;

	for(auto & ngbIt : ngbSet)
	{
		out_buff = HDO_out[ngbIt].get_g_ndNums();
		in_buff = HDO_in[ngbIt].get_comp_ndNums_buffer();
		count = HDO_out[ngbIt].get_num_items();

		MPI_Isend(static_cast<void *>(out_buff),count,MPI_INT,ngbIt,ngbIndex,comm,mpiOutRequest+ngbIndex);
		MPI_Irecv(static_cast<void *>(in_buff),count,MPI_INT,ngbIt,MPI_ANY_TAG,comm,mpiInRequest+ngbIndex);
		++ngbIndex;

	}

	// ensure MPI comms are complete
	MPI_Waitall(ngbSet.size(),mpiInRequest,mpiStatus);

}

void TLBM_Partition::swap_halo_speed_data()
{
	int ngbIndex = 0;
	int count;

	int * in_buff;
	int * out_buff;

	for(auto & ngbIt : ngbSet)
	{
		out_buff = HDO_out[ngbIt].get_spds();
		in_buff = HDO_in[ngbIt].get_comp_spds_buffer();
		count = HDO_out[ngbIt].get_num_items();

		MPI_Isend(static_cast<void *>(out_buff),count,MPI_INT,ngbIt,ngbIndex,comm,mpiOutRequest+ngbIndex);
		MPI_Irecv(static_cast<void *>(in_buff),count,MPI_INT,ngbIt,MPI_ANY_TAG,comm,mpiInRequest+ngbIndex);
		++ngbIndex;

	}

	// ensure MPI comms are complete
	MPI_Waitall(ngbSet.size(),mpiInRequest,mpiStatus);
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
	
	if (timeAvg == 1)
	{
   		update_time_avg();
	}	

	// ensure MPI comms are complete
	MPI_Waitall(ngbSet.size(),mpiInRequest,mpiStatus);


	// distribute incoming halo data
	insert_halo_data(fOut);

}

