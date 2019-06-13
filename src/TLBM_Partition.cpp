#include "TLBM_Partition.h"
#include <stdexcept>
#include <string>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <exception>


TLBM_Partition::TLBM_Partition(int r, int s):
rank(r),size(s)
{
  
  tlbm_initialize();

}

TLBM_Partition::~TLBM_Partition(){
  delete myLattice;
  delete adjMatrix;
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


  return 0; //<-- indicates no problem with initialization
}

void TLBM_Partition::load_parts(){
 // read parts.lbm and obtain information about lattice points in my partition as well
// as my neighbors.
   std::ifstream parts("parts.lbm");
   int p;
   int gNdInd = 0; // global node index of current lattice point
   int localNdInd = 0; // local node index counter
   partSizes = std::vector<int>(size,0);
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
	   gNdInd+=1;
   }
   parts.close(); // needed?

   // compute write offsets


//   // say something about nodes in each partition
//   if (rank == 0){
//	   std::cout << "part Sizes = ";
//	   for(auto i = partSizes.begin();i != partSizes.end(); ++i){
//		   std::cout << *i << ' ';
//	   }
//	   std::cout << std::endl;
//   }
//
//   printf("Rank %d size of localNdList = %d \n",rank,static_cast<int>(localNdList.size()));
    writeOffset = 0;
    for (auto i = partSizes.begin(); i < partSizes.begin()+rank;++i){

    	writeOffset += *i;
    }
//    printf("Rank %d write offset = %d \n",rank,writeOffset);
    numLnodes = partSizes[rank]; // just to make this easier
//    printf("Rank %d creating adjList with %d speeds \n",rank,myLattice->get_numSpd());


}

void TLBM_Partition::create_adj_matrix(){
  unsigned int numSpd = myLattice->get_numSpd();
  adjMatrix = new int[numLnodes*numSpd];
  const int * ex = myLattice->get_ex();
  const int * ey = myLattice->get_ey();
  const int * ez = myLattice->get_ez();
  int tgt, node;
  for (int nd = 0; nd < numLnodes;nd++){
	  node = localNdList[nd];// note this is a global node number of the local node
	  for(unsigned spd = 0; spd < numSpd; spd++){
		  tgt = get_tgt_index(node,ex[spd],ey[spd],ez[spd]);
		  adjMatrix[getIDx(numSpd,nd,spd)] = tgt;
	  }
  }
}

int TLBM_Partition::get_tgt_index(int gInd, int ex, int ey, int ez){

	LatticeIndex tgt; //ret.X = 0; ret.Y = 0; ret.Z = 0;

	LatticeIndex src =  get_xyz_index(gInd);

	tgt.X = (src.X + ex)%thisProblem.nx;
	tgt.Y = (src.Y + ey)%thisProblem.ny;
	tgt.Z = (src.Z + ez)%thisProblem.nz;

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

