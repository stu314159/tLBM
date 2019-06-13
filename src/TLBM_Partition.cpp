#include "TLBM_Partition.h"
#include <stdexcept>
#include <string>
#include <cstdio>

TLBM_Partition::TLBM_Partition(int r, int s):
rank(r),size(s)
{
  
  tlbm_initialize();

}

TLBM_Partition::~TLBM_Partition(){
  delete myLattice;
}

int TLBM_Partition::tlbm_initialize(){

  thisProblem.load_input();

  // construct appropriate lattice type
  if (thisProblem.latticeType == std::string("D3Q15"))
  {
	 myLattice = new D3Q15LatticeStructure;


  }else if(thisProblem.latticeType == std::string("D3Q19"))
  {
	  myLattice = new D3Q19LatticeStructure;


  }else {
	  throw std::invalid_argument("Invalid Lattice Structure!");
  }



  return 0; //<-- indicates no problem with initialization
}

void TLBM_Partition::load_parts(){
 // read parts.lbm and obtain information about lattice points in my partition as well
// as my neighbors.

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

