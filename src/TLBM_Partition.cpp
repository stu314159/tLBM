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

  if (thisProblem.latticeType == std::string("D3Q15"))
  {
	 myLattice = new D3Q15LatticeStructure;
	 printf("Partition %d creating a 15-speed lattice\n",rank);

  }else if(thisProblem.latticeType == std::string("D3Q19"))
  {
	  myLattice = new D3Q19LatticeStructure;
	  printf("Partition %d creating a 19-speed lattice \n",rank);

  }else {
	  throw std::invalid_argument("Invalid Lattice Structure!");
  }

//  // initialize myLattice pointer depending on the lattice type
//  switch (thisProblem.latticeType)
//  {
//  case std::string("D3Q15"):
//	  myLattice = new D3Q15LatticeStructure;
//	  break;
//  case std::string("D3Q19"):
//	  myLattice = new D3Q19LatticeStructure;
//	  break;
//  case default:
//	  throw std::invalid_argument("Invalid Lattice Structure");
//  }

  return 0; //<-- indicates no problem with initialization
}

