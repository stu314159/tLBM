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

