#ifndef TLBM_PARTITION_H
#define TLBM_PARTITION_H
#include "Problem.h"
#include "LatticeStructure.h"
#include "D3Q15LatticeStructure.h"
#include "D3Q19LatticeStructure.h"

struct LatticeIndex{
	// essentially a 3-tuple to hold integer indices into a 3D lattice
	// assumes Lattice points are ordered: first X-direction, then Y-direction
	// then Z-direction
	int X;
	int Y;
	int Z;
};

class TLBM_Partition{

  public:
    Problem thisProblem;
    TLBM_Partition(int rank, int size);
    ~TLBM_Partition();
    LatticeIndex get_xyz_index(int gInd);
    int get_gInd(int x, int y, int z);
    int get_gInd(LatticeIndex myXYZ);
 

  private:
    int rank;
    int size;
    int tlbm_initialize();
    LatticeStructure * myLattice;

};

#endif
