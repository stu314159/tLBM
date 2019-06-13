#ifndef TLBM_PARTITION_H
#define TLBM_PARTITION_H
#include "Problem.h"
#include "LatticeStructure.h"
#include "D3Q15LatticeStructure.h"
#include "D3Q19LatticeStructure.h"
#include "TLBM_definitions.h" // global constants and definitions.  Hack-ish?

#include <list>
#include <vector>
#include <map>



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
    LatticeStructure * myLattice;
    std::vector<int> localNdList; // my lattice points
    std::vector<int> partSizes; // number of LPs in each partition
    std::map<int,int> globalToLocal;
    std::map<int,int> localToGlobal;
    int numLnodes;
    int writeOffset;
    int * adjList = NULL;
    int tlbm_initialize();
    void load_parts();

};

#endif
