#ifndef TLBM_PARTITION_H
#define TLBM_PARTITION_H

#include "TLBM_definitions.h" // global constants and definitions.  Hack-ish?
#include "Problem.h"
#include "LatticeStructure.hpp"
#include "D3Q15LatticeStructure.hpp"
#include "D3Q19LatticeStructure.hpp"
#include "D3Q27LatticeStructure.hpp"
#include "HaloDataOrganizer.hpp"
#include "HaloDataObject.hpp"


#include <list>
#include <vector>
#include <map>
#include <set>



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
    int get_tgt_index(int gInd, int ex, int ey, int ez);
    int get_gInd(int x, int y, int z);
    int get_gInd(LatticeIndex myXYZ);


  private:
    int rank;
    int size;
    LatticeStructure<real> * myLattice;
    std::vector<int> localNdList; // my lattice points (global node numbers)
    std::vector<int> partSizes; // number of LPs in each partition
    std::vector<int> partsG; // partition assignment for each node by global node number
    std::set<int> boundaryNdList;
    std::map<int,int> globalToLocal;
    std::map<int,int> localToGlobal;
    std::map<int,HaloDataOrganizer<real>> HDO_out_dict;
    std::map<int,HaloDataOrganizer<real>> HDO_in_dict;
    std::set<int> ngbSet;
    int numLnodes;
    int writeOffset;
    int * adjMatrix = NULL;
    int tlbm_initialize();
    void load_parts();
    void create_adj_matrix();
    void compute_halo_data();

    static inline unsigned getIDx(int nSpd, int nIdx, int spd){
    	return nIdx*nSpd + spd;
    	// return spd*nnods + nIdx; // use this if it performs faster.
    }

};

#endif
