#ifndef TLBM_PARTITION_H
#define TLBM_PARTITION_H

#include <mpi.h>

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
    TLBM_Partition(int rank, int size, MPI_Comm comm);
    ~TLBM_Partition();
    LatticeIndex get_xyz_index(int gInd);
    int get_tgt_index(int gInd, int ex, int ey, int ez);
    int get_gInd(int x, int y, int z);
    int get_gInd(LatticeIndex myXYZ);
    int get_cut_size();
    int get_num_ts();
    int get_ts_rep_freq();
    int get_plot_freq();
    void take_LBM_time_step(bool isEven);


  private:
    int rank;
    int size;
    MPI_Comm comm;
    LatticeStructure<real> * myLattice;
    std::vector<int> localNdList; // my lattice points (global node numbers)
    std::vector<int> partSizes; // number of LPs in each partition
    std::vector<int> partsG; // partition assignment for each node by global node number
    std::set<int> boundaryNdList;
    std::set<int> haloNodes;
    std::map<int,int> globalToLocal;
    std::map<int,int> localToGlobal;
    HaloDataOrganizer<real> HDO_out;
    HaloDataOrganizer<real> HDO_in;
    std::set<int> ngbSet;
    int numLnodes;
    int numHaloNodes;
    int totalNodes;
    int writeOffset;
    int * adjMatrix = NULL;
    int * ndType = NULL;
    real * fEven = NULL;
    real * fOdd = NULL;
    real * fIn = NULL;
    real * fOut = NULL;
    real * ux = NULL;
    real * uy = NULL;
    real * uz = NULL;
    real * rho = NULL;
    int tlbm_initialize();
    void load_parts();
    void create_adj_matrix();
    void compute_halo_data();
    void make_adj_matrix_local();
    void allocate_arrays();
    void load_ndType();
    void initialize_data_arrays();
    void write_node_ordering();

    static inline unsigned getIDx(int nSpd, int nIdx, int spd){
    	return nIdx*nSpd + spd;
    	// return spd*nnods + nIdx; // use this if it performs faster.
    }

};

#endif
