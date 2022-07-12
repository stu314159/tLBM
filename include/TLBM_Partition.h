#ifndef TLBM_PARTITION_H
#define TLBM_PARTITION_H

#include <mpi.h>
#include "H5Cpp.h"

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
    int get_warmupTs();
    int get_num_global_nodes();
    int is_restart();
    int get_time_avg_flag();
    real get_data_member(const real * f, const int nd, const int spd);
    void set_data_member(real * f, const real val, const int nd, const int spd);
    void take_LBM_time_step(bool isEven);
    void write_data();
    void write_time_avg_data();
    void process_node_list(real * fOut, real * fIn, const std::set<int>& nodeList);
    static inline unsigned getIDx(int nnodes, int nIdx, int spd){
//    	return nIdx*nSpd + spd;
        return spd*nnodes + nIdx; // use this if it performs faster.
    }
    void print_adjacency(const int nIdx);


  private:
    int rank;
    int size;
    MPI_Comm comm;
    LatticeStructure<real> * myLattice;
    std::vector<int> localNdList; // my lattice points (global node numbers) - this includes Halo nodes
    std::vector<int> partSizes; // number of LPs in each partition
    std::vector<int> partsG; // partition assignment for each node by global node number
    std::set<int> boundaryNdList; // local node number set of boundary nodes
    std::set<int> interiorNdList; // local node number set of non-boundary nodes
    std::set<int> haloNodes;
    std::map<int,int> globalToLocal;
    std::map<int,int> localToGlobal;

    std::map<int, std::set<int> > forceCalcMap;

    HaloDataOrganizer<real> HDO_out;
    HaloDataOrganizer<real> HDO_in;
    std::set<int> ngbSet;
    int numLnodes;
    int numHaloNodes;
    int totalNodes;
    int writeOffset;
    int dataWriteNum;
    int * adjMatrix = NULL;
    int * ndType = NULL;
    real * fEven = NULL;
    real * fOdd = NULL;
    real * fEq = NULL;
    real * fIn = NULL;
    real * fOut = NULL;
    real * ux = NULL;
    real * uy = NULL;
    real * uz = NULL;
    real * rho = NULL;
    real * uAvg = NULL;
    real * vAvg = NULL;
    real * wAvg = NULL;
    real * rhoAvg = NULL;
    real * Fx = NULL;
    real * Fy = NULL;
    real * Fz = NULL;

//    // Device Arrays
//    real * d_fIn = NULL;
//    real * d_fOut = NULL;
//    real * d_fEven = NULL;
//    real * d_fOdd = NULL;
//    real * d_ux = NULL;
//    real * d_uy = NULL;
//    real * d_uz = NULL;
//    real * d_uAvg = NULL;
//    real * d_vAvg = NULL;
//    real * d_wAvg = NULL;
//    real * d_rhoAvg = NULL;
//    real * d_Fx = NULL;
//    real * d_Fy = NULL;
//    real * d_Fz = NULL;
//    int * d_INL = NULL;
//    int * d_BNL = NULL;
//    int * d_ndType = NULL;

    bool timeAvg;

    MPI_Request * mpiOutRequest = NULL;
    MPI_Request * mpiInRequest = NULL;
    MPI_Status * mpiStatus = NULL;
    int tlbm_initialize();
    void load_parts();
    void create_adj_matrix();
    void compute_halo_data();
    void make_adj_matrix_local();
    void allocate_arrays();
    void load_ndType();
    void initialize_data_arrays();
    void load_restart_data();
    void finalize_halo_data_arrays();
    void write_node_ordering();
    void make_interior_node_list();
    void stream_node_data(real * fOut, const real * fIn, const int nd);
    void extract_halo_data(real * fOut);
    void insert_halo_data(real * fOut);
    void initiate_data_exchange();
    void update_time_avg();

    void make_force_calc_map();
    void calc_force();

};

#endif
