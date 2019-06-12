#ifndef TLBM_PARTITION_H
#define TLBM_PARTITION_H
#include "Problem.h"
#include "LatticeStructure.h"
#include "D3Q15LatticeStructure.h"
#include "D3Q19LatticeStructure.h"

class TLBM_Partition{

  public:
    Problem thisProblem;
    TLBM_Partition(int rank, int size);
    ~TLBM_Partition();
 

  private:
    int rank;
    int size;
    int tlbm_initialize();
    LatticeStructure * myLattice;

};

#endif
