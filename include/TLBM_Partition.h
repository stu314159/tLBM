#ifndef TLBM_PARTITION_H
#define TLBM_PARTITION_H
#include "Problem.h"

class TLBM_Partition{

  public:
    Problem thisProblem;
    TLBM_Partition(int rank, int size);
    ~TLBM_Partition();
 

  private:
    int rank;
    int size;
    int tlbm_initialize();

};

#endif
