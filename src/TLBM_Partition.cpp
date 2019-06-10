#include "TLBM_Partition.h"

TLBM_Partition::TLBM_Partition(int r, int s):
rank(r),size(s)
{
  
  tlbm_initialize();

}

TLBM_Partition::~TLBM_Partition(){

}

int TLBM_Partition::tlbm_initialize(){

  thisProblem.load_input();

  return 0; //<-- indicates no problem with initialization
}

