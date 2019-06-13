#ifndef LATTICESTRUCTURE_H
#define LATTICESTRUCTURE_H

#include "TLBM_definitions.h"
class LatticeStructure
{
public:
  LatticeStructure();
  ~LatticeStructure();
  unsigned int get_numSpd(){return numSpd;}
  

protected:
  unsigned int numSpd;

};

#endif
