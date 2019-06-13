#ifndef D3Q15LATTICESTRUCTURE_H
#define D3Q15LATTICESTRUCTURE_H

#include "LatticeStructure.h"

class D3Q15LatticeStructure : public LatticeStructure
{

public:
  D3Q15LatticeStructure();
  ~D3Q15LatticeStructure();

private:

  const int ex[15] = {0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1};
  const int ey[15] = {0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1};
  const int ez[15] = {0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1};
  const float w[15] = {2.f/9.f,1.f/9.f,1.f/9,1.f/9.f,1.f/9.f,1.f/9.f,1.f/9.f,
		    1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f,
		1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f};
  const int bbSpd[15] = {0,2,1,4,3,6,5,14,13,12,11,10,9,8,7};

};


#endif
