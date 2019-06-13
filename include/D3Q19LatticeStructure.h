/*
 * D3Q19LatticeStructure.h
 *
 *  Created on: Jun 12, 2019
 *      Author: sblair
 */

#ifndef INCLUDE_D3Q19LATTICESTRUCTURE_H_
#define INCLUDE_D3Q19LATTICESTRUCTURE_H_

#include "LatticeStructure.h"

class D3Q19LatticeStructure : public LatticeStructure
{

public:
  D3Q19LatticeStructure();
  ~D3Q19LatticeStructure();

private:

  const int ex[19] = {0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0};
  const int ey[19] = {0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1};
  const int ez[19] = {0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1};
  const float w[19] = {3.f/9.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,
		    1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,
		1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f};
  const int bbSpd[19] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13,
		  12, 11, 18, 17, 16, 15};

};



#endif /* INCLUDE_D3Q19LATTICESTRUCTURE_H_ */
