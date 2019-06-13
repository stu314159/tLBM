/*
 * D3Q19LatticeStructure.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: sblair
 */

#ifndef INCLUDE_D3Q19LATTICESTRUCTURE_HPP_
#define INCLUDE_D3Q19LATTICESTRUCTURE_HPP_

#include "LatticeStructure.hpp"

template < class T>
class D3Q19LatticeStructure : public LatticeStructure<T>
{

public:
  D3Q19LatticeStructure();
  ~D3Q19LatticeStructure();

private:
  const int numSpd = 19;
  const int ex[19] = {0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0};
  const int ey[19] = {0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1};
  const int ez[19] = {0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1};
  const T w[19] = {3.f/9.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,
		    1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,
		1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f};
  const int bbSpd[19] = {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13,
		  12, 11, 18, 17, 16, 15};

};

template <class T>
D3Q19LatticeStructure<T>::D3Q19LatticeStructure():
LatticeStructure<T>()
{

}

template <class T>
D3Q19LatticeStructure<T>::~D3Q19LatticeStructure(){

}



#endif /* INCLUDE_D3Q19LATTICESTRUCTURE_HPP_ */
