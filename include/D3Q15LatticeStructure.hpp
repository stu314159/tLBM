/*
 * D3Q15LatticeStructure.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: sblair
 */

#ifndef INCLUDE_D3Q15LATTICESTRUCTURE_HPP_
#define INCLUDE_D3Q15LATTICESTRUCTURE_HPP_
#include "LatticeStructure.hpp"

template < class T >
class D3Q15LatticeStructure : public LatticeStructure<T>
{

public:
  D3Q15LatticeStructure();
  ~D3Q15LatticeStructure();

private:
  const int numSpd = 15;
  const int ex[15] = {0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1};
  const int ey[15] = {0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1};
  const int ez[15] = {0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1};
  const T w[15] = {2.f/9.f,1.f/9.f,1.f/9,1.f/9.f,1.f/9.f,1.f/9.f,1.f/9.f,
		    1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f,
		1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f};
  const int bbSpd[15] = {0,2,1,4,3,6,5,14,13,12,11,10,9,8,7};

};

template < class T >
D3Q15LatticeStructure<T>::D3Q15LatticeStructure():
LatticeStructure<T>()
{

}

template < class T >
D3Q15LatticeStructure<T>::~D3Q15LatticeStructure(){

}


#endif /* INCLUDE_D3Q15LATTICESTRUCTURE_HPP_ */
