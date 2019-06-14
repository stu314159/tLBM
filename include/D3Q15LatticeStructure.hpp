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
  static const int numSpd=15;
  int ex[numSpd];
  int ey[numSpd];
  int ez[numSpd];
  T w[numSpd];
  int bbSpd[numSpd];

};

template < class T >
D3Q15LatticeStructure<T>::D3Q15LatticeStructure():
LatticeStructure<T>(),
ex{0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1},
ey{0,0,0,1,-1,0,0,1,1,-1,-1,1,1,-1,-1},
ez{0,0,0,0,0,1,-1,1,1,1,1,-1,-1,-1,-1},
w{2.f/9.f,1.f/9.f,1.f/9,1.f/9.f,1.f/9.f,1.f/9.f,1.f/9.f,
    1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f,
1.f/72.f,1.f/72.f,1.f/72.f,1.f/72.f},
bbSpd{0,2,1,4,3,6,5,14,13,12,11,10,9,8,7}
{
// set base-class pointers
LatticeStructure<T>::setEx(ex);
LatticeStructure<T>::setEy(ey);
LatticeStructure<T>::setEz(ez);
LatticeStructure<T>::setW(w);
LatticeStructure<T>::setBB(bbSpd);
LatticeStructure<T>::set_numSpd(numSpd);

}

template < class T >
D3Q15LatticeStructure<T>::~D3Q15LatticeStructure(){

}


#endif /* INCLUDE_D3Q15LATTICESTRUCTURE_HPP_ */
