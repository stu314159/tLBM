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
  void set_inlet_bc_macro(T * fIn, T * uz, T * rho, const T u_bc, const int nd);

private:
  static const int numSpd = 19;
  int ex[numSpd];
  int ey[numSpd];
  int ez[numSpd];
  T w[numSpd];
  int bbSpd[numSpd];

};

template <class T>
D3Q19LatticeStructure<T>::D3Q19LatticeStructure():
LatticeStructure<T>(),
ex{0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0},
ey{0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1},
ez{0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1},
w{3.f/9.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,1.f/18.f,
    1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,
1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f},
bbSpd{0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13,
	  12, 11, 18, 17, 16, 15}
{
	// set base-class pointers
	// set base-class pointers
	LatticeStructure<T>::setEx(ex);
	LatticeStructure<T>::setEy(ey);
	LatticeStructure<T>::setEz(ez);
	LatticeStructure<T>::setW(w);
	LatticeStructure<T>::setBB(bbSpd);
	LatticeStructure<T>::set_numSpd(numSpd);
}

template <class T>
D3Q19LatticeStructure<T>::~D3Q19LatticeStructure(){

}

template <class T>
void D3Q19LatticeStructure<T>::set_inlet_bc_macro(T * fIn,T * uz, T * rho, const T u_bc, const int nd)
{

}


#endif /* INCLUDE_D3Q19LATTICESTRUCTURE_HPP_ */
