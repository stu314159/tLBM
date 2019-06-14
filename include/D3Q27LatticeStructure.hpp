/*
 * D3Q27LatticeStructure.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: sblair
 */

#ifndef INCLUDE_D3Q27LATTICESTRUCTURE_HPP_
#define INCLUDE_D3Q27LATTICESTRUCTURE_HPP_

#include "LatticeStructure.hpp"

template < class T>
class D3Q27LatticeStructure : public LatticeStructure<T>
{

public:
  D3Q27LatticeStructure();
  ~D3Q27LatticeStructure();

private:
  static const int numSpd = 27;
  int ex[numSpd];
  int ey[numSpd];
  int ez[numSpd];
  T w[numSpd];
  int bbSpd[numSpd];

};

template <class T>
D3Q27LatticeStructure<T>::D3Q27LatticeStructure():
LatticeStructure<T>(),
ex{0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1,
	  0,0,0,0,1,1,1,1,-1,-1,-1,-1},
ey{0,0,0,1,-1,0,0,1,-1,1,-1,0,0,0,0,1,
		  1,-1,-1,1,1,-1,-1,1,1,-1,-1},
ez{0,0,0,0,0,1,-1,0,0,0,0,1,-1,1,-1,1,
			  -1,1,-1,1,-1,1,-1,1,-1,1,-1},
w{8./27.,2./27.,2./27.,2./27.,2./27.,2./27.,2./27.,
				    1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
				    1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
				    1./216.,1./216.,1./216.,1./216.,
				1./216.,1./216.,1./216.,1./216.},
bbSpd{0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,
					  16,15,26,25,24,23,22,21,20,19}
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
D3Q27LatticeStructure<T>::~D3Q27LatticeStructure(){

}



#endif /* INCLUDE_D3Q27LATTICESTRUCTURE_HPP_ */
