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
  void set_inlet_bc_macro(const T * fIn,T* ux, T* uy, T * uz, T * rho,
		  const T u_bc, const int nd);

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

template < class T >
void D3Q15LatticeStructure<T>::set_inlet_bc_macro(const T * fIn, T* ux, T* uy,
		T * uz, T * rho, const T u_bc, const int nd)
{
//	f.uz = f.u_bc;
//		f.ux = 0.; f.uy = 0.;
//		f.rho = (1./(1. - f.uz))*(2.*
//				(f.f[6]+f.f[11]+f.f[12]+f.f[13]+f.f[14])+
//				(f.f[0]+f.f[1]+f.f[2]+f.f[3]+f.f[4]));
	T f6,f11,f12,f13,f14,f0,f1,f2,f3,f4;
	f6 = fIn[this->getIDx(numSpd,nd,6)];
	f11 = fIn[this->getIDx(numSpd,nd,11)];
	f12 = fIn[this->getIDx(numSpd,nd,12)];
	f13 = fIn[this->getIDx(numSpd,nd,13)];
	f14 = fIn[this->getIDx(numSpd,nd,14)];
	f0 = fIn[this->getIDx(numSpd,nd,0)];
	f1 = fIn[this->getIDx(numSpd,nd,1)];
	f2 = fIn[this->getIDx(numSpd,nd,2)];
	f3 = fIn[this->getIDx(numSpd,nd,3)];
	f4 = fIn[this->getIDx(numSpd,nd,4)];
	ux[nd] = 0; uy[nd] = 0; uz[nd] = u_bc;
	rho[nd] = (1./(1. - u_bc))*(2.*
			(f6 + f11 + f12 + f13 + f14) +
			(f0 + f1 + f2 + f3 + f4));


}



#endif /* INCLUDE_D3Q15LATTICESTRUCTURE_HPP_ */
