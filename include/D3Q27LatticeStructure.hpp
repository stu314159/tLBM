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
  void set_inlet_bc_macro(const T * fIn, T* ux, T* uy, T * uz, T * rho, const T u_bc, const int nd);

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

template <class T>
void D3Q27LatticeStructure<T>::set_inlet_bc_macro(const T * fIn, T* ux, T* uy,T * uz,
		T * rho, const T u_bc, const int nd)
{
//	f.uz = f.u_bc;
//		f.ux = 0; f.uy = 0.;
//		f.rho = (1./(1. - f.uz))*(2.*(f.f[6]+f.f[14]+f.f[12]+
//				f.f[18]+f.f[16]+f.f[26]+f.f[24]+f.f[22]+f.f[20])+
//				(f.f[0]+f.f[1]+f.f[2]+f.f[3]+f.f[4]+
//						f.f[7]+f.f[8]+f.f[9]+f.f[10]));

	T f6, f14, f12, f18, f16, f26, f24, f22, f20, f0, f1, f2, f3, f4, f7, f8, f9, f10;
	f0 = fIn[this->getIDx(numSpd,nd,0)];
	f1 = fIn[this->getIDx(numSpd,nd,1)];
	f2 = fIn[this->getIDx(numSpd,nd,2)];
	f3 = fIn[this->getIDx(numSpd,nd,3)];
	f4 = fIn[this->getIDx(numSpd,nd,4)];
	f6 = fIn[this->getIDx(numSpd,nd,6)];
	f7 = fIn[this->getIDx(numSpd,nd,7)];
	f8 = fIn[this->getIDx(numSpd,nd,8)];
	f9 = fIn[this->getIDx(numSpd,nd,9)];
	f10 = fIn[this->getIDx(numSpd,nd,10)];
	f12 = fIn[this->getIDx(numSpd,nd,12)];
	f14 = fIn[this->getIDx(numSpd,nd,14)];
	f16 = fIn[this->getIDx(numSpd,nd,16)];
	f18 = fIn[this->getIDx(numSpd,nd,18)];
	f20 = fIn[this->getIDx(numSpd,nd,20)];
	f22 = fIn[this->getIDx(numSpd,nd,22)];
	f24 = fIn[this->getIDx(numSpd,nd,24)];
	f26 = fIn[this->getIDx(numSpd,nd,26)];

	ux[nd] = 0.; uy[nd] = 0.; uz[nd] = u_bc;
	rho[nd] = (1./(1. - u_bc))*(2*(f6 + f14 + f12 + f18 + f16 + f26 + f24 + f22 + f20) +
			(f0 + f1 + f2 + f3 + f4 + f7 + f8 + f9 + f10));


}

#endif /* INCLUDE_D3Q27LATTICESTRUCTURE_HPP_ */
