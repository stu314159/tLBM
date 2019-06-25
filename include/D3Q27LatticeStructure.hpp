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
  void set_inlet_bc_macro(const T * fIn, T* ux, T* uy, T * uz,
		  T * rho, const T u_bc, const int nd);
  void set_inlet_bc_micro(T* fIn, const T* fEq, const int nd);
  void set_outlet_bc_macro(const T * fIn, T* uz, T* rho, const T rho_bc, const int nd);
  void set_outlet_bc_micro(T* fIn, const T* fEq, const int nd);

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
void D3Q27LatticeStructure<T>::set_inlet_bc_micro(T* fIn, const T* fEq, const int n)
{
	int sp[9]={5,11,13,15,17,19,21,23,25};
	int bbSp[9]={6,14,12,18,16,26,24,22,20};
	int numBB = 9;
	for(int s=0 ; s<numBB ; ++s)
	{
		fIn[this->getIDx(numSpd,n,sp[s])]= fEq[this->getIDx(numSpd,n,sp[s])]+
				fIn[this->getIDx(numSpd,n,bbSp[s])]-fEq[this->getIDx(numSpd,n,bbSp[s])];
	}

}

template <class T>
void D3Q27LatticeStructure<T>::set_outlet_bc_micro(T* fIn, const T* fEq, const int nd)
{
	int sp[9]={6,14,12,18,16,26,24,22,20};
	int bbSp[9]={5,11,13,15,17,19,21,23,25};
	int numBB = 9;
	for(int s=0;s<numBB;s++)
	{
		fIn[this->getIDx(numSpd,nd,sp[s])] = fEq[this->getIDx(numSpd,nd,sp[s])] +
				fIn[this->getIDx(numSpd,nd,bbSp[s])]-fEq[this->getIDx(numSpd,nd,bbSp[s])];
	}

}

template <class T>
void D3Q27LatticeStructure<T>::set_inlet_bc_macro(const T * fIn, T* ux, T* uy,T * uz,
		T * rho, const T u_bc, const int nd)
{

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

template <class T>
void D3Q27LatticeStructure<T>::set_outlet_bc_macro(const T* fIn, T* uz, T* rho, const T rho_bc, const int nd)
{
	rho[nd] = rho_bc;
	T f0,f1,f2,f3,f4,f5,f7,f8,f9,f10,f11,f13,f15,f17,f19,f21,f23,f25;
	f0 = fIn[this->getIDx(numSpd,nd,0)];
	f1 = fIn[this->getIDx(numSpd,nd,1)];
	f2 = fIn[this->getIDx(numSpd,nd,2)];
	f3 = fIn[this->getIDx(numSpd,nd,3)];
	f4 = fIn[this->getIDx(numSpd,nd,4)];
	f5 = fIn[this->getIDx(numSpd,nd,5)];
	f7 = fIn[this->getIDx(numSpd,nd,7)];
	f8 = fIn[this->getIDx(numSpd,nd,8)];
	f9 = fIn[this->getIDx(numSpd,nd,9)];
	f10 = fIn[this->getIDx(numSpd,nd,10)];
	f11 = fIn[this->getIDx(numSpd,nd,11)];
	f13 = fIn[this->getIDx(numSpd,nd,13)];
	f15 = fIn[this->getIDx(numSpd,nd,15)];
	f17 = fIn[this->getIDx(numSpd,nd,17)];
	f19 = fIn[this->getIDx(numSpd,nd,19)];
	f21 = fIn[this->getIDx(numSpd,nd,21)];
	f23 = fIn[this->getIDx(numSpd,nd,23)];
	f25 = fIn[this->getIDx(numSpd,nd,25)];
	uz[nd] = -1. + (1./rho_bc)*(2.0*(f5+f11+f13+f15+f17+f19+f21+f23+f25)+(f0+f1+f2+f3+f4+f7+f8+f9+f10));

}

#endif /* INCLUDE_D3Q27LATTICESTRUCTURE_HPP_ */
