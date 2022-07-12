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
  void set_inletW_bc_macro(const T * fIn, T* ux, T* uy, T * uz, T * rho,
		  const T u_bc, const int nd, const int nnodes);
  void bounce_back(T* fIn, const int nd, const int nnodes);
  void set_inletW_bc_micro(T* fIn, const T* fEq, const int nd, const int nnodes,
		  const int nlnodes);
  void set_outletE_bc_macro(const T* fIn, T* uz, T* rho, const T rho_bc, const int nd,
		  const int nnodes);
  void set_outletE_bc_micro(T* fIn, const T* fEq, const int nd, const int nnodes,
		  const int nlnodes);

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

template < class T >
void D3Q19LatticeStructure<T>::bounce_back(T* fIn, const int nd, const int nnodes){
	T tmp[numSpd];
	for(int s=0;s<numSpd;++s)
	{
		tmp[s] = fIn[this->getIDx(nnodes,nd,bbSpd[s])];
	}
	for(int s=0;s<numSpd;++s)
	{
		fIn[this->getIDx(nnodes,nd,s)]=tmp[s];
	}

}


template <class T>
void D3Q19LatticeStructure<T>::set_inletW_bc_micro(T* fIn, const T* fEq, const int n,
		const int nnodes, const int nlnodes)
{
	int sp[5]={5,11,12,15,16};
	int bbSp[5]={6,14,13,18,17};
	int numBB = 5;
	for(int s=0 ; s<numBB ; ++s)
	{
		fIn[this->getIDx(nnodes,n,sp[s])]= fEq[this->getIDx(nlnodes,n,sp[s])]+
				fIn[this->getIDx(nnodes,n,bbSp[s])]-fEq[this->getIDx(nlnodes,n,bbSp[s])];
	}
}

template <class T>
void D3Q19LatticeStructure<T>::set_outletE_bc_micro(T* fIn, const T* fEq, const int nd,
		const int nnodes, const int nlnodes)
{
	int sp[5]={6,14,13,18,17};
	int bbSp[5]={5,11,12,15,16};
	int numBB = 5;
	for(int s=0;s<numBB;s++)
	{
		fIn[this->getIDx(nnodes,nd,sp[s])] = fEq[this->getIDx(nlnodes,nd,sp[s])] +
				fIn[this->getIDx(nnodes,nd,bbSp[s])]-fEq[this->getIDx(nlnodes,nd,bbSp[s])];
	}

}


template <class T>
void D3Q19LatticeStructure<T>::set_inletW_bc_macro(const T * fIn,T* ux, T* uy,T * uz,
		T * rho, const T u_bc, const int nd, const int nnodes)
{


	T f0, f1, f2, f3, f4, f6, f7, f8, f9, f10, f13, f14, f17, f18;
	f0 = fIn[this->getIDx(nnodes,nd,0)];
	f1 = fIn[this->getIDx(nnodes,nd,1)];
	f2 = fIn[this->getIDx(nnodes,nd,2)];
	f3 = fIn[this->getIDx(nnodes,nd,3)];
	f4 = fIn[this->getIDx(nnodes,nd,4)];
	f6 = fIn[this->getIDx(nnodes,nd,6)];
	f7 = fIn[this->getIDx(nnodes,nd,7)];
	f8 = fIn[this->getIDx(nnodes,nd,8)];
	f9 = fIn[this->getIDx(nnodes,nd,9)];
	f10 = fIn[this->getIDx(nnodes,nd,10)];
	f13 = fIn[this->getIDx(nnodes,nd,13)];
	f14 = fIn[this->getIDx(nnodes,nd,14)];
	f17 = fIn[this->getIDx(nnodes,nd,17)];
	f18 = fIn[this->getIDx(nnodes,nd,18)];
    uz[nd] = u_bc; ux[nd] = 0; uy[nd] = 0;
    rho[nd] = (1./(1.-u_bc))*(2.*(f6 + f13 + f14 + f17 + f18) +
    		(f0 + f1 + f2 + f3 + f4 + f7 + f8 + f9 + f10));

}

template <class T>
void D3Q19LatticeStructure<T>::set_outletE_bc_macro(const T* fIn, T* uz, T* rho, const T rho_bc,
		const int nd, const int nnodes)
{
	rho[nd] = rho_bc;
	T f0,f1,f2,f3,f4,f5,f7,f8,f9,f10,f11,f12,f15,f16;
	f0 = fIn[this->getIDx(nnodes,nd,0)];
	f1 = fIn[this->getIDx(nnodes,nd,1)];
	f2 = fIn[this->getIDx(nnodes,nd,2)];
	f3 = fIn[this->getIDx(nnodes,nd,3)];
	f4 = fIn[this->getIDx(nnodes,nd,4)];
	f5 = fIn[this->getIDx(nnodes,nd,5)];
	f7 = fIn[this->getIDx(nnodes,nd,7)];
	f8 = fIn[this->getIDx(nnodes,nd,8)];
	f9 = fIn[this->getIDx(nnodes,nd,9)];
	f10 = fIn[this->getIDx(nnodes,nd,10)];
	f11 = fIn[this->getIDx(nnodes,nd,11)];
	f12 = fIn[this->getIDx(nnodes,nd,12)];
	f15 = fIn[this->getIDx(nnodes,nd,15)];
	f16 = fIn[this->getIDx(nnodes,nd,16)];

	uz[nd] = -1. + (1./rho_bc)*(2.*(f5+f11+f12+f15+f16)+(f0+f1+f2+f3+f4+f7+f8+f9+f10));

}


#endif /* INCLUDE_D3Q19LATTICESTRUCTURE_HPP_ */
