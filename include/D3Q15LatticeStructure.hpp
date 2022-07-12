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
  void set_inletW_bc_macro(const T * fIn,T* ux, T* uy, T * uz, T * rho,
		  const T u_bc, const int nd, const int nnodes);
  void set_inletW_bc_micro(T* fIn, const T* fEq, const int nd, const int nnodes,
		  const int nlnodes);
  void set_outletE_bc_macro(const T * fIn, T * uz, T * rho,const T rho_bc, const int nd,
		  const int nnodes);
  void set_outletE_bc_micro(T* fIn, const T* fEq, const int nd, const int nnodes,
		  const int nlnodes);
  void bounce_back(T* fIn, const int nd, const int nnodes);

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
void D3Q15LatticeStructure<T>::bounce_back(T* fIn, const int nd, const int nnodes){
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

template < class T >
D3Q15LatticeStructure<T>::~D3Q15LatticeStructure(){

}

template <class T>
void D3Q15LatticeStructure<T>::set_inletW_bc_micro(T* fIn, const T* fEq, const int n,
		const int nnodes, const int nlnodes)
{
//	int sp[5]={5,7,8,9,10};
//	  int bbSp[5]={6,11,12,13,14};
//	  int numBB = 5;
//	  for(int s=0;s<numBB;s++)
//	  {
//		  f.f[sp[s]]=f.fEq[sp[s]]+f.f[bbSp[s]]-f.fEq[bbSp[s]];
//	  }

	int sp[5] = {5,7,8,9,10};
	int bbSp[5] = {6,11,12,13,14};
	int numBB = 5;
	for(int s=0 ; s<numBB ; ++s)
	{
		fIn[this->getIDx(nnodes,n,sp[s])]= fEq[this->getIDx(nlnodes,n,sp[s])]+
				fIn[this->getIDx(nnodes,n,bbSp[s])]-fEq[this->getIDx(nlnodes,n,bbSp[s])];
	}

}

template < class T>
void D3Q15LatticeStructure<T>::set_outletE_bc_micro(T* fIn, const T* fEq, const int nd,
		const int nnodes, const int nlnodes)
{
	int sp[5]={6,11,12,13,14};
	int bbSp[5]={5,7,8,9,10};
	int numBB = 5;
	for(int s=0;s<numBB;s++)
	{
		fIn[this->getIDx(nnodes,nd,sp[s])] = fEq[this->getIDx(nlnodes,nd,sp[s])] +
				fIn[this->getIDx(nnodes,nd,bbSp[s])]-fEq[this->getIDx(nlnodes,nd,bbSp[s])];
	}

}

template < class T >
void D3Q15LatticeStructure<T>::set_inletW_bc_macro(const T * fIn, T* ux, T* uy,
		T * uz, T * rho, const T u_bc, const int nd, const int nnodes)
{

	T f6,f11,f12,f13,f14,f0,f1,f2,f3,f4;
	f6 = fIn[this->getIDx(nnodes,nd,6)];
	f11 = fIn[this->getIDx(nnodes,nd,11)];
	f12 = fIn[this->getIDx(nnodes,nd,12)];
	f13 = fIn[this->getIDx(nnodes,nd,13)];
	f14 = fIn[this->getIDx(nnodes,nd,14)];
	f0 = fIn[this->getIDx(nnodes,nd,0)];
	f1 = fIn[this->getIDx(nnodes,nd,1)];
	f2 = fIn[this->getIDx(nnodes,nd,2)];
	f3 = fIn[this->getIDx(nnodes,nd,3)];
	f4 = fIn[this->getIDx(nnodes,nd,4)];
	ux[nd] = 0; uy[nd] = 0; uz[nd] = u_bc;
	rho[nd] = (1./(1. - u_bc))*(2.*
			(f6 + f11 + f12 + f13 + f14) +
			(f0 + f1 + f2 + f3 + f4));


}

template < class T >
void D3Q15LatticeStructure<T>::set_outletE_bc_macro(const T * fIn, T* uz, T* rho, const T rho_bc,
		const int nd, const int nnodes)
{
	rho[nd] = rho_bc;
	T f0,f1,f2,f3,f4,f5,f7,f8,f9,f10;
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
	uz[nd] = -1.0 + (1.0/rho_bc)*(2.0*(f5+f7+f8+f9+f10)+(f0+f1+f2+f3+f4));

}



#endif /* INCLUDE_D3Q15LATTICESTRUCTURE_HPP_ */
