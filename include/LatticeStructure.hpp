/*
 * LatticeStructure.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: sblair
 */

#ifndef INCLUDE_LATTICESTRUCTURE_HPP_
#define INCLUDE_LATTICESTRUCTURE_HPP_

#include "TLBM_definitions.h"
#include <cmath>

template < class T >
class LatticeStructure
{
public:
  LatticeStructure();
  virtual ~LatticeStructure();
  int get_numSpd();
  const int * get_ex();
  const int * get_ey();
  const int * get_ez();
  const T * get_w();
  const int * get_bbSpd();
  void setEx(int * x){ex = x;}
  void setEy(int * y){ey = y;}
  void setEz(int * z){ez = z;}
  void setW(T * tw) {w = tw;}
  void setBB(int * bb){bbSpd = bb;}
  void set_numSpd(int ns){numSpd = ns;}
  inline unsigned getIDx(int nSpd, int nIdx, int spd){
      	return nIdx*nSpd + spd;
      	// return spd*nnods + nIdx; // use this if it performs faster.
      }

  void compute_macroscopic_data(T * ux, T * uy, T * uz, T * rho,
		  const T * fIn, const int nd);
  void compute_equilibrium(T * fEq,
		  const T* ux, const T* uy, const T* uz, const T* rho, const int nd);
  virtual void bounce_back(T * fIN, const int nd) = 0;
  virtual void set_inlet_bc_macro(const T * fIn, T* ux, T* uy, T * uz, T * rho,
		  const T u_bc, const int nd) = 0;
  virtual void set_inlet_bc_micro(T* fIn,const T* fEq, const int nd) = 0;
  virtual void set_outlet_bc_macro(const T * fIn, T * uz, T * rho, const T rho_bc,
		  const int nd) = 0;
  virtual void set_outlet_bc_micro(T* fIn, const T* fEq, const int nd) = 0;

  void set_uz_bc(T* fIn, const T* ux, const T* uy, T* uz, const T* rho,
		  const T u_bc, const int nd);

  void compute_strain_tensor(T* S, const T* fIn, const T* fEq, const int nd);
  T apply_turbulence_model(T omega, const T* S, const T cs);

  void relax(T* fIn, const T* fEq, const T omega, const int nd);

  void compute_piflat(T* piFlat,const T* fIn, const T* fEq, const int nd);


protected:
  int numSpd;
  int * ex;
  int * ey;
  int * ez;
  T * w;
  int * bbSpd;
};

template <class T>
void LatticeStructure<T>::compute_macroscopic_data(T * ux, T * uy, T * uz, T * rho,
		const T * fIn, const int nd){
	T ux_l = 0;
	T uy_l = 0;
	T uz_l = 0;
	T rho_l = 0;
	T f;
	for ( int spd = 0; spd < numSpd; ++spd)
	{
		f = fIn[getIDx(numSpd,nd,spd)];
		rho_l += f;
		ux_l += f*ex[spd];
		uy_l += f*ey[spd];
		uz_l += f*ez[spd];
	}
	ux_l /= rho_l;
	uy_l /= rho_l;
	uz_l /= rho_l;
	ux[nd] = ux_l; uy[nd] = uy_l; uz[nd] = uz_l; rho[nd] = rho_l;

}

template <class T>
void LatticeStructure<T>::compute_equilibrium(T* fEq,
		const T* ux, const T* uy, const T* uz, const T* rho, const int nd)
{
	T cu;
	for (int spd = 0; spd<numSpd; ++spd)
	{
		cu = 3.0*(ex[spd]*ux[nd]+ey[spd]*uy[nd]+ez[spd]*uz[nd]);
		fEq[getIDx(numSpd,nd,spd)]=w[spd]*rho[nd]*
				(1.0+cu+0.5*(cu*cu)-3./2.*(ux[nd]*ux[nd]+uy[nd]*uy[nd]+uz[nd]*uz[nd]));
	}
}

template <class T>
void LatticeStructure<T>::compute_strain_tensor(T* S,const T* fIn, const T* fEq, const int nd)
{
	T e[3];
	const int nDim = 3;
	for(int spd = 0; spd<numSpd;++spd)
	{
		e[0] = ex[spd]; e[1] = ey[spd]; e[2] = ez[spd];
		for(int i = 0; i<nDim; ++i)
		{
			for(int j = 0; j<nDim; ++j)
			{
				S[i*nDim+j]+=e[i]*e[j]*(fIn[this->getIDx(numSpd,nd,spd)] - fEq[this->getIDx(numSpd,nd,spd)]);
			}
		}
	}

}

template <class T>
void LatticeStructure<T>::compute_piflat(T* piFlat, const T* fIn, const T* fEq, const int nd)
{
	T fNeq;
	for(int spd = 0; spd<numSpd; ++spd)
	{
		int idx = getIDx(numSpd,nd,spd);
		fNeq = fIn[idx] - fEq[idx];
		piFlat[0]+=ex[spd]*ey[spd]*fNeq;
		piFlat[1]+=ey[spd]*ex[spd]*fNeq;
		piFlat[2]+=ez[spd]*ex[spd]*fNeq;

	}

}

template <class T>
T LatticeStructure<T>::apply_turbulence_model(T omega, const T* S, const T cs)
{
	T omega_r;
	T nu, nu_e;
	nu = ((1./omega)- 0.5)/3.0;
	T P;
	P = sqrt(S[0]*S[0]+S[4]*S[4]+S[8]*S[8] + 2.0*(S[1]*S[1]+S[2]*S[2] + S[5]*S[5]));
	P*=cs; P = sqrt(P+nu*nu)-nu;
	nu_e = P/6.;
	omega_r = 1./(3.*(nu+nu_e)+0.5);
	return omega_r;

}

template <class T>
void LatticeStructure<T>::relax(T* fIn, const T* fEq, const T omega, const int nd)
{
	for (int spd = 0; spd<numSpd; ++spd)
	{
		int idx = getIDx(numSpd,nd,spd);
		fIn[idx] = fIn[idx] - omega*(fIn[idx] - fEq[idx]);
	}
}

//template <class T>
//void LatticeStructure<T>::bounce_back(T * fIn, const int nd)
//{
////	for (int spd = 0; spd < numSpd; ++spd)
////	{
////		fOut[getIDx(numSpd,nd,bbSpd[spd])] = fIn[getIDx(numSpd,nd,spd)];
////	}
//}

template <class T>
void LatticeStructure<T>::set_uz_bc(T* fIn, const T* ux, const T* uy, T* uz,
		const T* rho, const T u_bc, const int nd)
{
	for ( int spd = 0; spd<numSpd; ++spd)
	{
		// pushes fIn[spd] values towards one such that uz = u_bc, ux = uy = 0.
		fIn[this->getIDx(numSpd,nd,spd)] += 3.*rho[nd]*w[spd]*
				(ez[spd]*(u_bc - uz[nd]) + ex[spd]*(0.-ux[nd]) + ey[spd]*(0. - uy[spd]));
	}

}

template <class T>
int LatticeStructure<T>::get_numSpd(){

	return numSpd;
}


template <class T>
const int * LatticeStructure<T>::get_ex(){
	return ex;
}

template <class T>
const int * LatticeStructure<T>::get_ey(){
	return ey;
}

template <class T>
const int * LatticeStructure<T>::get_ez(){
	return ez;
}

template <class T>
const T * LatticeStructure<T>::get_w(){
	return w;
}

template <class T>
const int * LatticeStructure<T>::get_bbSpd(){
	return bbSpd;
}

template <class T>
LatticeStructure<T>::LatticeStructure():
numSpd(0),ex(0),ey(0),ez(0),w(0),bbSpd(0)
{

}

template <class T>
LatticeStructure<T>::~LatticeStructure(){

}



#endif /* INCLUDE_LATTICESTRUCTURE_HPP_ */
