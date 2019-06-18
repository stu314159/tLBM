/*
 * LatticeStructure.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: sblair
 */

#ifndef INCLUDE_LATTICESTRUCTURE_HPP_
#define INCLUDE_LATTICESTRUCTURE_HPP_

#include "TLBM_definitions.h"

template < class T >
class LatticeStructure
{
public:
  LatticeStructure();
  ~LatticeStructure();
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
  static inline unsigned getIDx(int nSpd, int nIdx, int spd){
      	return nIdx*nSpd + spd;
      	// return spd*nnods + nIdx; // use this if it performs faster.
      }


protected:
  int numSpd;
  int * ex;
  int * ey;
  int * ez;
  T * w;
  int * bbSpd;
};

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
