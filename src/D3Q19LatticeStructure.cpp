/*
 * D3Q19LatticeStructure.cpp
 *
 *  Created on: Jun 12, 2019
 *      Author: sblair
 */

#include "LatticeStructure.h"
#include "D3Q19LatticeStructure.h"

template <class T>
D3Q19LatticeStructure<T>::D3Q19LatticeStructure():
LatticeStructure()
{
numSpd = 19;
}

template <class T>
D3Q19LatticeStructure<T>::~D3Q19LatticeStructure(){

}

void tempFun1(){
	D3Q19LatticeStructure<real> tempObj;
}

