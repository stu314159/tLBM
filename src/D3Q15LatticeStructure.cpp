#include "LatticeStructure.h"
#include "D3Q15LatticeStructure.h"

template < class T >
D3Q15LatticeStructure<T>::D3Q15LatticeStructure():
LatticeStructure()
{
numSpd = 15;
}

template < class T >
D3Q15LatticeStructure<T>::~D3Q15LatticeStructure(){

}


void tempFun2(){
	D3Q15LatticeStructure<real> tempObj;
}
