/*
 * HaloDataObject.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: stu
 */

#ifndef INCLUDE_HALODATAOBJECT_HPP_
#define INCLUDE_HALODATAOBJECT_HPP_

//#include "TLBM_definitions.h"
#include <map>
#include <vector>


template <class T>
class HaloDataObject
{
public:
	HaloDataObject();
	~HaloDataObject();
	void insert_item(int gnn,int spd);
	std::vector<int> & operator[](int gnd);

private:
	std::map<int,std::vector<int>> DataMap;
	int numItems;
	T* buffer;
};


template <class T>
HaloDataObject<T>::HaloDataObject():
numItems(0), buffer(0)
{

}

template <class T>
HaloDataObject<T>::~HaloDataObject()
{

}

template <class T>
std::vector<int> & HaloDataObject<T>::operator[](int gnn)
{
	return DataMap[gnn];
}
template <class T>
void HaloDataObject<T>::insert_item(int gnn,int spd)
{
	numItems+=1;
	DataMap[gnn].push_back(spd);

}


#endif /* INCLUDE_HALODATAOBJECT_HPP_ */
