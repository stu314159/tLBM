/*
 * HaloDataObject.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: stu
 */

#ifndef INCLUDE_HALODATAOBJECT_HPP_
#define INCLUDE_HALODATAOBJECT_HPP_

#include "TLBM_definitions.h"
#include <map>
#include <vector>


template <class T>
class HaloDataObject
{
public:
	HaloDataObject();
	~HaloDataObject();
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



#endif /* INCLUDE_HALODATAOBJECT_HPP_ */
