/*
 * HaloDataOrganizer.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: stu
 */

#ifndef INCLUDE_HALODATAORGANIZER_HPP_
#define INCLUDE_HALODATAORGANIZER_HPP_

#include <map>
#include "TLBM_definitions.h"
#include "HaloDataObject.hpp"

template <class T>
class HaloDataOrganizer
{
public:
	HaloDataOrganizer();
	~HaloDataOrganizer();

private:
	// map keys: neighboring partition rank.
	// map values: a Halo Data Ojbect that contains

	std::map<int,HaloDataObject<T>> Halo;

};

template <class T>
HaloDataOrganizer<T>::HaloDataOrganizer()
{

}

template <class T>
HaloDataOrganizer<T>::~HaloDataOrganizer()
{

}



#endif /* INCLUDE_HALODATAORGANIZER_HPP_ */