/*
 * HaloDataOrganizer.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: stu
 */

#ifndef INCLUDE_HALODATAORGANIZER_HPP_
#define INCLUDE_HALODATAORGANIZER_HPP_

#include <map>
//#include "TLBM_definitions.h"
#include "HaloDataObject.hpp"

template <class T>
class HaloDataOrganizer
{
public:
	HaloDataOrganizer();
	//HaloDataOrganizer(int ngb);
	~HaloDataOrganizer();
	void add_neighbor(int ngbNum);
	int get_num_neighbors();
	int get_cut_size();
	HaloDataObject<T>& operator[](int k);
	void print_halo();

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

template <class T>
void HaloDataOrganizer<T>::print_halo()
{
	for(const auto &ngb_it : Halo)
	{
		printf("Neighbor Partition %d: \n",ngb_it.first);
		ngb_it.second.print_halo();
	}

}

template <class T>
int HaloDataOrganizer<T>::get_cut_size()
{
	int numCuts = 0;
	for (const auto & keyIter : Halo)
	{
		numCuts += keyIter.second.get_num_items();
	}
	return numCuts; // fix this in a minute

}

template <class T>
HaloDataObject<T> & HaloDataOrganizer<T>::operator[](int k)
{
	return Halo[k];
}

template <class T>
int HaloDataOrganizer<T>::get_num_neighbors()
{
	return Halo.size();
}

template <class T>
void HaloDataOrganizer<T>::add_neighbor(int ngbNum)
{
	Halo[ngbNum] = HaloDataObject<T>();
}


#endif /* INCLUDE_HALODATAORGANIZER_HPP_ */
