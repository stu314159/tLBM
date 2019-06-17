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
#include <set>
#include <iostream>
#include <cstdio>
#include <string>


template <class T>
class HaloDataObject
{
public:
	HaloDataObject();
	~HaloDataObject();
	void insert_item(int gnn,int spd);
	std::set<int> & operator[](int gnd);
	void print_halo() const;
	int print_num_nodes() const;
	int get_num_items() const;
	int get_num_nodes()const;
	std::set<int> get_halo_nodes() const;

private:
	std::map<int,std::set<int>> DataMap;
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
std::set<int> HaloDataObject<T>::get_halo_nodes() const
{
	std::set<int> haloNodes;
	for (const auto & gnn : DataMap)
	{
		haloNodes.insert(gnn.first);
	}
	return haloNodes;

}


template <class T>
int HaloDataObject<T>::get_num_items() const
{
	return numItems;
}

template <class T>
int HaloDataObject<T>::get_num_nodes() const
{
	return DataMap.size();
}

template <class T>
int HaloDataObject<T>::print_num_nodes() const
{
	return DataMap.size(); // why, again, do I have this function?
}

template <class T>
void HaloDataObject<T>::print_halo() const
{
	// iterate over all keys (global nodes) for the object
	for(const auto &gnd_iter : DataMap)
	{
		printf("For global node %d: ",gnd_iter.first);
		std::string spdList ("");
		for(const auto & spd_iter : gnd_iter.second)
		{
			spdList += std::to_string(spd_iter) + " ";
		}
		printf("%s \n",spdList.c_str());

	}

}

template <class T>
std::set<int> & HaloDataObject<T>::operator[](int gnn)
{
	return DataMap[gnn];
}
template <class T>
void HaloDataObject<T>::insert_item(int gnn,int spd)
{
//	printf("For node %d, inserting speed %d \n",gnn,spd);
	numItems+=1;
	DataMap[gnn].insert(spd);


}


#endif /* INCLUDE_HALODATAOBJECT_HPP_ */
