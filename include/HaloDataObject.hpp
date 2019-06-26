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
	T* get_buffer() const;
	int * get_ndNums(){return ndNums;}
	int * get_spds(){return spds;}

	void allocate_arrays();
	void fill_nums_and_speeds(std::map<int,int> & globalToLocal);
	static inline unsigned getIDx(int nSpd, int nIdx, int spd){
		return nIdx*nSpd + spd;
		// return spd*nnods + nIdx; // use this if it performs faster.
	}
	void extract_halo_data(const T* fOut, const int numSpd);
	void insert_halo_data(T* fOut, const int numSpd);


private:
	std::map<int,std::set<int>> DataMap; // map between global node numbers and speeds
	int numItems;
	T* buffer;
	int * ndNums;
	int * spds;

};


template <class T>
void HaloDataObject<T>::allocate_arrays()
{
  buffer = new T[numItems];
  ndNums = new int[numItems];
  spds = new int[numItems];
}

template <class T>
HaloDataObject<T>::HaloDataObject():
numItems(0), buffer(0), ndNums(0), spds(0)
{

}

template <class T>
void HaloDataObject<T>::extract_halo_data(const T* fOut, const int numSpd)
{
	for(int i = 0; i<numItems; ++i)
	{
		buffer[i] = fOut[getIDx(numSpd,ndNums[i],spds[i])];
	}

}

template <class T>
void HaloDataObject<T>::insert_halo_data(T* fOut, const int numSpd)
{
	for(int i = 0; i<numItems; ++i)
	{
		fOut[getIDx(numSpd,ndNums[i],spds[i])] = buffer[i];
	}
}

template <class T>
void HaloDataObject<T>::fill_nums_and_speeds(std::map<int,int> & globalToLocal)
{
	int entry = 0;
	for(const auto & gNdIt : DataMap)
	{
		int gNd = gNdIt.first;
		int lNd = globalToLocal[gNd];
		for(const auto & spdIt : gNdIt.second)
		{
			ndNums[entry] = lNd;
			spds[entry] = spdIt;
			++entry;
		}

	}

}

template <class T>
HaloDataObject<T>::~HaloDataObject()
{
  delete [] buffer;
  delete [] ndNums;
  delete [] spds;
}

template <class T>
T* HaloDataObject<T>::get_buffer() const
{
	return buffer;
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
