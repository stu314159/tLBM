/*
 * HaloDataOrganizer.hpp
 *
 *  Created on: Jun 13, 2019
 *      Author: stu
 */

#ifndef INCLUDE_HALODATAORGANIZER_HPP_
#define INCLUDE_HALODATAORGANIZER_HPP_

#include <map>
#include <set>
#include <algorithm>
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
	int get_num_halo_nodes();
	HaloDataObject<T>& operator[](int k);
	void print_halo();
	std::set<int> get_halo_nodes() const;
	void allocate_halo_arrays();
	int check_ndNums_and_spds(std::map<int,int> & globalToLocal);
	void fill_arrays(std::map<int,int> & globalToLocal);
	static inline unsigned getIDx(int nnodes, int nIdx, int spd){
//		return nIdx*nSpd + spd;
		return spd*nnodes + nIdx; // use this if it performs faster.
	}
	void extract_halo_data(const T* fOut, const int nnodes);
	void insert_halo_data(T* fOut, const int numSpd);

private:
	// map keys: neighboring partition rank.
	// map values: a Halo Data Object that contains a data structure describing the data to be exchanged

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
void HaloDataOrganizer<T>::extract_halo_data(const T* fOut, const int nnodes)
{
	for(auto & ngbIT : Halo)
	{
		ngbIT.second.extract_halo_data(fOut,nnodes);
	}

}

template <class T>
void HaloDataOrganizer<T>::insert_halo_data(T* fOut, const int nnodes)
{
	for(auto & ngbIT: Halo)
	{
		ngbIT.second.insert_halo_data(fOut,nnodes);
	}
}

template <class T>
void HaloDataOrganizer<T>::allocate_halo_arrays()
{
	for(auto & haloIt : Halo)
	{
		haloIt.second.allocate_arrays();
	}
}

template <class T>
void HaloDataOrganizer<T>::fill_arrays(std::map<int,int> & globalToLocal)
{
	for(auto & haloIt : Halo)
	{
		haloIt.second.fill_nums_and_speeds(globalToLocal);
	}

}

template <class T>
int HaloDataOrganizer<T>::check_ndNums_and_spds(std::map<int,int> & globalToLocal)
{
	int errCt = 0;
	for(auto & haloIt : Halo)
	{
		errCt += haloIt.second.check_ndNums_and_spds(globalToLocal);
	}

	return errCt;
}


template <class T>
std::set<int> HaloDataOrganizer<T>::get_halo_nodes() const
{
	std::set<int> haloSet;
	std::set<int> tempSet;

	for(const auto & haloIt : Halo)
	{
		tempSet = haloIt.second.get_halo_nodes();
		for (const auto & tIt : tempSet)
		{
			haloSet.insert(tIt);
		}

	}

	return haloSet;

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
	return numCuts;

}

template <class T>
int HaloDataOrganizer<T>::get_num_halo_nodes()
{
	int numHalo = 0;
	for (const auto & keyIter : Halo)
	{
		numHalo += keyIter.second.get_num_nodes();
	}
	return numHalo;
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
