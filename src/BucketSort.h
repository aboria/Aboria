/*
 * BucketSort.h
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of RD_3D.
 *
 * RD_3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RD_3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with RD_3D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 22 Oct 2012
 *      Author: robinsonm
 */

#ifndef BUCKETSORT_H_
#define BUCKETSORT_H_

#include <boost/iterator/iterator_facade.hpp>
#include "Vector.h"
#include "Constants.h"
#include "Log.h"
#include <vector>
#include <iostream>
#include <set>


namespace Aboria {

const int CELL_EMPTY = -1;

template<typename T, typename F>
class BucketSort {
public:
	class const_iterator
//	  : public boost::iterator_facade<
//	        const_iterator
//	      , const std::tuple<const typename T::value_type&,const Vect3d&>
//	      , boost::forward_traversal_tag
//	    >
	{
	 public:
		  typedef const std::tuple<const typename T::value_type&,const Vect3d&>* pointer;
		  typedef std::forward_iterator_tag iterator_category;
		  typedef const std::tuple<const typename T::value_type&,const Vect3d&> value_type;
		  typedef const std::tuple<const typename T::value_type&,const Vect3d&> reference;
		  typedef std::ptrdiff_t difference_type;



		  const_iterator()
	      : m_node(),my_index(-1),self(false) {
			  cell_empty.push_back(CELL_EMPTY);
			  m_node = cell_empty.begin();
		  }

	    void go_to_next_candidate() {
	    	if (*m_node != CELL_EMPTY) {
	    		m_node = bucket_sort->linked_list.begin() + *m_node;
	    		//std::cout << "going to new particle *mnode = "<<*m_node<<std::endl;
	    	}
	    	while ((*m_node == CELL_EMPTY) && (surrounding_cell_offset_i != surrounding_cell_offset_end)) {
	    		if (self && (surrounding_cell_offset_i == surrounding_cell_offset_end-1)) {
	    			m_node = cell_i;
	    			while (*m_node != my_index) {
	    				m_node = bucket_sort->linked_list.begin() + *m_node;
	    			}
	    		} else {
		    		//std::cout << "going to new_cell with offset = "<<*surrounding_cell_offset_i<<std::endl;
	    			m_node = cell_i + *surrounding_cell_offset_i;
	    		}
	    		surrounding_cell_offset_i++;
	    	}
//	    	if (surrounding_cell_offset_i == surrounding_cell_offset_end) {
//	    		std::cout <<"finished all cells"<<std::endl;
//	    	} else {
//	    		std::cout <<"*m_node = "<<*m_node<<std::endl;
//	    	}

	    }

	    explicit const_iterator(const BucketSort* bucket_sort,
	    		const Vect3d centre,
	    		const int my_index = -1, const bool self = false)
	    : bucket_sort(bucket_sort),
	      my_index(my_index),self(self),
	      centre(centre) {
	    	cell_empty.push_back(CELL_EMPTY);
	    	m_node = cell_empty.begin();
//	    	if (((centre.array() < bucket_sort->low.array()).any()) ||
//	    			((centre.array() >= bucket_sort->high.array()).any())) {
//	    		return;
//	    	}
	    	cell_i = bucket_sort->cells.begin() + bucket_sort->find_cell_index(centre);
	    	surrounding_cell_offset_i = bucket_sort->surrounding_cell_offsets.begin();
	    	if (self) {
	    		surrounding_cell_offset_end = surrounding_cell_offset_i
	    						+ (bucket_sort->surrounding_cell_offsets.size()-1)/2 + 1;
	    	} else {
	    		//std::cout << "num offsets = "<<bucket_sort->surrounding_cell_offsets.size()<<std::endl;
	    		surrounding_cell_offset_end = bucket_sort->surrounding_cell_offsets.end();
	    	}

	    	increment();

	    }

	    reference operator *() {
	    	return dereference();
	    }
	    reference operator ->() {
	    	return dereference();
	    }
	    const_iterator& operator++() {
	    	increment();
	    	return *this;
	   }
	    const_iterator operator++(int) {
	    	const_iterator tmp(*this);
	    	operator++();
	    	return tmp;
	    }
	    inline bool operator==(const const_iterator& rhs) {
	    	return equal(rhs);
	    }
	    inline bool operator!=(const const_iterator& rhs){
	    	return !operator==(rhs);
	    }

	 private:
	    friend class boost::iterator_core_access;

	    bool equal(const_iterator const& other) const {
	    	//std::cout <<" testing equal *m_node = "<<*m_node<<" other.m_node = "<<*(other.m_node)<<std::endl;
	        return *m_node == *(other.m_node);
	    }

	    void increment() {
	    	//std::cout <<" increment "<<std::endl;
	    	go_to_next_candidate();
	    	while (*m_node != CELL_EMPTY) {
	    		const Vect3d p = bucket_sort->return_vect3d(bucket_sort->begin_iterator[*m_node]);

	    		dx = centre-bucket_sort->correct_position_for_periodicity(centre, p);
	    		//std::cout << "testing candidate with position "<<p<<" and dx = "<<dx<<std::endl;
	    		//if (dx.squaredNorm() <= radius2) {
	    		if ((abs(dx[0]) < bucket_sort->max_interaction_radius) &&
	    				(abs(dx[1]) < bucket_sort->max_interaction_radius) &&
	    				(abs(dx[2]) < bucket_sort->max_interaction_radius)) {

	    	    	//std::cout << "found candidate with position"<<p<<std::endl;
	    	    	break;
	    		} else {
	    			go_to_next_candidate();
	    		}

	    	}
	    }


	    reference dereference() const
	    { return std::tie(bucket_sort->begin_iterator[*m_node],dx); }

	    const BucketSort* bucket_sort;
	    std::vector<int>::const_iterator m_node;
	    std::vector<int>::const_iterator cell_i;
	    //Value* const linked_list;
	    const int my_index;
	    const bool self;
	    const Vect3d centre;
	    Vect3d dx;
	    std::vector<int> cell_empty;
	//    std::vector<Vect3d>::const_iterator positions;
	//    std::vector<Value>::const_iterator linked_list;
	    std::vector<int>::const_iterator surrounding_cell_offset_i,surrounding_cell_offset_end;
	};

	BucketSort(Vect3d low, Vect3d high, Vect3b periodic, F return_vect3d):
		return_vect3d(return_vect3d),
		low(low),high(high),domain_size(high-low),periodic(periodic),
		empty_cell(CELL_EMPTY) {
		LOG(2,"Creating bucketsort data structure with lower corner = "<<low<<" and upper corner = "<<high);
		const double dx = (high-low).maxCoeff()/10.0;
		reset(low, high, dx,periodic);
	}

	void reset(const Vect3d& low, const Vect3d& high, double _max_interaction_radius, const Vect3b& periodic);

	inline const Vect3d& get_low() {return low;}
	inline const Vect3d& get_high() {return high;}
	inline const Vect3b& get_periodic() {return periodic;}

	void embed_points(const T begin_iterator, const T end_iterator);
	const_iterator find_broadphase_neighbours(const Vect3d& r, const int my_index, const bool self) const;
	const_iterator end() const;
	Vect3d correct_position_for_periodicity(const Vect3d& source_r, const Vect3d& to_correct_r) const;
	Vect3d correct_position_for_periodicity(const Vect3d& to_correct_r) const;

private:
	inline int vect_to_index(const Vect3i& vect) const {
		return vect[0] * num_cells_along_yz + vect[1] * num_cells_along_axes[1] + vect[2];
	}
	inline int find_cell_index(const Vect3d &r) const {
		const Vect3i celli = ((r-low).cwiseProduct(inv_cell_size) + Vect3d(1.0,1.0,1.0)).cast<int>();
		ASSERT((celli[0] > 0) && (celli[0] < num_cells_along_axes[0]-1), "position is outside of x-range "<<r);
		ASSERT((celli[1] > 0) && (celli[1] < num_cells_along_axes[1]-1), "position is outside of y-range "<<r);
		ASSERT((celli[2] > 0) && (celli[2] < num_cells_along_axes[2]-1), "position is outside of z-range "<<r);
		//std::cout << " looking in cell " << celli <<" out of total cells " << num_cells_along_axes << " at position " << r<< std::endl;
		return vect_to_index(celli);
	}

	T begin_iterator,end_iterator;
	const F return_vect3d;
    std::vector<int> cells;
    std::vector<std::vector<int> > ghosting_indices_pb;
    std::vector<std::pair<int,int> > ghosting_indices_cb;
    std::vector<int> dirty_cells;
	std::vector<int> linked_list;
	std::vector<int> neighbr_list;
	Vect3d low,high,domain_size;
	Vect3b periodic;
	Vect3d cell_size,inv_cell_size;
	Vect3i num_cells_along_axes;
	int num_cells_along_yz;
	double max_interaction_radius;
	std::vector<int> surrounding_cell_offsets;
	const int empty_cell;
};

template<typename T, typename F>
void BucketSort<T,F>::embed_points(const T _begin_iterator, const T _end_iterator) {

	begin_iterator = _begin_iterator;
	end_iterator = _end_iterator;
	const unsigned int n = std::distance(begin_iterator,end_iterator);
	//std::cout <<"embedding "<<n<<" particles"<<std::endl;
	linked_list.assign(n, CELL_EMPTY);
	if (dirty_cells.size()>0) {
		for (int i: dirty_cells) {
			cells[i] = CELL_EMPTY;
		}
		dirty_cells.clear();
	} else {
		cells.assign(cells.size(), CELL_EMPTY);
	}
	//const bool particle_based = dirty_cells.size() < cells.size();
	const bool particle_based = true; //TODO: fix cell_based neighbour ghosting list
	const bool use_dirty = n < cells.size();
	if (use_dirty) {
		int i = 0;
		for (auto it = begin_iterator; it != end_iterator; ++it, ++i) {
			const int celli = find_cell_index(return_vect3d(*it));
			const int cell_entry = cells[celli];

			// Insert into own cell
			cells[celli] = i;
			dirty_cells.push_back(celli);
			linked_list[i] = cell_entry;

			// Insert into ghosted cells
			if (particle_based) {
				for (int j: ghosting_indices_pb[celli]) {
					const int cell_entry = cells[j];
					cells[j] = i;
					dirty_cells.push_back(j);
				}
			}
		}
	} else {
		int i = 0;
		for (auto it = begin_iterator; it != end_iterator; ++it, ++i) {
			const int celli = find_cell_index(return_vect3d(*it));
			const int cell_entry = cells[celli];

			// Insert into own cell
			cells[celli] = i;
			linked_list[i] = cell_entry;

			// Insert into ghosted cells
			if (particle_based) {
				for (int j: ghosting_indices_pb[celli]) {
					const int cell_entry = cells[j];
					cells[j] = i;
				}
			}
		}
	}


	if (!particle_based) {
		for (std::vector<std::pair<int,int> >::iterator index_pair = ghosting_indices_cb.begin(); index_pair != ghosting_indices_cb.end(); ++index_pair) {
			//BOOST_FOREACH(std::pair<int,int> index_pair, ghosting_indices) {
			cells[index_pair->first] = cells[index_pair->second];
		}
	}
}
template<typename T, typename F>
void BucketSort<T,F>::reset(const Vect3d& _low, const Vect3d& _high, double _max_interaction_radius, const Vect3b& _periodic) {
	LOG(2,"Resetting bucketsort data structure:");
	LOG(2,"\tMax interaction radius = "<<_max_interaction_radius);
	high = _high;
	low = _low;
	domain_size = high-low;
	periodic = _periodic;
	LOG(2,"\tPeriodic = "<<periodic);

	max_interaction_radius = _max_interaction_radius;
	Vect3i num_cells_without_ghost = ((high-low)/max_interaction_radius).cast<int>();
	Vect3d new_high = high;
	for (int i = 0; i < 3; ++i) {
		if (num_cells_without_ghost[i]==0) {
			LOG(2,"\tNote: Dimension "<<i<<" has no length, setting cell side equal to interaction radius.");
			new_high[i] = low[i] + max_interaction_radius;
			num_cells_without_ghost[i] = 1;
		}
	}
	num_cells_along_axes = num_cells_without_ghost + Vect3i(3,3,3);
	LOG(2,"\tNumber of cells along each axis = "<<num_cells_along_axes);
	cell_size = (new_high-low).cwiseQuotient(num_cells_without_ghost.cast<double>());
	LOG(2,"\tCell sizes along each axis = "<<cell_size);
	inv_cell_size = Vect3d(1,1,1).cwiseQuotient(cell_size);
	num_cells_along_yz = num_cells_along_axes[2]*num_cells_along_axes[1];
	const unsigned int num_cells = num_cells_along_axes.prod();
	cells.assign(num_cells, CELL_EMPTY);
	dirty_cells.clear();
	//TODO: assumed 3d
	surrounding_cell_offsets.clear();
	for (int i = -1; i < 2; ++i) {
		for (int j = -1; j < 2; ++j) {
			for (int k = -1; k < 2; ++k) {
				surrounding_cell_offsets.push_back(vect_to_index(Vect3i(i,j,k)));
			}
		}
	}


	ghosting_indices_pb.assign(num_cells, std::vector<int>());
	ghosting_indices_cb.clear();
	for (int i = 0; i < NDIM; ++i) {
		if (!periodic[i]) continue;
		int j,k;
		switch (i) {
			case 0:
				j = 1;
				k = 2;
				break;
			case 1:
				j = 0;
				k = 2;
				break;
			case 2:
				j = 0;
				k = 1;
				break;
			default:
				break;
		}

		Vect3i tmp;
		const int n = num_cells_along_axes[i];
		for (int jj = 0; jj < num_cells_along_axes[j]-1; ++jj) {
			tmp[j] = jj;
			for (int kk = 0; kk < num_cells_along_axes[k]-1; ++kk) {
				tmp[k] = kk;
				tmp[i] = n-3;
				const int index_from1 = vect_to_index(tmp);
				tmp[i] = 0;
				const int index_to1 = vect_to_index(tmp);
				ghosting_indices_pb[index_from1].push_back(index_to1);
				//ghosting_indices_cb.push_back(std::pair<int,int>(index_to1,index_from1));
				tmp[i] = 1;
				const int index_from2 = vect_to_index(tmp);
				tmp[i] = n-2;
				const int index_to2 = vect_to_index(tmp);
				ghosting_indices_pb[index_from2].push_back(index_to2);
				//ghosting_indices_cb.push_back(std::pair<int,int>(index_to2,index_from2));
			}
		}
	}
	/*
	 * collapse redirections
	 */
	for (int i = 0; i < num_cells; ++i) {
		std::set<int> ghosting_cells;
		for (int j: ghosting_indices_pb[i]) {
			ghosting_cells.insert(j);
			for (int k: ghosting_indices_pb[j]) {
				ghosting_cells.insert(k);
				for (int m: ghosting_indices_pb[k]) {
					ghosting_cells.insert(m);
				}
			}
		}
		ghosting_indices_pb[i].resize(ghosting_cells.size());
		std::copy(ghosting_cells.begin(),ghosting_cells.end(),ghosting_indices_pb[i].begin());
	}
}
template<typename T, typename F>
typename BucketSort<T,F>::const_iterator BucketSort<T,F>::find_broadphase_neighbours(const Vect3d& r,const int my_index, const bool self) const {
	return const_iterator(this,r,my_index,self);
}
template<typename T, typename F>
typename BucketSort<T,F>::const_iterator BucketSort<T,F>::end() const {
	return const_iterator();
}
template<typename T, typename F>
Vect3d BucketSort<T,F>::correct_position_for_periodicity(const Vect3d& source_r, const Vect3d& to_correct_r) const {
	Vect3d corrected_r = to_correct_r - source_r;
	for (int i = 0; i < NDIM; ++i) {
		if (!periodic[i]) continue;
		if (corrected_r[i] > domain_size[i]/2.0) corrected_r[i] -= domain_size[i];
		else if (corrected_r[i] < -domain_size[i]/2.0) corrected_r[i] += domain_size[i];
	}
	return corrected_r + source_r;
}

template<typename T, typename F>
Vect3d BucketSort<T,F>::correct_position_for_periodicity(const Vect3d& to_correct_r) const {
	Vect3d corrected_r = to_correct_r;
	for (int i = 0; i < NDIM; ++i) {
		if (!periodic[i]) continue;
		while (corrected_r[i] >= high[i]) corrected_r[i] -= domain_size[i];
		while (corrected_r[i] < low[i]) corrected_r[i] += domain_size[i];
	}
	return corrected_r;
}



}

#endif /* BUCKETSORT_H_ */
