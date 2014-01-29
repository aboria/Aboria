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

#include "Vector.h"
#include "Constants.h"
#include "Log.h"
#include <vector>
#include <iostream>

#include <boost/foreach.hpp>

namespace Tyche {

const int CELL_EMPTY = -1;

template <class Value>
class node_iter
  : public boost::iterator_facade<
        node_iter<Value>
      , Value
      , boost::forward_traversal_tag
    >
{
 public:
    node_iter()
      : m_node(0) {}

    explicit node_iter(Value* cell_i,
    		Value* linked_list,
    		const std::vector<int>& surrounding_cell_offsets,
    		const int my_index = -1, const bool self = false)
    : cell_i(cell_i),
      linked_list(linked_list),
      my_index(my_index),self(self) {
    	surrounding_cell_offset_i = surrounding_cell_offsets.begin();
    	if (self) {
    		surrounding_cell_offset_end = surrounding_cell_offset_i
    						+ (surrounding_cell_offsets.size()-1)/2;
    	} else {
    		surrounding_cell_offset_end = surrounding_cell_offsets.end();
    	}
    	m_node = cell_i + *surrounding_cell_offset_i;
    }

 private:
    friend class boost::iterator_core_access;

    bool equal(node_iter<Value> const& other) const {
        return *m_node == *(other.m_node);
    }

    void increment() {
    	m_node = linked_list[*m_node];
    	if (*m_node == CELL_EMPTY) surrounding_cell_offset_i++;
    	if (surrounding_cell_offset_i != surrounding_cell_offset_end) {
    		m_node = cell_i + *surrounding_cell_offset_i;
    	} else if (self) {
    		m_node = linked_list[my_index];
    	}
    }

    Value& dereference() const
    { return *m_node; }

    Value* m_node;
    Value* cell_i;
    Value* const linked_list;
    const int my_index;
    const bool self;
    std::vector<int>::iterator surrounding_cell_offset_i,surrounding_cell_offset_end;
};

class BucketSort {
public:
	typedef node_iter<int const> const_iterator;

	BucketSort(Vect3d low, Vect3d high, Vect3b periodic):
		low(low),high(high),domain_size(high-low),periodic(periodic),
		empty_cell(CELL_EMPTY) {
		LOG(2,"Creating bucketsort data structure with lower corner = "<<low<<" and upper corner = "<<high);
		const double dx = (high-low).maxCoeff()/10.0;
		reset(low, high, dx);
	}

	void reset(const Vect3d& low, const Vect3d& high, double _max_interaction_radius);
	inline const Vect3d& get_low() {return low;}
	inline const Vect3d& get_high() {return high;}

	template<typename T, typename F>
	void embed_points(const T& begin_iterator, const T& end_iterator, const F& return_vect3d);
	const_iterator find_broadphase_neighbours(const Vect3d& r, const int my_index, const bool self);
    const_iterator end();
	Vect3d correct_position_for_periodicity(const Vect3d& source_r, const Vect3d& to_correct_r);
	Vect3d correct_position_for_periodicity(const Vect3d& to_correct_r);

private:
	inline int vect_to_index(const Vect3i& vect) {
		return vect[0] * num_cells_along_yz + vect[1] * num_cells_along_axes[1] + vect[2];
	}
	inline int find_cell_index(const Vect3d &r) {
		const Vect3i celli = ((r-low).cwiseProduct(inv_cell_size) + Vect3d(1.0,1.0,1.0)).cast<int>();
		ASSERT((celli[0] >= 0) && (celli[0] < num_cells_along_axes[0]), "position is outside of x-range");
		ASSERT((celli[1] >= 0) && (celli[1] < num_cells_along_axes[1]), "position is outside of y-range");
		ASSERT((celli[2] >= 0) && (celli[2] < num_cells_along_axes[2]), "position is outside of z-range");
		//std::cout << " looking in cell " << celli <<" out of total cells " << num_cells_along_axes << " at position " << r<< std::endl;
		return vect_to_index(celli);
	}

    std::vector<int> cells;
    std::vector<std::vector<int> > ghosting_indices_pb;
    std::vector<std::pair<int,int> > ghosting_indices_cb;
    std::vector<int> dirty_cells;
	std::vector<int> linked_list;
	std::vector<int> neighbr_list;
	Vect3d low,high,domain_size;
	const Vect3b periodic;
	Vect3d cell_size,inv_cell_size;
	Vect3i num_cells_along_axes;
	int num_cells_along_yz;
	double max_interaction_radius;
	std::vector<int> surrounding_cell_offsets;
	const int empty_cell;
};

template<typename T, typename F>
void BucketSort::embed_points(const T& begin_iterator, const T& end_iterator, const F& return_vect3d) {
	const unsigned int n = std::distance(begin_iterator,end_iterator);
	linked_list.assign(n, CELL_EMPTY);
	const bool particle_based = dirty_cells.size() < cells.size();
	if (particle_based) {
		BOOST_FOREACH(int i, dirty_cells) {
			cells[i] = CELL_EMPTY;
		}
	} else {
		cells.assign(cells.size(), CELL_EMPTY);
	}

	dirty_cells.clear();
	int i = 0;
	for (auto it = begin_iterator; it != end_iterator; ++it, ++i) {
		const int celli = find_cell_index(return_vect3d(it));
		const int cell_entry = cells[celli];

		// Insert into own cell
		cells[celli] = i;
		dirty_cells.push_back(celli);
		linked_list[i] = cell_entry;

		// Insert into ghosted cells
		if (particle_based) {
			BOOST_FOREACH(int j, ghosting_indices_pb[celli]) {
				const int cell_entry = cells[j];
				cells[j] = i;
				dirty_cells.push_back(j);
				linked_list[i] = cell_entry;
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


}

#endif /* BUCKETSORT_H_ */
