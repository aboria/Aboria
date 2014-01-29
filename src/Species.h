/*
 * species.h
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
 *  Created on: 11 Oct 2012
 *      Author: robinsonm
 */

#ifndef SPECIES_H_
#define SPECIES_H_

#include <vector>
#include <boost/foreach.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include "MyRandom.h"
#include "Vector.h"
#include "StructuredGrid.h"
#include "BucketSort.h"

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>



namespace Tyche {

#define DATA_typename ParticleInfo
#define DATA_names   (r)(r0)(alive)(id)(saved_index)
#define DATA_types   (Vect3d)(Vect3d)(bool)(int)(int)
#include "Data.h"

template<int DataSize = 1, typename NeighbourSearch = BucketSort>
class Particles {
public:
	class get_pos_and_clear_dirty;
	class Value {
		friend class Particles;
		friend class get_pos_and_clear_dirty;
	public:
		Value(const Vect3d& r, const Vect3d& r0, const bool alive, const size_t id, const size_t saved_index):
			r(r),
			r0(r0),
			alive(alive),dirty(true),
			id(id),
			saved_index(saved_index)
		{}
		Value& operator=(const Value &rhs) {
			if (this != &rhs) {
				r = rhs.r;
				r0 = rhs.r0;
				alive = rhs.alive;
				id = rhs.id;
				saved_index = rhs.saved_index;
				data = rhs.data;
			}
			return *this;
		}
		const Vect3d& get_position() const {
			return r;
		}
		Vect3d& get_position_nonconst() {
			dirty = true;
			return r;
		}
		const Vect3d& get_old_position() const {
			return r0;
		}
		Vect3d& get_old_position_nonconst() {
			return r0;
		}
		const double* get_data() const {
			return data.data();
		}
		double* get_data(const size_t i) {
			return data.data();
		}
		bool is_alive() const {
			return alive;
		}
		const std::size_t& get_id() const {
			return id;
		}
		const std::size_t& get_saved_index() const {
			return saved_index;
		}
		void mark_for_deletion() {
			alive = false;
		}
		template<typename T>
		T::neighbour_iterator get_in_radius(const T& particles, const double radius) {
			return particles.get_in_geometry(
					particles.neighbour_search.find_broadphase_neighbours(get_position(), radius, index,false),
					particles.neighbour_search.end(),
					Sphere(get_position(),radius)
			);
		}
	private:
		Vect3d r,r0;
		bool alive,dirty;
		std::size_t id;
		std::size_t index,saved_index;
		std::vector<size_t> neighbour_indicies;
		boost::array<double,DataSize> data;
	};

	typedef typename std::vector<Value>::iterator iterator;
	class get_pos_and_clear_dirty {
		const Vect3d& operator()(iterator i) {
			i->dirty = false;
			return i->get_position();
		}
	};
	const int SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE = -1;


	Particles():
		next_id(0),
		neighbour_search(Vect3d(0,0,0),Vect3d(1,1,1),Vect3b(false,false,false))
	{}

	Value& operator[](std::size_t idx) {
		return data[idx];
	};
	const Value& operator[](std::size_t idx) const {
		return const_cast<Value>(*this)[idx];
	};
	typename std::vector<Value>::iterator begin() {
		return data.begin();
	}
	typename std::vector<Value>::iterator end() {
		return data.end();
	}
//	template <typename T>
//	typename neighbour_search_iterator<T> begin_neighbour_search(const T* particles_to_search) {
//		//return data.begin();
//	}
	void fill_uniform(const Vect3d low, const Vect3d high, const unsigned int N) {
		//TODO: assumes a 3d rectangular region
		boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, boost::uniform_real<>(0,1));
		const Vect3d dist = high-low;
		for(int i=0;i<N;i++) {
			add_particle(Vect3d(uni()*dist[0],uni()*dist[1],uni()*dist[2])+low);
		}
	}

	void delete_particles() {
		data.erase (std::remove_if( std::begin(data), std::end(data),
				[](Value& p) { return !p.is_alive(); }),
				std::end(data)
		);

	}
	void clear() {
		data.clear();
	}
	size_t size() const {
		return data.size();
	}
	void add_particle(const Vect3d& position) {
		data.push_back(Value(position,position,true, next_id++, SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE));
	}
	void add_particle(const Vect3d& position, const Vect3d& old_position) {
		data.push_back(Value(position,old_position,true, next_id++, SPECIES_SAVED_INDEX_FOR_NEW_PARTICLE));
	}
	const Vect3d& get_position(const size_t i) const {
		return data[i].get_position();
	}
	Vect3d& get_position_nonconst(const size_t i) {
		return data[i].get_position_nonconst();
	}
	const Vect3d& get_old_position(const size_t i) const {
		return data[i].get_old_position();
	}
	Vect3d& get_old_position_nonconst(const size_t i) {
		return data[i].get_old_position_nonconst();
	}
	const double* get_data(const size_t i) const {
		return data[i].get_data();
	}
	double* get_data(const size_t i) {
		return data[i].get_data_nonconst();
	}
	bool is_alive(const size_t i) const {
		return data[i].is_alive();
	}
	const size_t& get_id(const size_t i) const {
		return data[i].get_id();
	}

	void mark_for_deletion(const size_t i) {
		data[i].mark_for_deletion();
	}

	void save_indicies() {
		const int n = data.size();
		for (int i = 0; i < n; ++i) {
			data[i].saved_index = i;
		}
	}

	void init_neighbour_search(const Vect3d& low, const Vect3d& high, const double length_scale) {
		neighbour_search.reset(low,high,length_scale);
		neighbour_search.embed_points(begin(),end(),get_pos_and_clear_dirty());
	}
	void refresh_neighbour_search() {
		neighbour_search.embed_points(begin(),end(),get_pos_and_clear_dirty());
	}

	vtkSmartPointer<vtkUnstructuredGrid> get_vtk_grid() {
		vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkIntArray> newInt = vtkSmartPointer<vtkIntArray>::New();
		newInt->SetName("id");
		const vtkIdType n = size();
		newPts->SetNumberOfPoints(n);
		newInt->SetNumberOfValues(n);
		for (int i = 0; i < n; ++i) {
			//std::cout << "adding mol to vtk at position "<<mols.r[i]<<std::endl;
			newPts->SetPoint(i,get_position(i)[0],get_position(i)[1],get_position(i)[2]);
			newInt->SetValue(n,get_id(i));
		}
		newPts->ComputeBounds();

		grid->SetPoints(newPts);
		grid->GetPointData()->AddArray(newInt);

		return grid;
	}
private:
	std::vector<Value> data;
	NeighbourSearch neighbour_search;
	bool dirty;
	int next_id;
};


static const StructuredGrid empty_grid;

class Species {
public:
	Species(double D):D(D),grid(NULL) {
		id = species_count++;
		clear();
	}
	Species(double D, const StructuredGrid* grid):D(D),grid(grid)  {
		id = species_count++;
		clear();
	}
	static std::auto_ptr<Species> New(double D) {
		return std::auto_ptr<Species>(new Species(D));
	}
	void clear() {
		mols.clear();
		if (grid!=NULL) copy_numbers.assign(grid->size(),0);
	}
	void set_grid(const StructuredGrid* new_grid) {
		grid = new_grid;
		if (grid!=NULL) copy_numbers.assign(grid->size(),0);
	}
	//void fill_uniform(const Vect3d low, const Vect3d high, const unsigned int N);
	template<typename T>
	void fill_uniform(const T& geometry, const unsigned int N);
	template<typename T>
	unsigned int count_mols_in(const T& geometry);
	void get_concentrations(const StructuredGrid& calc_grid, std::vector<double>& mol_concentrations, std::vector<double>& compartment_concentrations) const;
	void get_concentration(const StructuredGrid& calc_grid, std::vector<double>& concentration) const;
	void get_concentration(const Vect3d low, const Vect3d high, const Vect3i n, std::vector<double>& concentration) const;
	vtkSmartPointer<vtkUnstructuredGrid> get_vtk();
	std::string get_status_string() const;
	friend std::ostream& operator<<( std::ostream& out, const Species& b ) {
		return out << b.get_status_string();
	}

	double D;
	Particles<1> mols;
	std::vector<int> copy_numbers;
	std::vector<int> mol_copy_numbers;
	const StructuredGrid* grid;
	int id;
	std::vector<double> tmpx,tmpy,tmpz;
private:
	static int species_count;
};

template<typename T>
void Species::fill_uniform(const T& geometry, const unsigned int N) {
	LOG(2,"Adding "<<N<<" molecules of Species ("<<id<<") within "<<geometry);
	for(int i=0;i<N;i++) {
		mols.add_particle(geometry.get_random_point_in());
	}
}

template<typename T>
unsigned int Species::count_mols_in(const T& geometry) {
	unsigned int count;
	const int n = mols.size();
	for (int i = 0; i < n; ++i) {
		if (geometry.is_in(mols.get_position(i))) count++;
	}
	return count;
}


extern Species null_species;
}


#endif /* SPECIES_H_ */
