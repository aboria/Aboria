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

#ifndef PARTICLES_H_
#define PARTICLES_H_

#include "BucketSort.h"
#include <vector>
#include <random>
#include <string>
//#include <boost/array.hpp>
//#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include "Zip.h"
#include "Vector.h"
//#include "MyRandom.h"

#ifndef HAVE_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#endif

namespace Aboria {


template<typename DataType>
class DataNames {
public:
	std::string static get(unsigned int n) {
		return "data_" + std::to_string(n);
	}
};

template<int I,typename ParticlesType>
struct Elem {
	typedef typename std::tuple_element<I,typename ParticlesType::data_type>::type type;
};

template<typename ParticlesType>
struct Elem<POSITION, ParticlesType> {
	typedef Vect3d type;
};

template<typename ParticlesType>
struct Elem<ID, ParticlesType> {
	typedef std::size_t type;
};

template<typename ParticlesType>
struct Elem<ALIVE, ParticlesType> {
	typedef std::size_t type;
};





template<typename DataType=std::tuple<double> >
class Particles {
	template<typename T>
	friend class Particles;
public:

	class value_type {
		friend class Particles;
	public:
		typedef std::mt19937 generator_type;
		value_type():uni(0,1),normal(0,1){}
		value_type(const value_type& rhs) {
			r = rhs.r;
			data = rhs.data;
			alive = rhs.alive;
			id = rhs.id;
			generator = rhs.generator;
		}
		value_type(const Vect3d &position):
			r(position),uni(0,1),normal(0,1) {}
		~value_type() {

		}
		value_type& operator=(const value_type &rhs) {
			if (this != &rhs) {
				data = rhs.dattypea;
			}
			return *this;
		}

		bool operator==(const value_type &rhs) const {
			return id == rhs.id;
		}

		void deep_copy(const value_type &rhs) {
			if (this != &rhs) {
				r = rhs.r;
				data = rhs.data;
				alive = rhs.alive;
				id = rhs.id;
				generator = rhs.generator;
			}
		}

//		Vect3d& get_position() {
//			return r;
//		}
		const Vect3d& get_position() const {
			return r;
		}
		void set_position(const Vect3d& arg) {
			r = arg;
		}
		const DataType& get_data() const {
			return data;
		}
		DataType& get_data() {
			return data;
		}
//		template<int N>
//		std::tuple_element<N,DataType>::type& get_data_elem() {
//			return std::get<N>(data);
//		}
		template<int I>
		const typename Elem<I, Particles<DataType> >::type& get_elem() const {
			return std::get<I>(data);
		}
		template<int I>
		typename Elem<I, Particles<DataType> >::type& get_elem() {
			return std::get<I>(data);
		}
		template<int I>
		void set_elem(const typename Elem<I, Particles<DataType> >::type& arg) {
			std::get<I>(data) = arg;
		}

		size_t get_index() {
			return index;
		}
		bool get_alive() const {
			return alive;
		}
		void set_alive(bool aliveIn) const {
			alive = aliveIn;
		}
		generator_type get_generator() {
			return generator;
		}
		double rand_uniform() {
			return uni(generator);
		}

		double rand_normal() {
			return normal(generator);
		}
		const std::size_t get_id() const {
			return id;
		}
		void mark_for_deletion() {
			alive = false;
		}

		template<typename T>
		boost::iterator_range<typename T::element_type::NeighbourSearch_type::const_iterator> get_neighbours(const T particles) {
			ASSERT(particles->searchable==true,"particles is not searchable");
			return boost::make_iterator_range(
			 particles->neighbour_search.find_broadphase_neighbours(get_position(), index,false),
			 particles->neighbour_search.end());
		}

		template<typename T>
		Vect3d correct_position_for_periodicity(const T particles, const Vect3d& position) {
			return particles->neighbour_search.correct_position_for_periodicity(r, position);
		}
	private:

		Vect3d r;
		bool alive;
		std::size_t id,index;
		DataType data;
		generator_type generator;
		std::uniform_real_distribution<double> uni;
		std::normal_distribution<double> normal;

	};
	typedef DataType data_type;
	typedef typename std::vector<value_type> vector_type;
	typedef typename vector_type::size_type size_type;
	typedef typename vector_type::size_type difference_type;
	typedef typename vector_type::iterator iterator;
	typedef typename vector_type::const_iterator const_iterator;
	struct get_pos {
		const Vect3d& operator()(const value_type& i) const {
			return i.get_position();
		}
	};
	typedef BucketSort<const_iterator,get_pos> NeighbourSearch_type;


	Particles():
		next_id(0),
		neighbour_search(Vect3d(0,0,0),Vect3d(1,1,1),Vect3b(false,false,false),get_pos()),
		searchable(false),track_ids(false),
		seed(time(NULL))
	{}

	Particles(const double seed):
		next_id(0),
		neighbour_search(Vect3d(0,0,0),Vect3d(1,1,1),Vect3b(false,false,false),get_pos()),
		searchable(false),track_ids(false),
		seed(seed)
	{}

	Particles(const Particles<DataType> &other):
			data(other.data),
			neighbour_search(other.neighbour_search),
			next_id(other.next_id),
			searchable(other.searchable),
			track_ids(other.track_ids),
			seed(other.seed),
			id_to_index(id_to_index)
	{}
	Particles(iterator first, iterator last):
		data(first,last),
		neighbour_search(Vect3d(0,0,0),Vect3d(1,1,1),Vect3b(false,false,false),get_pos()),
		searchable(false),
		track_ids(false),
		seed(0)
	{}

	static ptr<Particles<DataType> > New() {
		return ptr<Particles<DataType> >(new Particles<DataType> ());
	}

	value_type& operator[](std::size_t idx) {
		return data[idx];
	}
	const value_type& operator[](std::size_t idx) const {
		return const_cast<value_type>(*this)[idx];
	}
	iterator begin() {
		return data.begin();
	}
	iterator end() {
		return data.end();
	}


	void delete_particles() {
		const int n = data.size();
		for (int index = 0; index < n; ++index) {
			value_type& i = data[index];
			if (i.alive==false) {
				i.deep_copy(*(data.cend()-1));
				if (track_ids) id_to_index[i.id] = index;
				data.pop_back();
			}
		}
		if (searchable) neighbour_search.embed_points(data.cbegin(),data.cend());
	}
	void clear() {
		data.clear();
	}

	iterator erase (iterator i) {
		if (i != end()-1) {
			i->deep_copy(*(data.cend()-1));
			if (track_ids) id_to_index[i->id] = i-begin();
			if (searchable) neighbour_search.copy_points(i,end());

		}
		if (searchable) neighbour_search.delete_point(end());
		data.pop_back();
	}

	iterator erase (iterator first, iterator last) {
		for(iterator i=first;i!=last;i++) {
			erase(i);
		}
	}

	iterator insert (iterator position, const value_type& val) {
		data.insert(position,val);
	}

	void insert (iterator position, size_type n, const value_type& val) {
		data.insert(position,n,val);
	}

	template <class InputIterator>
	void insert (iterator position, InputIterator first, InputIterator last) {
		data.insert(position,first,last);
	}


	size_t size() const {
		return data.size();
	}

	void save_indicies() {
		const int n = data.size();
		for (int i = 0; i < n; ++i) {
			data[i].saved_index = i;
		}
	}

	void init_neighbour_search(const Vect3d& low, const Vect3d& high, const double length_scale, const Vect3b& periodic) {
		neighbour_search.reset(low,high,length_scale,periodic);
		neighbour_search.embed_points(data.cbegin(),data.cend());
		searchable = true;
	}
	
	void set_track_ids(bool set) {
		if (!track_ids && set) {
			const int n = data.size();
			for (int i = 0; i < n; ++i) {
				id_to_index[data[i].id] = i;
			}
		} else if (!set) {
			id_to_index.clear();
		}
	}

	value_type* find_id(const int id) {
		std::map<size_t,size_t>::iterator it = id_to_index.find(id);
		if (it==id_to_index.end()) {
			return NULL;
		} else {
			return &(data[*it]);
		}
	}

	template<typename F>
	void reset_neighbour_search(const double length_scale, F f) {
		std::for_each(begin(),end(),[&f](value_type& i) {
			i.r = f(i);
		});
		neighbour_search.reset(neighbour_search.get_low(),neighbour_search.get_high(),length_scale,neighbour_search.get_periodic());
		neighbour_search.embed_points(data.cbegin(),data.cend());
		searchable = true;
	}

	void reset_neighbour_search(const double length_scale) {
		neighbour_search.reset(neighbour_search.get_low(),neighbour_search.get_high(),length_scale,neighbour_search.get_periodic());
		neighbour_search.embed_points(data.cbegin(),data.cend());
		searchable = true;
	}

	const double get_lengthscale() const {
		return neighbour_search.get_lengthscale();
	}
	const Vect3d& get_low() const {
		return neighbour_search.get_low();
	}
	const Vect3d& get_high() const {
		return neighbour_search.get_high();
	}



	void push_back (const value_type& val) {
		data.push_back(val);
		if (searchable) neighbour_search.update_begin_and_end(data.cbegin(),data.cend());
		const int index = data.size();
		iterator i = end()-1;
		i->r = val.r;
		//if (std::get<1>(i->data).norm()>1) std::cout <<"adding bad orientation"<<std::endl;
		i->id = this->next_id++;
		i->generator.seed(i->id*seed);
		i->alive = true;
		i->index = index;
		if (track_ids) id_to_index[i->id] = index;
		if (searchable) neighbour_search.add_point(i);
	}

	void push_back(const Vect3d& position) {
		this->push_back(value_type(position));
	}

	void pop_back() {
		erase(end()-1);
	}

	template<typename F>
	void create_particles(const int n, F f) {
		const int old_size = data.size();
		data.resize(old_size+n);
		if (searchable) neighbour_search.update_begin_and_end(data.cbegin(),data.cend());
		int index = old_size;
		for (auto i=data.begin()+old_size; i!=data.end();i++,index++) {
			i->id = this->next_id++;
			i->generator.seed(i->id*seed);
			i->alive = true;
			i->index = index;
			i->r = f(*i);
			if (track_ids) id_to_index[i->id] = index;
			if (searchable) neighbour_search.add_point(i);
		}
	}

	boost::iterator_range<typename NeighbourSearch_type::const_iterator> get_neighbours(const Vect3d position) {
		return boost::make_iterator_range(
				neighbour_search.find_broadphase_neighbours(position, -1,false),
				neighbour_search.end());
	}

	template<typename F>
	void create_particles_cylinder(const Vect3d& min, const Vect3d& max, const int n, double dem_diameter, F f) {
		int shift=0;
		bool shiftz=false;
		double radius=(max[0]-min[0])/2-dem_diameter/2;
		const Vect3d origin(radius,radius,dem_diameter/2);
		int jmax=radius/(dem_diameter); //number of particles on radius
		for (int i = 0; i < n; ++i) { //z index
			radius=(max[0]-min[0])/2-dem_diameter/2;
			for (int j = 0; j < jmax; ++j) { //radius index
				radius=radius-dem_diameter;
				int kmax= radius*2*PI/dem_diameter; //number of particles on circumference

				double angle=shift*30+180*shiftz; // shift and shiftz avoid that vacuum is concentrated in certain point

				int index = data.size();
				data.resize(data.size()+kmax);
				if (searchable) neighbour_search.update_begin_and_end(data.cbegin(),data.cend());
				for (int k = 0; k < kmax; ++k) {
					double d_angle=dem_diameter/radius;
					angle=angle+d_angle;
					data[index].r[0] = origin[0] + radius*cos(angle);
					data[index].r[1] = origin[1] + radius*sin(angle);
					data[index].r[2] = origin[2] + dem_diameter*i;
					data[index].id = next_id++;
					data[index].index = index;

					data[index].generator.seed(data[index].id*seed);
					if (track_ids) id_to_index[data[index].id] = index;

					data[index].alive = true;
					f(data[index]);
					index++;
				}
				shift++;
			}
			shiftz=!shiftz;
		}
		if (searchable) neighbour_search.embed_points(data.cbegin(),data.cend());
	}


	template<typename F>
	void create_particles_grid(const Vect3d& min, const Vect3d& max, const Vect3i& n, F f) {
		const int nparticles = n.prod();
		int index = data.size();
		data.resize(data.size()+nparticles);
		if (searchable) neighbour_search.update_begin_and_end(data.cbegin(),data.cend());
		const Vect3d dx = (max-min)/n;
		for (int i = 0; i < n[0]; ++i) {
			for (int j = 0; j < n[1]; ++j) {
				for (int k = 0; k < n[2]; ++k) {
					data[index].r = min + Vect3d(i+0.5,j+0.5,k+0.5)*dx;
					data[index].id = next_id++;
					data[index].generator.seed(data[index].id*seed);
					data[index].alive = true;
					data[index].index = index;

					if (track_ids) id_to_index[data[index].id] = index;
					f(data[index]);

					index++;
				}
			}
		}
		if (searchable) neighbour_search.embed_points(data.cbegin(),data.cend());

	}

	template<typename F>
	void update_positions(iterator b, iterator e, F f) {
		std::for_each(b,e,[&f](value_type& i) {
			i.r = f(i);
		});
		if (searchable) {
			enforce_domain(neighbour_search.get_low(),neighbour_search.get_high(),neighbour_search.get_periodic());
			neighbour_search.embed_points(data.cbegin(),data.cend());
		}
	}
	template<typename F>
	void update_positions(F f) {
		std::for_each(begin(),end(),[&f](value_type& i) {
			i.r = f(i);
		});
		if (searchable) {
			enforce_domain(neighbour_search.get_low(),neighbour_search.get_high(),neighbour_search.get_periodic());
			neighbour_search.embed_points(data.cbegin(),data.cend());
		}
	}

	void update_positions() {
		if (searchable) {
			enforce_domain(neighbour_search.get_low(),neighbour_search.get_high(),neighbour_search.get_periodic());
			neighbour_search.embed_points(data.cbegin(),data.cend());
		}
	}

	template<typename F>
	void update_positions_sequential(iterator b, iterator e, F f) {
		const Vect3d low = neighbour_search.get_low();
		const Vect3d high = neighbour_search.get_high();
		const Vect3b periodic = neighbour_search.get_periodic();
		for(iterator i = b; i != e; i++) {
			i->r = f(*i);
			for (int d = 0; d < 3; ++d) {
				if (periodic[d]) {
					while (i->r[d]<low[d]) {
						i->r[d] += (high[d]-low[d]);
					}
					while (i->r[d]>=high[d]) {
						i->r[d] -= (high[d]-low[d]);
					}
				} else {
					if ((i->r[d]<low[d]) || (i->r[d]>=high[d])) {
						i->mark_for_deletion();
					}
				}
			}
			neighbour_search.update_point(i);
		}
	}


#ifndef HAVE_VTK
	struct save_elem {
		typedef int result_type;
		save_elem(int index, vtkSmartPointer<vtkFloatArray>* datas):
			index(index),datas(datas){}
		template<typename T>
		result_type operator()(result_type state, T& t) const {
			datas[state]->SetValue(index,t);
			return state + 1;
		}
		//template<>
		result_type operator()(result_type state, Vect3d& t) const {
			datas[state]->SetTuple3(index,t[0],t[1],t[2]);
			return state + 1;
		}
		result_type operator()(result_type state, const Vect3d& t) const {
			datas[state]->SetTuple3(index,t[0],t[1],t[2]);
			return state + 1;
		}
		int index;
		vtkSmartPointer<vtkFloatArray>* datas;
	};
	struct read_elem {
			typedef int result_type;
			read_elem(int index, vtkSmartPointer<vtkFloatArray>* datas):
				index(index),datas(datas){}
			template<typename T>
			result_type operator()(result_type state, T& t) const {
				t = datas[state]->GetValue(index);
				return state + 1;
			}
			//template<>
			result_type operator()(result_type state, Vect3d& t) const {
				double *data = datas[state]->GetTuple3(index);
				t[0] = data[0];
				t[1] = data[1];
				t[2] = data[2];
				return state + 1;
			}
			int index;
			vtkSmartPointer<vtkFloatArray>* datas;
	};
	struct setup_data {
		typedef int result_type;
		setup_data(vtkSmartPointer<vtkFloatArray>* datas):
			datas(datas){}

		template<typename T>
		result_type operator()(result_type state, T& t) const {
			datas[state]->SetNumberOfComponents(1);
			return state + 1;
		}
		result_type operator()(result_type state, const Vect3d& t) const {
			datas[state]->SetNumberOfComponents(3);
			return state + 1;
		}
		result_type operator()(result_type state, Vect3d& t) const {
			datas[state]->SetNumberOfComponents(3);
			return state + 1;
		}
		//template<>

		vtkSmartPointer<vtkFloatArray>* datas;
	};

	vtkSmartPointer<vtkUnstructuredGrid> get_grid() {
		if (!cache_grid) cache_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		copy_to_vtk_grid(cache_grid);
		return cache_grid;
	}

	void  copy_to_vtk_grid(vtkUnstructuredGrid *grid) {
		vtkSmartPointer<vtkPoints> points = grid->GetPoints();
		if (!points) {
			points = vtkSmartPointer<vtkPoints>::New();
			grid->SetPoints(points);
		}
		vtkSmartPointer<vtkCellArray> cells = grid->GetCells();
		if (!cells) {
			cells = vtkSmartPointer<vtkCellArray>::New();
			grid->SetCells(1,cells);
		}
		vtkSmartPointer<vtkUnsignedCharArray> cell_types = grid->GetCellTypesArray();
		vtkSmartPointer<vtkIntArray> ids = vtkIntArray::SafeDownCast(grid->GetPointData()->GetArray("id"));
		if (!ids) {
			ids = vtkSmartPointer<vtkIntArray>::New();
			ids->SetName("id");
			grid->GetPointData()->AddArray(ids);
		}
		constexpr size_t dn = std::tuple_size<DataType>::value;
		vtkSmartPointer<vtkFloatArray> datas[dn];
		for (int i = 0; i < dn; ++i) {
			std::string name = DataNames<DataType>::get(i);
			datas[i] = vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name.c_str()));
			if (!datas[i]) {
				datas[i] = vtkSmartPointer<vtkFloatArray>::New();
				datas[i]->SetName(name.c_str());
				grid->GetPointData()->AddArray(datas[i]);
			}
		}
		boost::fusion::fold(DataType(),0,setup_data(datas));

		const vtkIdType n = size();
		points->SetNumberOfPoints(n);
		cells->Reset();
		cell_types->Reset();
		ids->SetNumberOfTuples(n);
		for (int i = 0; i < dn; ++i) {
			datas[i]->SetNumberOfTuples(n);
		}
		int j = 0;

		for(auto& i: *this) {
			const int index = j++;
			//std::cout <<"copying point at "<<i.get_position()<<" with id = "<<i.get_id()<<std::endl;
			points->SetPoint(index,i.get_position()[0],i.get_position()[1],i.get_position()[2]);
			cells->InsertNextCell(1);
			cells->InsertCellPoint(index);
			cell_types->InsertNextTuple1(1);
			ids->SetValue(index,i.get_id());


			boost::fusion::fold(i.get_data(),0,save_elem(index,datas));
		}
	}


	void  copy_from_vtk_grid(vtkUnstructuredGrid *grid) {
			vtkSmartPointer<vtkPoints> points = grid->GetPoints();
			CHECK(points,"No points in vtkUnstructuredGrid");
			vtkSmartPointer<vtkCellArray> cells = grid->GetCells();
			CHECK(points,"No cells in vtkUnstructuredGrid");
			vtkSmartPointer<vtkUnsignedCharArray> cell_types = grid->GetCellTypesArray();
			vtkSmartPointer<vtkIntArray> ids = vtkIntArray::SafeDownCast(grid->GetPointData()->GetArray("id"));
			CHECK(ids,"No id array in vtkUnstructuredGrid");
			constexpr size_t dn = std::tuple_size<DataType>::value;

			vtkSmartPointer<vtkFloatArray> datas[dn];

			const vtkIdType n = ids->GetSize();

			for (int i = 0; i < dn; ++i) {
				std::string name = DataNames<DataType>::get(i);
				datas[i] = vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name.c_str()));
				CHECK(datas[i],"No data array "<<name<<" in vtkUnstructuredGrid");
				CHECK(datas[i]->GetNumberOfTuples()==n,"Data array "<<name<<" has size != id array. data size = "<<datas[i]->GetNumberOfTuples()<<". id size = "<<n);
			}

			this->clear();

			for (int j = 0; j < n; ++j) {
				value_type particle;
				particle.r  = points->GetPoint(j);
				particle.id = ids->GetValue(j);
				boost::fusion::fold(particle.get_data(),0,read_elem(j,datas));
				this->push_back(particle);
			}
		}
#endif

private:

	void enforce_domain(const Vect3d& low, const Vect3d& high, const Vect3b& periodic, const bool remove_deleted_particles = false) {
		std::for_each(begin(),end(),[low,high,periodic](value_type& i) {
			for (int d = 0; d < 3; ++d) {
				if (periodic[d]) {
					while (i.r[d]<low[d]) {
						i.r[d] += (high[d]-low[d]);
					}
					while (i.r[d]>=high[d]) {
						i.r[d] -= (high[d]-low[d]);
					}
				} else {
					if ((i.r[d]<low[d]) || (i.r[d]>=high[d])) {
						i.mark_for_deletion();
					}
				}
			}
		});
		if (remove_deleted_particles && ((periodic[0]==false)||(periodic[1]==false)||(periodic[2]==false))) {
			const int n = data.size();
			for (int index = 0; index < n; ++index) {
				value_type& i = data[index];
				if (i.alive==false) {
					i.deep_copy(*(data.cend()-1));
					if (track_ids) id_to_index[i.id] = index;
					data.pop_back();
				}
			}
		}
	}


	vector_type data;
	bool searchable,track_ids;
	int next_id;
	const double seed;
	NeighbourSearch_type neighbour_search;
	std::map<size_t,size_t> id_to_index;
#ifndef HAVE_VTK
	vtkSmartPointer<vtkUnstructuredGrid> cache_grid;
#endif
};



}


#endif /* SPECIES_H_ */
