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

#include <boost/mpl/vector.hpp>
#include <boost/mpl/find.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
namespace mpl = boost::mpl;

#include "Zip.h"
#include "Vector.h"
#include "Variable.h"
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





template<typename ... TYPES> 
class Particles {
//	template<typename T>
//	friend class Particles;
public:

    typedef typename std::tuple<position::value_type,id::value_type,alive::value_type,typename TYPES::value_type...> tuple_type;
	typedef typename mpl::vector<position,id,alive,TYPES...> type_vector;

    template<typename T>
    struct elem_by_type {
        typedef T type;
        typedef typename T::value_type value_type;
        typedef typename mpl::find<type_vector,T>::type iter;
        BOOST_MPL_ASSERT_NOT(( boost::is_same< typename mpl::end<type_vector>::type, typename iter::type > ));
        static const size_t index = iter::pos::value;
    };
    template<unsigned int I>
    struct elem_by_index {
        BOOST_MPL_ASSERT_RELATION( (mpl::size<type_vector>::type::value), >, I );
        typedef typename mpl::at<type_vector,mpl::int_<I> > type;
        typedef typename type::value_type value_type;
        static const size_t index = I;
    };

   	class value_type {
		friend class Particles;
	public:
		typedef std::mt19937 generator_type;
		value_type():m_uni(0,1),m_normal(0,1){}
		value_type(const value_type& rhs) {
			m_data = rhs.m_data;
			m_generator = rhs.m_generator;
		}
		value_type(const Vect3d &r):
			m_uni(0,1),m_normal(0,1) {
            this->set<position>(r);
        }
		~value_type() {

		}
		value_type& operator=(const value_type &rhs) {
			if (this != &rhs) {
				m_data = rhs.m_data;
			}
			return *this;
		}

		bool operator==(const value_type &rhs) const {
			return get<id>(*this) == get<id>(rhs);
		}

		void deep_copy(const value_type &rhs) {
			if (this != &rhs) {
				m_data = rhs.m_data;
				m_generator = rhs.m_generator;
			}
		}

		template<typename T>
		const typename elem_by_type<T>::value_type & get() const {
			return std::get<elem_by_type<T>::index>(m_data);
		}

		template<typename T>
		typename elem_by_type<T>::value_type & get() {
			return std::get<elem_by_type<T>::index>(m_data);
		}

		template<typename T>
		void set(const typename elem_by_type<T>::value_type & arg) {
			std::get<elem_by_type<T>::index>(m_data) = arg;
		}

		generator_type get_generator() {
			return m_generator;
		}
		double rand_uniform() {
			return m_uni(m_generator);
		}

		double rand_normal() {
			return m_normal(m_generator);
		}

	private:

        tuple_type m_data;
		generator_type m_generator;
		std::uniform_real_distribution<double> m_uni;
		std::normal_distribution<double> m_normal;

	};

	typedef typename std::vector<value_type> vector_type;
	typedef typename vector_type::size_type size_type;
	typedef typename vector_type::size_type difference_type;
	typedef typename vector_type::iterator iterator;
	typedef typename vector_type::const_iterator const_iterator;
	struct get_pos {
		const Vect3d& operator()(const value_type& i) const {
			return i.template get<position>();
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

	Particles(const Particles<TYPES...> &other):
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

	static ptr<Particles<TYPES...> > New() {
		return ptr<Particles<TYPES...> >(new Particles<TYPES...> ());
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
			while (i.template get<alive>() == false) {
				i.deep_copy(*(data.cend()-1));
				if (track_ids) id_to_index[i.template get<id>()] = index;
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
			if (track_ids) id_to_index[i->template get<id>()] = i-begin();
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
			i.template set<position>(f(i));
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
		i->template set<position>(val.template get<position>());
		i->template set<id>(this->next_id++);
		i->m_generator.seed(i->template get<id>()*seed);
		i->template set<alive>(true);
		if (track_ids) id_to_index[i->template get<id>()] = index;
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
			i->template set<id>(this->next_id++);
			i->m_generator.seed(i->template get<id>()*seed);
			i->template set<alive>(true);
			i->template set<position>(f(*i));
			if (track_ids) id_to_index[i->template get<id>()] = index;
			if (searchable) neighbour_search.add_point(i);
		}
	}

	boost::iterator_range<typename NeighbourSearch_type::const_iterator> get_neighbours(const Vect3d position) {
        ASSERT(searchable == true,"ERROR: using get_neighbours before initialising neighbour search. Please call the init_neighbour_search function before using get_neighbours");
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
			i.template set<position>(f(i));
		});
		if (searchable) {
			enforce_domain(neighbour_search.get_low(),neighbour_search.get_high(),neighbour_search.get_periodic());
			neighbour_search.embed_points(data.cbegin(),data.cend());
		}
	}
	template<typename F>
	void update_positions(F f) {
		std::for_each(begin(),end(),[&f](value_type& i) {
			i.template set<position>(f(i));
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
	struct write_from_tuple {
		write_from_tuple(tuple_type &write_from, int index, vtkSmartPointer<vtkFloatArray>* datas):
			write_from(write_from),index(index),datas(datas){}

        template< typename U > 
		typename boost::enable_if<boost::is_arithmetic<typename std::tuple_element<U::value,tuple_type>::type> >::type
        operator()(U i) {
            typedef typename std::tuple_element<U::value,tuple_type>::type data_type;
            data_type &write_from_elem = std::get<U::value>(write_from);
			datas[i]->SetValue(index,write_from_elem);
        }

        template< typename U >
		typename boost::enable_if<boost::is_same<typename std::tuple_element<U::value,tuple_type>::type,Vect3d> >::type
        operator()(U i) {
            typedef typename std::tuple_element<U::value,tuple_type>::type data_type;
            data_type &write_from_elem = std::get<U::value>(write_from);
            datas[i]->SetTuple3(index,write_from_elem[0],
                                      write_from_elem[1],
                                      write_from_elem[2]);
        }

        tuple_type &write_from;
		int index;
		vtkSmartPointer<vtkFloatArray>* datas;
	};

	struct read_into_tuple {
	    read_into_tuple(tuple_type &read_into, int index, vtkSmartPointer<vtkFloatArray>* datas):
		    read_into(read_into),index(index),datas(datas){}

        template< typename U > 
		typename boost::enable_if<boost::is_arithmetic<typename std::tuple_element<U::value,tuple_type>::type> >::type
        operator()(U i) {
            typedef typename std::tuple_element<U::value,tuple_type>::type data_type;
            data_type &read_into_elem = std::get<U::value>(read_into);
            read_into_elem = datas[i]->GetValue(index);
        }

        template< typename U >
		typename boost::enable_if<boost::is_same<typename std::tuple_element<U::value,tuple_type>::type,Vect3d> >::type
        operator()(U i) {
             typedef typename std::tuple_element<U::value,tuple_type>::type data_type;
            data_type &read_into_elem = std::get<U::value>(read_into);
			double *data = datas[i]->GetTuple3(index);
            read_into_elem[0] = data[0];
            read_into_elem[1] = data[1];
            read_into_elem[2] = data[2];
        }


        tuple_type &read_into;
		int index;
		vtkSmartPointer<vtkFloatArray>* datas;
	};
	
    struct setup_datas_for_writing {
		setup_datas_for_writing(size_t n, vtkSmartPointer<vtkFloatArray>* datas, vtkUnstructuredGrid *grid):
			n(n),datas(datas),grid(grid){}
        template< typename U > void operator()(U i) {
            typedef typename mpl::at<type_vector,U>::type variable_type;
            const char *name = variable_type().name;
            datas[i] = vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
			if (!datas[i]) {
				datas[i] = vtkSmartPointer<vtkFloatArray>::New();
				datas[i]->SetName(name);
				grid->GetPointData()->AddArray(datas[i]);
			}
			typedef typename mpl::at<type_vector,U>::type::value_type data_type;
            if (boost::is_same<data_type, Vect3d>::value) {
                datas[i]->SetNumberOfComponents(3);
            } else {
                datas[i]->SetNumberOfComponents(1);
            }
			datas[i]->SetNumberOfTuples(n);
        }

        size_t n;
		vtkSmartPointer<vtkFloatArray>* datas;
		vtkUnstructuredGrid *grid;
    };

    struct setup_datas_for_reading {
		setup_datas_for_reading(size_t n,vtkSmartPointer<vtkFloatArray>* datas,vtkUnstructuredGrid *grid):
			n(n),datas(datas),grid(grid){}
        template< typename U > void operator()(U i) {
            std::string name = mpl::at<type_vector,typename U::value>::type::name;
			datas[i] = vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name.c_str()));
			CHECK(datas[i],"No data array "<<name<<" in vtkUnstructuredGrid");
			CHECK(datas[i]->GetNumberOfTuples()==n,"Data array "<<name<<" has size != id array. data size = "<<datas[i]->GetNumberOfTuples()<<". id size = "<<n);
        }

        size_t n;
		vtkSmartPointer<vtkFloatArray>* datas;
		vtkUnstructuredGrid *grid;
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
		

		const vtkIdType n = size();

		constexpr size_t dn = std::tuple_size<tuple_type>::value;
		vtkSmartPointer<vtkFloatArray> datas[dn];
        mpl::for_each<mpl::range_c<int,1,mpl::size<type_vector>::type::value> > (setup_datas_for_writing(n,datas,grid));
		points->SetNumberOfPoints(n);
		cells->Reset();
		cell_types->Reset();
		int j = 0;

		for(auto& i: *this) {
			const int index = j++;
			//std::cout <<"copying point at "<<i.get_position()<<" with id = "<<i.get_id()<<std::endl;
            const Vect3d &r = i.template get<position>();
			points->SetPoint(index,r[0],r[1],r[2]);
			cells->InsertNextCell(1);
			cells->InsertCellPoint(index);
			cell_types->InsertNextTuple1(1);

            mpl::for_each<mpl::range_c<int,1,mpl::size<type_vector>::type::value> > (write_from_tuple(i.m_data,index,datas));
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
			constexpr size_t dn = std::tuple_size<tuple_type>::value;

			vtkSmartPointer<vtkFloatArray> datas[dn];

			const vtkIdType n = ids->GetSize();

            mpl::for_each<mpl::range_c<int,1,mpl::size<type_vector>::type::value> > (setup_datas_for_reading(n,datas,grid));

			this->clear();

			for (int j = 0; j < n; ++j) {
				value_type particle;
				particle.r  = points->GetPoint(j);
                mpl::for_each<mpl::range_c<int,1,mpl::size<type_vector>::type::value> > (read_into_tuple(particle.m_data,index,datas));
				this->push_back(particle);
			}
		}
#endif

private:


	void enforce_domain(const Vect3d& low, const Vect3d& high, const Vect3b& periodic, const bool remove_deleted_particles = false) {
		std::for_each(begin(),end(),[low,high,periodic](value_type& i) {
            Vect3d &r = i.template get<position>();
			for (int d = 0; d < 3; ++d) {
				if (periodic[d]) {
					while (r[d]<low[d]) {
						r[d] += (high[d]-low[d]);
					}
					while (r[d]>=high[d]) {
						r[d] -= (high[d]-low[d]);
					}
				} else {
					if ((r[d]<low[d]) || (r[d]>=high[d])) {
						i.template set<alive>(false);
                    }
				}
			}
		});
		if (remove_deleted_particles && ((periodic[0]==false)||(periodic[1]==false)||(periodic[2]==false))) {
			const int n = data.size();
			for (int index = 0; index < n; ++index) {
				value_type& i = data[index];
				if (i.template get<alive>()==false) {
					i.deep_copy(*(data.cend()-1));
					if (track_ids) id_to_index[i.template get<id>()] = index;
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


template<typename T, typename VT>
const typename T::value_type & get(const VT & arg) {
    return arg.template get<T>();
}
template<typename T, typename VT>
typename T::value_type & get(VT & arg) {
    return arg.template get<T>();
}
template<typename T, typename VT>
void set(VT & arg, const typename T::value_type & data) {
    arg.template set<T>(data);
}


}


#endif /* SPECIES_H_ */
