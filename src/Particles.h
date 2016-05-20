/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Aboria.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <vector>
#include <random>
#include <string>
//#include <boost/array.hpp>
//#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/adaptors.hpp>

#include "Zip.h"
#include "Get.h"
#include "Vector.h"
#include "Variable.h"
#include "Traits.h"
#include "BucketSearch2.h"
#include "OctTree.h"

#ifdef HAVE_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#endif

#ifdef HAVE_THRUST
#include <thrust/device_vector.h>
#endif


namespace Aboria {

/// \brief A STL-compatable container of particles in 3D space
///
///  Each particle has a 3D position and user-defined data-package 
///  (for other variables such as velocity, density etc) and is 
///  optionally embedded within a cuboidal spatial domain (for neighbourhood searches) 
///  that can be periodic or not. Each particle also has its own random number 
///  generator that is seeded via its own unique id.
///
///  For example, the following creates a set of particles which each have 
///  (along with the standard variables such as position, id etc) a 
///  data package consisting of one double variable type named scalar.
///
///  \code 
///     using namespace Aboria;
///     ABORIA_VARIABLE(scalar,double,"my scalar")
///     typedef Particles<scalar> MyParticles;
///     MyParticles particles();
///  \endcode
///
///  \param TYPES a list of one or more variable types
///
///  \see #ABORIA_VARIABLE
template<typename VAR=std::tuple<>, unsigned int D=3, template <typename,typename> class VECTOR=std::vector, typename TRAITS_USER=Traits<VECTOR> > 
class Particles {
public:

    typedef Particles<VAR,D,VECTOR,TRAITS_USER> particles_type;

    typedef TraitsCommon<VAR,D,TRAITS_USER> traits_type;

    /// a boost mpl vector type containing a vector of Variable 
    /// attached to the particles (includes position, id and 
    /// alive flag as well as all user-supplied variables)
    typedef typename traits_type::mpl_type_vector mpl_type_vector;

    /// 
    /// a tuple type containing a list of references to value_types for each Variable
    typedef typename traits_type::value_type value_type;

    typedef typename traits_type::reference reference;

   
    /// type used to hold data (tuple of vectors or similar)
    typedef typename traits_type::data_type data_type;
    /// type used to hold and return size information
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::size_type difference_type;
    /// non-const iterator type
    typedef typename traits_type::iterator iterator;
    /// const iterator type
    typedef typename traits_type::const_iterator const_iterator;

    UNPACK_TRAITS(traits_type)

    /// Contructs an empty container with no searching or id tracking enabled
    Particles():
        next_id(0),
        searchable(false),
        seed(time(NULL))
    {}

    /// Contructs an empty container with no searching or id tracking enabled
    /// \param seed initialises the base random seed for the container
    Particles(const uint32_t seed):
        next_id(0),
        searchable(false),
        seed(seed)
    {}

    /// copy-constructor. performs deep copying of all particles
    Particles(const particles_type &other):
            data(other.data),
            bucket_search(other.bucket_search),
            next_id(other.next_id),
            searchable(other.searchable),
            seed(other.seed),
            id_to_index(id_to_index)
    {}

    /// range-based copy-constructor. performs deep copying of all 
    /// particles from \p first to \p last
    Particles(iterator first, iterator last):
        data(first,last),
        searchable(false),
        seed(0)
    {
        if (searchable) embed_points(begin(),end());
    }


    
    //
    // STL Container
    //
    
    /// push the particle \p val to the back of the container
    void push_back (const value_type& val, bool update_neighbour_search=true) {
        traits_type::push_back(data,val);
        const int index = size();
        reference i = *(end()-1);
        Aboria::set<position>(i,Aboria::get<position>(val));
        Aboria::set<id>(i,this->next_id++);
        Aboria::get<random>(i).seed(seed + uint32_t(Aboria::get<id>(i)));
        Aboria::set<alive>(i,true);
        if (searchable && update_neighbour_search) bucket_search.add_points_at_end(begin(),end()-1,end());
    }

    /// push a new particle with position \p position
    /// to the back of the container
    void push_back(const double_d& pos) {
        value_type i;
        Aboria::set<position>(i,pos);
        this->push_back(i);
    }

    /// push the particles in \p particles to the back of the container
    void push_back (const particles_type& particles) {
        for (const value_type& i: particles) {
            this->push_back(i,false);
        }
        if (searchable) bucket_search.add_points_at_end(data.begin(),data.end()-particles.size(),data.end());
    }

    /// pop (delete) the particle at the end of the container 
    void pop_back() {
        erase(end()-1);
    }

    /// returns a reference to the particle at position \p idx
    reference operator[](std::size_t idx) {
        return traits_type::index(data, idx);
    }

    /*
    /// returns a const reference to the particle at position \p idx
    const_value_type operator[](std::size_t idx) const {
        return traits_type::index_const(data, idx);
    }
    */

    /// returns an iterator to the beginning of the container
    iterator begin() {
        return traits_type::begin(data);
    }

    /// returns an iterator to the end of the container
    iterator end() {
        return traits_type::end(data);
    }


    /// sets container to empty and deletes all particles
    void clear() {
        return traits_type::clear(data);
    }

    /// erase the particle pointed to by the iterator \p i.
    /// NOTE: This will potentially reorder the particles
    /// if neighbourhood searching is on, then this is updated
    iterator erase (iterator i, bool update_neighbour_search = true) {
        if (i != end()-1) {
            *i = *(end()-1);
            traits_type::pop_back(data);
            return i;
        } else {
            traits_type::pop_back(data);
            return end();
        }

        if (searchable && update_neighbour_search) bucket_search.embed_points(begin(),end());
    }


    /// erase the particles between the iterators \p first
    /// and \p last
    /// \see erase(iterator)
    iterator erase (iterator first, iterator last) {
        for(iterator i=first;i!=last-1;i++) {
            erase(i,false);
        }
        iterator return_iterator = erase(last-1,false);
        if (searchable) bucket_search.embed_points(begin(),end());
        return erase(last-1);
    }

    /// insert a particle \p val into the container at \p position
    iterator insert (iterator position, const value_type& val) {
        traits_type::insert(data,position,val);
    }

    /// insert a \p n copies of the particle \p val into the container at \p position
    void insert (iterator position, size_type n, const value_type& val) {
        traits_type::insert(data,position,n,val);
    }

    /// insert a range of particles pointed to by \p first and \p last at \p position 
    template <class InputIterator>
    void insert (iterator position, InputIterator first, InputIterator last) {
        traits_type::insert(data,position,first,last);
        data.insert(position,first,last);
    }

    /// return the total number of particles in the container
    size_t size() const {
        return this->get<position>().size();
    }


    //
    // Neighbourhood Searching
    //

    /// initialise the neighbourhood searching for the particle container.
    /// The neighbourhood searching is performed over a cuboidal domain from 
    /// \p low to \p high.
    ///
    /// \param low the lowest point in the search domain
    /// \param high the highest point in the search domain
    /// \param length_scale the maximum search radius required. particles closer 
    /// to each other than this length are considered neighbours
    /// \param periodic a boolean 3d vector indicating whether each dimension 
    /// is periodic (true) or not (false)
    void init_neighbour_search(const double_d& low, const double_d& high, const double length_scale, const bool_d& periodic) {
        bucket_search.set_domain(low,high,periodic,double_d(length_scale));
        enforce_domain(bucket_search.get_min(),bucket_search.get_max(),bucket_search.get_periodic());
        searchable = true;
    }

    /// get all particles within the neighbourhood length scale of a point. 
    /// NOTE: you must call init_neighbour_search() before using this function
    /// \param position the centre of the search region
    /// \see init_neighbour_search
    boost::iterator_range<typename BucketSearch<traits_type>::const_iterator> get_neighbours(const double_d position) {
        ASSERT(searchable == true,"ERROR: using get_neighbours before initialising neighbour search. Please call the init_neighbour_search function before using get_neighbours");
        return boost::make_iterator_range(
                bucket_search.find_broadphase_neighbours(position, -1,false),
                bucket_search.end());
    }

    /// set the length scale of the neighbourhood search to be equal to \p length_scale
    /// \see init_neighbour_search()
    void reset_neighbour_search(const double length_scale) {
        bucket_search.set_domain(bucket_search.get_min(),
                                    bucket_search.get_max(),
                                    bucket_search.get_periodic(),
                                    double_d(length_scale));
        bucket_search.embed_points(begin(),end());
        searchable = true;
    }

    /// return the length scale of the neighbourhood search
    /// \see init_neighbour_search()
    const double get_lengthscale() const {
        return bucket_search.get_bucket_side_length();
    }

    /// return the lower extent of the neighbourhood search
    /// \see init_neighbour_search()
    const double_d& get_min() const {
        return bucket_search.get_min();
    }

    /// return the upper extent of the neighbourhood search
    /// \see init_neighbour_search()
    const double_d& get_max() const {
        return bucket_search.get_max();
    }

    /// Update the neighbourhood search data. This function must be
    /// called after altering the particle positions (e.g. with 
    /// `set<position>(particle,new_position)`) in order for accurate
    /// neighbourhood searching
    /// \see get_neighbours()
    void update_positions() {
        if (searchable) {
            enforce_domain(bucket_search.get_min(),bucket_search.get_max(),bucket_search.get_periodic());
        }
    }


    //
    // Particle Creation/Deletion
    //
    
    /// deletes all particles with alive==false from the container
    /// NOTE: this will reorder the particles in the container, invalidating
    /// any iterators
    /// \param update_neighbour_search updates neighbourhood search
    /// information if true (default=true)
    void delete_particles(const bool update_neighbour_search = true) {
        for (int index = 0; index < size(); ++index) {
            iterator i = begin() + index;
            while (Aboria::get<alive>(*i) == false) {
                if ((index < size()-1) && (size() > 1)) {
                    //std::swap(*i,*(end()-1));
                    *i = *(end()-1);
                    pop_back();
                    i = begin() + index;
                } else {
                    pop_back();
                    break;
                }
            }
        }
        if (searchable && update_neighbour_search) bucket_search.embed_points(begin(),end());
    }


    //
    // Vector getters/setters
    //
    template<typename T>
    const typename TRAITS_USER::template vector_type<typename T::value_type>::type & get() const {
        return data.get<T>();
    }

private:
    struct set_variable {
        set_variable(data_type &data_to_set, data_type &internal_data, int n):
            data_to_set(data_to_set),internal_data(internal_data),n(n) {}
        template< typename U > void operator()(U i) {
            typedef typename mpl::at<mpl_type_vector,U>::type variable_type;
            const char *name = variable_type().name;
            CHECK(n == data_to_set.get<U::value>().size(), "error, data vectors not the same size")
            internal_data.set<U::value>(data_to_set.get<U::value>());
        }

        data_type &data_to_set;
        data_type &internal_data;
        int n;
    };
public:

    void set(data_type &data_to_set) {
        mpl::for_each<mpl::range_c<int,1,mpl::size<mpl_type_vector>::type::value> > (set_variable(data_to_set,data,data_to_set.get<0>().size()));
    }


#ifdef HAVE_VTK
    
    /// get a vtk unstructured grid version of the particle container
    /// This grid is cached internally. The first time this function is 
    /// called the grid is created and updated with the particle data.
    /// If the particle data is updated subsequently you need to set \p refresh=true
    /// for the grid to be updated
    vtkSmartPointer<vtkUnstructuredGrid> get_grid(bool refresh = false) {
        if (!cache_grid) {
            cache_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
            copy_to_vtk_grid(cache_grid);
        } else if (refresh == true) {
            copy_to_vtk_grid(cache_grid);
        }
        return cache_grid;
    }

    ///  copy the particle data to a VTK unstructured grid
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
        mpl::for_each<mpl::range_c<int,1,mpl::size<mpl_type_vector>::type::value> > (setup_datas_for_writing(n,datas,grid));
        points->SetNumberOfPoints(n);
        cells->Reset();
        cell_types->Reset();
        int j = 0;

        for(auto& i: *this) {
            const int index = j++;
            //std::cout <<"copying point at "<<i.get_position()<<" with id = "<<i.get_id()<<std::endl;
            const double_d &r = get<position>(i);
            points->SetPoint(index,r[0],r[1],r[2]);
            cells->InsertNextCell(1);
            cells->InsertCellPoint(index);
            cell_types->InsertNextTuple1(1);

            mpl::for_each<mpl::range_c<int,1,mpl::size<mpl_type_vector>::type::value> > (write_from_tuple(i.m_data,index,datas));
        }
    }


    ///  update the particle data according to a supplied VTK unstructured grid
    ///  it is assumed that \p grid has been generated from copy_to_vtk_grid() or
    ///  get_grid()
    ///  \see get_grid()
    ///  \see copy_to_vtk_grid()
    void  copy_from_vtk_grid(vtkUnstructuredGrid *grid) {
            vtkSmartPointer<vtkPoints> points = grid->GetPoints();
            CHECK(points,"No points in vtkUnstructuredGrid");
            vtkSmartPointer<vtkCellArray> cells = grid->GetCells();
            CHECK(points,"No cells in vtkUnstructuredGrid");
            vtkSmartPointer<vtkUnsignedCharArray> cell_types = grid->GetCellTypesArray();
            constexpr size_t dn = std::tuple_size<tuple_type>::value;

            vtkSmartPointer<vtkFloatArray> datas[dn];

            const vtkIdType n = points->GetNumberOfPoints();

            mpl::for_each<mpl::range_c<int,1,mpl::size<mpl_type_vector>::type::value> > (setup_datas_for_reading(n,datas,grid));

            this->clear();

            for (int j = 0; j < n; ++j) {
                value_type particle;
                const double *r = points->GetPoint(j);
                set<position>(particle,double_d(r[0],r[1],r[2]));
                mpl::for_each<mpl::range_c<int,1,mpl::size<mpl_type_vector>::type::value> > (read_into_tuple(particle.m_data,j,datas));
                this->push_back(particle);
            }
        }
#endif

private:

    /// enforce a cuboidal domain. Any particles outside this domain for 
    /// non-periodic dimensions will have alive set to false. For periodic dimensions
    /// the particle will be placed accordingly back within the domain
    /// \param low lower extent of the domain
    /// \param high upper extent of the domain
    /// \param periodic boolean 3d vector setting each dimension to be periodic (true)
    /// or non-periodic (false)
    /// \param remove_deleted_particles if true, removes particles with alive==false from 
    /// the container (default = true)
    void enforce_domain(const double_d& low, const double_d& high, const bool_d& periodic, const bool remove_deleted_particles = true) {
        std::for_each(begin(),end(),[low,high,periodic](reference i) {
            double_d &r = Aboria::get<position>(i);
            for (unsigned int d = 0; d < dimension; ++d) {
                if (periodic[d]) {
                    while (r[d]<low[d]) {
                        r[d] += (high[d]-low[d]);
                    }
                    while (r[d]>=high[d]) {
                        r[d] -= (high[d]-low[d]);
                    }
                } else {
                    if ((r[d]<low[d]) || (r[d]>=high[d])) {
                    Aboria::set<alive>(i,false);
                    }
                }
            }
        });
        if (remove_deleted_particles && ((periodic[0]==false)||(periodic[1]==false)||(periodic[2]==false))) {
            delete_particles();
        } else {
            bucket_search.embed_points(begin(),end());
        }
    }


    data_type data;
    bool searchable;
    int next_id;
    const uint32_t seed;
    std::map<size_t,size_t> id_to_index;
    BucketSearch<traits_type> bucket_search;


#ifdef HAVE_VTK
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
        typename boost::enable_if<boost::is_same<typename std::tuple_element<U::value,tuple_type>::type,double_d> >::type
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
        typename boost::enable_if<boost::is_same<typename std::tuple_element<U::value,tuple_type>::type,double_d> >::type
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
            typedef typename mpl::at<mpl_type_vector,U>::type variable_type;
            const char *name = variable_type().name;
            datas[i] = vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
            if (!datas[i]) {
                datas[i] = vtkSmartPointer<vtkFloatArray>::New();
                datas[i]->SetName(name);
                grid->GetPointData()->AddArray(datas[i]);
            }
            typedef typename mpl::at<mpl_type_vector,U>::type::value_type data_type;
            if (boost::is_same<data_type, double_d>::value) {
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
            typedef typename mpl::at<mpl_type_vector,U>::type variable_type;
            const char *name = variable_type().name;
            datas[i] = vtkFloatArray::SafeDownCast(grid->GetPointData()->GetArray(name));
            CHECK(datas[i],"No data array "<<name<<" in vtkUnstructuredGrid");
            CHECK(datas[i]->GetNumberOfTuples()==n,"Data array "<<name<<" has size != id array. data size = "<<datas[i]->GetNumberOfTuples()<<". id size = "<<n);
        }

        size_t n;
        vtkSmartPointer<vtkFloatArray>* datas;
        vtkUnstructuredGrid *grid;
    };

    vtkSmartPointer<vtkUnstructuredGrid> cache_grid;
#endif
};



}


#endif /* SPECIES_H_ */
