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

#ifdef HAVE_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#endif


namespace Aboria {

/// 
/// helper class to find an element of type_vector from a Variable type T
template<typename T, typename type_vector>
struct elem_by_type {
    typedef T type;
    typedef typename T::value_type value_type;

    /// 
    /// iter is a boost mpl iterator to the found element in Variable T
    typedef typename mpl::find<type_vector,T>::type iter;
    BOOST_MPL_ASSERT_NOT(( boost::is_same< typename mpl::end<type_vector>::type, typename iter::type > ));

    /// 
    /// index contains the index of the found element
    static const size_t index = iter::pos::value;
};


/// helper class to find an element of type_vector from an
/// unsigned int index I
template<unsigned int I, typename type_vector>
struct elem_by_index {
    BOOST_MPL_ASSERT_RELATION( (mpl::size<type_vector>::type::value), >, I );
    typedef typename mpl::at<type_vector,mpl::int_<I> > type;

    /// 
    /// value_type is the variable's value_type at index I 
    typedef typename type::value_type value_type;
    static const size_t index = I;
};


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
template<typename ... TYPES> 
class Particles {
public:
    /// a boost mpl vector type containing a vector of Variable 
    /// attached to the particles (includes position, id and 
    /// alive flag as well as all user-supplied variables)
    typedef typename mpl::vector<position,id,alive,TYPES...> type_vector;

    /// 
    /// a tuple type containing a list of value_types for each Variable
    typedef typename std::tuple<position::value_type,id::value_type,alive::value_type,typename TYPES::value_type...> tuple_type;

    
    /// \brief class used to store information on each particle. 
    ///
    /// Each particle has a number of attached variables, 
    /// including by default a position, id and alive flag. 
    /// The data for each variable is contained in value_type::m_data with type
    /// Particles::tuple_type
    ///
    /// Each particle also has a random number generator with a seed based on
    /// its id. setting this seed is the responsibility of the container class
    class Particle { 
        friend class Particles; 
    public: 
        /// default random number generator
        typedef std::mt19937 generator_type; 

        /// default constructor
        Particle(): m_uni(0,1), m_normal(0,1) {} 

        /// copy-contructor for particles
        Particle(const Particle& rhs): 
            m_uni(rhs.m_uni), 
            m_normal(rhs.m_normal), 
            m_data(rhs.m_data), 
            m_generator(rhs.m_generator) {} 

        /// create particle with a given position
        explicit Particle(const Vect3d &r): 
            m_uni(0,1), 
            m_normal(0,1) { 
                this->set<position>(r); 
            } 

        ~Particle() {}

        /// assignment operator only copies all Variable data.
        /// random number generator is not copied
        Particle& operator=(const Particle &rhs) {
            if (this != &rhs) {
                m_data = rhs.m_data;
            }
            return *this;
        }

        /// particles are equal if they have the same id Variable
        bool operator==(const Particle &rhs) const {
            return get<id>(*this) == get<id>(rhs);
        }

        /// perform a deep particle copy, random number generator is 
        /// also copied as well as Variable data
        void deep_copy(const Particle &rhs) {
            if (this != &rhs) {
                m_data = rhs.m_data;
                m_generator = rhs.m_generator;
            }
        }

        /// get the value_type of a stored Variable. get() is templated using 
        /// the Variable type itself.
        /// \return const reference to a T::value_type
        template<typename T>
        const typename elem_by_type<T,type_vector>::value_type& get() const {
            return std::get<elem_by_type<T,type_vector>::index>(m_data);
        }

        /// non-const version of get() const
        /// \return reference to a T::value_type
        template<typename T>
        typename elem_by_type<T,type_vector>::value_type& get() {
            return std::get<elem_by_type<T,type_vector>::index>(m_data);
        }

        /// set the value_type of a stored Variable
        /// \param T a Variable type
        /// \param arg a type compatible with T::value_type
        template<typename T>
        void set(const typename elem_by_type<T,type_vector>::value_type& arg) {
            std::get<elem_by_type<T,type_vector>::index>(m_data) = arg;
        }

        /// getter for internal random number generator
        generator_type get_generator() {
            return m_generator;
        }

        /// generate a uniform random number between 0 and 1
        double rand_uniform() {
            std::uniform_real_distribution<double> uni(0,1);
            return uni(m_generator);
        }

        /// generate a normally distributed random number
        double rand_normal() {
            std::normal_distribution<double> normal(0,1);
            return normal(m_generator);
        }

    private:

        tuple_type m_data;
        generator_type m_generator;
        std::uniform_real_distribution<double> m_uni;
        std::normal_distribution<double> m_normal;

    };

    /// typedef Particle to value_type
    typedef Particle value_type;
    /// vector type used to hold a collection of value_type
    typedef typename std::vector<value_type> vector_type;
    /// type used to hold and return size information
    typedef typename vector_type::size_type size_type;
    typedef typename vector_type::size_type difference_type;
    /// non-const iterator type
    typedef typename vector_type::iterator iterator;
    /// const iterator type
    typedef typename vector_type::const_iterator const_iterator;
    struct get_pos {
        const Vect3d& operator()(const value_type& i) const {
            return i.template get<position>();
        }
    };
    /// external type used to implement neighbourhood searching
    typedef BucketSort<const_iterator,get_pos> NeighbourSearch_type;


    /// Contructs an empty container with no searching or id tracking enabled
    Particles():
        next_id(0),
        neighbour_search(Vect3d(0,0,0),Vect3d(1,1,1),Vect3b(false,false,false),get_pos()),
        searchable(false),track_ids(false),
        seed(time(NULL))
    {}

    /// Contructs an empty container with no searching or id tracking enabled
    /// \param seed initialises the base random seed for the container
    Particles(const double seed):
        next_id(0),
        neighbour_search(Vect3d(0,0,0),Vect3d(1,1,1),Vect3b(false,false,false),get_pos()),
        searchable(false),track_ids(false),
        seed(seed)
    {}

    /// copy-constructor. performs deep copying of all particles
    Particles(const Particles<TYPES...> &other):
            data(other.data),
            neighbour_search(other.neighbour_search),
            next_id(other.next_id),
            searchable(other.searchable),
            track_ids(other.track_ids),
            seed(other.seed),
            id_to_index(id_to_index)
    {}

    /// range-based copy-constructor. performs deep copying of all 
    /// particles from \p first to \p last
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

    
    //
    // STL Container
    //
    
    /// push the particle \p val to the back of the container
    void push_back (const value_type& val) {
        data.push_back(val);
        if (searchable) neighbour_search.update_begin_and_end(data.cbegin(),data.cend());
        const int index = data.size();
        iterator i = end()-1;
        i->template set<position>(val.template get<position>());
        i->template set<id>(this->next_id++);
        i->m_generator.seed(i->template get<id>()+seed);
        i->template set<alive>(true);
        if (track_ids) id_to_index[i->template get<id>()] = index;
        if (searchable) neighbour_search.add_point(i);
    }

    /// push a new particle with position \p position
    /// to the back of the container
    void push_back(const Vect3d& position) {
        this->push_back(value_type(position));
    }

    /// pop (delete) the particle at the end of the container 
    void pop_back() {
        erase(end()-1);
    }

    /// returns a reference to the particle at position \p idx
    value_type& operator[](std::size_t idx) {
        return data[idx];
    }

    /// returns a const reference to the particle at position \p idx
    const value_type& operator[](std::size_t idx) const {
        return const_cast<value_type>(*this)[idx];
    }

    /// returns an iterator to the beginning of the container
    iterator begin() {
        return data.begin();
    }

    /// returns an iterator to the end of the container
    iterator end() {
        return data.end();
    }

    

    /// sets container to empty and deletes all particles
    void clear() {
        data.clear();
    }

    /// erase the particle pointed to by the iterator \p i.
    /// NOTE: This will potentially reorder the particles
    /// if neighbourhood searching is on, then this is updated
    iterator erase (iterator i) {
        if (i != end()-1) {
            i->deep_copy(*(data.cend()-1));
            if (track_ids) id_to_index[i->template get<id>()] = i-begin();
            if (searchable) neighbour_search.copy_points(i,end()-1);
            if (searchable) neighbour_search.delete_point(end()-1);
            data.pop_back();
            return i;
        } else {
            if (searchable) neighbour_search.delete_point(end()-1);
            data.pop_back();
            return end();
        }
    }


    /// erase the particles between the iterators \p first
    /// and \p last
    /// \see erase(iterator)
    iterator erase (iterator first, iterator last) {
        for(iterator i=first;i!=last-1;i++) {
            erase(i);
        }
        return erase(last-1);
    }

    /// insert a particle \p val into the container at \p position
    iterator insert (iterator position, const value_type& val) {
        data.insert(position,val);
    }

    /// insert a \p n copies of the particle \p val into the container at \p position
    void insert (iterator position, size_type n, const value_type& val) {
        data.insert(position,n,val);
    }

    /// insert a range of particles pointed to by \p first and \p last at \p position 
    template <class InputIterator>
    void insert (iterator position, InputIterator first, InputIterator last) {
        data.insert(position,first,last);
    }

    /// return the total number of particles in the container
    size_t size() const {
        return data.size();
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
    void init_neighbour_search(const Vect3d& low, const Vect3d& high, const double length_scale, const Vect3b& periodic) {
        neighbour_search.reset(low,high,length_scale,periodic);
        neighbour_search.embed_points(data.cbegin(),data.cend());
        searchable = true;
    }

    /// get all particles within the neighbourhood length scale of a point. 
    /// NOTE: you must call init_neighbour_search() before using this function
    /// \param position the centre of the search region
    /// \see init_neighbour_search
    boost::iterator_range<typename NeighbourSearch_type::const_iterator> get_neighbours(const Vect3d position) {
        ASSERT(searchable == true,"ERROR: using get_neighbours before initialising neighbour search. Please call the init_neighbour_search function before using get_neighbours");
        return boost::make_iterator_range(
                neighbour_search.find_broadphase_neighbours(position, -1,false),
                neighbour_search.end());
    }

    /// set the length scale of the neighbourhood search to be equal to \p length_scale
    /// \see init_neighbour_search()
    void reset_neighbour_search(const double length_scale) {
        neighbour_search.reset(neighbour_search.get_low(),neighbour_search.get_high(),length_scale,neighbour_search.get_periodic());
        neighbour_search.embed_points(data.cbegin(),data.cend());
        searchable = true;
    }

    /// return the length scale of the neighbourhood search
    /// \see init_neighbour_search()
    const double get_lengthscale() const {
        return neighbour_search.get_lengthscale();
    }

    /// return the lower extent of the neighbourhood search
    /// \see init_neighbour_search()
    const Vect3d& get_low() const {
        return neighbour_search.get_low();
    }

    /// return the upper extent of the neighbourhood search
    /// \see init_neighbour_search()
    const Vect3d& get_high() const {
        return neighbour_search.get_high();
    }

    /// Update the neighbourhood search data. This function must be
    /// called after altering the particle positions (e.g. with 
    /// `set<position>(particle,new_position)`) in order for accurate
    /// neighbourhood searching
    /// \see get_neighbours()
    void update_positions() {
        if (searchable) {
            enforce_domain(neighbour_search.get_low(),neighbour_search.get_high(),neighbour_search.get_periodic());
        }
    }

    /// Update the neighbourhood search data continuously while altering the particle
    /// positions from \p b to \p e. This function iterates through the container from 
    /// \p b to \p e and updates the position of each particle using \p f(). It then
    /// refeshes the neighbourhood searching before moving onto the next particle.
    /// \param f a functor that takes a reference to a particle value_type and returns
    /// a new position for that particle
    /// \see get_neighbours()
    template<typename F>
    void update_positions_sequential(iterator b, iterator e, F f) {
        const Vect3d low = neighbour_search.get_low();
        const Vect3d high = neighbour_search.get_high();
        const Vect3b periodic = neighbour_search.get_periodic();
        for(iterator i = b; i != e; i++) {
            Vect3d& r = i->template get<position>();
            r = f(*i);
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
                        i->template set<alive>(false);
                    }
                }
            }
            neighbour_search.update_point(i);
        }
    }

    //
    // ID Tracking
    //

    /// turn on or off id tracking
    /// \see find_id(const int)
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

    /// return a pointer to the particle with id==\p id
    /// \return if the particle with id==\p id is within the container 
    /// it returns a pointer to it, otherwise NULL
    value_type* find_id(const int id) {
        std::map<size_t,size_t>::iterator it = id_to_index.find(id);
        if (it==id_to_index.end()) {
            return NULL;
        } else {
            return &(data[*it]);
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
        for (int index = 0; index < data.size(); ++index) {
            typename vector_type::iterator i = data.begin() + index;
            while (i->template get<alive>() == false) {
                if ((index < data.size()-1) && (data.size() > 1)) {
                    i->deep_copy(*(data.cend()-1));
                    if (track_ids) id_to_index[i->template get<id>()] = index;
                    data.pop_back();
                    i = data.begin() + index;
                } else {
                    data.pop_back();
                    break;
                }
            }
        }
        if (searchable && update_neighbour_search) neighbour_search.embed_points(data.cbegin(),data.cend());
    }

    /// Create \p n particles with a user-defined function that sets each 
    /// particles position
    /// \param n number of particles to create
    /// \param f a functor that returns a Vect3d position
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

            mpl::for_each<mpl::range_c<int,1,mpl::size<type_vector>::type::value> > (setup_datas_for_reading(n,datas,grid));

            this->clear();

            for (int j = 0; j < n; ++j) {
                value_type particle;
                const double *r = points->GetPoint(j);
                particle.template set<position>(Vect3d(r[0],r[1],r[2]));
                mpl::for_each<mpl::range_c<int,1,mpl::size<type_vector>::type::value> > (read_into_tuple(particle.m_data,j,datas));
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
    void enforce_domain(const Vect3d& low, const Vect3d& high, const Vect3b& periodic, const bool remove_deleted_particles = true) {
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
            delete_particles();
        } else {
            neighbour_search.embed_points(data.cbegin(),data.cend());
        }
    }


    vector_type data;
    bool searchable,track_ids;
    int next_id;
    const double seed;
    NeighbourSearch_type neighbour_search;
    std::map<size_t,size_t> id_to_index;

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
            typedef typename mpl::at<type_vector,U>::type variable_type;
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


/// get a variable from a particle \p arg
/// \param T Variable type
/// \param VT Particle type (optional)
/// \param arg the particle
/// \return a const reference to a T::value_type holding the variable data
template<typename T, typename VT>
const typename T::value_type & get(const VT & arg) {
    return arg.template get<T>();
}

/// get a variable from a particle \p arg
/// \param T Variable type
/// \param VT Particle type (optional)
/// \param arg the particle
/// \return a reference to a T::value_type holding the variable data
template<typename T, typename VT>
typename T::value_type & get(VT & arg) {
    return arg.template get<T>();
}

/// set a variable attached to a particle \p arg
/// \param T Variable type
/// \param VT Particle type (optional)
/// \param arg the particle
/// \param data a reference to a T::value_type. This will be copied to
/// the variable attached to particle \p arg
template<typename T, typename VT>
void set(VT & arg, const typename T::value_type & data) {
    arg.template set<T>(data);
}


}


#endif /* SPECIES_H_ */
