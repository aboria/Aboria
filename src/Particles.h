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
#include "BucketSearchSerial.h"
//#include "OctTree.h"
#include "CudaInclude.h"

#ifdef HAVE_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#endif

#include "detail/Particles.h"


namespace Aboria {
namespace detail {

template <typename Reference>
struct resize_lambda {
    uint32_t seed;
    int next_id;
    const size_t *start_id_pointer;

    resize_lambda(const uint32_t &seed, const int &next_id, const size_t *start_id_pointer):
            seed(seed),next_id(next_id),start_id_pointer(start_id_pointer) {}

    CUDA_HOST_DEVICE
    void operator()(Reference i) const {
        Aboria::get<alive>(i) = true;

        const size_t index = &Aboria::get<id>(i)-start_id_pointer;
        Aboria::get<id>(i) = index + next_id;

        generator_type& gen = Aboria::get<generator>(i);
        gen.seed(seed + uint32_t(Aboria::get<id>(i)));

    }
};

template <typename Reference>
struct set_seed_lambda {
    uint32_t seed;

    set_seed_lambda(const uint32_t &seed):
            seed(seed) {}

    CUDA_HOST_DEVICE
    void operator()(Reference i) const {
        generator_type& gen = Aboria::get<generator>(i);
        gen.seed(seed + uint32_t(Aboria::get<id>(i)));
    }
};

template <typename ConstReference>
struct is_alive {
    CUDA_HOST_DEVICE
    bool operator()(ConstReference i) const {
        return Aboria::get<alive>(i);
    }
};



}
    /*
#ifdef HAVE_THUST
template <typename T, typename Alloc> 
using default_vector = thrust::host_vector<T,Alloc>;
#else
template <typename T, typename Alloc> 
using default_vector = std::vector<T,Alloc>;
#endif
*/

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
template<typename VAR=std::tuple<>, unsigned int D=3, template <typename,typename> class VECTOR=std::vector, template <typename> class SearchMethod=bucket_search_serial, typename TRAITS_USER=Traits<VECTOR> > 
class Particles {
public:

    ///
    /// This type
    typedef Particles<VAR,D,VECTOR,SearchMethod,TRAITS_USER> particles_type;

    ///
    /// The traits type used to build up the Particle container.
    /// Contains Level 0 vector class and dimension information
    typedef TraitsCommon<VAR,D,TRAITS_USER> traits_type;

    /// 
    /// a tuple type containing value_types for each Variable
    typedef typename traits_type::value_type value_type;

    /// 
    /// a tuple type containing references to value_types for each Variable
    typedef typename traits_type::reference reference;

    /// 
    /// a tuple type containing const_references to value_types for each Variable
    typedef typename traits_type::const_reference const_reference;

    /// 
    /// a tuple type containing raw references to value_types for each Variable
    typedef typename traits_type::raw_pointer raw_pointer;

    /// 
    /// a tuple type containing raw references to value_types for each Variable
    typedef typename traits_type::raw_reference raw_reference;

    /// 
    /// a tuple type containing raw references to value_types for each Variable
    typedef typename traits_type::raw_const_reference raw_const_reference;
   
    /// 
    /// type used to hold data (a tuple of vectors)
    typedef typename traits_type::data_type data_type;

    /// 
    /// type used to hold and return size information
    typedef typename traits_type::size_type size_type;

    /// 
    /// type used to hold and return the difference of iterator, or distance between
    typedef typename traits_type::size_type difference_type;

    /// 
    /// non-const iterator type
    /// \see zip_iterator
    typedef typename traits_type::iterator iterator;

    /// 
    /// const iterator type
    /// \see zip_iterator
    typedef typename traits_type::const_iterator const_iterator;

    ///
    /// the neighbourhood data structure that the particles will be embedded into
    typedef SearchMethod<traits_type> search_type;

    ///
    /// the query class that is associated with search_type
    typedef typename search_type::query_type query_type;

    /// a boost mpl vector type containing a vector of Variable 
    /// attached to the particles (includes position, id and 
    /// alive flag as well as all user-supplied variables)
    typedef typename traits_type::mpl_type_vector mpl_type_vector;

    template <typename T>
    using elem_by_type = detail::get_elem_by_type<T,mpl_type_vector>;
    template <typename T>
    using return_type = typename detail::getter_helper<typename data_type::tuple_type>::template return_type<elem_by_type<T>::index>;

    ///
    /// the number of spatial dimensions 
    static const unsigned int dimension = D;

    ///
    /// a type to store a vector of doubles with given dimension
    typedef Vector<double,dimension> double_d;

    ///
    /// a type to store a vector of doubles with given dimension
    typedef Vector<double,dimension> int_d;

    ///
    /// a type to store a vector of bool with given dimension
    typedef Vector<bool,dimension> bool_d;

    ///
    /// the tag type for the default position variable
    typedef typename traits_type::position position;


    /// Contructs an empty container with no searching or id tracking enabled
    Particles():
        next_id(0),
        searchable(false),
        seed(time(NULL))
    {}

    /// Constructs a container with `size` particles
    Particles(const size_t size):
        next_id(0),
        searchable(false),
        seed(time(NULL))
    {
        resize(size);
    }

    /// copy-constructor. performs deep copying of all particles
    Particles(const particles_type &other):
            data(other.data),
            search(other.search),
            next_id(other.next_id),
            searchable(other.searchable),
            seed(other.seed)
    {}

    /// range-based copy-constructor. performs deep copying of all 
    /// particles from \p first to \p last
    Particles(iterator first, iterator last):
        data(traits_type::construct(first,last)),
        searchable(false),
        seed(0)
    {}

    
    //
    // STL Container
    //
    

    /// resize the continer. Note that if new particles are created, they
    /// are NOT added to the neighbour search structure, and might be outside
    /// the domain
    void resize(size_type n) {
        size_t old_n = this->size();
        traits_type::resize(data,n);         
        if (n > old_n) {
            const size_t *start_id_pointer = 
                iterator_to_raw_pointer(get<id>(data).begin() + old_n); 
            detail::for_each(begin()+old_n, end(), 
                detail::resize_lambda<raw_reference>(seed,next_id,start_id_pointer));
            next_id += n-old_n;
        }
    }
    
    /// push the particle \p val to the back of the container (if its within
    /// the searchable domain)
    void push_back (const value_type& val, bool update_neighbour_search=true) {
        // add val to container
        traits_type::push_back(data,val);

        // overwrite id, alive and random generator
        reference i = *(end()-1);
        Aboria::get<id>(i) = this->next_id++;
        Aboria::get<generator>(i) = generator_type((seed + uint32_t(Aboria::get<id>(i))));
        Aboria::get<alive>(i) = true;

        if (update_neighbour_search) {
            update_positions(end()-1,end());
        }
    }

    /// set the base seed of the container. Note that the random number generator for
    /// each particle is set to \p value plus the particle's id
    void set_seed(const uint32_t value) {
        seed = value;
        detail::for_each(begin(),end(),
                detail::set_seed_lambda<raw_reference>(seed));
        
    }

    /// push a new particle with position \p position
    /// to the back of the container
    void push_back(const double_d& pos) {
        value_type i;
        Aboria::get<position>(i) = pos;
        this->push_back(i);
    }

    /// push the particles in \p particles to the back of the container
    void push_back (const particles_type& particles) {
        for (const value_type& i: particles) {
            this->push_back(i,false);
        }
        update_positions(end()-particles.size(),end());
    }

    /// pop (delete) the particle at the end of the container 
    void pop_back(bool update_neighbour_search=true) {
        erase(end()-1,update_neighbour_search);
    }

    /// returns a reference to the particle at position \p idx
    reference operator[](std::size_t idx) {
        return traits_type::index(data, idx);
    }

    /// returns a const_reference to the particle at position \p idx
    const_reference operator[](std::size_t idx) const {
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

    /// returns an iterator to the beginning of the container
    const_iterator begin() const {
        return traits_type::cbegin(data);
    }

    /// returns an iterator to the end of the container
    const_iterator end() const {
        return traits_type::cend(data);
    }

    /// returns a const_iterator to the beginning of the container
    const_iterator cbegin() const {
        return traits_type::cbegin(data);
    }

    /// returns an iterator to the end of the container
    const_iterator cend() const {
        return traits_type::cend(data);
    }

    /// sets container to empty and deletes all particles
    void clear() {
        return traits_type::clear(data);
    }

    /// erase the particle pointed to by the iterator \p i.
    /// NOTE: This will potentially reorder the particles
    /// if neighbourhood searching is on, then this is updated
    iterator erase (iterator i, const bool update_neighbour_search = true) {
        const size_t i_position = i-begin();
        *get<alive>(i) = false;
        update_positions(i,end());
        return begin()+i_position;

        /*
        if (update_neighbour_search) {
            search.before_delete_particles_range(i_position,1);
        }

        if (i_position == size()-1) {
            // just pop off back element
            traits_type::pop_back(data);
            i = end();
        } else {
            if (search.ordered()) {
                // have to maintain order and move everthing...
                traits_type::erase(data,i);
            } else {
                // copy end element to i and pop off end element
                *i = *(end()-1);
                traits_type::pop_back(data);
            }
        }
            
        if (update_neighbour_search) {
            if (search.after_delete_particles_range(begin(),end(),i_position,1)) {
                reorder(search.get_order().begin(),search.get_order().end());
            }
        }
        */
    }


    /// erase the particles between the iterators \p first
    /// and \p last
    /// \see erase(iterator)
    iterator erase (iterator first, iterator last, const bool update_neighbour_search = true) {
        const size_t index_end = last-begin();
        detail::fill(get<alive>(first),get<alive>(last),false);
        update_positions(first,end());
        return begin() + index_end;


        /*
        const size_t n_before_range = first-begin();
        const size_t n = last-first;
        const size_t n_after_range = end()-last;

        if (update_neighbour_search) {
            search.before_delete_particles_range(n_before_range,n);
        }

        if (n_after_range > n && !search.ordered()) {
            // move elements at end to deleted region
            detail::copy(end()-n,end(),first);
            traits_type::resize(data,size()-n);         
        } else {
            // order is maintained
            traits_type::erase(data,first,last);
        }
        
        if (update_neighbour_search) {
            if (search.after_delete_particles_range(begin(),end(),n_before_range,n)) {
                reorder(search.get_order().begin(),search.get_order().end());
            }
        }
        
        return end()-n_after_range; 
        */
    }

    /// insert a particle \p val into the container at \p position
    iterator insert (iterator position, const value_type& val) {
        return traits_type::insert(data,position,val);
    }

    /// insert a \p n copies of the particle \p val into the container at \p position
    void insert (iterator position, size_type n, const value_type& val) {
        traits_type::insert(data,position,n,val);
    }

    /// insert a range of particles pointed to by \p first and \p last at \p position 
    template <class InputIterator>
    iterator insert (iterator position, InputIterator first, InputIterator last) {
        return insert_dispatch(position, first, last, std::integral_constant<bool,
                            std::is_same<InputIterator,iterator>::value 
                            || std::is_same<InputIterator,const_iterator>::value>());
    }

    /// return the total number of particles in the container
    size_t size() const {
        return Aboria::get<position>(data).size();
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
    void init_neighbour_search(const double_d& low, const double_d& high, const bool_d& periodic,
                                const unsigned int n_particles_in_leaf=10) {
        LOG(2, "Particles:init_neighbour_search: low = "<<low<<" high = "<<high<<" periodic = "<<periodic<<" n_particles_in_leaf = "<<n_particles_in_leaf);

        search.set_domain(low,high,periodic,n_particles_in_leaf);
        update_positions(begin(),end());

        searchable = true;
    }

    void init_id_search() {
        LOG(2, "Particles:init_id_search");
        search.init_id_map();
        update_positions(begin(),end());
        searchable = true;
    }

    const query_type& get_query() const {
        ASSERT(searchable,"init_neighbour_search not called on this particle set");
        return search.get_query();
    }


    double_d correct_dx_for_periodicity(const double_d& uncorrected_dx) const {
        double_d dx = uncorrected_dx;
        //double_d domain_width = get_max()-get_min();
        const bool_d& periodic = get_periodic();
        for (size_t d=0; d<traits_type::dimension; ++d) {
            if (periodic[d]) {
                const double domain_width = get_max()[d]-get_min()[d];
                while (dx[d] > domain_width/2) dx[d] -= domain_width;
                while (dx[d] <= -domain_width/2) dx[d] += domain_width;
            }
        }
        return dx;
    }

    /// return the length scale of the neighbourhood search
    /// \see init_neighbour_search()
    double_d get_lengthscale() const {
        return search.get_min_bucket_size();
    }

    /// return the lower extent of the neighbourhood search
    /// \see init_neighbour_search()
    const double_d& get_min() const {
        return search.get_min();
    }

    /// return the upper extent of the neighbourhood search
    /// \see init_neighbour_search()
    const double_d& get_max() const {
        return search.get_max();
    }

    /// return the periodicty of the neighbourhood search
    /// \see init_neighbour_search()
    const bool_d& get_periodic() const {
        return search.get_periodic();
    }

    /// Update the neighbourhood search data. This function must be
    /// called after altering the particle positions (e.g. with 
    /// `set<position>(particle,new_position)`) in order for accurate
    /// neighbourhood searching
    /// \see get_neighbours()
    void update_positions(iterator update_begin, iterator update_end) {
        if (search.update_positions(begin(),end(),update_begin,update_end)) {
            reorder(update_begin,update_end,
                    search.get_alive_indicies().begin(),
                    search.get_alive_indicies().end());
        }
    }

    void update_positions() {
        update_positions(begin(),end());
    }

    
    // Need to be mark as device to enable get functions being device/host
    CUDA_HOST_DEVICE
    const typename data_type::tuple_type & get_tuple() const { 
        return data.get_tuple(); 
    }

    // Need to be mark as device to enable get functions being device/host
    CUDA_HOST_DEVICE
    typename data_type::tuple_type & get_tuple() { 
        return data.get_tuple(); 
    }

    /*
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
    */

    friend std::ostream& operator<< (std::ostream& stream, const Particles& particles) {
        traits_type::header_to_stream(stream);
        stream << '\n';
        for (const_iterator i=particles.cbegin(); i!=particles.cend(); ++i) { 
            traits_type::to_stream(i,stream);
            stream << '\n';
        }
        return stream;
    }

    friend std::istream& operator>> (std::istream& stream, Particles& particles) {
        stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        for (iterator i=particles.begin(); i!=particles.end(); ++i) { 
            traits_type::from_stream(i);
        }
        return stream;
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        traits_type::serialize(data,ar,version);
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

        constexpr size_t dn = mpl::size<mpl_type_vector>::type::value;
        vtkSmartPointer<vtkFloatArray> datas[dn];
        mpl::for_each<mpl::range_c<int,1,dn> > (
                detail::setup_datas_for_writing<reference>(n,datas,grid)
                );
        points->SetNumberOfPoints(n);
        cells->Reset();
        cell_types->Reset();
        int j = 0;

        double write_point[3];
        const unsigned int max_d = std::min(3u,D);
        for(auto i: *this) {
            const int index = j++;
            //std::cout <<"copying point at "<<i.get_position()<<" with id = "<<i.get_id()<<std::endl;
            const double_d &r = get<position>(i);
            for (int d=0; d<max_d; ++d) {
                write_point[d] = r[d];
            }
            points->SetPoint(index,write_point);
            cells->InsertNextCell(1);
            cells->InsertCellPoint(index);
            cell_types->InsertNextTuple1(1);

            mpl::for_each<mpl::range_c<int,1,dn> > (
                    detail::write_from_tuple<reference>(
                        i.get_tuple(),
                        index,
                        datas,
                        seed + uint32_t(Aboria::get<id>(i))
                        )
                    );
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
            constexpr size_t dn = mpl::size<mpl_type_vector>::type::value;

            vtkSmartPointer<vtkFloatArray> datas[dn];

            const vtkIdType n = points->GetNumberOfPoints();

            mpl::for_each<mpl::range_c<int,1,dn> > (
                    detail::setup_datas_for_reading<reference>(n,datas,grid)
                    );

            this->clear();

            const unsigned int max_d = std::min(3u,traits_type::dimension);
            for (int j = 0; j < n; ++j) {
                value_type particle;
                const double *r = points->GetPoint(j);
                for (int d=0; d<max_d; ++d) {
                    get<position>(particle)[d] = r[d];
                }
                mpl::for_each<mpl::range_c<int,1,dn> > (
                        detail::read_into_tuple<reference>(particle.get_tuple(),j,datas)
                        );
                this->push_back(particle);
            }
        }
#endif

private:
    typedef typename traits_type::vector_unsigned_int vector_unsigned_int;
    typedef typename traits_type::vector_int vector_int;


    // gather the given update range using the order given 
    void reorder(iterator update_begin, iterator update_end, 
                 const typename vector_int::const_iterator& order_start,
                 const typename vector_int::const_iterator& order_end) {
        LOG(2,"Particles: reordering particles");
        ASSERT(update_end==end(),"if triggering a reorder, should be updating the end");
        const size_t n_update = update_end-update_begin;
        const size_t n_alive = order_end-order_start;
        const size_t old_n = size();
        const size_t new_n = old_n-(n_update-n_alive);
        if (n_alive > old_n/2) {
            traits_type::resize(other_data,new_n);
            // copy non-update region to other data buffer
            std::copy(begin(),update_begin,traits_type::begin(other_data));
            // gather update_region according to order to other data buffer
            detail::gather(order_start,order_end,
                    traits_type::begin(data),
                    traits_type::begin(other_data)+(update_begin-begin()));
            // swap to using other data buffer
            data.swap(other_data);
            search.update_iterators(begin(),end());
        } else {
            traits_type::resize(other_data,n_alive);
            // gather update_region to other buffer
            detail::gather(order_start,order_end,
                    traits_type::begin(data),
                    traits_type::begin(other_data));
            traits_type::resize(data,new_n);         
            // copy other buffer back to current data buffer
            detail::copy(traits_type::begin(other_data),
                    traits_type::end(other_data),
                    update_begin);
            search.update_iterators(begin(),end());
        }

 
            /*
            if (n_update > n/2) {
                traits_type::resize(other_data,new_n);         
                detail::scatter_if(begin(),end(),
                        transform_iterator{my_index - alive_sum}
                        get<alive>(begin()),
                        traits_type::begin(other_data));
                data.swap(other_data);
                search.update_iterators(begin(),end());
            } else {
                const size_t i_update = alive_sum_start-alive_sum_begin;
                if (other_data.size() < n_update) other_data.resize(n_update);
                detail::scatter_if(begin()+i_update,end(),
                        transform_iterator{my_index - alive_sum},
                        get<alive>(begin()),
                        traits_type::begin(other_data));
                traits_type::resize(data,new_n);         
                detail::copy(traits_type::begin(other_data),
                             traits_type::begin(data)+i_update);
                search.update_iterators(begin(),end());
            }
            */
    }

    template <class InputIterator>
    iterator insert_dispatch (iterator position, InputIterator first, InputIterator last, 
            std::false_type) {
        //TODO: this is very inefficient, but need to handle normal iterators
        //to value_type. Could improve....
        if (first==last) return position;
        for (InputIterator it = last; it != first; --it) {
            position = insert(position, *it);
        }
        return insert(position, *first);
    }

    template <class InputIterator>
    iterator insert_dispatch (iterator position, InputIterator first, InputIterator last, 
            std::true_type) {
        return traits_type::insert(data,position,first,last);
    }


    data_type data;
    data_type other_data;
    int next_id;
    bool searchable;
    uint32_t seed;
    search_type search;
    vector_int m_delete_indicies;


#ifdef HAVE_VTK
    vtkSmartPointer<vtkUnstructuredGrid> cache_grid;
#endif
};



}


#endif /* SPECIES_H_ */
