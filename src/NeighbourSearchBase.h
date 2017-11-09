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

//
// Acknowledgement: This source was modified from the Thrust example bucket_sort2d.cu
//


#ifndef NEIGHBOUR_SEARCH_BASE_H_
#define NEIGHBOUR_SEARCH_BASE_H_

#include "detail/Algorithms.h"
#include "detail/SpatialUtil.h"
#include "detail/Distance.h"
#include "Traits.h"
#include "CudaInclude.h"
#include "Vector.h"
#include "Get.h"
#include "Log.h"
#include <stack>

namespace Aboria {


template <typename IteratorType>
struct iterator_range_with_transpose {
    typedef IteratorType iterator;
    typedef typename iterator::traits_type traits_type;
    typedef typename traits_type::double_d double_d;
    IteratorType m_begin;
    IteratorType m_end;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range_with_transpose()
    {}
    /*
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range_with_transpose(IteratorType&& begin, IteratorType&& end, const double_d& transpose):
        m_begin(std::move(begin)),m_end(std::move(end)),m_transpose(transpose) 
    {}
    */
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range_with_transpose(const IteratorType& begin, const IteratorType& end, const double_d &transpose):
        m_begin(begin),m_end(end),m_transpose(transpose) 
    {}

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range_with_transpose(const IteratorType& begin, const IteratorType& end):
        m_begin(begin),m_end(end),m_transpose(0) 
    {}

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const IteratorType &begin() const { return m_begin; }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const IteratorType &end() const { return m_end; }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    IteratorType &begin() { return m_begin; }
    
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    IteratorType &end() { return m_end; }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const double_d& get_transpose() { return m_transpose; }
    double_d m_transpose;

};

template <typename IteratorType>
struct iterator_range{
    typedef IteratorType iterator;
    IteratorType m_begin;
    IteratorType m_end;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range()
    {}

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range(const iterator_range& arg) = default;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range(iterator_range&& arg) = default;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range& operator=(const iterator_range& ) = default;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator_range(const IteratorType& begin, const IteratorType& end):
        m_begin(begin),m_end(end)
    {}

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const IteratorType &begin() const { return m_begin; }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    const IteratorType &end() const { return m_end; }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    IteratorType &begin() { return m_begin; }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    IteratorType &end() { return m_end; }
};

template <typename IteratorType>
CUDA_HOST_DEVICE
iterator_range<IteratorType> make_iterator_range(IteratorType&& begin, IteratorType&& end) {
    return iterator_range<IteratorType>(begin,end);
}

template <typename Derived, typename Traits, typename QueryType>
class neighbour_search_base {
public:

    typedef QueryType query_type;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::iterator iterator;
    typedef typename Traits::vector_unsigned_int vector_unsigned_int;
    typedef typename Traits::vector_int vector_int;
    typedef typename Traits::reference reference;
    typedef typename Traits::raw_reference raw_reference;

    const Derived& cast() const { return static_cast<const Derived&>(*this); }
    Derived& cast() { return static_cast<Derived&>(*this); }

    neighbour_search_base():m_id_map(false) {
        LOG_CUDA(2,"neighbour_search_base: constructor, setting default domain");
        const double min = std::numeric_limits<double>::min();
        const double max = std::numeric_limits<double>::max();
        set_domain(double_d(min/3.0),double_d(max/3.0),bool_d(false),10,false); 
    };

    static constexpr bool ordered() {
        return true;
    }

    template <unsigned int D, typename Reference>
    struct enforce_domain_lambda {
        typedef Vector<double,D> double_d;
        typedef Vector<bool,D> bool_d;
        typedef position_d<D> position;
        static const unsigned int dimension = D;
        const double_d low,high;
        const bool_d periodic;

        enforce_domain_lambda(const double_d &low, const double_d &high, const bool_d &periodic):
            low(low),high(high),periodic(periodic) {}

        CUDA_HOST_DEVICE
        void operator()(Reference i) const {
            double_d r = Aboria::get<position>(i);
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
#ifdef __CUDA_ARCH__
                        LOG_CUDA(2,"removing particle");
#else
                        LOG(2,"removing particle with r = "<<r);
#endif
                        Aboria::get<alive>(i) = uint8_t(false);
                    }
                }
            }
            Aboria::get<position>(i) = r;
        }
    };


    /// resets the domain extents, periodicity and bucket size
    /// \param low the lower extent of the search domain
    /// \param high the upper extent of the search domain
    /// \param _max_interaction_radius the side length of each bucket
    /// \param periodic a boolean vector indicating wether each dimension
    void set_domain(const double_d &min_in, const double_d &max_in, const bool_d& periodic_in, const double n_particles_in_leaf=10, const bool not_in_constructor=true) {
        LOG(2,"neighbour_search_base: set_domain:");
        m_domain_has_been_set = not_in_constructor;
        m_bounds.bmin = min_in;
        m_bounds.bmax = max_in;
        m_periodic = periodic_in;
        m_n_particles_in_leaf = n_particles_in_leaf; 
        if (not_in_constructor) {
            cast().set_domain_impl();
        }
        LOG(2,"\tbounds = "<<m_bounds);
	    LOG(2,"\tparticles_in_leaf = "<<m_n_particles_in_leaf);
        LOG(2,"\tperiodic = "<<m_periodic);
    }

    size_t find_id_map(const size_t id) const {
        const size_t n = m_particles_end-m_particles_begin;
        return detail::lower_bound(m_id_map_key.begin(),
                                   m_id_map_key.begin()+n,
                                  id) 
                                - m_id_map_key.begin();
    }

    void init_id_map() {
        m_id_map = true;
        m_id_map_key.clear();
        m_id_map_value.clear();
    }

    /*
    void init_id_map(const bool reorder=false) {
        m_id_map = true;
	    LOG(2,"neighbour_search_base: init_id_map");
        const size_t n = m_particles_end - m_particles_begin;
        m_id_map_key.resize(n);
        m_id_map_value.resize(n);
        
        if (reorder) {
            detail::gather(m_alive_indices.begin(),m_alive_indices.end(),
                           get<id>(m_particles_begin),
                           m_id_map_key.begin());
        } else {
            detail::copy(get<id>(m_particles_begin),get<id>(m_particles_end),
                     m_id_map_key.begin());

        }
        detail::sequence(m_id_map_value.begin(),m_id_map_value.end());
        detail::sort_by_key(m_id_map_key.begin(),m_id_map_key.end(),
                            m_id_map_value.begin());

        if (ABORIA_LOG_LEVEL >= 4) { 
            print_id_map();
        }
    }

    void update_id_map(const size_t dist) {
	    LOG(2,"neighbour_search_base: update_id_map:"<<dist);
        const size_t n = m_particles_end - m_particles_begin;
        m_id_map_key.resize(n);
        m_id_map_value.resize(n);

        auto map_key_start = m_id_map_key.end()-dist;
        auto map_value_start = m_id_map_value.end()-dist;
        auto particles_start = m_particles_end-dist;

        detail::copy(get<id>(particles_start),get<id>(m_particles_end),map_key_start);
        detail::sequence(map_value_start,m_id_map_value.end(),n-dist);
        detail::sort_by_key(map_key_start,m_id_map_key.end(),map_value_start);
    }

    void copy_id_map(const size_t from_id, const size_t to_id, const int to_index) {
        const size_t from_map_index = find_id_map(from_id);
        const size_t to_map_index = find_id_map(to_id);

	    LOG(2,"neighbour_search_base: copy_id_map: from_id="<<from_id<<" to_id="<<to_id<<" to_index="<<to_index<<" from_map_index="<<from_map_index<<" to_map_index="<<to_map_index);

        ASSERT(from_map_index < (m_particles_end-m_particles_begin), "id not found");
        ASSERT(to_map_index < (m_particles_end-m_particles_begin), "id not found");
        ASSERT(m_id_map_value[to_map_index] == to_index,"found index is incorrect");

        // copy_from id gets a new index
        m_id_map_value[from_map_index] = to_index;
        // mark to id for deletion
        m_id_map_value[to_map_index] = -1;
    }

    void delete_marked_id_map() {
	    LOG(2,"neighbour_search_base: delete_marked_id_map");
        // now delete them
        // TODO: would be usefule to have an aboria `make_zip_iterator`
        typedef zip_iterator<typename Traits::template tuple<
                                typename vector_int::iterator, 
                                typename vector_int::iterator>,
                             mpl::vector<>> pair_zip_type;

        auto start_zip = pair_zip_type(m_id_map_key.begin(), m_id_map_value.begin());
        auto end_zip = pair_zip_type(m_id_map_key.end(), m_id_map_value.end());

        size_t first_dead_index = detail::stable_partition(start_zip,end_zip,
                CUDA_DEVICE
                [&](typename pair_zip_type::reference const& i) {
                return detail::get_impl<1>(i.get_tuple()) >= 0;
                }) - start_zip;

        m_id_map_key.erase(m_id_map_key.begin()+first_dead_index,
                m_id_map_key.end());
        m_id_map_value.erase(m_id_map_value.begin()+first_dead_index,
                m_id_map_value.end());
    }


    void delete_id_map_in_range(const size_t i, const size_t n) {
        // TODO: would be usefule to have an aboria `make_zip_iterator`
        typedef zip_iterator<typename Traits::template tuple<
                                typename vector_int::iterator, 
                                typename vector_int::iterator>,
                             mpl::vector<>> pair_zip_type;

        
        const size_t np = m_particles_end-m_particles_begin; 
        if (n*std::log(np) < np) {
            if (n == 1) {
                const size_t is_map_index = find_id_map(get<
            detail::for_each(start_zip,end_zip,
                CUDA_DEVICE
                [&](typename pair_zip_type::reference const& p) {
                    size_t& index = detail::get_impl<1>(p.get_tuple());
                    if (index < i+n) {
                        // to be deleted
                        index = -1;
                    } else {
                        // to be shifted
                        index -= n;
                    }
                });

        auto start_zip = pair_zip_type(m_id_map_key.begin()+i, m_id_map_value.begin()+i);

        if (

        if (cast().ordered()) {
            // particles will have kept ordering when deleting, so need to shift end indices
#if defined(__CUDACC__)
            typedef typename thrust::detail::iterator_category_to_system<
                typename vector_int::iterator::iterator_category
                >::type system;
            detail::counting_iterator<unsigned int,system> count(i);
#else
            detail::counting_iterator<unsigned int> count(i);
#endif
            const query_type& query = get_query();
            const bool ordered = cast().ordered();
            detail::for_each(count,count+n,
                CUDA_DEVICE
                [&](const unsigned int index) {
                    const size_t id_map_index = query.find(index)-query.get_particles_begin();
                    if (index < i+n) {
                        // to be deleted
                        index = -1;
                    } else {
                        // to be shifted
                        index -= n;
                    }
                });
        } else {
            // just delete
            auto end_zip = pair_zip_type(m_id_map_key.end(), m_id_map_value.end());
            auto end_zip = pair_zip_type(m_id_map_key.begin()+(i+n), m_id_map_value.begin()+(i+n));
            detail::for_each(start_zip,end_zip,
                CUDA_DEVICE
                [&](typename pair_zip_type::reference const& p) {
                    size_t& index = detail::get_impl<1>(p.get_tuple());
                    index = -1;
                });
        }

        delete_marked_id_map();

    }
        */
    void print_id_map() {
        std::cout << "particle ids:\n";
        for (auto i = m_particles_begin; i != m_particles_end; ++i) {
            std::cout << *get<id>(i) << ',';
        }
        std::cout << std::endl;

        std::cout << "alive indices:\n";
        for (auto i = m_alive_indices.begin(); i != m_alive_indices.end(); ++i) {
            std::cout << *i << ',';
        }
        std::cout << std::endl;

        std::cout << "id map (id,index):\n";
        for (int i = 0; i < m_id_map_key.size(); ++i) {
            std::cout << "(" << m_id_map_key[i] << "," << m_id_map_value[i] << ")\n";
        }
        std::cout << std::endl;
    }


    // can have added new particles (so iterators might be invalid)
    // can't have less particles (erase will update search)
    // can have alive==false particles 
    // only update id_map for new particles if derived class sets reorder, or
    // particles have been added
    bool
    update_positions(iterator begin, iterator end, 
                     iterator update_begin, iterator update_end, 
                     const bool delete_dead_particles=true) {

        //if (!m_domain_has_been_set && !m_id_map) return false;

	    LOG(2,"neighbour_search_base: update_positions: updating "<<update_end-update_begin<<" points");

        const size_t previous_n = m_particles_end-m_particles_begin;
        m_particles_begin = begin;
        m_particles_end = end;
        const size_t dead_and_alive_n = end-begin;
        CHECK(dead_and_alive_n >= previous_n,"error, particles got deleted somehow");

        CHECK(!cast().ordered() || (update_begin==begin && update_end==end), 
                "ordered search data structure can only update the entire particle set");

        const size_t new_n = dead_and_alive_n-previous_n;

        // make sure update range covers new particles
        CHECK(new_n == 0 || (update_end==end && update_end-update_begin>=new_n),
                "detected "<<new_n<<" new particles, which are not covered by update range");
        
        const size_t update_start_index = update_begin-begin;
        const size_t update_end_index = update_end-begin;
        const bool update_end_point = update_end==end;
        const size_t update_n = update_end_index-update_start_index;
        if (update_n == 0) return false;

        // enforce domain
        if (m_domain_has_been_set) {
            detail::for_each(update_begin, update_end,
                enforce_domain_lambda<Traits::dimension,raw_reference>(
                    get_min(),get_max(),get_periodic()));
        }

        // m_alive_sum will hold a cummulative sum of the living
        // num_dead holds total number of the dead
        // num_alive_new holds total number of the living that are new particles
        int num_dead = 0;
        m_alive_sum.resize(update_n);
        if (delete_dead_particles)  {
            detail::exclusive_scan(
                    get<alive>(update_begin),get<alive>(update_end),
                    m_alive_sum.begin(),
                    0);
            const int num_alive = m_alive_sum.back() +        
                                  static_cast<int>(*get<alive>(update_end-1));
            num_dead = update_n-num_alive;
            /*
            if (update_n > new_n) {
                const int num_alive_old = m_alive_sum[update_n-new_n+1];
                num_alive_new = num_alive-num_alive_old;
            } else {
                num_alive_new = num_alive;
            }
            */
        } else {
            detail::sequence(m_alive_sum.begin(),m_alive_sum.end());
        }

        CHECK(update_end==end || num_dead==0, 
                "cannot delete dead points if not updating the end of the vector"); 

	    LOG(2,"neighbour_search_base: update_positions: found "<<num_dead<<" dead points, and "<<new_n<<" new points");

        // m_alive_indices holds particle set indicies that are alive
        m_alive_indices.resize(update_n-num_dead);

#if defined(__CUDACC__)
        typedef typename thrust::detail::iterator_category_to_system<
            typename vector_int::iterator::iterator_category
            >::type system;
        detail::counting_iterator<unsigned int,system> count_start(update_start_index);
        detail::counting_iterator<unsigned int,system> count_end(update_end_index);
#else
        detail::counting_iterator<unsigned int> count_start(update_start_index);
        detail::counting_iterator<unsigned int> count_end(update_end_index);
#endif

        // scatter alive indicies to m_alive_indicies
        detail::scatter_if(count_start,count_end, //items to scatter
                           m_alive_sum.begin(), // map
                           get<alive>(update_begin), // scattered if true
                           m_alive_indices.begin()
                );


        if (m_domain_has_been_set) {
            LOG(2,"neighbour_search_base: update_positions_impl:");
            cast().update_positions_impl(update_begin,update_end,new_n);
        }
        if (m_id_map) {
            LOG(2,"neighbour_search_base: update_id_map:");
            // if no new particles, no dead, no reorder, or no init than can assume that
            // previous id map is correct
            if (cast().ordered() || new_n > 0 || num_dead > 0 || m_id_map_key.size()==0) {
                m_id_map_key.resize(dead_and_alive_n-num_dead);
                m_id_map_value.resize(dead_and_alive_n-num_dead);

                // before update range
                if (update_start_index > 0) {
                    detail::sequence(m_id_map_value.begin(),
                                     m_id_map_value.begin()+update_start_index);
                    detail::copy(get<id>(begin),
                                 get<id>(begin)+update_start_index,
                                 m_id_map_key.begin());
                }


                // after update range
                ASSERT(update_end_index == dead_and_alive_n, "if not updateing last particle then should not get here");

                // update range
                /*
                detail::transform(m_alive_indices.begin(),m_alive_indices.end(),
                             m_id_map_value.begin()+update_start_index,
                             [&](const int index) { 
                                const int index_into_update = index - update_start_index;
                                const int num_dead_before_index = index_into_update - 
                                            m_alive_sum[index_into_update];
                                return index - num_dead_before_index; 
                             });
                             */
                detail::sequence(m_id_map_value.begin()+update_start_index,
                                 m_id_map_value.end(),
                                 update_start_index);
                auto raw_id = iterator_to_raw_pointer(get<id>(begin));
                detail::transform(m_alive_indices.begin(),m_alive_indices.end(),
                             m_id_map_key.begin()+update_start_index,
                             [=] CUDA_HOST_DEVICE (const int index) { 
                                return raw_id[index]; 
                             });
                detail::sort_by_key(m_id_map_key.begin(),
                                    m_id_map_key.end(),
                                    m_id_map_value.begin());
#ifndef __CUDA_ARCH__
                if (4 <= ABORIA_LOG_LEVEL) { 
                    print_id_map();
                }
#endif
            }
        }

        query_type& query = cast().get_query_impl();
        query.m_id_map_key = iterator_to_raw_pointer(m_id_map_key.begin());
        query.m_id_map_value = iterator_to_raw_pointer(m_id_map_value.begin());
        query.m_particles_begin = iterator_to_raw_pointer(m_particles_begin);
        query.m_particles_end = iterator_to_raw_pointer(m_particles_end);

        return cast().ordered() || num_dead > 0;
            
    }

    void update_iterators(iterator begin, iterator end) {
        m_particles_begin = begin;
        m_particles_end = end;
        query_type& query = cast().get_query_impl();
        query.m_particles_begin = iterator_to_raw_pointer(m_particles_begin);
        query.m_particles_end = iterator_to_raw_pointer(m_particles_end);
	    LOG(2,"neighbour_search_base: update iterators");
        cast().update_iterator_impl();
    }

    
    const query_type& get_query() const {
        return cast().get_query_impl();
    }

    const vector_int& get_alive_indicies() const {
        return m_alive_indices;
    }

    const vector_int& get_id_map() const {
        return m_id_map;
    }

    const double_d& get_min() const { return m_bounds.bmin; }
    const double_d& get_max() const { return m_bounds.bmax; }
    const bool_d& get_periodic() const { return m_periodic; }
    bool domain_has_been_set() const { return m_domain_has_been_set; }

protected:
    iterator m_particles_begin;
    iterator m_particles_end;
    vector_int m_alive_sum;
    vector_int m_alive_indices;
    vector_int m_id_map_key;
    vector_int m_id_map_value;
    bool m_id_map;
    bool_d m_periodic;
    bool m_domain_has_been_set;
    detail::bbox<Traits::dimension> m_bounds;
    double m_n_particles_in_leaf; 
};

// assume that these iterators, and query functions, are only called from device code
template <typename Traits>
class ranges_iterator {
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::value_type p_value_type;
    typedef typename Traits::raw_reference p_reference;
    typedef typename Traits::raw_pointer p_pointer;

public:
    typedef Traits traits_type;
    typedef const p_pointer pointer;
	typedef std::random_access_iterator_tag iterator_category;
    typedef const p_reference reference;
    typedef const p_reference value_type;
	typedef std::ptrdiff_t difference_type;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    ranges_iterator() {}

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    ranges_iterator(const p_pointer& begin):
        m_current_p(begin)
    {}

        
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    reference operator ->() const {
        return dereference();
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    ranges_iterator& operator++() {
        increment();
        return *this;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    ranges_iterator operator++(int) {
        ranges_iterator tmp(*this);
        operator++();
        return tmp;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    ranges_iterator operator+(int n) {
        ranges_iterator tmp(*this);
        tmp.increment(n);
        return tmp;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    size_t operator-(ranges_iterator start) const {
        return get_by_index<0>(m_current_p) 
                - get_by_index<0>(start.m_current_p);
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    inline bool operator==(const ranges_iterator& rhs) const {
        return equal(rhs);
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    inline bool operator!=(const ranges_iterator& rhs) const {
        return !operator==(rhs);
    }

private:
    friend class boost::iterator_core_access;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    bool equal(ranges_iterator const& other) const {
        return m_current_p == other.m_current_p;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    CUDA_HOST_DEVICE
    reference dereference() const { 
        return *m_current_p; 
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    void increment() {
        ++m_current_p;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    void increment(const int n) {
        m_current_p += n;
    }

    p_pointer m_current_p;
};

/// A const iterator to a set of neighbouring points. This iterator implements
/// a STL forward iterator type
// assume that these iterators, and query functions, are only called from device code
template <typename Traits>
class linked_list_iterator {
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::value_type p_value_type;
    typedef typename Traits::raw_reference p_reference;
    typedef typename Traits::raw_pointer p_pointer;

public:
    typedef Traits traits_type;
    typedef const p_pointer pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const p_reference reference;
    typedef const p_reference value_type;
	typedef std::ptrdiff_t difference_type;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    linked_list_iterator(): 
        m_current_index(detail::get_empty_id()) {
#if defined(__CUDA_ARCH__)
            CHECK_CUDA((!std::is_same<typename Traits::template vector<double>,
                                      std::vector<double>>::value),
                       "Cannot use std::vector in device code");
#endif
    }

   
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    linked_list_iterator(
            const int index,
            const p_pointer& particles_begin,
            int* const linked_list_begin):
        m_current_index(index),
        m_particles_begin(particles_begin),
        m_linked_list_begin(linked_list_begin)
    {}

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    linked_list_iterator(const linked_list_iterator& other):
        m_current_index(other.m_current_index),
        m_particles_begin(other.m_particles_begin),
        m_linked_list_begin(other.m_linked_list_begin)
    {}

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    void operator=(const linked_list_iterator& other) {
        m_current_index = other.m_current_index;
        if (get_by_index<0>(m_particles_begin) != 
            get_by_index<0>(other.m_particles_begin)) {
            m_particles_begin = other.m_particles_begin;
        }
        m_linked_list_begin = other.m_linked_list_begin;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    reference operator ->() {
        return dereference();
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    linked_list_iterator& operator++() {
        increment();
        return *this;
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    linked_list_iterator operator++(int) {
        linked_list_iterator tmp(*this);
        operator++();
        return tmp;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    linked_list_iterator operator+(int n) {
        linked_list_iterator tmp(*this);
        for (int i = 0; i < n; ++i) {
            tmp.increment();
        }
        return tmp;
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    size_t operator-(linked_list_iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    inline bool operator==(const linked_list_iterator& rhs) const {
        return equal(rhs);
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    inline bool operator!=(const linked_list_iterator& rhs) const {
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    bool increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (linked_list_iterator):"); 
#endif
        if (m_current_index != detail::get_empty_id()) {
            m_current_index = m_linked_list_begin[m_current_index];
#ifndef __CUDA_ARCH__
            LOG(4,"\tgoing to new particle m_current_index = "<<m_current_index);
#endif
        } 

#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (linked_list_iterator): m_current_index = "<<m_current_index); 
#endif
        if (m_current_index == detail::get_empty_id()) {
            return false;
        } else {
            return true;
        }
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    bool equal(linked_list_iterator const& other) const {
        return m_current_index == other.m_current_index;
    }


    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    reference dereference() const
    { return *(m_particles_begin + m_current_index); }


    int m_current_index;
    p_pointer m_particles_begin;
    int* m_linked_list_begin;

};


template <typename Traits, typename Iterator>
class index_vector_iterator {
    typedef index_vector_iterator<Traits,Iterator> iterator;
    typedef typename Traits::double_d double_d;
    typedef typename Traits::bool_d bool_d;
    typedef typename Traits::value_type p_value_type;
    typedef typename Traits::raw_reference p_reference;
    typedef typename Traits::raw_pointer p_pointer;

public:
    typedef Traits traits_type;
    typedef const p_pointer pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef const p_reference reference;
    typedef const p_reference value_type;
	typedef std::ptrdiff_t difference_type;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    index_vector_iterator() 
    {}
       
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    index_vector_iterator(
            Iterator begin,
            const p_pointer& particles_begin):
        m_current_index(begin),
        m_particles_begin(particles_begin)
    {}


    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    reference operator ->() {
        return dereference();
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator& operator++() {
        increment();
        return *this;
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    size_t operator-(iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }
    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    void increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (index_vector_iterator):"); 
#endif
        ++m_current_index;
#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (index_vector_iterator): m_current_index = "<<m_current_index); 
#endif
    }

    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        return m_current_index == other.m_current_index;
    }


    ABORIA_HOST_DEVICE_IGNORE_WARN
    CUDA_HOST_DEVICE
    reference dereference() const
    { return *(m_particles_begin + *m_current_index); }


    Iterator m_current_index;
    p_pointer m_particles_begin;
};

template<class T, std::size_t N>
class static_vector
{
    // properly aligned uninitialized storage for N T's
    typename std::aligned_storage<sizeof(T), alignof(T)>::type data[N];
    size_t m_size = 0;
 
public:
    CUDA_HOST_DEVICE
    void push_back (const T& val) {
        ASSERT_CUDA(m_size < N);
        new(data+m_size) T(val);
        ++m_size;
    }

    CUDA_HOST_DEVICE
    void push_back (T&& val) {
        ASSERT_CUDA(m_size < N);
        new(data+m_size) T(val);
        ++m_size;
    }

    // Create an object in aligned storage
    CUDA_HOST_DEVICE
    template<typename ...Args> void emplace_back(Args&&... args) {
        ASSERT_CUDA(m_size < N);
        new(data+m_size) T(std::forward<Args>(args)...);
        ++m_size;
    }

    CUDA_HOST_DEVICE
    void resize (size_t n) {
        m_size = n;
    }

    CUDA_HOST_DEVICE
    size_t size () const {
        return m_size;
    }

    CUDA_HOST_DEVICE
    bool empty() const {
        return m_size==0;
    }

    CUDA_HOST_DEVICE
    void pop_back() {
        --m_size;
        reinterpret_cast<T*>(data+m_size)->~T();
    }

    CUDA_HOST_DEVICE
    T& back() {
        return *reinterpret_cast<T*>(data+m_size-1);
    }

    CUDA_HOST_DEVICE
    const T& back() const {
        return *reinterpret_cast<const T*>(data+m_size-1);
    }

    // Access an object in aligned storage
    CUDA_HOST_DEVICE
    const T& operator[](std::size_t pos) const {
        return *reinterpret_cast<const T*>(data+pos);
    }

    CUDA_HOST_DEVICE
    T& operator[](std::size_t pos) {
        return *reinterpret_cast<T*>(data+pos);
    }
 
    // Delete objects from aligned storage
    CUDA_HOST_DEVICE
    ~static_vector() {
        for(std::size_t pos = 0; pos < m_size; ++pos) {
            reinterpret_cast<T*>(data+pos)->~T();
        }
    }

};

template <typename Query>
class depth_first_iterator {
    typedef depth_first_iterator<Query> iterator;
    typedef typename Query::child_iterator child_iterator;
    static const unsigned int dimension = Query::dimension;
    typedef Vector<double,dimension> double_d;
    typedef Vector<int,dimension> int_d;
    typedef detail::bbox<dimension> box_type;

public:
    typedef typename child_iterator::value_type value_type;
    typedef child_iterator pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef typename child_iterator::reference reference;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    depth_first_iterator() {}
       
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    CUDA_HOST_DEVICE
    depth_first_iterator(const child_iterator& start_node,
                        const unsigned tree_depth,
                        const Query *query
                  ):
        m_query(query)
    {
        ASSERT_CUDA(tree_depth <= m_stack_max_size);
        if (start_node != false) {
            m_stack.push_back(start_node);
        } else {
#ifndef __CUDA_ARCH__
            LOG(3,"\tdepth_first_iterator (constructor): start is false, no children to search.");
#endif
        }
    }

    
    CUDA_HOST_DEVICE
    depth_first_iterator(const iterator& copy):
        m_query(copy.m_query)
    {
        for (int i = 0; i < copy.m_stack.size(); ++i) {
            m_stack.push_back(copy.m_stack[i]);
        }
    }

    ~depth_first_iterator() {
    }

    CUDA_HOST_DEVICE
    iterator& operator=(const iterator& copy) {
        m_stack.resize(copy.m_stack.size());
        for (int i = 0; i < copy.m_stack.size(); ++i) {
            m_stack[i] = copy.m_stack[i];
        }

        m_query = copy.m_query;
        return *this;
    }


    box_type get_bounds() const {
        return m_stack.back().get_bounds();
    }

    const child_iterator& get_child_iterator() const {
        return m_stack.back();
    }

    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }
    CUDA_HOST_DEVICE
    reference operator ->() {
        return dereference();
    }
    CUDA_HOST_DEVICE
    iterator& operator++() {
        increment();
        return *this;
    }
    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }
    CUDA_HOST_DEVICE
    size_t operator-(iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }
    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }
    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;

    

    CUDA_HOST_DEVICE
    void increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (depth_first_iterator): depth = "<<m_stack.size()<<" top child number "<<m_stack.back().get_child_number()); 
#endif
        if (m_query->is_leaf_node(*m_stack.back())) {
            ++m_stack.back();
            while (!m_stack.empty() && m_stack.back() == false) {
                LOG(4,"\tpop stack (depth_first_iterator):"); 
                m_stack.pop_back();
            }
        } else {
            LOG(4,"\tpush stack (depth_first_iterator):"); 
            m_stack.push_back(m_query->get_children(m_stack.back()++));
        }

#ifndef __CUDA_ARCH__
        LOG(4,"\tend increment (depth_first_iterator):"); 
#endif
    }

    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        if (m_stack.empty() || other.m_stack.empty()) {
            return m_stack.empty() == other.m_stack.empty();
        } else {
            return m_stack.back() == other.m_stack.back();
        }
    }

    CUDA_HOST_DEVICE
    reference dereference() const
    { return *m_stack.back(); }

    
    const static unsigned m_stack_max_size = Query::m_max_tree_depth;
    static_vector<child_iterator,m_stack_max_size> m_stack;
    const Query *m_query;
};

template <typename Query, int LNormNumber>
class tree_query_iterator {
    typedef tree_query_iterator<Query,LNormNumber> iterator;
    typedef typename Query::child_iterator child_iterator;
    static const unsigned int dimension = Query::dimension;
    typedef Vector<double,dimension> double_d;
    typedef Vector<int,dimension> int_d;
    typedef detail::bbox<dimension> box_type;

public:
    typedef typename child_iterator::value_type value_type;
    typedef typename child_iterator::pointer pointer;
	typedef std::forward_iterator_tag iterator_category;
    typedef typename child_iterator::reference reference;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    tree_query_iterator() {}
       
    /// this constructor is used to start the iterator at the head of a bucket 
    /// list
    CUDA_HOST_DEVICE
    tree_query_iterator(const child_iterator& start,
                  const double_d& query_point,
                  const double_d& max_distance,
                  const unsigned tree_depth,
                  const Query *query,
                  const bool ordered=false
                  ):
        m_query_point(query_point),
        m_inv_max_distance(1.0/max_distance),
        m_query(query)
    {
        ASSERT_CUDA(tree_depth <= m_stack_max_size);
        if (start != false) {
            m_stack.push_back(start);
            go_to_next_leaf();
        } else {
#ifndef __CUDA_ARCH__
            LOG(3,"\ttree_query_iterator (constructor) with query pt = "<<m_query_point<<"): start is false, no children to search.");
#endif
        }

        if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
            LOG(3,"\ttree_query_iterator (constructor) with query pt = "<<m_query_point<<"): search region outside domain or no children to search.");
#endif
        } else {
#ifndef __CUDA_ARCH__
            LOG(3,"\tocttree_query_iterator (constructor) with query pt = "<<m_query_point<<"):  found bbox = "<<m_query->get_bounds(m_stack.back()));
#endif
        }
    }

    CUDA_HOST_DEVICE
    tree_query_iterator(const iterator& copy):
        m_query_point(copy.m_query_point),
        m_inv_max_distance(copy.m_inv_max_distance),
        m_query(copy.m_query)

    {
        for (int i = 0; i < copy.m_stack.size(); ++i) {
            m_stack.push_back(copy.m_stack[i]);
        }
    }

    CUDA_HOST_DEVICE
    ~tree_query_iterator() {
    }

    CUDA_HOST_DEVICE
    box_type get_bounds() const {
        return m_stack.back().get_bounds();
    }

    CUDA_HOST_DEVICE
    iterator& operator=(const iterator& copy) {
        m_query_point = copy.m_query_point;
        m_inv_max_distance = copy.m_inv_max_distance;

        m_stack.resize(copy.m_stack.size());
        for (int i = 0; i < copy.m_stack.size(); ++i) {
            m_stack[i] = copy.m_stack[i];
        }

        m_query = copy.m_query;
        return *this;
    }

    /*
    iterator& operator=(const octtree_depth_first_iterator<Query>& copy) {
        m_stack = copy.m_stack;
        go_to_next_leaf();
        return *this;
    }
    */

    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }

    CUDA_HOST_DEVICE
    reference operator ->() {
        return dereference();
    }

    CUDA_HOST_DEVICE
    iterator& operator++() {
        increment();
        return *this;
    }

    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }
    
    CUDA_HOST_DEVICE
    size_t operator-(iterator start) const {
        size_t count = 0;
        while (start != *this) {
            start++;
            count++;
        }
        return count;
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

 private:
    friend class boost::iterator_core_access;


    CUDA_HOST_DEVICE
    bool child_is_within_query(const child_iterator& node) {
        const box_type& bounds = m_query->get_bounds(node);
        double accum = 0;
        for (int j = 0; j < dimension; j++) {
            const bool less_than_bmin = m_query_point[j]<bounds.bmin[j];
            const bool more_than_bmax = m_query_point[j]>bounds.bmax[j];

            // dist 0 if between min/max, or distance to min/max if not
            const double dist = (less_than_bmin^more_than_bmax)
                                *(less_than_bmin?
                                        (bounds.bmin[j]-m_query_point[j]):
                                        (m_query_point[j]-bounds.bmax[j]));

            accum = detail::distance_helper<LNormNumber>::
                        accumulate_norm(accum,dist*m_inv_max_distance[j]); 
        }
        return (accum < 1.0);
    }

    CUDA_HOST_DEVICE
    void increment_stack() {
        while (!m_stack.empty()) {
            ++m_stack.back();
#ifndef __CUDA_ARCH__
            LOG(4,"\tincrement stack with child "<<m_stack.back().get_child_number());
#endif
            if (m_stack.back() == false) {
#ifndef __CUDA_ARCH__
                LOG(4,"\tincrement_stack: pop");
#endif
                m_stack.pop_back();
            } else {
                break;
            }
        } 
    }


    CUDA_HOST_DEVICE
    void go_to_next_leaf() {
        bool exit = m_stack.empty();
        while (!exit) {
            child_iterator& node = m_stack.back();
#ifndef __CUDA_ARCH__
            LOG(3,"\tgo_to_next_leaf with child "<<node.get_child_number()<<" with bounds "<<node.get_bounds());
#endif
            if (child_is_within_query(node)) { // could be in this child
#ifndef __CUDA_ARCH__
                LOG(4,"\tthis child is within query, so going to next child");
#endif
                if (m_query->is_leaf_node(*node)) {
                    exit = true;
                } else {
#ifndef __CUDA_ARCH__
                    LOG(4,"\tdive down");
#endif
                    m_stack.push_back(m_query->get_children(node));
                }
            } else { // not in this one, so go to next child, or go up if no more children
#ifndef __CUDA_ARCH__
                LOG(4,"\tthis child is NOT within query, so going to next child");
#endif
                increment_stack();
                exit = m_stack.empty();
            }
        }
#ifndef __CUDA_ARCH__
        if (4 <= ABORIA_LOG_LEVEL) { 
            if (m_stack.empty()) {
                LOG(4,"\tgo_to_next_leaf: stack empty, finishing");
            } else {
                LOG(4,"\tgo_to_next_leaf: found leaf, finishing");
            }
        }
#endif

    }

    CUDA_HOST_DEVICE
    void increment() {
#ifndef __CUDA_ARCH__
        LOG(4,"\tincrement (octtree_query_iterator):"); 
#endif
        increment_stack();
        go_to_next_leaf();

        if (m_stack.empty()) {
#ifndef __CUDA_ARCH__
            LOG(3,"\tend increment (octree_query_iterator): no more nodes"); 
#endif
        } else {
#ifndef __CUDA_ARCH__
            LOG(3,"\tend increment (octree_query_iterator): looking in bbox "<<m_query->get_bounds(m_stack.back())); 
#endif
        }
    }
    
    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        if (m_stack.empty() || other.m_stack.empty()) {
            return m_stack.empty() == other.m_stack.empty();
        } else {
            return m_stack.back() == other.m_stack.back();
        }
    }


    CUDA_HOST_DEVICE
    reference dereference() const
    { return *m_stack.back(); }


    //unsigned m_stack_size;
    const static unsigned m_stack_max_size = Query::m_max_tree_depth;
    static_vector<child_iterator,m_stack_max_size> m_stack;
    double_d m_query_point;
    double_d m_inv_max_distance;
    const Query *m_query;
};


template <unsigned int D>
class lattice_iterator {
    typedef lattice_iterator<D> iterator;
    typedef Vector<double,D> double_d;
    typedef Vector<int,D> int_d;

    // make a proxy int_d in case you ever
    // want to get a pointer object to the
    // reference (which are both of the
    // same type)
    struct proxy_int_d: public int_d {
        CUDA_HOST_DEVICE
        proxy_int_d():
            int_d() 
        {}

        CUDA_HOST_DEVICE
        proxy_int_d(const int_d& arg):
            int_d(arg) 
        {}

        CUDA_HOST_DEVICE
        proxy_int_d& operator&() {
            return *this;
        }
        
        CUDA_HOST_DEVICE
        const proxy_int_d& operator&() const {
            return *this;
        }

        CUDA_HOST_DEVICE
        const proxy_int_d& operator*() const {
            return *this;
        }

        CUDA_HOST_DEVICE
        proxy_int_d& operator*() {
            return *this;
        }

        CUDA_HOST_DEVICE
        const proxy_int_d* operator->() const {
            return this;
        }

        CUDA_HOST_DEVICE
        proxy_int_d* operator->() {
            return this;
        }
    };

   
    int_d m_min;
    int_d m_max;
    proxy_int_d m_index;
    int_d m_size;
    bool m_valid;
public:
    typedef proxy_int_d pointer;
	typedef std::random_access_iterator_tag iterator_category;
    typedef const proxy_int_d& reference;
    typedef proxy_int_d value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    lattice_iterator():
        m_valid(false)
    {}

    CUDA_HOST_DEVICE
    lattice_iterator(const int_d &min, 
                     const int_d &max):
        m_min(min),
        m_max(max),
        m_index(min),
        m_size(minus(max,min)),
        m_valid(true)
    {}

    CUDA_HOST_DEVICE
    lattice_iterator(const int_d &min, 
                     const int_d &max,
                     const int_d &index):
        m_min(min),
        m_max(max),
        m_index(index),
        m_size(minus(max,min)),
        m_valid(true)
    {}

    /*
    CUDA_HOST_DEVICE
    iterator& operator=(const iterator& copy) {
        m_index = copy.m_index;
        m_valid = copy.m_valid;
        ASSERT_CUDA(m_valid?(m_index >= m_min).all()&&(m_index < m_max).all():true);
        return *this;
    }
    */

    CUDA_HOST_DEVICE
    iterator& operator=(const int_d& copy) {
        m_index = copy;
        m_valid = (m_index >= m_min).all()&&(m_index < m_max).all();
        return *this;
    }

    CUDA_HOST_DEVICE
    explicit operator size_t() const {
        return collapse_index_vector(m_index);
    }

    CUDA_HOST_DEVICE
    const lattice_iterator& get_child_iterator() const {
        return *this;
    }

    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }

    CUDA_HOST_DEVICE
    reference operator ->() const {
        return dereference();
    }

    CUDA_HOST_DEVICE
    iterator& operator++() {
        increment();
        return *this;
    }

    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }

    CUDA_HOST_DEVICE
    iterator operator+(const int n) {
        iterator tmp(*this);
        tmp.increment(n);
        return tmp;
    }

    CUDA_HOST_DEVICE
    iterator& operator+=(const int n) {
        increment(n);
        return *this;
    }

    CUDA_HOST_DEVICE
    iterator& operator-=(const int n) {
        increment(-n);
        return *this;
    }

    CUDA_HOST_DEVICE
    iterator operator-(const int n) {
        iterator tmp(*this);
        tmp.increment(-n);
        return tmp;
    }

    CUDA_HOST_DEVICE
    size_t operator-(const iterator& start) const {
        int distance;
        if (!m_valid) {
            const int_d test = minus(minus(start.m_max,1),start.m_index);
            distance = start.collapse_index_vector(minus(minus(start.m_max,1),start.m_index))+1;
        } else if (!start.m_valid) {
            distance = collapse_index_vector(minus(m_index,m_min));
        } else {
            distance = collapse_index_vector(minus(m_index,start.m_index));
        }
        assert(distance > 0);
        return distance;
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const bool rhs) const {
        return equal(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator!=(const bool rhs) const {
        return !operator==(rhs);
    }

private:

    static inline 
    CUDA_HOST_DEVICE
    int_d minus(const int_d& arg1, const int_d& arg2) {
        int_d ret;
        for (size_t i = 0; i < D; ++i) {
            ret[i] = arg1[i]-arg2[i];
        }
        return ret;
    }

    static inline 
    CUDA_HOST_DEVICE
    int_d minus(const int_d& arg1, const int arg2) {
        int_d ret;
        for (size_t i = 0; i < D; ++i) {
            ret[i] = arg1[i]-arg2;
        }
        return ret;
    }

    CUDA_HOST_DEVICE
    int collapse_index_vector(const int_d &vindex) const {
        int index = 0;
        unsigned int multiplier = 1.0;
        for (int i = D-1; i>=0; --i) {
            if (i != D-1) {
                multiplier *= m_size[i+1];
            }
            index += multiplier*vindex[i];
        }
        return index;
    }

    CUDA_HOST_DEVICE
    int_d reassemble_index_vector(const int index) const {
        int_d vindex;
        int i = index;
        for (int d = D-1; d>=0; --d) {
            double div = (double)i / m_size[d];
            vindex[d] = std::round((div-std::floor(div)) * m_size[d]);
            i = std::floor(div);
        }
        return vindex;
    }

    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        if (!other.m_valid) return !m_valid;
        if (!m_valid) return !other.m_valid;
        for (size_t i = 0; i < D; ++i) {
           if (m_index[i] != other.m_index[i]) {
               return false;
           }
        }
        return true;
    }

    CUDA_HOST_DEVICE
    bool equal(const bool other) const {
        return m_valid==other;
    }

    CUDA_HOST_DEVICE
    reference dereference() const { 
        return m_index; 
    }

    CUDA_HOST_DEVICE
    void increment() {
        for (int i=D-1; i>=0; --i) {
            ++m_index[i];
            if (m_index[i] < m_max[i]) break;
            if (i != 0) {
                m_index[i] = m_min[i];
            } else {
                m_valid = false;
            }
        }
    }

    CUDA_HOST_DEVICE
    void increment(const int n) {
        if (n==1) {
            increment();
        } else {
            int collapsed_index = collapse_index_vector(m_index);
            m_index = reassemble_index_vector(collapsed_index += n);
        }
    }
};



template <typename Query, int LNormNumber>
class lattice_iterator_within_distance {
    typedef lattice_iterator_within_distance<Query,LNormNumber> iterator;
    static const unsigned int dimension = Query::dimension;
    typedef Vector<double,dimension> double_d;
    typedef Vector<int,dimension> int_d;

    // make a proxy int_d in case you ever
    // want to get a pointer object to the
    // reference (which are both of the
    // same type)
    struct proxy_int_d: public int_d {
        CUDA_HOST_DEVICE
        proxy_int_d():
            int_d() 
        {}

        CUDA_HOST_DEVICE
        proxy_int_d(const int_d& arg):
            int_d(arg) 
        {}

        CUDA_HOST_DEVICE
        proxy_int_d& operator&() {
            return *this;
        }
        
        CUDA_HOST_DEVICE
        const proxy_int_d& operator&() const {
            return *this;
        }

        CUDA_HOST_DEVICE
        const proxy_int_d& operator*() const {
            return *this;
        }

        CUDA_HOST_DEVICE
        proxy_int_d& operator*() {
            return *this;
        }

        CUDA_HOST_DEVICE
        const proxy_int_d* operator->() const {
            return this;
        }

        CUDA_HOST_DEVICE
        proxy_int_d* operator->() {
            return this;
        }
    };

    double_d m_query_point;
    double_d m_inv_max_distance;
    int m_quadrant; 
    const Query *m_query;
    bool m_valid;
    int_d m_min;
    proxy_int_d m_index;
    int_d m_base_index;
public:
    typedef proxy_int_d pointer;
	typedef std::random_access_iterator_tag iterator_category;
    typedef const proxy_int_d& reference;
    typedef proxy_int_d value_type;
	typedef std::ptrdiff_t difference_type;

    CUDA_HOST_DEVICE
    lattice_iterator_within_distance():
        m_valid(false)
    {}

    CUDA_HOST_DEVICE
    lattice_iterator_within_distance(const double_d& query_point,
                     const double_d& max_distance,
                     const Query* query):
        m_query_point(query_point),
        m_inv_max_distance(1.0/max_distance),
        m_quadrant(0),
        m_query(query),
        m_valid(true)
    {
        if (outside_domain(query_point,max_distance)) {
            m_valid = false;
        } else {
            reset_min_and_index(); 
        }
    }

    CUDA_HOST_DEVICE
    explicit operator size_t() const {
        return m_query->m_point_to_bucket_index.collapse_index_vector(m_index);
    }


    CUDA_HOST_DEVICE
    reference operator *() const {
        return dereference();
    }

    CUDA_HOST_DEVICE
    reference operator ->() const {
        return dereference();
    }

    CUDA_HOST_DEVICE
    iterator& operator++() {
        increment();
        return *this;
    }

    CUDA_HOST_DEVICE
    iterator operator++(int) {
        iterator tmp(*this);
        operator++();
        return tmp;
    }

    CUDA_HOST_DEVICE
    size_t operator-(const iterator& start) const {
        int distance = 0;
        iterator tmp = start;
        while (tmp != *this) {
            ++distance;
            ++tmp;
        }
        return distance;
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const iterator& rhs) const {
        return equal(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator==(const bool rhs) const {
        return equal(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator!=(const iterator& rhs) const {
        return !operator==(rhs);
    }

    CUDA_HOST_DEVICE
    inline bool operator!=(const bool rhs) const {
        return !operator==(rhs);
    }

private:
    CUDA_HOST_DEVICE
    bool equal(iterator const& other) const {
        if (!other.m_valid) return !m_valid;
        if (!m_valid) return !other.m_valid;
        for (size_t i = 0; i < dimension; ++i) {
           if (m_index[i] != other.m_index[i]) {
               return false;
           }
        }
        return true;
    }

    CUDA_HOST_DEVICE
    bool equal(const bool other) const {
        return m_valid==other;
    }

    CUDA_HOST_DEVICE
    reference dereference() const { 
        return m_index; 
    }


    CUDA_HOST_DEVICE
    bool ith_quadrant_bit(const int i) const {
        return (1 == ((m_quadrant >> i) & 1));
    }

    CUDA_HOST_DEVICE
    void reset_min_and_index() {
        bool no_buckets = true;

        LOG_CUDA(3,"lattice_iterator_within_distance: reset_min_and_index:begin");
        while (m_valid && no_buckets) {
            for (int i = 0; i < dimension; ++i) {
                m_min[i] = m_query->m_point_to_bucket_index.get_min_index_by_quadrant(
                        m_query_point[i],i,ith_quadrant_bit(i));
            }

            // check distance
            double accum = 0;
            for (int j = 0; j < dimension; ++j) {
                const double dist = 
                    m_query->m_point_to_bucket_index.get_dist_to_bucket(
                            m_query_point[j],m_base_index[j],m_min[j],j);
                ASSERT_CUDA(dist >= 0);
                //std::cout <<"dist= "<<dist<< " "<<dist*m_inv_max_distance[j]<< std::endl;
                accum = detail::distance_helper<LNormNumber>::
                    accumulate_norm(accum,dist*m_inv_max_distance[j]); 
            }
            //std::cout <<"accum = "<<accum<< std::endl;

            no_buckets = accum > 1.0;

            //std::cout <<" m_min = "<<m_min<<" m_quadrant = "<<m_quadrant << std::endl;

            // if good, check that this quadrant is within domain
            if (!no_buckets) {
                for (int i=0; i<dimension; i++) {
                    if (ith_quadrant_bit(i)) {
                        if (m_min[i] < 0) {
                            m_min[i] = 0;
                        } else if (m_min[i] > m_query->m_end_bucket[i]) {
                            no_buckets = true;
                            m_min[i] = m_query->m_end_bucket[i];
                        }
                    } else {
                        if (m_min[i] < 0) {
                            no_buckets = true;
                            m_min[i] = 0;
                        } else if (m_min[i] > m_query->m_end_bucket[i]) {
                            m_min[i] = m_query->m_end_bucket[i];
                        }
                    }
                }
            }
            if (no_buckets) {
                // if no buckets, move onto next quadrant
                ++m_quadrant;
                if (m_quadrant >= (1<<dimension)) {
                    m_valid = false;
                }
            } else {
                // we got a non empty quadrent, lets go!
                m_index = m_min;
            }
        }
        //std::cout <<"m_valid = "<<m_valid<<" m_min = "<<m_min<< "m_index = "<<m_index<<" m_quadrant = "<<m_quadrant << std::endl;
        LOG_CUDA(3,"lattice_iterator_within_distance: reset_min_and_index:end");
    }

    CUDA_HOST_DEVICE
    bool outside_domain(const double_d& position,const double_d& max_distance) {
        m_base_index = m_query->m_point_to_bucket_index.find_bucket_index_vector(position);
        int_d start = m_query->m_point_to_bucket_index.find_bucket_index_vector(position-max_distance);
        int_d end = m_query->m_point_to_bucket_index.find_bucket_index_vector(position+max_distance);

        bool no_buckets = false;
        for (int i=0; i<dimension; i++) {
            if (start[i] > m_query->m_end_bucket[i]) {
                no_buckets = true;
            }
            if (end[i] < 0) {
                no_buckets = true;
            }
        }
        return no_buckets;
    } 


    CUDA_HOST_DEVICE
    void increment() {
        LOG_CUDA(3,"lattice_iterator_within_distance: increment :begin");
        for (int i=dimension-1; i>=0; --i) {
            double distance = 0;

            // increment or decrement index depending on the current
            // quadrant
            bool potential_bucket = true;
            if (ith_quadrant_bit(i)) {
                ++m_index[i];
                potential_bucket = m_index[i] <= m_query->m_end_bucket[i];
            } else {
                --m_index[i];
                potential_bucket = m_index[i] >= 0;
            }

            //std::cout <<"m_min = "<<m_min<< "m_index = "<<m_index<<" m_quadrant = "<<m_quadrant << std::endl;

            // if index is outside domain don't bother calcing
            // distance
            if (potential_bucket) { 
                double accum = 0;
                for (int j = 0; j < dimension; ++j) {
                    const double dist = 
                        m_query->m_point_to_bucket_index.get_dist_to_bucket(
                                m_query_point[j],m_base_index[j],m_index[j],j);
                    ASSERT_CUDA(dist >= 0);
                    //std::cout <<"dist= "<<dist<< " "<<dist*m_inv_max_distance[j]<< std::endl;
                    accum = detail::distance_helper<LNormNumber>::
                        accumulate_norm(accum,dist*m_inv_max_distance[j]); 
                }
                //std::cout <<"accum = "<<accum<< std::endl;

                potential_bucket = accum <= 1.0;
            }

            // if bucket is still good break out of loop
            if (potential_bucket) break;

            // must be outside distance or outside domain, so reset index back to min
            m_index[i] = m_min[i];

            // if gone through all dimensions, move to next quadrant
            if (i == 0) {
                ++m_quadrant;
                if (m_quadrant < (1<<dimension)) {
                    reset_min_and_index();
                } else {
                    // if gone through all quadrants, iterator
                    // is now invalid
                    m_valid = false;
                }
            }
        }
        LOG_CUDA(3,"lattice_iterator_within_distance: increment :end");
    }

};



}


#endif /* BUCKETSEARCH_H_ */
