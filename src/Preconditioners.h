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


#ifndef PRECONDITIONERS_H_
#define PRECONDITIONERS_H_

#include <type_traits>

#ifdef HAVE_EIGEN
#include "detail/Operators.h"

namespace Aboria {


template <typename _Scalar>
class RASMPreconditioner
{
    typedef _Scalar Scalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef typename traits_type::template vector_type<
        typename std::remove_const<pointer>::type
        >::type bucket_pointer_vector_type;
    typedef typename traits_type::template vector_type<bucket_pointer_vector_type>::type connectivity_type;
    typedef typename traits_type::template vector_type<size_t>::type connectivity_type;
    typedef typename traits_type::template vector_type<size_t>::type particles_type;
    connectivity_type m_row_buffer;
    connectivity_type m_row_buckets;
    connectivity_type m_col_buffer;
    connectivity_type m_col_buckets;


  public:
    typedef typename Vector::StorageIndex StorageIndex;
    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };

    RASMPreconditioner() : m_isInitialized(false) {}

    template<typename MatType>
    explicit RASMPreconditioner(const MatType& mat) : m_invdiag(mat.cols())
    {
      compute(mat);
    }

    Index rows() const { return m_invdiag.size(); }
    Index cols() const { return m_invdiag.size(); }

    void store_domain_indicies(double_d low, double_d high) {
        const size_t domain_index = m_number_of_domains++;
        m_row_buffer.resize(m_number_of_domains);
        m_row_indicies.resize(m_number_of_domains);
        m_row_particles.push_back(a_particles);
        m_col_buffer.resize(m_number_of_domains);
        m_col_indicies.resize(m_number_of_domains);
        m_col_particles.push_back(b_particles);

        const query_range a_range = a_query.get_buckets_near_point(
                                        0.5*(high-low)+low,
                                        0.5*(high-low)+m_buffer);
        
        for(reference_a a_bucket: a_range) {
            auto a_particles = a_query.get_bucket_particles(a_bucket);
            for (particle_a_reference a_particle: a_particles) {
                const size_t index = ;
                if ((get<position>(a_particle) < low).any()
                        || (get<position>(a_particle) > high).any()) {
                    m_row_buffer[domain_index].push_back(index);
                } else {
                    m_row_indicies[domain_index].push_back(index);
                }
            }
        }

        const query_range b_range = b_query.get_buckets_near_point(
                                        0.5*(high-low)+low,
                                        0.5*(high-low)+m_buffer);

        for(reference_b b_bucket: b_range) {
            auto b_particles = b_query.get_bucket_particles(b_bucket);
            for (particle_b_reference b_particle: b_particles) {
                const size_t index = ;
                if ((get<position>(b_particle) < low).any()
                        || (get<position>(b_particle) > high).any()) {
                    m_col_buffer[domain_index].push_back(index);
                } else {
                    m_col_indicies[domain_index].push_back(index);
                }
            }
        }
    }

    template <typename query_type_a>
    size_t analyze_dive(reference_a bucket) {
        size_t count;
        if (a_query.is_leaf_node(bucket)) {
            auto particles = a_query.get_bucket_particles(bucket);
            count = std::distance(particles.begin(),particles.end());
            
        } else {
            count = factorize_dive(a_query.get_child1(&bucket))
                    + factorize_dive(a_query.get_child2(&bucket));
        }
        if (count < goal) {
            return count;
        } else {
            store_domain_indicies(
                    a_query.get_bounds_low(bucket),
                    a_query.get_bounds_high(bucket));
        }
    }

    // for a tree, use data structure to find a good division
    template <typename block_type>
    void analyze_impl_block(const Index i, const block_type& block) const {
        particle_type_a a;
        query_type_a a_query;
        typedef query_type_a::reference reference_a;
        for (reference root: a.get_root_buckets()) {
            if (a_query.is_leaf_node(root)) {
                // no tree!, just accept this one I guess
                good_bucket(root);
            } else {
                // dive tree, find bucket with given number of particles
                factorize_dive(root);
            }
        }
    }

    // for a regular grid, assume particles are evenly distributed
    template <typename block_type>
    void analyze_impl_block(const Index i, const block_type& block) const {
        particle_type_a a;
        query_type_a a_query;
        typedef query_type_a::reference reference_a;

        const double_d& low = a_query.get_bounds_low();
        const double_d& high = a_query.get_bounds_high();
        const size_t N = a.size();
        const double total_volume = (high-low).prod();
        const double box_volume = double(m_goal_n)/double(N)*total_volume;
        const double box_side_length = std::pow(box_volume,1.0/dimension);
        const unsigned_int_d size = 
                    floor((high-low)/box_side_length)
                    .template cast<unsigned int>();
        for (int i=0; i<dimension; ++i) {
            if (size[i] == 0) {
                size[i] = 1;
            }
        }
        side_length = (high-low)/size;
        iterator_range<lattice_iterator<dimension>> range(
                lattice_iterator<dimension>(int_d(0),size,int_d(0)),
                ++lattice_iterator<dimension>(int_d(0),size,size));
        for (lattice_iterator<dimension>::reference box: range) {
            store_domain_indicies(box*side_length,(box+1)*side_length);
        }
    }

    template<std::size_t... I>
    void analyze_impl(detail::index_sequence<I...>) const {
        factorize_impl_block(start_row<I>,tuple_ns::get<I*NJ+I>(m_blocks));
    }
    
    template<typename MatType>
    RASMPreconditioner& analyzePattern(const MatType& )
    {
        analyze_impl(detail::make_index_sequence<NI>());
        return *this;
    }

    void factorize_domain(double_d low, double_d high) {
        const size_t domain_index = m_number_of_domains++;
        m_domain_indicies.push_back(connectivity_type::value_type());
        m_domain_buffer.push_back(connectivity_type::value_type());
        m_domain_mask.push_back(connectivity_type::value_type());
        connectivity_type::reference buffer = m_domain_buffer[domain_index];
        connectivity_type::reference indicies = m_domain_indicies[domain_index];
        connectivity_type::reference mask = m_domain_mask[domain_index];
            
        const query_range a_range = a_query.get_buckets_near_point(
                                        0.5*(high-low)+low,
                                        0.5*(high-low)+m_buffer);
        
        for(reference_a a_bucket: a_range) {
            auto a_particles = a_query.get_bucket_particles(a_bucket);
            for (particle_a_reference a_particle: a_particles) {
                const size_t index = ;
                if ((get<position>(a_particle) < low).any()
                        || (get<position>(a_particle) > high).any()) {
                    buffer.push_back(index);
                } else {
                    mask.push_back(indicies.size()+buffer.size());
                    indicies.push_back(index);
                }
            }
        }

        matrix_type domain_matrix(indicies.size()+buffer.size());



        
    }

    template <typename query_type_a>
    size_t analyze_dive(reference_a bucket) {
        size_t count;
        if (a_query.is_leaf_node(bucket)) {
            auto particles = a_query.get_bucket_particles(bucket);
            count = std::distance(particles.begin(),particles.end());
            
        } else {
            count = factorize_dive(a_query.get_child1(&bucket))
                    + factorize_dive(a_query.get_child2(&bucket));
        }
        if (count < goal) {
            return count;
        } else {
            store_domain_indicies(
                    a_query.get_bounds_low(bucket),
                    a_query.get_bounds_high(bucket));
        }
    }



    // for a regular grid, assume particles are evenly distributed
    template <typename block_type>
    void factorize_impl_block(const Index i, const block_type& block) const {
        particle_type_a a;
        query_type_a a_query;
        typedef query_type_a::reference reference_a;

        const double_d& low = a_query.get_bounds_low();
        const double_d& high = a_query.get_bounds_high();
        const size_t N = a.size();
        const double total_volume = (high-low).prod();
        const double box_volume = double(m_goal_n)/double(N)*total_volume;
        const double box_side_length = std::pow(box_volume,1.0/dimension);
        const unsigned_int_d size = 
                    floor((high-low)/box_side_length)
                    .template cast<unsigned int>();
        for (int i=0; i<dimension; ++i) {
            if (size[i] == 0) {
                size[i] = 1;
            }
        }
        side_length = (high-low)/size;
        iterator_range<lattice_iterator<dimension>> range(
                lattice_iterator<dimension>(int_d(0),size,int_d(0)),
                ++lattice_iterator<dimension>(int_d(0),size,size));
        for (lattice_iterator<dimension>::reference box: range) {
            factorize_domain(box*side_length,(box+1)*side_length);
        }
    }

    template<std::size_t... I>
    void factorize_impl(detail::index_sequence<I...>) const {
        factorize_impl_block(start_row<I>,tuple_ns::get<I*NJ+I>(m_blocks));
    }
    
    template<typename MatType>
    RASMPreconditioner& factorize(const MatType& mat)
    {

        for (int i = 0; i < m_number_of_domains; ++i) {
            
        }
        factorize_impl(detail::make_index_sequence<NI>());
        // loop over leaf buckets and factorize relevent sub-matricies
        // of mat
        
      m_invdiag.resize(mat.cols());
      for(int j=0; j<mat.outerSize(); ++j)
      {
        typename MatType::InnerIterator it(mat,j);
        while(it && it.index()!=j) ++it;
        if(it && it.index()==j && it.value()!=Scalar(0))
          m_invdiag(j) = Scalar(1)/it.value();
        else
          m_invdiag(j) = Scalar(1);
      }
      m_isInitialized = true;
      return *this;
    }
    
    template<typename MatType>
    RASMPreconditioner& compute(const MatType& mat)
    {
      return factorize(mat);
    }

    /** \internal */
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
        // loop over domains and invert relevent sub-matricies in
        // mat
        Vector domain_x;
        Vector domain_b;
        for (int i = 0; i < m_number_of_domains; ++i) {
            connectivity_type::reference buffer = m_domain_buffer[i];
            connectivity_type::reference indicies = m_domain_indicies[i];
            connectivity_type::reference mask = m_domain_mask[i];
            
            const size_t nb = indicies.size()+buffer.size();
            domain_x.resize(nb);
            domain_b.resize(nb);

            // copy x values from big vector
            size_t sub_index = 0;
            for (int j = 0; j < indicies.size(); ++i) {
                domain_x[sub_index++] = x[indicies[j]];
            }
            for (int j = 0; j < buffer.size(); ++i) {
                domain_x[sub_index++] = x[buffer[j]];
            }

            // invert submatrix and mask off the buffer indicies
            domain_b = m_domain_factorized_matrix[i].solve(domain_x);

            // copy b values to big vector
            for (int j = 0; j < indicies.size(); ++i) {
                b[indicies[j]] = domain_b[mask[j]];
            }
        }
    }
    

    template<typename Rhs> inline const Solve<RASMPreconditioner, Rhs>
    solve(const MatrixBase<Rhs>& b) const
    {
      eigen_assert(m_isInitialized && "RASMPreconditioner is not initialized.");
      eigen_assert(m_invdiag.size()==b.rows()
                && "RASMPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
      return Solve<RASMPreconditioner, Rhs>(*this, b.derived());
    }
    
    ComputationInfo info() { return Success; }

  protected:
    Vector m_invdiag;
    bool m_isInitialized;
};
