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


template <template <typename,typename> class Vector=std::vector>
class RASMPreconditioner
{
    typedef double Scalar;
    typedef Matrix<Scalar,Dynamic,1> Vector;
    typedef Vector<Vector<size_t>> connectivity_type;
    connectivity_type m_domain_indicies;
    connectivity_type m_domain_buffer;
    connectivity_type m_domain_mask;
    Scalar m_buffer;
    size_t m_goal;

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

    void set_buffer_size(double size) {
        m_buffer = size;
    }

    void set_number_of_particles_per_domain(size_t n) {
        m_goal = n;
    }

    template<typename MatType>
    RASMPreconditioner& analyzePattern(const MatType& )
    {
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

    // dive tree accumulating a count of particles in each bucket
    // if the count >= goal, then factorize the bucket
    // if a child is factorized, then factorize the other children
    template <typename Query, 
             typename Reference=typename Query::reference,
             typename Pair=std::pair<bool,size_t>>
    Pair factorize_dive(Query& query, Reference bucket) {
        size_t count = 0;
        bool factorized = false;
        if (query.is_leaf_node(bucket)) {
            auto particles = query.get_bucket_particles(bucket);
            count = std::distance(particles.begin(),particles.end());
            if (count >= m_goal) {
                factorize_domain(
                    a_query.get_bounds_low(bucket),
                    a_query.get_bounds_high(bucket));
                factorized = true;
            }
        } else {
            Pair child1 = factorize_dive(a_query.get_child1(&bucket));
            Pair child2 = factorize_dive(a_query.get_child2(&bucket));
            count = child1.second + child2.second;
            if (child1.first && child2.first) {
                factorized = true;
            if (!child1.first && child2.first) {
                // if one branch factorised, factorise others
                factorize_domain(
                    a_query.get_bounds_low(query.get_child1(&bucket)),
                    a_query.get_bounds_high(query.get_child1(&bucket)));
                factorized = true;
            } else if (child1.first && !child2.first) {
                // if one branch factorised, factorise others
                factorize_domain(
                    a_query.get_bounds_low(query.get_child2(&bucket)),
                    a_query.get_bounds_high(query.get_child2(&bucket)));
                factorized = true;
            } else if (count >= m_goal) {
                // if this branch needs factorizing, do it
                factorize_domain(
                    a_query.get_bounds_low(bucket),
                    a_query.get_bounds_high(bucket));
                factorized = true;
            }
        }
        return Pair(count,factorized);
    }

    template <typename Block>
    void factorize_impl_block(const Index i, const Block& block) const {
        typedef Block::row_particles_type row_particles_type;
        typedef row_particles_type::query_type row_query_type;
        row_particles_type& a = block.get_row_particles();
        row_query_type& query = a.get_query();
        typedef query_type::reference reference;
        if (query.is_tree()) {
            // for a tree, use data structure to find a good division
            for (reference root: a.get_root_buckets()) {
                if (query.is_leaf_node(root)) {
                    // no tree!, just accept this one I guess
                    factorize_domain(
                            query.get_bounds_low(root),
                            query.get_bounds_high(root));
                } else {
                    // dive tree, find bucket with given number of particles
                    factorize_dive(query,root);
                }
            }
        } else {
            // for a regular grid, assume particles are evenly distributed
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
    }

    template <typename RowParticles, typename ColParticles, typename F>
    template <typename Block, 
             typename RowParticles=Block::row_particles_type,
             typename RowSearch=RowParticles::search_type,
template <typename> class SearchMethod
             typename = typename std::enable_if<>
    void factorize_impl_block(const Index i, const KernelDense<RowParticles,ColParticles,F>& block) {
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

    template<unsigned int NI, unsigned int NJ, typename Blocks, std::size_t... I>
    void factorize_impl(const MatrixReplacement<NI,NJ,Blocks>& mat, detail::index_sequence<I...>) {
        int dummy[] = { 0, factorize_impl_block(start_row<I>,tuple_ns::get<I*NJ+I>(mat.m_blocks))... };
        static_cast<void>(dummy);
    }
    
    template<unsigned int NI, unsigned int NJ, typename Blocks>
    RASMPreconditioner& factorize(const MatrixReplacement<NI,NJ,Blocks>& mat)
    {
        factorize_impl(mat, detail::make_index_sequence<NI>());
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
