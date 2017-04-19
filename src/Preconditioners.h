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

#ifdef HAVE_EIGEN

namespace Aboria {

template <template<typename> class Solver=Eigen::HouseholderQR>
class RASMPreconditioner {
    typedef double Scalar;
    typedef size_t Index;
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> matrix_type;
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> vector_type;
    typedef Solver<matrix_type> solver_type;
    typedef std::vector<size_t> storage_vector_type;
    typedef std::vector<storage_vector_type> connectivity_type;
    connectivity_type m_domain_indicies;
    connectivity_type m_domain_buffer;
    connectivity_type m_domain_mask;
    std::vector<solver_type> m_domain_factorized_matrix;
    Scalar m_buffer;
    size_t m_goal;
    Index m_rows;
    Index m_cols;

  public:
    typedef typename vector_type::StorageIndex StorageIndex;
    enum {
      ColsAtCompileTime = Eigen::Dynamic,
      MaxColsAtCompileTime = Eigen::Dynamic
    };

    RASMPreconditioner() : 
        m_isInitialized(false),
        m_buffer(1.0),
        m_goal(10)
    {}

    template<typename MatType>
    explicit RASMPreconditioner(const MatType& mat) {
      compute(mat);
    }

    Index rows() const { return m_rows; }
    Index cols() const { return m_cols; }

    void set_buffer_size(double size) {
        m_buffer = size;
    }

    void set_number_of_particles_per_domain(size_t n) {
        m_goal = n;
    }

    template<typename Kernel, typename Query,
        unsigned int D = Query::dimension>
    void analyze_domain(const Index start_row, 
                        const Kernel& kernel, 
                        const Query& query, 
                        const Vector<double,D>& low, 
                        const Vector<double,D>& high) {
        const size_t domain_index = m_domain_indicies.size();
        m_domain_indicies.push_back(connectivity_type::value_type());
        m_domain_buffer.push_back(connectivity_type::value_type());
        m_domain_mask.push_back(connectivity_type::value_type());
        storage_vector_type& buffer = m_domain_buffer[domain_index];
        storage_vector_type& indicies = m_domain_indicies[domain_index];
        storage_vector_type& mask = m_domain_mask[domain_index];
            
        auto range = query.get_buckets_near_point(
                              0.5*(high-low)+low,
                              0.5*(high-low)+m_buffer);
        
        typedef typename Query::reference reference;
        typedef typename Query::traits_type::position position;
        typedef typename Query::particle_iterator::reference particle_reference;
        for(reference bucket: range) {
            auto particles = query.get_bucket_particles(bucket);
            for (particle_reference particle: particles) {
                const size_t index = &(get<position>(particle))
                                        -get<position>(query.get_particles_begin());
                if ((get<position>(particle) < low).any()
                        || (get<position>(particle) > high).any()) {
                    buffer.push_back(start_row + index);
                } else {
                    mask.push_back(indicies.size()+buffer.size());
                    indicies.push_back(start_row + index);
                }
            }
        }
    }

    // dive tree accumulating a count of particles in each bucket
    // if the count >= goal, then factorize the bucket
    // if a child is factorized, then factorize the other children
    template <typename Kernel,
             typename Query, 
             typename Reference=typename Query::reference,
             typename Pair=std::pair<bool,size_t>>
    Pair analyze_dive(const Index start_row, 
                      const Kernel& kernel, 
                      const Query& query, 
                      Reference bucket) {
        size_t count = 0;
        bool done = false;
        if (query.is_leaf_node(bucket)) {
            auto particles = query.get_bucket_particles(bucket);
            count = std::distance(particles.begin(),particles.end());
            if (count >= m_goal) {
                analyze_domain(start_row, kernel, query,
                    query.get_bucket_bounds_low(bucket),
                    query.get_bucket_bounds_high(bucket));
                done = true;
            }
        } else {
            Pair child1 = analyze_dive(start_row, kernel, query, query.get_child1(&bucket));
            Pair child2 = analyze_dive(start_row, kernel, query, query.get_child2(&bucket));
            count = child1.second + child2.second;
            if (child1.first && child2.first) {
                done = true;
            } else if (!child1.first && child2.first) {
                // if one branch factorised, factorise others
                analyze_domain(start_row, kernel, query,
                    query.get_bucket_bounds_low(query.get_child1(&bucket)),
                    query.get_bucket_bounds_high(query.get_child1(&bucket)));
                done = true;
            } else if (child1.first && !child2.first) {
                // if one branch factorised, factorise others
                analyze_domain(start_row, kernel, query,
                    query.get_bucket_bounds_low(query.get_child2(&bucket)),
                    query.get_bucket_bounds_high(query.get_child2(&bucket)));
                done = true;
            } else if (count >= m_goal) {
                // if this branch needs factorizing, do it
                analyze_domain(start_row, kernel, query,
                    query.get_bucket_bounds_low(bucket),
                    query.get_bucket_bounds_high(bucket));
                done = true;
            }
        }
        return Pair(count,done);
    }

    template <typename Kernel>
    void analyze_impl_block(const Index start_row, const Kernel& kernel) {
        typedef typename Kernel::row_particles_type row_particles_type;
        typedef typename Kernel::col_particles_type col_particles_type;
        typedef typename row_particles_type::query_type query_type;
        typedef typename query_type::reference reference;
        static const unsigned int dimension = query_type::dimension;
        typedef Vector<double,dimension> double_d;
        typedef Vector<unsigned int,dimension> unsigned_int_d;
        typedef Vector<int,dimension> int_d;

        static_assert(std::is_same<row_particles_type,col_particles_type>::value,
           "RASM preconditioner restricted to identical row and col particle sets");
        const row_particles_type& a = kernel.get_row_particles();
        CHECK(&a == &(kernel.get_col_particles()),
           "RASM preconditioner restricted to identical row and col particle sets");
        const query_type& query = a.get_query();
        if (query.is_tree()) {
            // for a tree, use data structure to find a good division
            for (reference root: query.get_root_buckets()) {
                if (query.is_leaf_node(root)) {
                    // no tree!, just accept this one I guess
                    analyze_domain(start_row, kernel, query,
                            query.get_bucket_bounds_low(root),
                            query.get_bucket_bounds_high(root));
                } else {
                    // dive tree, find bucket with given number of particles
                    analyze_dive(start_row,kernel,query,root);
                }
            }
        } else {
            // for a regular grid, assume particles are evenly distributed
            const double_d& low = query.get_bounds_low();
            const double_d& high = query.get_bounds_high();
            const size_t N = a.size();
            const double total_volume = (high-low).prod();
            const double box_volume = double(m_goal)/double(N)*total_volume;
            const double box_side_length = std::pow(box_volume,1.0/dimension);
            unsigned_int_d size = 
                floor((high-low)/box_side_length)
                .template cast<unsigned int>();
            for (int i=0; i<dimension; ++i) {
                if (size[i] == 0) {
                    size[i] = 1;
                }
            }
            const double_d side_length = (high-low)/size;
            iterator_range<lattice_iterator<dimension>> range(
                    lattice_iterator<dimension>(int_d(0),size,int_d(0)),
                    ++lattice_iterator<dimension>(int_d(0),size,size));
            for (typename lattice_iterator<dimension>::reference box: range) {
                analyze_domain(start_row, kernel, query,
                        box*side_length,(box+1)*side_length);
            }
        }
    }

    template <typename RowParticles, typename ColParticles>
    void analyze_impl_block(
            const Index start_row, 
            const KernelZero<RowParticles,ColParticles>& kernel) {
    }
     
    template<unsigned int NI, unsigned int NJ, typename Blocks, std::size_t... I>
    void analyze_impl(const MatrixReplacement<NI,NJ,Blocks>& mat, 
                        detail::index_sequence<I...>) {
        int dummy[] = { 0, 
          (analyze_impl_block(mat.template start_row<I>(),tuple_ns::get<I*NJ+I>(mat.m_blocks)),0)... 
            };
        static_cast<void>(dummy);
    }

    template<unsigned int NI, unsigned int NJ, typename Blocks>
    RASMPreconditioner& analyzePattern(const MatrixReplacement<NI,NJ,Blocks>& mat)
    {
        m_rows = mat.rows();
        m_cols = mat.cols();
        analyze_impl(mat, detail::make_index_sequence<NI>());
        return *this;
    }

    template <int _Options, typename _StorageIndex>
    RASMPreconditioner& analyzePattern(const Eigen::SparseMatrix<Scalar,_Options,_StorageIndex>& mat) {
        CHECK(m_domain_indicies.size()>0, "RASMPreconditioner::analyzePattern(): cannot analyze sparse matrix, call analyzePattern using a Aboria MatrixReplacement class first");
    }
    template<typename Derived>
    RASMPreconditioner& analyzePattern(const Eigen::DenseBase<Derived>& mat) {
        CHECK(m_domain_indicies.size()>0, "RASMPreconditioner::analyzePattern(): cannot analyze dense matrix, call analyzePattern need to pass a Aboria MatrixReplacement class first");
    }

    template<typename MatType>
    RASMPreconditioner& factorize(const MatType& mat)
    {
        eigen_assert(m_rows==mat.rows()
                && "RASMPreconditioner::solve(): invalid number of rows of mat");
        eigen_assert(m_cols==mat.cols()
                && "RASMPreconditioner::solve(): invalid number of rows of mat");

        m_domain_factorized_matrix.resize(m_domain_indicies.size());

        matrix_type domain_matrix;

        for (int domain_index = 0; domain_index < m_domain_factorized_matrix.size(); ++domain_index) {
            const storage_vector_type& buffer = m_domain_buffer[domain_index];
            const storage_vector_type& indicies = m_domain_indicies[domain_index];
            const storage_vector_type& mask = m_domain_mask[domain_index];
            solver_type& solver = m_domain_factorized_matrix[domain_index];

            domain_matrix.resize(indicies.size()+buffer.size(),
                                 indicies.size()+buffer.size());

            size_t i = 0;
            for (const size_t& big_index_i: indicies) {
                size_t j = 0;
                for (const size_t& big_index_j: indicies) {
                    domain_matrix(i,j++) = mat.coeff(big_index_i,big_index_j);
                }
                for (const size_t& big_index_j: buffer) {
                    domain_matrix(i,j++) = mat.coeff(big_index_i,big_index_j);
                }
                ++i;
            }
            for (const size_t& big_index_i: buffer) {
                size_t j = 0;
                for (const size_t& big_index_j: indicies) {
                    domain_matrix(i,j++) = mat.coeff(big_index_i,big_index_j);
                }
                for (const size_t& big_index_j: buffer) {
                    domain_matrix(i,j++) = mat.coeff(big_index_i,big_index_j);
                }
                ++i;
            }

            solver.compute(domain_matrix);
        }

        m_isInitialized = true;

        return *this;
    }
    
    template<typename MatType>
    RASMPreconditioner& compute(const MatType& mat)
    {
        analyzePattern(mat);
        return factorize(mat);
    }

    /** \internal */
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const
    {
        // loop over domains and invert relevent sub-matricies in
        // mat
        vector_type domain_x;
        vector_type domain_b;
        for (int i = 0; i < m_domain_indicies.size(); ++i) {
            const storage_vector_type& buffer = m_domain_buffer[i];
            const storage_vector_type& indicies = m_domain_indicies[i];
            const storage_vector_type& mask = m_domain_mask[i];
            
            const size_t nb = indicies.size()+buffer.size();
            domain_x.resize(nb);
            domain_b.resize(nb);

            // copy x values from big vector
            size_t sub_index = 0;
            for (int j = 0; j < indicies.size(); ++j) {
                domain_b[sub_index++] = b[indicies[j]];
            }
            for (int j = 0; j < buffer.size(); ++j) {
                domain_b[sub_index++] = b[buffer[j]];
            }

            // invert submatrix and mask off the buffer indicies
            domain_x = m_domain_factorized_matrix[i].solve(domain_b);

            // copy b values to big vector
            for (int j = 0; j < indicies.size(); ++j) {
                x[indicies[j]] = domain_x[mask[j]];
            }
        }
    }
    

    template<typename Rhs> 
    inline const Eigen::Solve<RASMPreconditioner, Rhs>
    solve(const Eigen::MatrixBase<Rhs>& b) const {
        eigen_assert(m_rows==b.rows()
                && "RASMPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
        eigen_assert(m_isInitialized 
                && "RASMPreconditioner is not initialized.");
        return Eigen::Solve<RASMPreconditioner, Rhs>(*this, b.derived());
    }
    
    Eigen::ComputationInfo info() { return Eigen::Success; }

  protected:
    bool m_isInitialized;
};

}

#endif //HAVE_EIGEN
#endif 
