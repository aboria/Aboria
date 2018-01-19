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


#include <fstream>

#ifdef HAVE_EIGEN
#include <unsupported/Eigen/SparseExtra>

namespace Aboria {

namespace detail {
    template<typename Function, 
             typename Dest, unsigned int NI, unsigned int NJ, typename Blocks, typename Rhs>
    void apply_function_to_diagonal_blocks(Function&& f, 
                             Dest& y, 
                             const MatrixReplacement<NI,NJ,Blocks>& mat, 
                             const Rhs& rhs,
                             std::integral_constant<unsigned int,NI>) {}

    template<typename Function,
             unsigned int NI, unsigned int NJ, typename Blocks>
    void apply_function_to_diagonal_blocks(Function&& f, 
                             const MatrixReplacement<NI,NJ,Blocks>& mat,
                             std::integral_constant<unsigned int,NI>) {}


    template<typename Function, 
             typename Dest, 
             unsigned int NI, 
             unsigned int NJ, 
             typename Blocks, 
             typename Rhs, 
             unsigned int I>
    void apply_function_to_diagonal_blocks(Function&& f,
                                        Dest& x, 
                                        const MatrixReplacement<NI,NJ,Blocks>& mat, 
                                        const Rhs& b, std::integral_constant<unsigned int,I>) {
        f(x.segment(mat.template start_row<I>(),mat.template size_row<I>()),
          b.segment(mat.template start_col<I>(),mat.template size_col<I>()),
          std::get<I*NJ+I>(mat.m_blocks));
        apply_function_to_diagonal_blocks(std::forward<Function>(f),x,mat,b,std::integral_constant<unsigned int,I+1>());
    }

    template<typename Function, 
             unsigned int NI, 
             unsigned int NJ, 
             typename Blocks, 
             unsigned int I>
    void apply_function_to_diagonal_blocks(Function&& f,
                                        const MatrixReplacement<NI,NJ,Blocks>& mat, 
                                        std::integral_constant<unsigned int,I>) {
        f(std::get<I*NJ+I>(mat.m_blocks));
        apply_function_to_diagonal_blocks(std::forward<Function>(f),mat,std::integral_constant<unsigned int,I+1>());
    }

    template<typename Function, typename Dest,
             unsigned int NI, unsigned int NJ, typename Blocks, typename Rhs>
    void apply_function_to_diagonal_blocks(Function&& function,
                                        Dest& x, 
                                        const MatrixReplacement<NI,NJ,Blocks>& mat, 
                                        const Rhs& b) {
        apply_function_to_diagonal_blocks(std::forward<Function>(function),x,mat,b,std::integral_constant<unsigned int,0>());
    }

    template<typename Function, unsigned int NI, unsigned int NJ, typename Blocks>
    void apply_function_to_diagonal_blocks(Function&& function,
                                        const MatrixReplacement<NI,NJ,Blocks>& mat) {
        apply_function_to_diagonal_blocks(std::forward<Function>(function),mat,std::integral_constant<unsigned int,0>());
    }
}

template <typename Solver>
class ReducedOrderPreconditioner {
    typedef Solver solver_type;
    typedef size_t Index;
    typedef double Scalar;
    typedef H2LibMatrix h2_matrix_type;
    Index m_rows;
    Index m_cols;
    size_t m_order;
    double m_tol;
    std::vector<size_t> m_col_sizes;
    std::vector<size_t> m_row_sizes;
    std::vector<std::shared_ptr<solver_type>> m_solvers;
    std::vector<std::shared_ptr<h2_matrix_type>> m_h2mats;

  public:
    typedef size_t StorageIndex;
    enum {
      ColsAtCompileTime = Eigen::Dynamic,
      MaxColsAtCompileTime = Eigen::Dynamic
    };

    ReducedOrderPreconditioner():m_order(4),m_tol(1e-5)
    {}

    template<typename MatType>
    explicit ReducedOrderPreconditioner(const MatType& mat):m_order(4),m_tol(1e-5) {
      compute(mat);
    }

    Index rows() const { return m_rows; }
    Index cols() const { return m_cols; }

    void set_order(const size_t order) {
        m_order = order;
    }

    void set_tolerance(const double tol) {
        m_tol = tol;
    }

    template<unsigned int NI, unsigned int NJ, typename Blocks>
    ReducedOrderPreconditioner& analyzePattern(const MatrixReplacement<NI,NJ,Blocks>& mat)
    {

        LOG(2,"ExtMatrixPreconditioner: analyze pattern");
        m_rows = mat.rows();
        m_cols = mat.cols();

        return *this;
    }
    

    struct factorize_block {
        std::vector<size_t> &m_col_sizes; 
        std::vector<size_t> &m_row_sizes; 
        std::vector<std::shared_ptr<solver_type>> &m_solvers; 
        std::vector<std::shared_ptr<h2_matrix_type>> &m_h2mats; 
        size_t m_order;
        double m_tol;
        int i;

        factorize_block(std::vector<size_t>& col_sizes,
                        std::vector<size_t>& row_sizes,
                        std::vector<std::shared_ptr<solver_type>>& solvers,
                        std::vector<std::shared_ptr<h2_matrix_type>>& h2mats,
                        size_t order,
                        double tol):
            m_col_sizes(col_sizes),m_row_sizes(row_sizes),m_solvers(solvers),m_h2mats(h2mats),m_order(order),m_tol(tol),i(0) {}


        template <typename Block>
        void operator()(const Block& block) {
            LOG(2,"ReducedOrderPreconditioner: block "<<i<<": non h2 block");
            m_solvers[i] = nullptr;
            m_col_sizes[i] = block.cols();
            m_row_sizes[i] = block.rows();
            ++i;
        }

        template <typename RowParticles, typename ColParticles, typename PositionF, typename F>
        void operator()(const KernelH2<RowParticles,ColParticles,PositionF,F>& kernel) {
            static const unsigned int dimension = RowParticles::dimension;

            LOG(2,"ReducedOrderPreconditioner: block "<<i<<": generate h2 matrix with order "<<m_order);
            m_h2mats[i] = std::make_shared<h2_matrix_type>(
                                    make_h2lib_matrix(kernel.get_row_elements(),
                                        kernel.get_col_elements(),
                       make_h2lib_black_box_expansion<dimension>(
                                            m_order,
                                            kernel.get_position_function()),
                                        kernel.get_kernel_function()));
            m_col_sizes[i] = kernel.cols();
            m_row_sizes[i] = kernel.rows();

            LOG(2,"ReducedOrderPreconditioner: block "<<i<<": factorise h2 matrix with tolerance "<<m_tol);
            m_solvers[i] = std::make_shared<solver_type>(m_h2mats[i]->get_ph2matrix(),m_h2mats[i]->get_pblock(),m_tol);

            //m_solvers[i]->setMaxIterations(m_inner_iterations);
            //LOG(2,"ExtMatrixPreconditioner: block "<<i<<": set precon");
            //m_solvers[i]->preconditioner().setDroptol(0.1);
            //m_solvers[i]->preconditioner().compute(m_str_ext_matrices[i]);
            LOG(2,"ReducedOrderPreconditioner: block "<<i<<": factorization complete");
            ++i;
        }
    };


    template <unsigned int NI, unsigned int NJ, typename Blocks>
    ReducedOrderPreconditioner& factorize(const MatrixReplacement<NI,NJ,Blocks>& mat)
    {
        LOG(2,"ReducedOrderPreconditioner: factorizing domain");

        m_rows = mat.rows();
        m_cols = mat.cols();
        m_solvers.resize(NI);
        m_h2mats.resize(NI);
        m_col_sizes.resize(NI);
        m_row_sizes.resize(NI);
        detail::apply_function_to_diagonal_blocks(
                factorize_block(m_col_sizes,m_row_sizes,m_solvers,m_h2mats,m_order,m_tol), 
                mat);
        m_isInitialized = true;

        return *this;
    }
    
    template<typename MatType>
    ReducedOrderPreconditioner& compute(const MatType& mat)
    {
        analyzePattern(mat);
        return factorize(mat);
    }

    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const {
        size_t row = 0;
        size_t col = 0;
        std::vector<double> buffer;

        for (int i = 0; i < m_solvers.size(); ++i) {
            auto b_segment = b.segment(col,m_col_sizes[i]);
            auto x_segment = x.segment(row,m_row_sizes[i]);
            if (m_solvers[i] != nullptr) { // solver only exists for h2 blocks
                LOG(2,"ReducedOrderPreconditioner: block "<<i<<" solve");
                typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> map_type; 
                buffer.resize(std::max(m_col_sizes[i],m_row_sizes[i]));
                map_type b_map(buffer.data(),m_row_sizes[i]);
                map_type x_map(buffer.data(),m_col_sizes[i]);
                b_map = b_segment;
                m_solvers[i]->solve(buffer);
                x_segment = x_map;
                m_h2mats[i]->matrix_vector_multiply(b_map,1,false,x_map);
                std::cout << "accuracy = "<<(b_map-b_segment).norm()<< std::endl;
            } else {
                x_segment = b_segment;
            }
            row += m_row_sizes[i];
            col += m_col_sizes[i];
            
        }
        LOG(2,"ReducedOrderPreconditioner: done solve_impl");
    }
    

    template<typename Rhs> 
    inline const Eigen::Solve<ReducedOrderPreconditioner, Rhs>
    solve(const Eigen::MatrixBase<Rhs>& b) const {
        eigen_assert(m_rows==b.rows()
                && "ReducedOrderPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
        eigen_assert(m_isInitialized 
                && "ReducedOrderPreconditioneris not initialized.");
        return Eigen::Solve<ReducedOrderPreconditioner, Rhs>(*this, b.derived());
    }
    
    Eigen::ComputationInfo info() { return Eigen::Success; }

  protected:
    bool m_isInitialized;
};


#if 0
template <unsigned int ReducedOrder, 
          typename InnerPreconditioner=Eigen::IncompleteLUT<double>,
          typename IterativeSolver=Eigen::DGMRES<Eigen::SparseMatrix<double>,
                                                InnerPreconditioner>>
class ExtMatrixPreconditioner {
    typedef double Scalar;
    typedef size_t Index;
    typedef Eigen::SparseMatrix<Scalar> sparse_matrix_type;
    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> vector_type;
    typedef Eigen::Matrix<int,Eigen::Dynamic,1> index_vector_type;
    typedef InnerPreconditioner solver_type;
    Index m_rows;
    Index m_cols;
    size_t m_inner_iterations;
    std::vector<index_vector_type> m_col_maps;
    std::vector<size_t> m_col_sizes;
    std::vector<size_t> m_row_sizes;
    std::vector<index_vector_type> m_row_maps;
    std::vector<std::shared_ptr<solver_type>> m_solvers;
    std::vector<sparse_matrix_type> m_ext_matrices;
    std::vector<sparse_matrix_type> m_str_ext_matrices;

  public:
    typedef typename vector_type::StorageIndex StorageIndex;
    enum {
      ColsAtCompileTime = Eigen::Dynamic,
      MaxColsAtCompileTime = Eigen::Dynamic
    };

    ExtMatrixPreconditioner():m_inner_iterations(10)
    {}

    template<typename MatType>
    explicit ExtMatrixPreconditioner(const MatType& mat):m_inner_iterations(10) {
      compute(mat);
    }

    Index rows() const { return m_rows; }
    Index cols() const { return m_cols; }

    template<unsigned int NI, unsigned int NJ, typename Blocks>
    ExtMatrixPreconditioner& analyzePattern(const MatrixReplacement<NI,NJ,Blocks>& mat)
    {

        LOG(2,"ExtMatrixPreconditioner: analyze pattern");
        m_rows = mat.rows();
        m_cols = mat.cols();

        return *this;
    }
    

    struct factorize_block {
        std::vector<size_t> &m_col_sizes; 
        std::vector<size_t> &m_row_sizes; 
        std::vector<index_vector_type> &m_col_maps; 
        std::vector<index_vector_type> &m_row_maps; 
        std::vector<std::shared_ptr<solver_type>> &m_solvers; 
        std::vector<sparse_matrix_type> &m_ext_matrices; 
        std::vector<sparse_matrix_type> &m_str_ext_matrices; 
        size_t m_inner_iterations;
        int i;

        factorize_block(std::vector<size_t>& col_sizes,
                        std::vector<size_t>& row_sizes,
                        std::vector<index_vector_type>& col_maps,
                        std::vector<index_vector_type>& row_maps,
                        std::vector<std::shared_ptr<solver_type>>& solvers,
                        std::vector<sparse_matrix_type>& ext_matrices,
                        std::vector<sparse_matrix_type>& str_ext_matrices,
                        size_t inner_iterations):
            m_col_sizes(col_sizes),m_row_sizes(row_sizes),m_col_maps(col_maps),m_row_maps(row_maps),m_solvers(solvers),m_ext_matrices(ext_matrices),m_str_ext_matrices(str_ext_matrices),m_inner_iterations(inner_iterations),i(0) {}


        template <typename Block>
        void operator()(const Block& block) {
            LOG(2,"ExtMatrixPreconditioner: block "<<i<<": non h2 block");
            m_solvers[i] = nullptr;
            m_col_sizes[i] = block.cols();
            m_row_sizes[i] = block.rows();
            ++i;
        }

        template <typename RowParticles, typename ColParticles, typename PositionF>
        void operator()(const KernelH2<RowParticles,ColParticles,PositionF>& kernel) {
            static const unsigned int dimension = RowParticles::dimension;

            auto h2 = make_h2_matrix(kernel.get_row_particles(),
                                     kernel.get_col_particles(),
                       make_black_box_expansion<dimension,ReducedOrder>(
                                    kernel.get_position_function()));
            LOG(2,"ExtMatrixPreconditioner: block "<<i<<": generate extended matrix");
            m_ext_matrices[i] = h2.gen_extended_matrix();
            m_str_ext_matrices[i] = h2.gen_stripped_extended_matrix();
            /*
            Eigen::saveMarket(m_ext_matrices[i],"ext_matrix.mat");
            Eigen::saveMarket(m_str_ext_matrices[i],"str_ext_matrix.mat");
            std::ofstream myfile;
            myfile.open ("ext_matrix.csv");
            myfile << Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(m_ext_matrices[i]);
            myfile.close();
            myfile.open ("str_ext_matrix.csv");
            myfile << Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>(m_str_ext_matrices[i]);
            myfile.close();
            */

            m_col_sizes[i] = kernel.cols();
            m_row_sizes[i] = kernel.rows();
            m_col_maps[i] = h2.gen_column_map();
            m_row_maps[i] = h2.gen_row_map();
            LOG(2,"ExtMatrixPreconditioner: block "<<i<<": set precon");
            //LOG(2,"ExtMatrixPreconditioner: block "<<i<<": create solver");
            //m_solvers[i] = std::make_shared<solver_type>(m_ext_matrices[i]);
            m_solvers[i] = std::make_shared<solver_type>();
            m_solvers[i]->setDroptol(0.1);
            m_solvers[i]->compute(m_str_ext_matrices[i]);
            if (m_solvers[i]->info() != Eigen::Success) {
                ERROR("ExtMatrixPreconditioner inner preconditioner could not factorize");
            }
            //m_solvers[i]->setMaxIterations(m_inner_iterations);
            //LOG(2,"ExtMatrixPreconditioner: block "<<i<<": set precon");
            //m_solvers[i]->preconditioner().setDroptol(0.1);
            //m_solvers[i]->preconditioner().compute(m_str_ext_matrices[i]);
            LOG(2,"ExtMatrixPreconditioner: block "<<i<<": factorization complete");
            ++i;
        }
    };


    template <unsigned int NI, unsigned int NJ, typename Blocks>
    ExtMatrixPreconditioner& factorize(const MatrixReplacement<NI,NJ,Blocks>& mat)
    {
        LOG(2,"ExtMatrixPreconditioner: factorizing domain");

        m_rows = mat.rows();
        m_cols = mat.cols();
        m_solvers.resize(NI);
        m_ext_matrices.resize(NI);
        m_str_ext_matrices.resize(NI);
        m_col_sizes.resize(NI);
        m_row_sizes.resize(NI);
        m_col_maps.resize(NI);
        m_row_maps.resize(NI);
        detail::apply_function_to_diagonal_blocks(
                factorize_block(m_col_sizes,m_row_sizes,m_col_maps,m_row_maps,
                                m_solvers,m_ext_matrices,m_str_ext_matrices,m_inner_iterations), 
                mat);
        m_isInitialized = true;

        return *this;
    }
    
    template<typename MatType>
    ExtMatrixPreconditioner& compute(const MatType& mat)
    {
        analyzePattern(mat);
        return factorize(mat);
    }

    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const {
        vector_type m_ext_b;
        vector_type m_ext_x;
        size_t row = 0;
        size_t col = 0;

        for (int i = 0; i < m_solvers.size(); ++i) {
            if (m_solvers[i] != nullptr) { // solver only exists for h2 blocks
                LOG(2,"ExtMatrixPreconditioner: block "<<i<<" solve");
                // construct b
                m_ext_b.resize(m_ext_matrices[i].rows());
                for (int j = 0; j < m_row_maps[i].size(); ++j) {
                    m_ext_b[m_row_maps[i][j]] = b[row + j];
                }
                for (int j = m_row_maps[i].size(); j < m_ext_matrices[i].rows(); ++j) {
                    m_ext_b[j] = 0;
                }
                m_ext_x = m_ext_b;
                //for (int j = 0; j < m_ext_matrices[i].cols(); ++j) {
                //    m_ext_x[j] = 0;
                //}
                //for (int j = 0; j < m_ext_matrices[i].rows(); ++j) {
                //    std::cout << "ext_b["<<j<<"] = "<< m_ext_b[j] << std::endl;
                //}

                //TODO: solve with guess?
                //m_ext_x = m_solvers[i]->solveWithGuess(m_ext_b,m_ext_x);
                //m_ext_x = m_solvers[i]->solve(m_ext_b);
                Eigen::Index iters = m_inner_iterations; 
                double tol_error = 1e-10;
                std::cout << "m_ext_b norm = "<<m_ext_b.norm() << std::endl;
                std::cout << "m_ext_x norm = "<<m_ext_x.norm() << std::endl;
                std::cout << "residual norm = "<<(m_ext_matrices[i]*m_ext_x-m_ext_b).norm() << std::endl;
                Eigen::internal::gmres(m_ext_matrices[i],m_ext_b,m_ext_x,*m_solvers[i],
                        iters,2*m_inner_iterations,tol_error);
                //LOG(2,"ExtMatrixPreconditioner: solve complete: #iterations: " << m_solvers[i]->iterations() << ", estimated error: " << m_solvers[i]->error() << " true error = "<<(m_ext_matrices[i]*m_ext_x-m_ext_b).norm());
                LOG(2,"ExtMatrixPreconditioner: solve complete: #iterations: " << iters << " true error = "<<(m_ext_matrices[i]*m_ext_x-m_ext_b).norm());

                // filter to x
                for (int j = 0; j < m_col_maps[i].size(); ++j) {
                    x[col + j] = m_ext_x[m_col_maps[i][j]];
                }

                // increment row/col by number of particles (NOT size of ext vectors)
                row += m_row_maps[i].size();
                col += m_col_maps[i].size();
            } else {
                LOG(2,"ExtMatrixPreconditioner: block "<<i<<" non h2 block");
                for (int j = 0; j < m_col_sizes[i]; ++j) {
                    x[col + j] = b[row + j];
                }
                // increment row/col by the size of the block
                row += m_row_sizes[i];
                col += m_col_sizes[i];
            }
            
        }
        LOG(2,"ExtMatrixPreconditioner: done solve_impl");
    }
    

    template<typename Rhs> 
    inline const Eigen::Solve<ExtMatrixPreconditioner, Rhs>
    solve(const Eigen::MatrixBase<Rhs>& b) const {
        eigen_assert(m_rows==b.rows()
                && "ExtMatrixPreconditioner::solve(): invalid number of rows of the right hand side matrix b");
        eigen_assert(m_isInitialized 
                && "ExtMatrixPreconditioner is not initialized.");
        return Eigen::Solve<ExtMatrixPreconditioner, Rhs>(*this, b.derived());
    }
    
    Eigen::ComputationInfo info() { return Eigen::Success; }

  protected:
    bool m_isInitialized;
};
#endif


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
    connectivity_type m_domain_random;
    std::vector<solver_type> m_domain_factorized_matrix;
    Scalar m_buffer;
    size_t m_goal;
    size_t m_random;
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
        m_buffer(0),
        m_goal(10),
        m_random(0)
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

    void set_number_of_random_particles(size_t n) {
        ERROR("RASM Preconditioner random support not working")
        m_random = n;
    }

    template<typename Kernel, typename Query,
        unsigned int D = Query::dimension,
        typename box_type = detail::bbox<D>
        >
    void analyze_domain(const Index start_row, 
                        const Kernel& kernel, 
                        const Query& query, 
                        const box_type& bounds) {
        const size_t domain_index = m_domain_indicies.size();
        m_domain_indicies.push_back(connectivity_type::value_type());
        m_domain_buffer.push_back(connectivity_type::value_type());
        storage_vector_type& buffer = m_domain_buffer[domain_index];
        storage_vector_type& indicies = m_domain_indicies[domain_index];
            
        //TODO: cant use chebyshev search function as does not yet
        //support anisotropic distance measures, need to add this then 
        //remove hack below
        auto range = query.template get_buckets_near_point<-1>(
                              0.5*(bounds.bmax-bounds.bmin)+bounds.bmin,
                              0.5*(bounds.bmax-bounds.bmin)+m_buffer);
        
        typedef typename Query::reference reference;
        typedef typename Query::traits_type::position position;
        typedef typename Query::particle_iterator::reference particle_reference;
        for(reference bucket: range) {
            auto particles = query.get_bucket_particles(bucket);
            for (particle_reference particle: particles) {
                const size_t index = &(get<position>(particle))
                                        -get<position>(query.get_particles_begin());
                if ((get<position>(particle) >= bounds.bmin-m_buffer).all()
                        && (get<position>(particle) < bounds.bmax+m_buffer).all()) {
                    if ((get<position>(particle) < bounds.bmin).any()
                            || (get<position>(particle) >= bounds.bmax).any()) {
                        buffer.push_back(start_row + index);
                    } else {
                        indicies.push_back(start_row + index);
                    }
                }
            }
        }
        ASSERT(buffer.size() > 0, "no particles in buffer");
        ASSERT(indicies.size() > 0, "no particles in domain");
        //std::cout << "analyze_domain with bounds = "<<bounds<<" indicies = "<<indicies.size()<<" buffer = "<<buffer.size() << " volume = "<<(bounds.bmax-bounds.bmin).prod()<<std::endl;
    }


    // dive tree accumulating a count of particles in each bucket
    // if the count >= goal, then factorize the bucket
    // if a child is factorized, then factorize the other children
    template <typename Kernel,
             typename Query, 
             typename Reference=typename Query::reference,
             typename Pair=std::pair<bool,size_t>,
             unsigned int D = Query::dimension,
             typename child_iterator = typename Query::child_iterator,
             typename box_type = detail::bbox<D>
             >
    Pair analyze_dive(const Index start_row, 
                      const Kernel& kernel, 
                      const Query& query, 
                      const child_iterator& ci) {
        size_t count = 0;
        bool done = false;
        if (query.is_leaf_node(*ci)) {
            auto particles = query.get_bucket_particles(*ci);
            count = std::distance(particles.begin(),particles.end());
        } else {
            std::vector<child_iterator> not_done;
            for (child_iterator cj = query.get_children(ci); cj != false; ++cj) {
                Pair child = analyze_dive(start_row, kernel, query, cj);
                done |= child.first;
                if (!child.first) {
                    not_done.push_back(cj);
                }
                count += child.second;
                
            }
            // if any child is done, need to analyze domain for siblings
            if (done) {
                for (child_iterator& cj: not_done) {
                    analyze_domain(start_row, kernel, query,
                            query.get_bounds(cj));
                }
            }
        }
        if (count >= m_goal && !done) {
            analyze_domain(start_row, kernel, query,
                    query.get_bounds(ci));
            done = true;
        }
        //std::cout << "analyze_dive bounds "<< query.get_bounds(ci) << std::endl;
        //std::cout << "leaf = "<<query.is_leaf_node(*ci)<<" count = "<<count<<" done = "<<done << std::endl;
        return Pair(done,count);
    }

    template <typename Kernel>
    void analyze_impl_block(const Index start_row, const Kernel& kernel) {
        typedef typename Kernel::row_elements_type row_elements_type;
        typedef typename Kernel::col_elements_type col_elements_type;
        typedef typename row_elements_type::query_type query_type;
        typedef typename query_type::reference reference;
        static const unsigned int dimension = query_type::dimension;
        typedef Vector<double,dimension> double_d;
        typedef Vector<unsigned int,dimension> unsigned_int_d;
        typedef Vector<int,dimension> int_d;
        typedef detail::bbox<dimension> box_type;
        typedef typename query_type::child_iterator child_iterator;

        static_assert(std::is_same<row_elements_type,col_elements_type>::value,
           "RASM preconditioner restricted to identical row and col particle sets");
        const row_elements_type& a = kernel.get_row_elements();
        CHECK(&a == &(kernel.get_col_elements()),
           "RASM preconditioner restricted to identical row and col particle sets");
        const query_type& query = a.get_query();
        if (query.is_tree()) {
            // for a tree, use data structure to find a good division
            for (child_iterator ci = query.get_children(); ci != false; ++ci) {
                if (query.is_leaf_node(*ci)) {
                    // no tree!, just accept this one I guess
                    analyze_domain(start_row, kernel, query,
                            query.get_bounds(ci));
                } else {
                    // dive tree, find bucket with given number of particles
                    analyze_dive(start_row,kernel,query,ci);
                }
            }
        } else {
            // for a regular grid, assume particles are evenly distributed
            const box_type& bounds = query.get_bounds();
            const size_t N = a.size();
            const double total_volume = (bounds.bmax-bounds.bmin).prod();
            const double box_volume = double(m_goal)/double(N)*total_volume;
            const double box_side_length = std::pow(box_volume,1.0/dimension);
            unsigned_int_d size = 
                floor((bounds.bmax-bounds.bmin)/box_side_length)
                .template cast<unsigned int>();
            for (int i=0; i<dimension; ++i) {
                if (size[i] == 0) {
                    size[i] = 1;
                }
            }
            const double_d side_length = (bounds.bmax-bounds.bmin)/size;
            iterator_range<lattice_iterator<dimension>> range(
                    lattice_iterator<dimension>(int_d::Zero(),size),
                    lattice_iterator<dimension>());
            for (typename lattice_iterator<dimension>::reference box: range) {
                analyze_domain(start_row, kernel, query,
                        box_type(box*side_length,(box+1)*side_length));
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
          (analyze_impl_block(mat.template start_row<I>(),std::get<I*NJ+I>(mat.m_blocks)),0)... 
            };
        static_cast<void>(dummy);
    }

    template<unsigned int NI, unsigned int NJ, typename Blocks>
    RASMPreconditioner& analyzePattern(const MatrixReplacement<NI,NJ,Blocks>& mat)
    {

        LOG(2,"RASMPreconditioner: analyze pattern");
        m_rows = mat.rows();
        m_cols = mat.cols();
        analyze_impl(mat, detail::make_index_sequence<NI>());

        

        const size_t n = m_domain_indicies.size();
        generator_type generator(time(NULL));
        std::uniform_int_distribution<int> uniform_domain(0,n-1);

        m_domain_random.resize(n);
        for (int domain_index = 0; domain_index < n; ++domain_index) {
            storage_vector_type& random = m_domain_random[domain_index];
            random.resize(m_random);
            for (int d = 0; d < m_random; ++d) {
                do {
                    int other_domain = uniform_domain(generator);
                    while (other_domain == domain_index) other_domain = uniform_domain(generator);
                    std::uniform_int_distribution<int> uniform_index(0,m_domain_indicies[other_domain].size()-1);
                    random[d] = m_domain_indicies[other_domain][uniform_index(generator)];
                } while ((random.begin()+d)!=std::find(random.begin(),random.begin()+d,random[d]));
            }
        }

        int count = 0;
        int minsize = 1000;
        int maxsize = 0;
        for (int domain_index = 0; domain_index < m_domain_indicies.size(); ++domain_index) {
            const int size = m_domain_indicies[domain_index].size();
            count += size;
            if (size < minsize) minsize = size;
            if (size > maxsize) maxsize = size;
        }
        LOG(2,"RASMPreconditioner: finished analysis, found "<<m_domain_indicies.size()<<" domains, with "<<minsize<<"--"<<maxsize<<" particles ("<<count<<" total).");
        return *this;
    }

    template <int _Options, typename _StorageIndex>
    RASMPreconditioner& analyzePattern(const Eigen::SparseMatrix<Scalar,_Options,_StorageIndex>& mat) {
        CHECK(m_domain_indicies.size()>0, "RASMPreconditioner::analyzePattern(): cannot analyze sparse matrix, call analyzePattern using a Aboria MatrixReplacement class first");
        return *this;
    }

    template <int _Options, typename _StorageIndex,  int RefOptions, typename RefStrideType>
    RASMPreconditioner& analyzePattern(const Eigen::Ref<const Eigen::SparseMatrix<Scalar,_Options,_StorageIndex>,RefOptions,RefStrideType>& mat) {
        CHECK(m_domain_indicies.size()>0, "RASMPreconditioner::analyzePattern(): cannot analyze sparse matrix, call analyzePattern using a Aboria MatrixReplacement class first");
        return *this;
    }

    template<typename Derived>
    RASMPreconditioner& analyzePattern(const Eigen::DenseBase<Derived>& mat) {
        CHECK(m_domain_indicies.size()>0, "RASMPreconditioner::analyzePattern(): cannot analyze dense matrix, call analyzePattern need to pass a Aboria MatrixReplacement class first");
        return *this;
    }

    template<typename MatType>
    RASMPreconditioner& factorize(const MatType& mat)
    {
        LOG(2,"RASMPreconditioner: factorizing domain");
        eigen_assert(m_rows==mat.rows()
                && "RASMPreconditioner::solve(): invalid number of rows of mat");
        eigen_assert(m_cols==mat.cols()
                && "RASMPreconditioner::solve(): invalid number of rows of mat");

        m_domain_factorized_matrix.resize(m_domain_indicies.size());

        matrix_type domain_matrix;

        for (int domain_index = 0; domain_index < m_domain_factorized_matrix.size(); ++domain_index) {
            const storage_vector_type& buffer = m_domain_buffer[domain_index];
            const storage_vector_type& indicies = m_domain_indicies[domain_index];
            const storage_vector_type& random = m_domain_random[domain_index];
            solver_type& solver = m_domain_factorized_matrix[domain_index];

            const size_t size = indicies.size()+buffer.size()+random.size();
            //std::cout << "domain "<<domain_index<<"indicies = "<<indicies.size()<<" buffer =  "<<buffer.size()<<" random = "<<random.size()<<std::endl;

            domain_matrix.resize(size,size);

            size_t i = 0;
            for (const size_t& big_index_i: indicies) {
                size_t j = 0;
                for (const size_t& big_index_j: indicies) {
                    domain_matrix(i,j++) = mat.coeff(big_index_i,big_index_j);
                }
                for (const size_t& big_index_j: buffer) {
                    domain_matrix(i,j++) = mat.coeff(big_index_i,big_index_j);
                }
                for (const size_t& big_index_j: random) {
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
                for (const size_t& big_index_j: random) {
                    domain_matrix(i,j++) = mat.coeff(big_index_i,big_index_j);
                }
                ++i;
            }
            for (const size_t& big_index_i: random) {
                size_t j = 0;
                for (const size_t& big_index_j: indicies) {
                    domain_matrix(i,j++) = mat.coeff(big_index_i,big_index_j);
                }
                for (const size_t& big_index_j: buffer) {
                    domain_matrix(i,j++) = mat.coeff(big_index_i,big_index_j);
                }
                for (const size_t& big_index_j: random) {
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
        x = b;
        for (int i = 0; i < m_domain_indicies.size(); ++i) {
            if (m_domain_indicies.size() == 0) continue;

            const storage_vector_type& buffer = m_domain_buffer[i];
            const storage_vector_type& indicies = m_domain_indicies[i];
            const storage_vector_type& random = m_domain_random[i];
            
            const size_t nb = indicies.size()+buffer.size()+random.size();
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
            for (int j = 0; j < random.size(); ++j) {
                domain_b[sub_index++] = b[random[j]];
            }

            // solve domain 
            domain_x = m_domain_factorized_matrix[i].solve(domain_b);

            // copy accumulate b values to big vector
            for (int j = 0; j < indicies.size(); ++j) {
                x[indicies[j]] = domain_x[j];
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
