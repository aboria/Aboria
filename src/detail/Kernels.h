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


#ifndef DETAIL_KERNELS_H_
#define DETAIL_KERNELS_H_

#include "Elements.h"

namespace Aboria {

namespace detail {

    template<typename RowElements, typename ColElements, typename F>
    struct kernel_helper {
        static const unsigned int dimension = RowElements::dimension;
        typedef Vector<double,dimension> double_d;
        typedef double_d const & const_position_reference;
        typedef typename RowElements::const_reference const_row_reference;
        typedef typename ColElements::const_reference const_col_reference;

        typedef typename std::result_of<F(const_position_reference, 
                                  const_row_reference, 
                                  const_col_reference)>::type FunctionReturn;

        typedef typename std::conditional<
                            std::is_arithmetic<FunctionReturn>::value,
                            Eigen::Matrix<FunctionReturn,1,1>,
                            FunctionReturn>::type Block;

        static_assert(Block::RowsAtCompileTime >= 0,"element type rows must be fixed");
        static_assert(Block::ColsAtCompileTime >= 0,"element type cols must be fixed");
    };

    
    template <typename RowElements, typename ColElements, typename Kernel,
              size_t QuadratureOrder = 8,
              typename = std::enable_if_t<
                        is_particles<RowElements>::value &&
                        is_particles<ColElements>::value>
                        >
    struct integrate_over_element {
        typedef typename detail::kernel_helper<RowElements,ColElements,Kernel> kernel_helper;
        static const unsigned int dimension = kernel_helper::dimension;
        typedef Vector<double,dimension> double_d;
        typedef typename kernel_helper::Block Block;
        typedef typename RowElements::position position;
        typedef typename kernel_helper::const_row_reference const_row_reference;
        typedef typename kernel_helper::const_col_reference const_col_reference;


        const Kernel& m_kernel;
        const RowElements& m_row; 
        const ColElements& m_col; 
        const bool m_periodic;
        integrate_over_element(const RowElements& row, 
                               const ColElements& col,
                               const Kernel& k):
            m_kernel(k),
            m_row(row),
            m_col(col),
            m_periodic(!col.get_periodic().any())
        {}

        Block operator()(const_row_reference a,
                         const_col_reference b) const {
            double_d dx; 
            if (m_periodic) { 
                dx = m_col.correct_dx_for_periodicity(get<position>(b)-get<position>(a));
            } else {
                dx = get<position>(b)-get<position>(a);
            }
            return m_kernel(dx,a,b);
        }

        Block operator()(const double_d& dx,
                         const_row_reference a,
                         const_col_reference b) const {
            return m_kernel(dx,a,b);
        }
    };

    template <typename RowElements, typename ColElements, typename Kernel,
              size_t QuadratureOrder = 8,
              typename = std::enable_if_t<
                        is_particles<RowElements>::value &&
                        is_elements<2,ColElements>::value>
                        >
    struct integrate_over_element {
        typedef typename detail::kernel_helper<RowElements,ColElements,F> kernel_helper;
        static const unsigned int dimension = kernel_helper::dimension;
        typedef kernel_helper::Block Block;
        typedef typename kernel_helper::const_row_reference const_row_reference;
        typedef typename kernel_helper::const_col_reference const_col_reference;
        typedef typename ColParticles::variable_type variable_type;
        typedef typename RowElements::position position;
        typedef typename detail::GaussLegendre<QuadratureOrder> quadrature_type;

        const Kernel& m_kernel;
        const RowElements& m_row; 
        const ColElements& m_col; 
        const ColParticles& m_colp; 
        const bool m_periodic;
        integrate_over_element(const RowElements& row, 
                               const ColElements& col,
                               const Kernel& k):
            m_kernel(k),
            m_row(row),
            m_col(col),
            m_colp(col.get_particles()),
            m_periodic(!col.get_periodic().any())
        {}

        Block operator()(const_row_reference a,
                          const_col_reference b) const {
            Block result = Block::Zeros();
            double_d a_b0,a_b1; 
            auto b0 = m_colp.get_query().find(get<variable_type>(b)[0]);
            auto b1 = m_colp.get_query().find(get<variable_type>(b)[1]);
            ASSERT(b0 != iterator_to_raw_pointer(m_colp.end()),"cannot find b0");
            ASSERT(b1 != iterator_to_raw_pointer(m_colp.end()),"cannot find b1");
            if (m_periodic) { 
                a_b0 = col.correct_dx_for_periodicity(*get<position>(b0)-get<position>(a));
                a_b1 = col.correct_dx_for_periodicity(*get<position>(b1)-get<position>(a));
            } else {
                a_b0 = *get<position>(b0)-get<position>(a);
                a_b1 = *get<position>(b1)-get<position>(a);
            }
            const double_d scale = 0.5*(a_b1-a_b0);
            const double_d offset = 0.5*(a_b0+a_b1);
            for (size_t i = 0; i < QuadratureOrder; ++i) {
                const double_d mapped_node = scale*quadrature_type::nodes[i] + offset;
                result += quadrature_type::weights[i] * m_kernel(mapped_node,a,b)
                
            }
            return result;
        }
    };

    template <typename Elements, typename Kernel,
              size_t QuadratureOrder = 8,
              typename = std::enable_if_t<
                        is_particles<Elements>::value
                        >
    struct integrate_over_element {
        typedef typename detail::kernel_helper<RowElements,ColElements,F> kernel_helper;
        const unsigned int dimension = kernel_helper::dimension;
        typedef kernel_helper::Block Block;
        typedef typename kernel_helper::const_row_reference const_row_reference;
        typedef typename kernel_helper::const_col_reference const_col_reference;


        const Kernel& m_kernel;
        const Elements& m_elements; 
        integrate_over_element(const RowElements& row, 
                               const ColElements& col,
                               const Kernel& k):
            m_kernel(k),
            m_row(row),
            m_col(col),
            m_periodic(!col.get_periodic().any())
        {}

        Block operator()(const_reference a) const {
            return m_kernel(get<position>(a));
        }
    };

    template <typename RowElements, typename ColElements, typename Kernel,
              size_t QuadratureOrder = 8,
              typename = std::enable_if_t<
                        is_particles<RowElements>::value &&
                        is_elements<2,ColElements>::value>
                        >
    struct integrate_over_element {
        typedef typename detail::kernel_helper<RowElements,ColElements,F> kernel_helper;
        const unsigned int dimension = kernel_helper::dimension;
        typedef kernel_helper::Block Block;
        typedef typename kernel_helper::const_row_reference const_row_reference;
        typedef typename kernel_helper::const_col_reference const_col_reference;
        typedef typename ColParticles::variable_type variable_type;
        typedef typename detail::GaussLegendre<QuadratureOrder> quadrature_type;

        const Kernel& m_kernel;
        const RowElements& m_row; 
        const ColElements& m_col; 
        const ColParticles& m_colp; 
        const bool m_periodic;
        integrate_over_element(const RowElements& row, 
                               const ColElements& col,
                               const Kernel& k):
            m_kernel(k),
            m_row(row),
            m_col(col),
            m_colp(col.get_particles()),
            m_periodic(!col.get_periodic().any())
        {}

        Block operator()(const_row_reference a,
                          const_col_reference b) const {
            Block result = Block::Zeros();
            double_d a_b0,a_b1; 
            auto b0 = m_colp.get_query().find(get<variable_type>(b)[0]);
            auto b1 = m_colp.get_query().find(get<variable_type>(b)[1]);
            ASSERT(b0 != iterator_to_raw_pointer(m_colp.end()),"cannot find b0");
            ASSERT(b1 != iterator_to_raw_pointer(m_colp.end()),"cannot find b1");
            if (m_periodic) { 
                a_b0 = col.correct_dx_for_periodicity(*get<position>(b0)-get<position>(a));
                a_b1 = col.correct_dx_for_periodicity(*get<position>(b1)-get<position>(a));
            } else {
                a_b0 = *get<position>(b0)-get<position>(a);
                a_b1 = *get<position>(b1)-get<position>(a);
            }
            const double_d scale = 0.5*(a_b1-a_b0);
            const double_d offset = 0.5*(a_b0+a_b1);
            for (size_t i = 0; i < QuadratureOrder; ++i) {
                const double_d mapped_node = scale*quadrature_type::nodes[i] + offset;
                result += quadrature_type::weights[i] * m_kernel(mapped_node,a,b)
                
            }
            return result;
        }
    };





    template <typename RowElements, typename ColElements, typename F>
    struct position_kernel {
        const unsigned int dimension = RowElements::dimension;
        typedef Vector<double,dimension> double_d;
        typedef double_d const & const_position_reference;
        typedef typename RowElements::const_reference const_row_reference;
        typedef typename ColElements::const_reference const_col_reference;


        F m_f;
        position_kernel(const F f):m_(f) {}
        double operator()(const double_d& dx,
                        const_row_reference a,
                        const_col_reference b) const {
                            return m_f(dx,get<position>(a),get<position>(b));
        }
    };

    template <typename Elements, size_t Repeats,
               typename = std::enable_if_t<
                        is_particles<Elements>::value
                        >
    struct integrate_chebyshev {
        const unsigned int dimension = Elements::dimension;
        typedef Vector<double,dimension> double_d;
        typedef double_d const & const_position_reference;
        typedef typename Elements::const_reference const_reference;
        typedef detail::bbox<dimension> box_type;
        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix;

        size_t m_order;
        size_t m_ncheb;
        box_type m_box;
        detail::ChebyshevRn<dimension> m_cheb;
        const Elements& m_elements;

        chebyshev_kernel(const Elements& elements, const size_t order, const box_type box):
            m_order(order),
            m_ncheb(std::pow(order,dimension)),
            m_box(box),
            m_cheb(order,box),
            m_elements(elements) 
        {}

        void operator()(eigen_matrix& result) const {
            for (int i = 0; i < elements.size(); ++i) {
                m_cheb.set_position(get<position>(elements)[i]);
                lattice_iterator<dimension> mj(int_d(0),int_d(m_order));
                for (int j=0; j<m_ncheb; ++j,++mj) {
                    result.block<Repeats,Repeats>(i*Repeats,j*Repeats) = m_cheb(*mj); 
                }
            }
        }
    };

    template <typename Elements,
               typename = std::enable_if_t<
                        is_elements<2,Elements>::value
                        >
    struct integrate_chebyshev {
        const unsigned int dimension = Elements::dimension;
        typedef typename Elements::particles_type particles_type;
        typedef typename Elements::variable_type variable_type;
        typedef typename particles_type::query_type query_type;
        typedef Vector<double,dimension> double_d;
        typedef double_d const & const_position_reference;
        typedef typename Elements::const_reference const_reference;
        typedef detail::bbox<dimension> box_type;
        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> eigen_matrix;

        size_t m_order;
        size_t m_ncheb;
        box_type m_box;
        detail::ChebyshevRn<dimension> m_cheb;
        const Elements& m_elements;

        chebyshev_kernel(const Elements& elements, const size_t order, const box_type box):
            m_order(order),
            m_ncheb(std::pow(order,dimension)),
            m_box(box),
            m_cheb(order,box),
            m_elements(elements)
        {}

        void operator()(eigen_matrix& result) const {
            result = eigen_matrix::Zeros(elements.size()*Repeats,m_ncheb);
            const auto& query = m_elements.get_particles().get_query();
            for (int i = 0; i < elements.size(); ++i) {
                auto pa = query.find(get<variable_type>(elements)[i][0]);
                auto pb = query.find(get<variable_type>(elements)[i][1]);
                ASSERT(pa != query.get_particles_begin()+m_query.number_of_particles(),"cannot find a");
                ASSERT(pb != query.get_particles_begin()+m_query.number_of_particles(),"cannot find b");
                const double_d& a = *get<position>(pa);
                const double_d& b = *get<position>(pb);
                const double_d scale = 0.5*(b-a);
                const double_d offset = 0.5*(a+b);
                for (size_t q = 0; q < QuadratureOrder; ++q) {
                    const double_d mapped_node = scale*quadrature_type::nodes[q] + offset;
                    m_cheb.set_position(mapped_node);
                    lattice_iterator<dimension> mj(int_d(0),int_d(m_order));
                    for (int j=0; j<m_ncheb; ++j,++mj) {
                        result.block<Repeats,Repeats>(i*Repeats,j*Repeats) += quadrature_type::weights[q]*cheb_rn(*mj);
                    }
                }
            }
        }
    };



    

} //namespace detail
} //namespace Aboria

#endif

