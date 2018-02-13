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
#include "detail/GaussLegendre.h"
#include "detail/Particles.h"

#ifdef HAVE_EIGEN

namespace Aboria {

namespace detail {
template <typename RowRef, typename ColRef, typename F>
struct kernel_helper_ref {

  typedef typename std::result_of<F(RowRef, ColRef)>::type FunctionReturn;

  typedef typename std::conditional<std::is_arithmetic<FunctionReturn>::value,
                                    Eigen::Matrix<FunctionReturn, 1, 1>,
                                    FunctionReturn>::type Block;

  static const int block_rows = Block::RowsAtCompileTime;
  static const int block_cols = Block::ColsAtCompileTime;

  static_assert(block_rows > 0,
                "kernel function must return fixed size matrix");
  static_assert(block_cols > 0,
                "kernel function must return fixed size matrix");
};

template <typename RowElements, typename ColElements, typename F>
struct kernel_helper {
  static const unsigned int dimension = RowElements::dimension;
  typedef Vector<double, dimension> double_d;
  typedef double_d const &const_position_reference;
  typedef typename RowElements::const_reference const_row_reference;
  typedef typename ColElements::const_reference const_col_reference;

  typedef
      typename std::result_of<F(const_row_reference, const_col_reference)>::type
          FunctionReturn;

  typedef typename std::conditional<std::is_arithmetic<FunctionReturn>::value,
                                    Eigen::Matrix<FunctionReturn, 1, 1>,
                                    FunctionReturn>::type Block;

  static const int BlockRows = Block::RowsAtCompileTime;
  static const int BlockCols = Block::ColsAtCompileTime;
  typedef typename Block::Scalar Scalar;
  typedef typename Eigen::Matrix<Scalar, BlockRows, 1> BlockLHSVector;
  typedef typename Eigen::Matrix<Scalar, BlockCols, 1> BlockRHSVector;

  static_assert(BlockRows >= 0, "element type rows must be fixed");
  static_assert(BlockCols >= 0, "element type cols must be fixed");
};

template <size_t D, typename F> struct position_kernel_helper {
  static const unsigned int dimension = D;
  typedef Vector<double, dimension> double_d;
  typedef double_d const &const_position_reference;

  typedef typename std::result_of<F(
      const_position_reference, const_position_reference)>::type FunctionReturn;

  typedef typename std::conditional<std::is_arithmetic<FunctionReturn>::value,
                                    Eigen::Matrix<FunctionReturn, 1, 1>,
                                    FunctionReturn>::type Block;

  static const int BlockRows = Block::RowsAtCompileTime;
  static const int BlockCols = Block::ColsAtCompileTime;
  typedef typename Block::Scalar Scalar;
  typedef typename Eigen::Matrix<Scalar, BlockRows, 1> BlockLHSVector;
  typedef typename Eigen::Matrix<Scalar, BlockCols, 1> BlockRHSVector;

  static_assert(Block::RowsAtCompileTime >= 0,
                "element type rows must be fixed");
  static_assert(Block::ColsAtCompileTime >= 0,
                "element type cols must be fixed");
};

/*
template<typename Elements, typename F>
struct single_kernel_helper {
    static const unsigned int dimension = Elements::dimension;
    typedef Vector<double,dimension> double_d;
    typedef double_d const & const_position_reference;
    typedef typename Elements::const_reference const_reference;

    typedef typename std::result_of<F(const_position_reference,
                              const_reference,
                              const_col_reference)>::type FunctionReturn;

    typedef typename std::conditional<
                        std::is_arithmetic<FunctionReturn>::value,
                        Eigen::Matrix<FunctionReturn,1,1>,
                        FunctionReturn>::type Block;

    static_assert(Block::RowsAtCompileTime >= 0,"element type rows must be
fixed"); static_assert(Block::ColsAtCompileTime >= 0,"element type cols must be
fixed");
};
*/

template <typename RowElements, typename ColElements, typename Kernel,
          size_t QuadratureOrder = 8, typename = void>
struct integrate_over_element {};

template <typename RowElements, typename ColElements, typename Kernel,
          size_t QuadratureOrder>
struct integrate_over_element<
    RowElements, ColElements, Kernel, QuadratureOrder,
    typename std::enable_if<is_particles<RowElements>::value &&
                            is_particles<ColElements>::value>::type> {

  static const unsigned int dimension = RowElements::dimension;
  typedef typename detail::position_kernel_helper<dimension, Kernel>
      kernel_helper;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef typename kernel_helper::Block Block;
  typedef typename RowElements::position position;
  typedef typename kernel_helper::const_row_reference const_row_reference;
  typedef typename kernel_helper::const_col_reference const_col_reference;
  typedef double_d const &const_position_reference;

  const Kernel &m_kernel;
  const RowElements &m_row;
  const ColElements &m_col;
  integrate_over_element(const RowElements &row, const ColElements &col,
                         const Kernel &k)
      : m_kernel(k), m_row(row), m_col(col) {}

  Block operator()(const_row_reference a, const_col_reference b) const {
    return m_kernel(get<position>(a), get<position>(b));
  }
};

template <typename RowElements, typename ColElements, typename Kernel,
          size_t QuadratureOrder>
struct integrate_over_element<
    RowElements, ColElements, Kernel, QuadratureOrder,
    typename std::enable_if<is_particles<RowElements>::value &&
                            is_elements<2, ColElements>::value>::type> {

  static const unsigned int dimension = RowElements::dimension;
  typedef typename detail::position_kernel_helper<dimension, Kernel>
      kernel_helper;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef typename kernel_helper::Block Block;
  typedef typename kernel_helper::const_row_reference const_row_reference;
  typedef typename kernel_helper::const_col_reference const_col_reference;
  typedef typename ColElements::variable_type variable_type;
  typedef typename ColElements::particles_type col_particles_type;
  typedef typename RowElements::position position;
  typedef typename detail::GaussLegendre<QuadratureOrder> quadrature_type;

  const Kernel &m_kernel;
  const RowElements &m_row;
  const ColElements &m_col;
  const col_particles_type &m_colp;
  integrate_over_element(const RowElements &row, const ColElements &col,
                         const Kernel &k)
      : m_kernel(k), m_row(row), m_col(col), m_colp(col.get_particles()) {}

  Block operator()(const_row_reference a, const_col_reference b) const {
    Block result = Block::Zeros();
    auto b0 = m_colp.get_query().find(get<variable_type>(b)[0]);
    auto b1 = m_colp.get_query().find(get<variable_type>(b)[1]);
    ASSERT(b0 != iterator_to_raw_pointer(m_colp.end()), "cannot find b0");
    ASSERT(b1 != iterator_to_raw_pointer(m_colp.end()), "cannot find b1");

    const double_d node_scale = 0.5 * (*get<position>(b1) - *get<position>(b0));
    const double_d node_offset =
        0.5 * (*get<position>(b0) + *get<position>(b1));
    for (size_t i = 0; i < QuadratureOrder; ++i) {
      const double_d mapped_node =
          node_scale * quadrature_type::nodes[i] + node_offset;
      result += quadrature_type::weights[i] * m_kernel(a, mapped_node);
    }
    return result;
  }
};

/*
template <typename Elements, typename Kernel,
          size_t QuadratureOrder = 8,
          typename = void>
struct single_integrate_over_element {};


template <typename Elements, typename Kernel, size_t QuadratureOrder>
struct single_integrate_over_element<Elements, Kernel, QuadratureOrder,
          typename std::enable_if_t<is_particles<Elements>::value>> {

    typedef typename detail::kernel_helper<Elements,Kernel> kernel_helper;
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

template <typename Elements, typename Kernel, size_t QuadratureOrder>
struct single_integrate_over_element<Elements, Kernel, QuadratureOrder,
          typename std::enable_if_t<is_elements<2,Elements>::value>> {


    typedef typename detail::kernel_helper<Elements,Kernel> kernel_helper;
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
            a_b0 =
col.correct_dx_for_periodicity(*get<position>(b0)-get<position>(a)); a_b1 =
col.correct_dx_for_periodicity(*get<position>(b1)-get<position>(a)); } else {
            a_b0 = *get<position>(b0)-get<position>(a);
            a_b1 = *get<position>(b1)-get<position>(a);
        }
        const double_d scale = 0.5*(a_b1-a_b0);
        const double_d offset = 0.5*(a_b0+a_b1);
        for (size_t i = 0; i < QuadratureOrder; ++i) {
            const double_d mapped_node = scale*quadrature_type::nodes[i] +
offset; result += quadrature_type::weights[i] * m_kernel(mapped_node,a,b)

        }
        return result;
    }
};
*/

template <typename RowElements, typename ColElements> struct zero_kernel {
  typedef typename RowElements::const_reference const_row_reference;
  typedef typename ColElements::const_reference const_col_reference;

  double operator()(const_row_reference a, const_col_reference b) const {
    return 0.0;
  }
};

template <typename RowElements, typename ColElements, typename F>
struct position_kernel {
  const static unsigned int dimension = RowElements::dimension;
  typedef Vector<double, dimension> double_d;
  typedef position_d<dimension> position;
  typedef double_d const &const_position_reference;
  typedef typename RowElements::const_reference const_row_reference;
  typedef typename ColElements::const_reference const_col_reference;
  typedef typename std::result_of<F(
      const_position_reference, const_position_reference)>::type FunctionReturn;
  F m_f;
  position_kernel(const F f) : m_f(f) {}
  FunctionReturn operator()(const_row_reference a,
                            const_col_reference b) const {
    return m_f(get<position>(a), get<position>(b));
  }
};

template <typename RowElements, typename ColElements, typename FRadius,
          typename F>
struct sparse_kernel {
  const static unsigned int dimension = RowElements::dimension;
  typedef Vector<double, dimension> double_d;
  typedef double_d const &const_position_reference;
  typedef position_d<dimension> position;
  typedef typename RowElements::const_reference const_row_reference;
  typedef typename ColElements::const_reference const_col_reference;
  typedef
      typename std::result_of<F(const_position_reference, const_row_reference,
                                const_col_reference)>::type FunctionReturn;
  typedef typename std::conditional<std::is_arithmetic<FunctionReturn>::value,
                                    Eigen::Matrix<FunctionReturn, 1, 1>,
                                    FunctionReturn>::type Block;

  F m_f;
  FRadius m_fradius;
  const ColElements &m_col;

  sparse_kernel(const ColElements &col, const FRadius fradius, const F f)
      : m_f(f), m_fradius(fradius), m_col(col) {}
  Block operator()(const_row_reference a, const_col_reference b) const {

    const double_d dx =
        m_col.correct_dx_for_periodicity(get<position>(b) - get<position>(a));
    if (dx.squaredNorm() < std::pow(m_fradius(a), 2)) {
      return Block(m_f(dx, a, b));
    } else {
      return Block::Zero();
    }
  }
};

template <typename RowElements,
          typename const_row_reference = typename RowElements::const_reference>
struct constant_radius {
  const double m_radius;
  constant_radius(const double radius) : m_radius(radius) {}
  double operator()(const_row_reference a) const { return m_radius; }
};

template <typename Elements, size_t Repeats, size_t QuadratureOrder,
          typename = void>
struct integrate_chebyshev {};

template <typename Elements, size_t Repeats, size_t QuadratureOrder>
struct integrate_chebyshev<
    Elements, Repeats, QuadratureOrder,
    typename std::enable_if<is_particles<Elements>::value>::type> {

  static const unsigned int dimension = Elements::dimension;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef double_d const &const_position_reference;
  typedef typename Elements::const_reference const_reference;
  typedef typename Elements::position position;
  typedef bbox<dimension> box_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigen_matrix;
  typedef Eigen::Matrix<double, Repeats, Repeats> Block;

  size_t m_order;
  size_t m_ncheb;
  detail::ChebyshevRn<dimension> m_cheb;
  const Elements &m_elements;

  integrate_chebyshev(const Elements &elements, const size_t order,
                      const box_type &box)
      : m_order(order), m_ncheb(std::pow(order, dimension)), m_cheb(order, box),
        m_elements(elements) {}

  template <typename Derived>
  void operator()(const Eigen::DenseBase<Derived> &result) {
    for (size_t i = 0; i < m_elements.size(); ++i) {
      m_cheb.set_position(get<position>(m_elements)[i]);
      lattice_iterator<dimension> mj(int_d::Constant(0),
                                     int_d::Constant(m_order));
      for (size_t j = 0; j < m_ncheb; ++j, ++mj) {
        const_cast<Eigen::DenseBase<Derived> &>(result)
            .template block<Repeats, Repeats>(i * Repeats, j * Repeats) =
            m_cheb(*mj) * Block::Identity();
      }
    }
  }
};

template <typename Elements, size_t Repeats, size_t QuadratureOrder>
struct integrate_chebyshev<
    Elements, Repeats, QuadratureOrder,
    typename std::enable_if<is_elements<2, Elements>::value>::type> {

  static const unsigned int dimension = Elements::dimension;
  typedef typename Elements::particles_type particles_type;
  typedef typename Elements::variable_type variable_type;
  typedef typename particles_type::query_type query_type;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef double_d const &const_position_reference;
  typedef typename Elements::const_reference const_reference;
  typedef typename particles_type::position position;
  typedef bbox<dimension> box_type;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigen_matrix;
  typedef typename detail::GaussLegendre<QuadratureOrder> quadrature_type;
  typedef Eigen::Matrix<double, Repeats, Repeats> Block;

  size_t m_order;
  size_t m_ncheb;
  detail::ChebyshevRn<dimension> m_cheb;
  const Elements &m_elements;

  integrate_chebyshev(const Elements &elements, const size_t order,
                      const box_type &box)
      : m_order(order), m_ncheb(std::pow(order, dimension)), m_cheb(order, box),
        m_elements(elements) {}

  template <typename Derived>
  void operator()(const Eigen::DenseBase<Derived> &result) const {
    const_cast<Eigen::DenseBase<Derived> &>(result).setZero();
    const auto &query = m_elements.get_particles().get_query();
    for (size_t i = 0; i < m_elements.size(); ++i) {
      auto pa = query.find(get<variable_type>(m_elements)[i][0]);
      auto pb = query.find(get<variable_type>(m_elements)[i][1]);
      ASSERT(pa != query.get_particles_begin() + query.number_of_particles(),
             "cannot find a");
      ASSERT(pb != query.get_particles_begin() + query.number_of_particles(),
             "cannot find b");
      const double_d &a = *get<position>(pa);
      const double_d &b = *get<position>(pb);
      const double_d scale = 0.5 * (b - a);
      const double_d offset = 0.5 * (a + b);
      for (size_t q = 0; q < QuadratureOrder; ++q) {
        const double_d mapped_node = scale * quadrature_type::nodes[q] + offset;
        m_cheb.set_position(mapped_node);
        lattice_iterator<dimension> mj(int_d(0), int_d(m_order));
        for (size_t j = 0; j < m_ncheb; ++j, ++mj) {
          const_cast<Eigen::DenseBase<Derived> &>(result)
              .template block<Repeats, Repeats>(i * Repeats, j * Repeats) +=
              quadrature_type::weights[q] * cheb_rn(*mj) * Block::Identity();
        }
      }
    }
  }
};

} // namespace detail
} // namespace Aboria

#endif // HAVE_EIGEN

#endif
