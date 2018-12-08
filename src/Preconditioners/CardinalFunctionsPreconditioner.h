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

#ifndef CARDINAL_FUNCTIONS_PRECONDITIONER_H_
#define CARDINAL_FUNCTIONS_PRECONDITIONER_H_

#include "detail/Preconditioners.h"

#ifdef HAVE_CAIRO
#include <cairo-svg.h>
#endif

namespace Aboria {

class CardinalFunctionsPreconditioner {
  typedef double Scalar;
  typedef size_t Index;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> vector_type;
  typedef std::vector<size_t> storage_vector_type;
  typedef std::vector<storage_vector_type> connectivity_type;

protected:
  bool m_isInitialized;

private:
  size_t m_random;
  double m_sigma;
  double m_M;

  connectivity_type m_domain_buffer;
  std::vector<vector_type> m_weights;
  Index m_rows;
  Index m_cols;

public:
  typedef typename vector_type::StorageIndex StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  CardinalFunctionsPreconditioner()
      : m_isInitialized(false), m_random(0), m_sigma(-1), m_M(1.0) {}

  template <typename MatType>
  explicit CardinalFunctionsPreconditioner(const MatType &mat) {
    compute(mat);
  }

  Index rows() const { return m_rows; }
  Index cols() const { return m_cols; }

  void set_number_of_random_particles(size_t n) { m_random = n; }
  void set_sigma(double value) { m_sigma = value; }
  void set_rejection_sampling_scale(double value) { m_M = value; }

  template <typename Kernel>
  void analyze_impl_block(const Index start_row, const Kernel &kernel) {
    typedef typename Kernel::row_elements_type row_elements_type;
    typedef typename Kernel::col_elements_type col_elements_type;
    typedef typename row_elements_type::query_type query_type;
    typedef typename query_type::traits_type traits_type;
    typedef typename query_type::child_iterator child_iterator;
    typedef typename traits_type::double_d double_d;
    typedef typename traits_type::int_d int_d;
    typedef typename traits_type::position position;

    static_assert(
        std::is_same<row_elements_type, col_elements_type>::value,
        "Cardinal Functions preconditioner restricted to identical row and col "
        "particle sets");
    const row_elements_type &a = kernel.get_row_elements();
    CHECK(
        &a == &(kernel.get_col_elements()),
        "Cardinal Functions preconditioner restricted to identical row and col "
        "particle "
        "sets");
    const query_type &query = a.get_query();
    std::default_random_engine generator;
    m_domain_buffer.resize(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
      storage_vector_type &buffer = m_domain_buffer[i];

      // add buffer particles through random sampling
      int nspecial = std::pow(3, query_type::dimension);
      // const int nspecial = 0;
      buffer.resize(m_random + nspecial);
      std::vector<child_iterator> buckets(m_random + nspecial);
      std::uniform_real_distribution<double> uniform(0, 1);
      std::normal_distribution<double> normal(0, m_sigma);

      // add special points
      // updates nspecial with actual number of special points
      lattice_iterator<query_type::dimension> special_it(
          int_d::Constant(0), int_d::Constant(3), int_d::Constant(0));
      for (nspecial = 0; special_it != false; ++special_it, ++nspecial) {
        const double_d &bmin = query.get_bounds().bmin;
        const double_d &bmax = query.get_bounds().bmax;
        const double_d pos =
            (*special_it) *
                (0.5 * (bmax - bmin) - std::numeric_limits<double>::epsilon()) +
            bmin;
        // std::cout <<"adding special point at pos = "<<pos<<std::endl;
        buckets[nspecial] = query.get_bucket(pos);
      }

      const double_d middle = get<position>(a)[i];
      if (m_sigma > 0) {
        const double scale2 = 1.0 / std::pow(m_sigma, 2);
        auto gaussianf = [&](const double_d &x) {
          return std::exp(-(x - middle).squaredNorm() * scale2);
        };
        std::generate(buckets.begin() + nspecial, buckets.end(), [&]() {
          double_d sp;
          bool accepted;
          do {
            for (size_t i = 0; i < query_type::dimension; i++) {
              sp[i] = normal(generator) + middle[i];
            }
            const bool in_domain =
                (sp < a.get_max()).all() && (sp >= a.get_min()).all();
            accepted =
                in_domain && uniform(generator) <
                                 kernel.get_position_function()(middle, sp) /
                                     (gaussianf(sp) * m_M);
          } while (!accepted);
          return query.get_bucket(sp);
        });
      } else {
        const double volume = (a.get_max() - a.get_min()).prod();
        std::generate(buckets.begin() + nspecial, buckets.end(), [&]() {
          double_d sp;
          bool accepted;
          do {
            for (size_t i = 0; i < query_type::dimension; i++) {
              sp[i] =
                  0.5 * (a.get_max()[i] - a.get_min()[i]) * uniform(generator) +
                  a.get_min()[i];
            }
            const bool in_domain =
                (sp < a.get_max()).all() && (sp >= a.get_min()).all();
            accepted =
                in_domain &&
                uniform(generator) <
                    kernel.get_position_function()(middle, sp) * volume / m_M;
          } while (!accepted);
          return query.get_bucket(sp);
        });
      }

      std::unordered_map<size_t, std::pair<child_iterator, size_t>> counts;
      for (int i = 0; i < buckets.size(); ++i) {
        auto bucket_index = query.get_bucket_index(*(buckets[i]));
        auto it = counts.find(bucket_index);
        if (it != counts.end()) {
          it->second.second++;
        } else {
          counts[bucket_index] = std::make_pair(buckets[i], 1);
        }
      }
      // for (auto i : counts) {
      // std::cout << "bucket index " << i.first << " with bounds "
      //<< query.get_bounds(i.second.first) << " has " << i.second.second
      //<< " counts" << std::endl;
      //}

      int out_index = 0;
      std::for_each(counts.begin(), counts.end(), [&](auto i) {
        auto ci = i.second.first;
        size_t count = i.second.second;
        auto pit = query.get_bucket_particles(*ci);
        auto num_particles = pit.distance_to_end();
        std::vector<int> bucket_indices(num_particles);
        std::iota(bucket_indices.begin(), bucket_indices.end(), 0);
        std::shuffle(bucket_indices.begin(), bucket_indices.end(), generator);
        const int trunc_count = std::min(count, bucket_indices.size());
        std::transform(
            bucket_indices.begin(), bucket_indices.begin() + trunc_count,
            buffer.begin() + out_index, [&](const int i) {
              return (&get<position>(*(pit + i)) - &get<position>(a)[0]) +
                     start_row;
            });
        // std::cout << "looking for " << count
        //<< " samples in buffer. Found at indicies ";
        // for (size_t i = out_index; i < out_index + trunc_count; ++i) {
        //  std::cout << buffer[i] << " ";
        //}
        // std::cout << std::endl;
        out_index += trunc_count;
      });
      buffer.resize(out_index);

      // ensure that cardinal index isn't in buffer
      std::remove_if(buffer.begin(), buffer.end(),
                     [&i](const int j) { return j == i; });

#ifdef HAVE_CAIRO
      const int image_size = 512;
      cairo_surface_t *surface = cairo_svg_surface_create(
          ("sampler" + std::to_string(i) + ".svg").c_str(), image_size,
          image_size);
      cairo_svg_surface_restrict_to_version(surface, CAIRO_SVG_VERSION_1_2);
      cairo_t *cr = cairo_create(surface);
      const double lw = 0.01;

      cairo_scale(cr, image_size, image_size);
      cairo_set_line_width(cr, lw);
      const double PI = boost::math::constants::pi<double>();
      cairo_set_source_rgba(cr, 0.5, 0, 0, 1.0);

      auto &pos = get<position>(a)[i];
      cairo_arc(cr, pos[0], pos[1], lw, 0, 2 * PI);
      cairo_fill(cr);

      cairo_set_source_rgba(cr, 0, 0, 0.5, 1.0);
      for (auto i : buffer) {
        auto &pos = get<position>(a)[i];
        cairo_arc(cr, pos[0], pos[1], lw, 0, 2 * PI);
        cairo_fill(cr);
      }
      cairo_destroy(cr);
      cairo_surface_destroy(surface);
#endif

      ASSERT(buffer.size() > 0, "no particles in buffer");
    }
  }

  template <typename RowParticles, typename ColParticles>
  void
  analyze_impl_block(const Index start_row,
                     const KernelZero<RowParticles, ColParticles> &kernel) {}

  template <unsigned int NI, unsigned int NJ, typename Blocks, std::size_t... I>
  void analyze_impl(const MatrixReplacement<NI, NJ, Blocks> &mat,
                    detail::index_sequence<I...>) {
    int dummy[] = {0, (analyze_impl_block(mat.template start_row<I>(),
                                          std::get<I * NJ + I>(mat.m_blocks)),
                       0)...};
    static_cast<void>(dummy);
  }

  template <unsigned int NI, unsigned int NJ, typename Blocks>
  CardinalFunctionsPreconditioner &
  analyzePattern(const MatrixReplacement<NI, NJ, Blocks> &mat) {

    LOG(2, "CardinalFunctionsPreconditioner: analyze pattern");
    m_rows = mat.rows();
    m_cols = mat.cols();
    analyze_impl(mat, detail::make_index_sequence<NI>());

    int minsize_buffer = 1000;
    int maxsize_buffer = 0;
    for (size_t domain_index = 0; domain_index < m_domain_buffer.size();
         ++domain_index) {
      const int size_buffer = m_domain_buffer[domain_index].size();
      if (size_buffer < minsize_buffer)
        minsize_buffer = size_buffer;
      if (size_buffer > maxsize_buffer)
        maxsize_buffer = size_buffer;
    }
    LOG(2, "CardinalFunctionsPreconditioner: finished analysis, found "
               << m_domain_buffer.size() << " domains, with " << minsize_buffer
               << "--" << maxsize_buffer << " buffer particles")
    return *this;
  }

  template <int _Options, typename _StorageIndex>
  CardinalFunctionsPreconditioner &analyzePattern(
      const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex> &mat) {
    CHECK(m_domain_buffer.size() > 0,
          "CardinalFunctionsPreconditioner::analyzePattern(): cannot analyze "
          "sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <int _Options, typename _StorageIndex, int RefOptions,
            typename RefStrideType>
  CardinalFunctionsPreconditioner &
  analyzePattern(const Eigen::Ref<
                 const Eigen::SparseMatrix<Scalar, _Options, _StorageIndex>,
                 RefOptions, RefStrideType> &mat) {
    CHECK(m_domain_buffer.size() > 0,
          "CardinalFunctionsPreconditioner::analyzePattern(): cannot analyze "
          "sparse "
          "matrix, "
          "call analyzePattern using a Aboria MatrixReplacement class first");
    return *this;
  }

  template <typename Derived>
  CardinalFunctionsPreconditioner &
  analyzePattern(const Eigen::DenseBase<Derived> &mat) {
    CHECK(m_domain_buffer.size() > 0,
          "CardinalFunctionsPreconditioner::analyzePattern(): cannot analyze "
          "dense "
          "matrix, "
          "call analyzePattern need to pass a Aboria MatrixReplacement class "
          "first");
    return *this;
  }

  template <typename MatType>
  CardinalFunctionsPreconditioner &factorize(const MatType &mat) {
    LOG(2, "CardinalFunctionsPreconditioner: factorizing domain");
    eigen_assert(static_cast<typename MatType::Index>(m_rows) == mat.rows() &&
                 "CardinalFunctionsPreconditioner::solve(): invalid number of "
                 "rows of mat");
    eigen_assert(static_cast<typename MatType::Index>(m_cols) == mat.cols() &&
                 "CardinalFunctionsPreconditioner::solve(): invalid number of "
                 "rows of mat");

    matrix_type domain_matrix;
    const size_t N = mat.rows();
    m_weights.resize(m_domain_buffer.size());
    for (size_t domain_index = 0; domain_index < m_domain_buffer.size();
         ++domain_index) {
      const storage_vector_type &buffer = m_domain_buffer[domain_index];
      vector_type &weights = m_weights[domain_index];

      const size_t size = 1 + buffer.size();
      // std::cout << "domain "<<domain_index<<"indicies =
      // "<<indicies.size()<<" buffer =  "<<buffer.size()<<" random =
      // "<<random.size()<<std::endl;

      domain_matrix.resize(size, N);
      weights.resize(size);

      size_t i = 0;
      for (const size_t &big_index_i : {domain_index}) {
        for (size_t j = 0; j < N; ++j) {
          domain_matrix(i, j) = mat.coeff(big_index_i, j);
        }
        ++i;
      }
      for (const size_t &big_index_i : buffer) {
        for (size_t j = 0; j < N; ++j) {
          domain_matrix(i, j) = mat.coeff(big_index_i, j);
        }
        ++i;
      }

      vector_type b = vector_type::Zero(N);
      b[domain_index] = 1;

      weights = domain_matrix.transpose().colPivHouseholderQr().solve(b);
      // weights = domain_matrix.transpose()
      //              .bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
      //              .solve(b);
      double relative_error =
          (domain_matrix.transpose() * weights - b).norm() / b.norm();
      if (relative_error > 1e-3 || std::isnan(relative_error)) {
        std::cout << "domain index = " << domain_index
                  << ": relative error = " << relative_error << std::endl;
      }
    }

    m_isInitialized = true;

    return *this;
  }

  template <typename MatType>
  CardinalFunctionsPreconditioner &compute(const MatType &mat) {
    analyzePattern(mat);
    return factorize(mat);
  }

  /** \internal */
  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {
    for (size_t i = 0; i < m_domain_buffer.size(); ++i) {
      const storage_vector_type &buffer = m_domain_buffer[i];
      const vector_type &weights = m_weights[i];

      // x = W * b
      x[i] = weights[0] * b[i];
      for (size_t j = 0; j < buffer.size(); ++j) {
        x[i] += weights[j] * b[buffer[j]];
      }
    }
  }

  template <typename Rhs>
  inline const Eigen::Solve<CardinalFunctionsPreconditioner, Rhs>
  solve(const Eigen::MatrixBase<Rhs> &b) const {
    eigen_assert(static_cast<typename Rhs::Index>(m_rows) == b.rows() &&
                 "CardinalFunctionsPreconditioner::solve(): invalid number of "
                 "rows of the "
                 "right hand side matrix b");
    eigen_assert(m_isInitialized &&
                 "CardinalFunctionsPreconditioner is not initialized.");
    return Eigen::Solve<CardinalFunctionsPreconditioner, Rhs>(*this,
                                                              b.derived());
  }

  Eigen::ComputationInfo info() { return Eigen::Success; }
};

} // namespace Aboria

#endif
