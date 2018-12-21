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

#ifndef DUAL_BREADTH_FIRST_SEARCH_BASE_H_
#define DUAL_BREADTH_FIRST_SEARCH_BASE_H_

namespace Aboria {

template <typename RowQuery, typename ColQuery>
class DualBreadthFirstSearchIterator {
  using iterator = DualBreadthFirstSearchIterator<RowQuery, ColQuery>;
  typedef typename RowQuery::child_iterator row_child_iterator;
  typedef typename ColQuery::child_iterator col_child_iterator;
  typedef typename RowQuery::traits_type traits_t;

  struct row_col_pair {
    row_child_iterator row;
    col_child_iterator col;
  };

  typedef typename traits_t::template vector_type<row_col_pair>::type level_t;
  typedef typename traits_t::template vector_type<int>::type num_children_t;
  const RowQuery *m_row_query;
  const ColQuery *m_col_query;
  level_t m_level;
  num_children_t m_counts;
  level_t m_next_level;
  int m_level_num;
  bool m_finished;
  bool m_filtered;

public:
  typedef level_t const value_type;
  typedef level_t const *pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef level_t const &reference;
  typedef std::ptrdiff_t difference_type;

  DualBreadthFirstSearchIterator() : m_finished(true), m_filtered(false) {}

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  DualBreadthFirstSearchIterator(const row_child_iterator &start_row_ci,
                                 const RowQuery *row_query,
                                 const col_child_iterator &start_col_ci,
                                 const ColQuery *col_query)
      : m_row_query(row_query), m_col_query(col_query), m_level_num(1),
        m_finished(false), m_filtered(false) {
    if (start_row_ci != false && start_col_ci != false) {
      m_level.resize(start_row_ci.distance_to_end() *
                     start_col_ci.distance_to_end());
      int i = 0;
      for (auto ci = start_row_ci; ci != false; ++ci) {
        for (auto cj = start_col_ci; cj != false; ++cj) {
          m_level[i++] = row_col_pair{ci, cj};
        }
      }
    } else {
      m_finished = true;
#ifndef __CUDA_ARCH__
      LOG(3, "\tbf_iterator (constructor): start is false, no "
             "children to search.");
#endif
    }
  }

  reference previous() { return m_next_level; }

  CUDA_HOST_DEVICE
  reference operator*() { return dereference(); }
  CUDA_HOST_DEVICE
  pointer operator->() { return &dereference(); }
  CUDA_HOST_DEVICE
  iterator &operator++() {
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
  inline bool operator==(const iterator &rhs) const { return equal(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator!=(const iterator &rhs) const { return !operator==(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator==(const bool rhs) const { return equal(rhs); }

  CUDA_HOST_DEVICE
  inline bool operator!=(const bool rhs) const { return !operator==(rhs); }

  struct count_children {
    const RowQuery m_row_query;
    const ColQuery m_col_query;

    CUDA_HOST_DEVICE
    int operator()(const row_col_pair &rc_pair) {

      const int nrow =
          m_row_query.is_leaf_node(*rc_pair.row)
              ? 1
              : m_row_query.get_children(rc_pair.row).distance_to_end();
      const int ncol =
          m_col_query.is_leaf_node(*rc_pair.col)
              ? 1
              : m_row_query.get_children(rc_pair.col).distance_to_end();

      return nrow * ncol;
    }
  };

  template <typename T> struct count_filtered {
    count_children m_cc;
    T m_filter;

    CUDA_HOST_DEVICE
    int operator()(const row_col_pair &rc_pair) {
      return m_filter(rc_pair) ? 0 : m_cc(rc_pair);
    }
  };

  struct copy_children_and_leafs {
    int *m_counts;
    const RowQuery m_row_query;
    const ColQuery m_col_query;
    row_col_pair *m_next_level;
    row_col_pair *m_level;

    CUDA_HOST_DEVICE
    void work(const int i) {
      auto &rc_pair = m_level[i];
      int next_level_index = m_counts[i];
      if (m_row_query.is_leaf_node(*rc_pair.row) &&
          m_col_query.is_leaf_node(*rc_pair.col)) {
        m_next_level[next_level_index] = rc_pair;
      } else if (m_row_query.is_leaf_node(*rc_pair.row)) {
        for (auto col_ci = m_col_query.get_children(rc_pair.col);
             col_ci != false; ++col_ci) {
          m_next_level[next_level_index++] = row_col_pair{rc_pair.row, col_ci};
        }
      } else if (m_col_query.is_leaf_node(*rc_pair.col)) {
        for (auto row_ci = m_row_query.get_children(rc_pair.row);
             row_ci != false; ++row_ci) {
          m_next_level[next_level_index++] = row_col_pair{row_ci, rc_pair.col};
        }
      } else {
        for (auto row_ci = m_row_query.get_children(rc_pair.row);
             row_ci != false; ++row_ci) {
          for (auto col_ci = m_col_query.get_children(rc_pair.col);
               col_ci != false; ++col_ci) {
            m_next_level[next_level_index++] = row_col_pair{row_ci, col_ci};
          }
        }
      }
    }

    CUDA_HOST_DEVICE
    void operator()(const int i) {
      // if filtered do nothing
      if (m_counts[i + 1] == m_counts[i])
        return;

      work(i);
    }
  };

  template <typename T> void filter(const T &f) {
    m_counts.resize(m_level.size());
    auto count_f =
        count_filtered<T>{count_children{*m_row_query, *m_col_query}, f};
    detail::transform_exclusive_scan(m_level.begin(), m_level.end(),
                                     m_counts.begin(), count_f, 0,
                                     detail::plus());
    const size_t nchildren =
        static_cast<int>(m_counts.back()) + count_f(m_level.back());

    m_next_level.resize(nchildren);

    m_filtered = true;
  }

  template <typename T> void filter_with_gather(const T &f, level_t &store) {
    m_counts.resize(m_level.size());
    detail::transform(
        m_level.begin(), m_level.end(), m_counts.begin(),
        count_filtered<T>{count_children{*m_row_query, *m_col_query}, f});

    store.resize(m_level.size());
    auto new_end =
        detail::copy_if(m_level.begin(), m_level.end(), m_counts.begin(),
                        store.begin(), [](const int i) { return i == 0; });

    store.erase(new_end, store.end());

    size_t nchildren = m_counts.back();
    detail::exclusive_scan(m_counts.begin(), m_counts.end(), m_counts.begin(),
                           0);

    nchildren += m_counts.back();

    m_next_level.resize(nchildren);

    m_filtered = true;
  }

private:
  CUDA_HOST_DEVICE
  void increment() {
#ifndef __CUDA_ARCH__
    LOG(4, "\tincrement (bf_iterator): m_level size = " << m_level.size());
#endif
    // m_level [ci0, ci1, ci2, ...]
    // exclusive scan for # children + # leafs: n_child
    // (std::vector<Vector<int,2>> = [{0,0}, {3,0}, {3,1}, ..., {N-#child
    // cin,NL-#child==0}] resize m_next_level(N) resize m_leafs(NL)

    if (m_filtered) {
      m_filtered = false;
    } else {
      m_counts.resize(m_level.size());
      detail::transform_exclusive_scan(
          m_level.begin(), m_level.end(), m_counts.begin(),
          count_children{*m_row_query, *m_col_query}, 0, detail::plus());

      // resize for new children and leafs
      const int nchildren =
          static_cast<int>(m_counts.back()) +
          count_children{*m_row_query, *m_col_query}(m_level.back());

      m_next_level.resize(nchildren);
    }

    m_finished = m_next_level.size() == 0;

    // don't bother doing next level if finished
    if (m_finished) {
      m_level.swap(m_next_level);
      m_level_num++;
      return;
    }

// tabulate m_level to copy children to m_next_level, or leafs to
// m_leafs
#if defined(__CUDACC__)
    typedef typename thrust::detail::iterator_category_to_system<
        typename traits_t::vector_int::iterator::iterator_category>::type
        system;
    thrust::counting_iterator<int, system> count(0);
#else
    auto count = traits_t::make_counting_iterator(0);
#endif

    auto copy_f = copy_children_and_leafs{
        iterator_to_raw_pointer(m_counts.begin()), *m_row_query, *m_col_query,
        iterator_to_raw_pointer(m_next_level.begin()),
        iterator_to_raw_pointer(m_level.begin())};
    detail::for_each(count, count + m_level.size() - 1, copy_f);

    if (m_counts.back() != static_cast<int>(m_next_level.size())) {
      // last ci is not filtered
      copy_f.work(m_level.size() - 1);
    }

    // swap level back to m_level and increment level count
    m_level.swap(m_next_level);
    m_level_num++;
#ifndef __CUDA_ARCH__
    LOG(4, "\tend increment (bf_iterator):");
#endif
  }

  CUDA_HOST_DEVICE
  bool equal(iterator const &other) const {
    return m_row_query == other.m_row_query &&
           m_col_query == other.m_col_query && m_level_num == other.m_level_num;
  }

  CUDA_HOST_DEVICE bool equal(const bool other) const {
    return m_finished != other;
  }

  CUDA_HOST_DEVICE
  reference dereference() const { return m_level; }
};

template <typename RowQuery, typename ColQuery>
auto create_dual_breadth_first_iterator(const RowQuery &row_query,
                                        const ColQuery &col_query) {
  return DualBreadthFirstSearchIterator<RowQuery, ColQuery>(
      row_query.get_children(), &row_query, col_query.get_children(),
      &col_query);
}

} // namespace Aboria

#endif
