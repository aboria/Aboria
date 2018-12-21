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

#ifndef BREADTH_FIRST_SEARCH_BASE_H_
#define BREADTH_FIRST_SEARCH_BASE_H_

#include "Vector.h"
#include "detail/Algorithms.h"

namespace Aboria {

template <typename Query> class BreadthFirstSearchIterator {
  typedef BreadthFirstSearchIterator<Query> iterator;
  typedef typename Query::child_iterator child_iterator;
  static const unsigned int dimension = Query::dimension;
  typedef Vector<double, dimension> double_d;
  typedef Vector<int, dimension> int_d;
  typedef typename Query::traits_type traits_type;
  typedef typename traits_type::template vector<child_iterator> vector_ci;
  typedef typename traits_type::template vector<int> vector_int;

  bool m_finished;
  size_t m_level_num;
  vector_ci m_level;
  vector_ci m_next_level;
  vector_int m_counts;
  const Query *m_query;
  bool m_filtered;

public:
  typedef const vector_ci value_type;
  typedef vector_ci const *pointer;
  typedef vector_ci const &reference;
  typedef std::forward_iterator_tag iterator_category;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  BreadthFirstSearchIterator() : m_finished(true), m_filtered(false) {}

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  CUDA_HOST_DEVICE
  BreadthFirstSearchIterator(const child_iterator &start_node,
                             const Query *query)
      : m_finished(false), m_level_num(1), m_query(query), m_filtered(false) {
    if (start_node != false) {
      m_level.resize(start_node.distance_to_end());
      int i = 0;
      for (auto ci = start_node; ci != false; ++ci, ++i) {
        m_level[i] = ci;
      }
    } else {
      m_finished = true;
#ifndef __CUDA_ARCH__
      LOG(3, "\tBreadthFirstSearchIterator (constructor): start is false, no "
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
    const Query m_query;

    CUDA_HOST_DEVICE
    int operator()(const child_iterator &ci) {
      const int nchild = m_query.is_leaf_node(*ci)
                             ? 1
                             : m_query.get_children(ci).distance_to_end();
      return nchild;
    }
  };

  template <typename T> struct count_filtered {
    count_children m_cc;
    T m_filter;

    CUDA_HOST_DEVICE
    int operator()(const child_iterator &ci) {
      return m_filter(ci) ? 0 : m_cc(ci);
    }
  };

  struct copy_children_and_leafs {
    int *m_counts;
    const Query m_query;
    child_iterator *m_next_level;
    child_iterator *m_level;
    int m_nchildren;

    void work(const int i) {
      auto ci = m_level[i];

      if (m_query.is_leaf_node(*ci)) {
        m_next_level[m_counts[i]] = ci;
      } else {
        auto child = m_query.get_children(ci);
        for (int next_level_index = m_counts[i]; child != false;
             ++child, ++next_level_index) {
          m_next_level[next_level_index] = child;
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
    auto count_f = count_filtered<T>{count_children{*m_query}, f};
    detail::transform_exclusive_scan(m_level.begin(), m_level.end(),
                                     m_counts.begin(), count_f, 0,
                                     detail::plus());

    const size_t nchildren =
        static_cast<int>(m_counts.back()) + count_f(m_level.back());

    m_next_level.resize(nchildren);

    m_filtered = true;
  }

  template <typename T> void filter_with_gather(const T &f, vector_ci &store) {
    m_counts.resize(m_level.size());
    detail::transform(m_level.begin(), m_level.end(), m_counts.begin(),
                      count_filtered<T>{count_children{*m_query}, f});

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
    LOG(4, "\tincrement (BreadthFirstSearchIterator): m_level size = "
               << m_level.size());
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
          count_children{*m_query}, 0, detail::plus());

      const size_t nchildren =
          m_counts.back() + count_children{*m_query}(m_level.back());

      m_next_level.resize(nchildren);
    }

    // resize for new children and leafs
    m_finished = m_next_level.size() == 0;

    // don't bother doing next level if they are all leafs
    if (m_finished) {
      m_level.swap(m_next_level);
      m_level_num++;
      return;
    }

// tabulate m_level to copy children to m_next_level, or leafs to
// m_leafs
#if defined(__CUDACC__)
    typedef typename thrust::detail::iterator_category_to_system<
        typename Query::traits_type::vector_int::iterator::iterator_category>::
        type system;
    thrust::counting_iterator<int, system> count(0);
#else
    auto count = Query::traits_type::make_counting_iterator(0);
#endif
    auto copy_f = copy_children_and_leafs{
        iterator_to_raw_pointer(m_counts.begin()), *m_query,
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
    LOG(4, "\tend increment (BreadthFirstSearchIterator):");
#endif
  }

  CUDA_HOST_DEVICE
  bool equal(iterator const &other) const {
    return m_query == other.m_query && m_level_num == other.m_level_num;
  }

  CUDA_HOST_DEVICE bool equal(const bool other) const {
    return m_finished != other;
  }

  CUDA_HOST_DEVICE
  reference dereference() const { return m_level; }
};

template <typename Query> auto make_flat_tree(const Query &query) {
  typedef typename Query::child_iterator child_iterator;
  typedef typename Query::traits_type traits_t;
  typedef typename traits_t::template vector_type<child_iterator>::type level_t;
  typedef typename traits_t::template vector_type<level_t>::type tree_t;

  tree_t tree;
  for (auto bf_it =
           BreadthFirstSearchIterator<Query>(query.get_children(), &query);
       bf_it != false; ++bf_it) {
    tree.push_back(*bf_it);
    bf_it.filter([&](auto ci) { return query.is_leaf_node(*ci); });
  }
  return tree;
}

} // namespace Aboria
#endif
