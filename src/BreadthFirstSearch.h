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

namespace Aboria {

template <typename Query> class BreadthFirstSearchIterator {
  typedef typename Query::child_iterator child_iterator;
  typedef typename Query::Traits traits_t;

  typedef
      typename traits_t::template vector_type<child_iterator>::type m_level_t;
  typedef typename traits_t::template vector_type<int>::type m_num_children_t;
  const Query *m_query;
  m_level_t m_current_level;
  m_num_children_t m_num_children;
  m_level_t m_next_level;
  int m_depth;

public:
  typedef const m_level_t value_type;
  typedef m_level_t *const pointer;
  typedef std::forward_iterator_tag iterator_category;
  typedef m_level_t &const reference;
  typedef std::ptrdiff_t difference_type;

  CUDA_HOST_DEVICE
  BreadthFirstSearchIterator() {}

  /// this constructor is used to start the iterator at the head of a bucket
  /// list
  CUDA_HOST_DEVICE
  BreadthFirstSearchIterator(const Query *query)
      : m_query(query), m_current_level(query->get_child_iterator()),
        m_depth(0) {}

  CUDA_HOST_DEVICE
  reference operator*() const { return dereference(); }
  CUDA_HOST_DEVICE
  reference operator->() { return dereference(); }
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

private:
  CUDA_HOST_DEVICE
  void increment() {
#ifndef __CUDA_ARCH__
    LOG(4, "\tincrement (breadth_first_iterator) m_depth = " << m_depth);
#endif
    // count number of children
    auto count_children = [](child_iterator ci) {
      return ci.number_of_children();
    };
    m_num_children.resize(m_current_level.size());
    detail::transform_exclusive_scan(
        m_current_level.begin(), m_current_level.end(), m_num_children.begin(),
        count_children, std::plus<int>());
    int num_children = *(m_num_children.end() - 1) +
                       count_children(*(m_current_level.end() - 1));

    // create new level
    m_next_level.resize(num_children);
    detail::for_each(
        traits_t::make_zip_iterator(traits_t::make_tuple(
            m_current_level.begin(), m_num_children.begin())),
        traits_t::make_zip_iterator(
            traits_t::make_tuple(m_current_level.end(), m_num_children.end())),
        m_next_level.end(),
        [_m_query = m_query, _m_next_level = iterator_to_raw_pointer(
                                 m_next_level.begin())](auto i) {
          auto parent = i.template <0>();
          for (auto ci = _m_query->get_child_iterator(parent),
                    int next_index = i.template get<1>();
               ci != false; ++ci, ++next_index) {
            _m_next_level[next_index] = ci;
          }
        });

    // swap to current level
    m_current_level.swap(m_next_level);

#ifndef __CUDA_ARCH__
    LOG(4, "\tend increment (breadth_first_iterator):");
#endif
  }

  CUDA_HOST_DEVICE
  bool equal(iterator const &other) const {
    return m_depth == other.m_depth && m_query == other.m_query;
  }

  CUDA_HOST_DEVICE
  bool equal(const bool other) const {
    return m_current_level.empty() != other;
  }

  CUDA_HOST_DEVICE
  reference dereference() const { return m_current_level; }
};

template <typename Query> auto make_flat_tree(const Query &query) {
  typedef typename Query::child_iterator child_iterator;
  typedef typename Query::Traits traits_t;
  typedef typename traits_t::template vector_type<child_iterator>::type level_t;
  typedef typename traits_t::template vector_type<m_level_t>::type tree_t;

  tree_t tree;
  for (auto bf_it = BreadthFirstSearchIterator(query); bf_it != false;
       ++bf_it) {
    tree.push_back(*bf_it)
  }
  return tree;
}

} // namespace Aboria
