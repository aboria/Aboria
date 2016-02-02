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


#ifndef PARTICLESARRAY_H_
#define PARTICLESARRAY_H_

#include "Particles.h"

namespace Aboria {

template<class C>
    auto join(C&& c)
    -> decltype(boost::make_iterator_range(std::begin(c),std::end(c))) {
        return boost::make_iterator_range(std::begin(c),std::end(c));
    }

template<class C, class D, class... Args>
    auto join(C&& c, D&& d, Args&&... args)
    -> decltype(boost::join(boost::join(boost::make_iterator_range(std::begin(c),std::end(c)),
                                                     boost::make_iterator_range(std::begin(d),std::end(d))),
                                     join(std::forward<Args>(args)...))) {
          return boost::join(boost::join(boost::make_iterator_range(std::begin(c),std::end(c)),
                                                       boost::make_iterator_range(std::begin(d),std::end(d))),
                                       join(std::forward<Args>(args)...));
    }

template <typename T>
class ParticlesArray {
public:
    typedef typename T::value_type value_type;
    typedef typename T::iterator set_iterator;
    typedef typename T::const_iterator const_set_iterator;

    ParticlesArray(std::initializer_list<T> sets): sets(sets) {}


    template <class Value>
    class node_iter
        : public boost::iterator_facade<
            node_iter<Value>
            , Value
            , boost::forward_traversal_tag
          >
    {
    public:
        node_iter()
            : m_node(0),set_index(0) {}

        explicit node_iter(Value* p)
            : m_node(p) {}

        template <class OtherValue>
        node_iter(
            node_iter<OtherValue> const& other
            , typename boost::enable_if<
                boost::is_convertible<OtherValue*,Value*>
                , enabler
              >::type = enabler()
        )
            : m_node(other.m_node) {}

    private:
        friend class boost::iterator_core_access;
        template <class> friend class node_iter;

        template <class OtherValue>
        bool equal(node_iter<OtherValue> const& other) const
        {
            return this->m_node == other.m_node;
        }

        void increment()
        { m_node = m_node->next(); }

        Value& dereference() const
        { return *m_node; }

        Value m_set_iterator;
        int set_index;
    };
    typedef impl::node_iterator<set_iterator> iterator;
    typedef impl::node_iterator<const_set_iterator> const_iterator;


    boost::iterator_range<typename NeighbourSearch_type::const_iterator> get_neighbours(const Vect3d position) {
        ASSERT(searchable == true,"ERROR: using get_neighbours before initialising neighbour search. Please call the init_neighbour_search function before using get_neighbours");
        return boost::make_iterator_range(
                neighbour_search.find_broadphase_neighbours(position, -1,false),
                neighbour_search.end());

    iterator begin() {
        return data.begin();
    }
    iterator end() {
        return data.end();
    }


    }


private:

    std::vector<T> sets;
}


}

#endif PARTICLESARRAY_H_


