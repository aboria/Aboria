/*
 * ParticlesArray.h
 *
 * Copyright 2012 Martin Robinson
 *
  * This file is part of RD_3D.
 *
 * RD_3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RD_3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with RD_3D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 11 Oct 2012
 *      Author: robinsonm
 */

#ifndef PARTICLESARRAY_H_
#define PARTICLESARRAY_H_

#include "Particles.h"

namespace Aboria {

template <typename T>
class ParticlesArray {
public:
    typedef typename T::value_type value_type;

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
            : m_node(0) {}

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

        Value* m_node;
    };
    typedef impl::node_iterator<value_type> iterator;
    typedef impl::node_iterator<value_type const> const_iterator;


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


