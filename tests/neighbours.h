/*
 * neighbours.h
 *
 * Copyright 2014 Martin Robinson
 *
 * This file is part of Aboria.
 *
 * Aboria is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Aboria is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Aboria.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 27 Nov 2014
 *      Author: robinsonm
 */

#ifndef NEIGHBOURS_H_
#define NEIGHBOURS_H_

#include <cxxtest/TestSuite.h>

#define LOG_LEVEL 1
#include "Aboria.h"

using namespace Aboria;


class NeighboursTest : public CxxTest::TestSuite {
public:
    void test_single_particle(void) {
    	typedef Particles<std::tuple<double> > Test_type;
    	Test_type test;
    	Vect3d min(-1,-1,-1);
    	Vect3d max(1,1,1);
    	Vect3d periodic(true,true,true);
    	double diameter = 0.1;
    	test.init_neighbour_search(min,max,diameter,periodic);
    	Test_type::value_type p;

    	p.set_position(Vect3d(0,0,0));
    	test.push_back(p);

    	int count = 0;
    	for (auto tpl: test.get_neighbours(Vect3d(diameter/2,diameter/2,0))) {
    		count++;
    	}
    	TS_ASSERT_EQUALS(count,1);

    	auto tpl = test.get_neighbours(Vect3d(diameter/2,diameter/2,0));
    	TS_ASSERT_EQUALS(tpl.size(),1);

    	tpl = test.get_neighbours(Vect3d(2*diameter,0,0));
    	TS_ASSERT_EQUALS(tpl.size(),0);
    }

    void test_two_particles(void) {
    	typedef Particles<std::tuple<double> > Test_type;
    	Test_type test;
    	Vect3d min(-1,-1,-1);
    	Vect3d max(1,1,1);
    	Vect3d periodic(true,true,true);
    	double diameter = 0.1;
    	test.init_neighbour_search(min,max,diameter,periodic);
    	Test_type::value_type p;

    	p.set_position(Vect3d(0,0,0));
    	test.push_back(p);

    	p.set_position(Vect3d(diameter/2,0,0));
    	test.push_back(p);

    	auto tpl = test.get_neighbours(Vect3d(1.1*diameter,0,0));
    	TS_ASSERT_EQUALS(tpl.size(),1);
    	const Test_type::value_type &pfound = std::get<0>(*tpl.begin());
    	TS_ASSERT_EQUALS(pfound.get_id(),test[1].get_id());

    	tpl = test.get_neighbours(Vect3d(0.9*diameter,0,0));
    	TS_ASSERT_EQUALS(tpl.size(),2);

    	tpl = test.get_neighbours(Vect3d(1.6*diameter,0,0));
    	TS_ASSERT_EQUALS(tpl.size(),0);
    }

    void test_create_particles(void) {
    	typedef Particles<std::tuple<double> > Test_type;
    	Test_type test;
    	Vect3d min(-1,-1,-1);
    	Vect3d max(1,1,1);
    	Vect3d periodic(true,true,true);
    	double diameter = 0.1;
    	test.init_neighbour_search(min,max,diameter,periodic);
    	int count = 0;
    	test.create_particles(2,[&count,&test,&diameter](Test_type::value_type& i) {
    		auto tpl = test.get_neighbours(Vect3d(diameter/2,0,0));
    		count++;
    		if (count == 1) {
    			TS_ASSERT_EQUALS(tpl.size(),0);
    			return Vect3d(0,0,0);
    		} else {
    			TS_ASSERT_EQUALS(tpl.size(),1);
    			return Vect3d(diameter/2,0,0);
    		}

    	});

    	auto tpl = test.get_neighbours(Vect3d(1.1*diameter,0,0));
    	TS_ASSERT_EQUALS(tpl.size(),1);
    	const Test_type::value_type &pfound = std::get<0>(*tpl.begin());
    	TS_ASSERT_EQUALS(pfound.get_id(),test[1].get_id());

    	tpl = test.get_neighbours(Vect3d(0.9*diameter,0,0));
    	TS_ASSERT_EQUALS(tpl.size(),2);

    	tpl = test.get_neighbours(Vect3d(1.6*diameter,0,0));
    	TS_ASSERT_EQUALS(tpl.size(),0);
    }
};



#endif /* NEIGHBOURS_H_ */
