/*
 * particle_container.h
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

#ifndef PARTICLE_CONTAINER_H_
#define PARTICLE_CONTAINER_H_

#include <cxxtest/TestSuite.h>

#define LOG_LEVEL 1
#include "Aboria.h"

using namespace Aboria;


class ParticleContainerTest : public CxxTest::TestSuite {
public:
    void test_add_particle1(void) {
    	typedef Particles<> Test_type;
    	Test_type test;
    	Test_type::value_type p;
    	test.push_back(p);
    	TS_ASSERT_EQUALS(test.size(),1);
    }

    void test_add_particle2(void) {
    	typedef Particles<std::tuple<double> > Test_type;
    	Test_type test;
    	Test_type::value_type p;
    	test.push_back(p);
    	TS_ASSERT_EQUALS(test.size(),1);
    }

    void test_add_delete_particle(void) {
    	typedef Particles<std::tuple<double> > Test_type;
    	Test_type test;
    	Test_type::value_type p;
    	test.push_back(p);
    	test.pop_back();
    	TS_ASSERT_EQUALS(test.size(),0);
    	test.push_back(p);
    	test.push_back(p);
    	test.push_back(p);
    	TS_ASSERT_EQUALS(test.size(),3);
    	test.erase(test.begin());
    	TS_ASSERT_EQUALS(test.size(),2);
    	test.erase(test.begin(),test.end()-1);
    	TS_ASSERT_EQUALS(test.size(),0);
    }

};




#endif /* PARTICLE_CONTAINER_H_ */
