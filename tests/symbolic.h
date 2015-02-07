/*
 * symbolic.h
 * 
 * Copyright 2015 Martin Robinson
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
 *  Created on: 5 Feb 2015
 *      Author: robinsonm
 */

#ifndef SYMBOLICTEST_H_
#define SYMBOLICTEST_H_


#include <cxxtest/TestSuite.h>

#define LOG_LEVEL 1
#include "Aboria.h"

using namespace Aboria;


class SymbolicTest : public CxxTest::TestSuite {
public:
    void test_create_double_vector(void) {
    	typedef Particles<std::tuple<double> > ParticlesType;
    	ParticlesType particles;
    	auto theDouble = get_vector<0>(particles);
    }

    void test_create_default_vectors(void) {
    	typedef Particles<std::tuple<double> > ParticlesType;
    	ParticlesType particles;
    	auto position = get_vector<POSITION>(particles);
    	auto id = get_vector<ID>(particles);
    	auto alive = get_vector<ALIVE>(particles);
    }

    void test_transform(void) {
    	typedef Particles<std::tuple<double> > ParticlesType;
    	ParticlesType particles;
    	auto theDouble = get_vector<0>(particles);
    	auto position = get_vector<POSITION>(particles);
    	auto id = get_vector<ID>(particles);
    	auto alive = get_vector<ALIVE>(particles);

    	theDouble = 0;

    	ParticlesType::value_type p;
    	particles.push_back(p);
    	particles.push_back(p);

    	theDouble = 0;

    	TS_ASSERT_EQUALS(particles[0].get_elem<0>(),0);
    	TS_ASSERT_EQUALS(particles[1].get_elem<0>(),0);

    	theDouble = 1;

    	TS_ASSERT_EQUALS(particles[0].get_elem<0>(),1);
    	TS_ASSERT_EQUALS(particles[1].get_elem<0>(),1);

    	theDouble = theDouble + 1;

    	TS_ASSERT_EQUALS(particles[0].get_elem<0>(),2);
    	TS_ASSERT_EQUALS(particles[1].get_elem<0>(),2);

    	position = Vect3d(1,2,3);

    	TS_ASSERT_EQUALS(particles[0].get_position()[0],1);
    	TS_ASSERT_EQUALS(particles[0].get_position()[1],2);
    	TS_ASSERT_EQUALS(particles[0].get_position()[2],3);

    	position = position * theDouble;

    	TS_ASSERT_EQUALS(particles[0].get_position()[0],2);
    	TS_ASSERT_EQUALS(particles[0].get_position()[1],4);
    	TS_ASSERT_EQUALS(particles[0].get_position()[2],6);

       	position = if_else(id == 0, Vect3d(0,0,0), Vect3d(3,2,1));

       	TS_ASSERT_EQUALS(particles[0].get_position()[0],0);
       	TS_ASSERT_EQUALS(particles[0].get_position()[1],0);
       	TS_ASSERT_EQUALS(particles[0].get_position()[2],0);

       	TS_ASSERT_EQUALS(particles[1].get_position()[0],3);
       	TS_ASSERT_EQUALS(particles[1].get_position()[1],2);
       	TS_ASSERT_EQUALS(particles[1].get_position()[2],1);


    }

    void test_neighbours(void) {
       	typedef Particles<std::tuple<double> > ParticlesType;
       	ParticlesType particles;

//       	Vect3d min(-1,-1,-1);
//       	Vect3d max(1,1,1);
//       	Vect3d periodic(true,true,true);
       	double diameter = 0.1;
//       	particles.init_neighbour_search(min,max,diameter,periodic);

       	auto theDouble = get_vector<0>(particles);
       	auto position = get_vector<POSITION>(particles);
       	auto id = get_vector<ID>(particles);
       	auto alive = get_vector<ALIVE>(particles);

       	particles.push_back(Vect3d(0,0,0));
       	particles.push_back(Vect3d(diameter*2,0,0));

       	theDouble = sum(particles, radius < diameter, 1);

       	TS_ASSERT_EQUALS(particles[0].get_elem<0>(),1);
       	TS_ASSERT_EQUALS(particles[1].get_elem<0>(),1);

       	position = if_else(id == 0, Vect3d(0,0,0), Vect3d(diameter/2.0,0,0));

       	theDouble = sum(particles, radius < diameter, 1);

       	TS_ASSERT_EQUALS(particles[0].get_elem<0>(),2);
       	TS_ASSERT_EQUALS(particles[1].get_elem<0>(),2);

       }

};

#endif /* SYMBOLICTEST_H_ */
