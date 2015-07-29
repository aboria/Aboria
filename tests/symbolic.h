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
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<scalar> ParticlesType;
    	ParticlesType particles;
    	auto scalar_ = get_vector<scalar>(particles);
    }

    void test_create_default_vectors(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<scalar> ParticlesType;
    	ParticlesType particles;
    	auto position_ = get_vector<position>(particles);
    	auto id_ = get_vector<id>(particles);
    	auto alive_ = get_vector<alive>(particles);
    }

    void test_transform(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<scalar> ParticlesType;
    	ParticlesType particles;
    	auto scalar_ = get_vector<scalar>(particles);
    	auto position_ = get_vector<position>(particles);
    	auto id_ = get_vector<id>(particles);
    	auto alive_ = get_vector<alive>(particles);

    	scalar_ = 0;

    	ParticlesType::value_type p;
    	particles.push_back(p);
    	particles.push_back(p);

    	scalar_ = 0;

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),0);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),0);

    	scalar_= 1;

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),1);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),1);

    	scalar_ = scalar_ + 1;

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),2);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),2);

    	position_ = Vect3d(1,2,3);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],1);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],2);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],3);

        position_ += Vect3d(1,2,3);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],2);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],4);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],6);

    	position_ = position_ * scalar_;

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],4);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],8);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],12);

       	position_ = if_else(id_ == 0, Vect3d(0,0,0), Vect3d(3,2,1));

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0);
        
        TS_ASSERT_EQUALS(get<position>(particles[1])[0],3);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],2);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],1);
    }

    void test_neighbours(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<scalar> ParticlesType;
       	ParticlesType particles;

        Vect3d min(-1,-1,-1);
        Vect3d max(1,1,1);
        Vect3d periodic(true,true,true);
       	double diameter = 0.1;
        particles.init_neighbour_search(min,max,diameter,periodic);

       	auto scalar_ = get_vector<scalar>(particles);
       	auto position_ = get_vector<position>(particles);
       	auto id_ = get_vector<id>(particles);
       	auto alive_ = get_vector<alive>(particles);

       	particles.push_back(Vect3d(0,0,0));
       	particles.push_back(Vect3d(diameter*2,0,0));

       	Label<0> a;
       	Label<1> b;
        Dx dx;
        Accumulate<std::plus<double> > sum;

       	scalar_ = sum(b=particles, norm(dx) < diameter, 1,0);

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),1);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),1);

       	position_ = if_else(id_ == 0, Vect3d(diameter/2.0,0,0), Vect3d(0,0,0));

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],diameter/2.0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],0);

    	position_ = 0.5*(position_ + scalar_);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0.5*(diameter/2.0 + 1));
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0.5);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0.5);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],0.5);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],0.5);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],0.5);

//       	old_posiiton   = position;
//       	for() {
//          ellipsoids = ellipsoid(position)
//          new_position = constrain(position + 
//       	position = position + sqrt(2*D*dt)*normal() + dt*interpolate(position, drift);
//       	position = any(particles, norm(dx()) < diameter, reflect(position,new_position,ellipsoids))
//       	}
//       	dx =position-old_position;
//       	msd = dot(dx,dx);
//       	mean = mean(dx);
//       	var = var(dx);

       	scalar_ = sum(b=particles, norm(dx) < diameter, 1);

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),2);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),2);

       	position_ = if_else(id_ == 0, Vect3d(0,0,0), Vect3d(diameter/2.0,diameter/2.0,diameter/2.0));
       	scalar_ = if_else(id_ == 0, 0, 1);

       	TS_ASSERT_EQUALS(get<scalar>(particles[0]),0);
       	TS_ASSERT_EQUALS(get<scalar>(particles[1]),1);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],diameter/2.0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],diameter/2.0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],diameter/2.0);

        Accumulate<std::plus<Vect3d> > sumVect;

    	position_ = sumVect(b=particles, norm(dx) < diameter, Vect3d(0,0,0) + 0.5*(scalar_[a]/2.0 + scalar_[b]/10.0),0);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0.05);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0.05);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0.05);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],0.55);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],0.55);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],0.55);

    }

};

#endif /* SYMBOLICTEST_H_ */
