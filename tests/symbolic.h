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

        Symbol<scalar> s;
        Label<0,ParticlesType> a(particles);
    }

    void test_create_default_vectors(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<scalar> ParticlesType;
    	ParticlesType particles;

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<alive> alive_;

    }

    void test_transform(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<scalar> ParticlesType;
    	ParticlesType particles;

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<alive> alive_;
        Symbol<scalar> s;
        Label<0,ParticlesType> a(particles);

    	s[a] = 0;

    	ParticlesType::value_type particle;
    	particles.push_back(particle);
    	particles.push_back(particle);

    	s[a] = 0;

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),0);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),0);

    	s[a] = 1;

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),1);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),1);

    	s[a] = s + 1;

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),2);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),2);

        s[a] = s[a] + 1;

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),3);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),3);

    	p[a] = Vect3d(1,2,3);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],1);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],2);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],3);

        p[a] += Vect3d(1,2,3);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],2);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],4);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],6);

    	p[a] = p * s;

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],6);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],12);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],18);

       	p[a] = if_else(id_ == 0, Vect3d(0,0,0), Vect3d(3,2,1));

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

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<alive> alive_;
        Symbol<scalar> s;
        Label<0,ParticlesType> a(particles);
        Label<1,ParticlesType> b(particles);
        Dx dx;
        Accumulate<std::plus<double> > sum;

       	particles.push_back(Vect3d(0,0,0));
       	particles.push_back(Vect3d(diameter*2,0,0));

       	s[a] = sum(b, norm(dx) < diameter, 1);

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),1);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),1);

       	p[a] = if_else(id_ == 0, Vect3d(diameter/2.0,0,0), Vect3d(0,0,0));

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],diameter/2.0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],0);

    	p[a] = 0.5*(p + s);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0.5*(diameter/2.0 + 1));
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0.5);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0.5);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],0.5);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],0.5);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],0.5);

       	s[a] = sum(b, norm(dx) < diameter, 1);

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),2);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),2);

       	p[a] = if_else(id_ == 0, Vect3d(0,0,0), Vect3d(diameter/2.0,diameter/2.0,diameter/2.0));
       	s[a] = if_else(id_ == 0, 0, 1);

       	TS_ASSERT_EQUALS(get<scalar>(particles[0]),0);
       	TS_ASSERT_EQUALS(get<scalar>(particles[1]),1);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],diameter/2.0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],diameter/2.0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],diameter/2.0);

        Accumulate<std::plus<Vect3d> > sumVect;

    	p[a] = sumVect(b, norm(dx) < diameter, Vect3d(0,0,0) + 0.5*(s[a]/2.0 + s[b]/10.0));

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0.05);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0.05);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0.05);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],0.55);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],0.55);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],0.55);

    }


    void test_level0_expressions(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")

    	typedef Particles<scalar> ParticlesType;
       	ParticlesType particles;

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<alive> alive_;
        Symbol<scalar> s;
        Label<0,ParticlesType> a(particles);
        Label<1,ParticlesType> b(particles);
        Dx dx;

        Accumulate<std::plus<double> > sum;
        Accumulate<Aboria::max<double> > max;
        max.set_init(-1);
        Accumulate<Aboria::min<double> > min;
        min.set_init(1000);

       	particles.push_back(Vect3d(0,0,0));
       	particles.push_back(Vect3d(2,0,0));

        double result = eval(sum(a, p[0]<1, 1));
    	TS_ASSERT_EQUALS(result,1);
       	result = eval(sum(a, true, 1));
    	TS_ASSERT_EQUALS(result,2);
       	result = eval(sum(a, true, p[0]));
    	TS_ASSERT_EQUALS(result,2);
        int result2 = eval(max(a, true, id_));
    	TS_ASSERT_EQUALS(result2,1);
        result2 = eval(min(a, true, id_));
    	TS_ASSERT_EQUALS(result2,0);
       	particles.push_back(Vect3d(0,0,0));
        result2 = eval(max(a, true, id_));
    	TS_ASSERT_EQUALS(result2,2);
    }


};

#endif /* SYMBOLICTEST_H_ */
