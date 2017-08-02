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


#ifndef SYMBOLICTEST_H_
#define SYMBOLICTEST_H_


#include <cxxtest/TestSuite.h>

#include "Aboria.h"

using namespace Aboria;


class SymbolicTest : public CxxTest::TestSuite {
public:
    void helper_create_double_vector(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<std::tuple<scalar>> ParticlesType;
    	ParticlesType particles;

        Symbol<scalar> s;
        Label<0,ParticlesType> a(particles);
    }

    void helper_create_default_vectors(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<std::tuple<scalar>> ParticlesType;
        typedef position_d<3> position;
    	ParticlesType particles;

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<alive> alive_;

    }

    void helper_transform(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<std::tuple<scalar>> ParticlesType;
        typedef position_d<3> position;
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

    	s[a] = s[a] + 1;

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),2);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),2);

        s[a] = s[a] + 1;

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),3);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),3);

    	p[a] = vdouble3(1,2,3);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],1);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],2);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],3);

        p[a] += vdouble3(1,2,3);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],2);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],4);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],6);

    	p[a] = p[a] * s[a];

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],6);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],12);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],18);

       	p[a] = if_else(id_[a] == 0, vdouble3(0,0,0), vdouble3(3,2,1));

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0);
        
        TS_ASSERT_EQUALS(get<position>(particles[1])[0],3);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],2);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],1);
    }

    void helper_neighbours(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")

    	typedef Particles<std::tuple<scalar>> ParticlesType;
        typedef position_d<3> position;
       	ParticlesType particles;

        vdouble3 min(-1,-1,-1);
        vdouble3 max(1,1,1);
        vdouble3 periodic(true,true,true);
       	double diameter = 0.1;
        particles.init_neighbour_search(min,max,periodic);

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<alive> alive_;
        Symbol<scalar> s;
        Label<0,ParticlesType> a(particles);
        Label<1,ParticlesType> b(particles);
        auto dx = create_dx(a,b);
        AccumulateWithinDistance<std::plus<double> > sum(diameter);

       	particles.push_back(vdouble3(0,0,0));
       	particles.push_back(vdouble3(diameter*2,0,0));

       	s[a] = sum(b, 1);

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),1);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),1);

       	p[a] = if_else(id_[a] == 0, vdouble3(diameter/2.0,0,0), vdouble3(0,0,0));

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],diameter/2.0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],0);

    	p[a] = 0.5*(p[a] + s[a]);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0.5*(diameter/2.0 + 1));
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0.5);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0.5);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],0.5);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],0.5);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],0.5);

       	s[a] = sum(b, 1);

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),2);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),2);

       	p[a] = if_else(id_[a] == 0, vdouble3(0,0,0), vdouble3(diameter/2.0,diameter/2.0,diameter/2.0));
       	s[a] = if_else(id_[a] == 0, 0, 1);

       	TS_ASSERT_EQUALS(get<scalar>(particles[0]),0);
       	TS_ASSERT_EQUALS(get<scalar>(particles[1]),1);

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],diameter/2.0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],diameter/2.0);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],diameter/2.0);

        AccumulateWithinDistance<std::plus<vdouble3> > sumVect(diameter);

    	p[a] = sumVect(b, vdouble3(0,0,0) + 0.5*(s[a]/2.0 + s[b]/10.0));

    	TS_ASSERT_EQUALS(get<position>(particles[0])[0],0.05);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[1],0.05);
    	TS_ASSERT_EQUALS(get<position>(particles[0])[2],0.05);

    	TS_ASSERT_EQUALS(get<position>(particles[1])[0],0.55);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[1],0.55);
    	TS_ASSERT_EQUALS(get<position>(particles[1])[2],0.55);

        //
        // test inf norm range sum
        //
        AccumulateWithinDistance<std::plus<double>, -1> box_sum(diameter);
        get<position>(particles)[0] = vdouble3(0,0,0);
        get<position>(particles)[1] = 0.99*vdouble3(diameter,diameter,diameter);
        particles.update_positions();

        s[a] = sum(b, 1);

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),1);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),1);

        s[a] = box_sum(b, 1);

    	TS_ASSERT_EQUALS(get<scalar>(particles[0]),2);
    	TS_ASSERT_EQUALS(get<scalar>(particles[1]),2);
    }


    void helper_level0_expressions(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")

    	typedef Particles<std::tuple<scalar>> ParticlesType;
        typedef position_d<3> position;
       	ParticlesType particles;

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<alive> alive_;
        Symbol<scalar> s;
        Label<0,ParticlesType> a(particles);
        Label<1,ParticlesType> b(particles);
        auto dx = create_dx(a,b);

        Accumulate<std::plus<double> > sum;
        Accumulate<Aboria::max<double> > max;
        max.set_init(-1);
        Accumulate<Aboria::min<double> > min;
        min.set_init(1000);

       	particles.push_back(vdouble3(0,0,0));
       	particles.push_back(vdouble3(2,0,0));

        double result = eval(sum(a, if_else(p[a][0]<1, 1, 0)));
    	TS_ASSERT_EQUALS(result,1);
       	result = eval(sum(a, 1));
    	TS_ASSERT_EQUALS(result,2);
       	result = eval(sum(a, p[a][0]));
    	TS_ASSERT_EQUALS(result,2);
        int result2 = eval(max(a, id_[a]));
    	TS_ASSERT_EQUALS(result2,1);
        result2 = eval(min(a, id_[a]));
    	TS_ASSERT_EQUALS(result2,0);
       	particles.push_back(vdouble3(0,0,0));
        result2 = eval(max(a, id_[a]));
    	TS_ASSERT_EQUALS(result2,2);
    }

    void test_default() {
        helper_create_default_vectors();
        helper_create_double_vector();
        helper_transform();
        helper_neighbours();
        helper_level0_expressions();
    }

};

#endif /* SYMBOLICTEST_H_ */
