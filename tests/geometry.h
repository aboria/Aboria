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


#ifndef GEOMETRY_TEST_H_
#define GEOMETRY_TEST_H_

#include <cxxtest/TestSuite.h>

#include <random>

#include "Aboria.h"

using namespace Aboria;


class GeometryTest : public CxxTest::TestSuite {
public:

    typedef std::mt19937 generator_type;
    generator_type generator;
	void test_diffusion_reflection(void) {
		const double tol = 1e-9;
		Sphere sp_centre(double3(0,0,0),0.5,false);
		double3 p1(0,0,0);
		double3 p2;
		std::uniform_real_distribution<double> uni(0,0.1);
		for (int i = 0; i < 100; ++i) {
			p2 = p1 + double3(uni(generator),uni(generator),uni(generator));
			reflect_once(p1,p2,sp_centre);
			TS_ASSERT_EQUALS(sp_centre.is_in(p2),false);
			p1 = p2;

		}
	}
	void testReflectOnce(void) {
		const double tol = GEOMETRY_TOLERANCE;

		Sphere sp_centre(double3(0,0,0),0.5,false);
		double3 p1(0.4,0,0);
		double3 p2(0.6,0,0);
		TS_ASSERT_EQUALS(reflect_once(p1,p2,sp_centre),true);
		TS_ASSERT_DELTA(p2[0],0.4,tol);
		TS_ASSERT_DELTA(p2[1],0,tol);
		TS_ASSERT_DELTA(p2[2],0,tol);

		double3 p3(0.6,0,0);
		double3 p4(0.4,0,0);
		TS_ASSERT_EQUALS(reflect_once(p3,p4,sp_centre),false);

		double3 p5(0.6,0,0);
		double3 p6(0.4,100,0);
		TS_ASSERT_EQUALS(reflect_once(p5,p6,sp_centre),false);

		double3 p7(0,0.4,0);
		double3 p8(0,0.7,0);
		TS_ASSERT_EQUALS(reflect_once(p7,p8,sp_centre),true);
		TS_ASSERT_DELTA(p8[0],0,tol);
		TS_ASSERT_DELTA(p8[1],0.3,tol);
		TS_ASSERT_DELTA(p8[2],0,tol);
	}

    void testReflectOnceSymbolic(void) {
		const double tol = GEOMETRY_TOLERANCE;

        ABORIA_VARIABLE(radius,double,"radius")
    	typedef Particles<std::tuple<radius>> ParticlesType;
        typedef ParticlesType::position position;
       	ParticlesType particles;
        Symbol<position> p;
        Symbol<radius> r;
        Label<0,ParticlesType> a(particles);

       	particles.push_back(double3(0,0,0));

		Sphere sp_centre(double3(0,0,0),0.5,false);
		double3 p1(0.4,0,0);
		double3 p2(0.6,0,0);
        p[a] =  reflect_(p1,p2,sp_centre);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0.4,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);

		GeometrySymbolic<Sphere> sphere_;		

        p[a] = double3(0,-0.4,0);
        p[a] = double3(0,-0.3,0) | sphere_(p,0.5,false);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],-0.3,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);

        p[a] = double3(0,-0.4,0);
        p[a] = double3(0,-0.2,0) | sphere_(p[a],0.3,true);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);

        p[a] = double3(0,-0.4,0);
        r[a] = 0.5;
        p[a] = double3(0,-0.3,0) | sphere_(p,r,false);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],-0.3,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);

        p[a] = double3(0,-0.4,0);
        r[a] = 0.3;
        p[a] = double3(0,-0.2,0) | sphere_(p[a],r[a],true);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);


       	ParticlesType sphere_particles;
        Label<0,ParticlesType> b(sphere_particles);
        Label<1,ParticlesType> c(sphere_particles);
        
       	sphere_particles.push_back(double3(0,0,0));
        r[b] = 0.3;

    	sphere_particles.init_neighbour_search(double3(-4,-4,-4),double3(4,4,4),0.3,bool3(false,false,false));

        p[a] = double3(0,0.4,0);

		GeometriesSymbolic<Sphere> spheres_;		

        p[a] = double3(0,-0.2,0) | spheres_(c,r[b]);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);


	}

	void testSphereIsIn(void) {
		Sphere sp_centre(double3(0,0,0),0.5,true);
		Sphere sp_off_centre(double3(1.0,0,0),0.5,false);
		TS_ASSERT_EQUALS(sp_centre.is_in(double3(0,0,0)),true);
		TS_ASSERT_EQUALS(sp_centre.is_in(double3(0.4,0,0)),true);
		TS_ASSERT_EQUALS(sp_centre.is_in(double3(0.6,0,0)),false);
		TS_ASSERT_EQUALS(sp_off_centre.is_in(double3(1.0,0,0)),false);
		TS_ASSERT_EQUALS(sp_off_centre.is_in(double3(1.4,0,0)),false);
		TS_ASSERT_EQUALS(sp_off_centre.is_in(double3(1.6,0,0)),true);
	}
	void testSphereLineXSurface(void) {
		Sphere sp_centre(double3(0,0,0),0.5,true);
		Sphere sp_off_centre(double3(1.0,0,0),0.5,false);
		const double tol = GEOMETRY_TOLERANCE;
		TS_ASSERT_EQUALS(sp_centre.lineXsurface(double3(0,0,0),double3(0,0,0.1)).first,-1)

#define CHECK_RESULT(a,b,c,d) \
		TS_ASSERT_DELTA(result.first,a,tol); \
		TS_ASSERT_DELTA(result.second[0],b,tol); \
		TS_ASSERT_DELTA(result.second[1],c,tol); \
		TS_ASSERT_DELTA(result.second[2],d,tol); \

		std::pair<double,double3> result = sp_centre.lineXsurface(double3(0,0,0),double3(0,0,1.0));
		CHECK_RESULT(0.5,0,0,1)

		result = sp_centre.lineXsurface(double3(0.6,0,0),double3(0.4,0,0));
		CHECK_RESULT(0.1,1,0,0)

		result = sp_centre.lineXsurface(double3(0.6,0,0),double3(-0.6,0,0));
		CHECK_RESULT(0.1,1,0,0)

		result = sp_off_centre.lineXsurface(double3(1,0,0),double3(1,0,1));
		CHECK_RESULT(0.5,0,0,-1)

		result = sp_off_centre.lineXsurface(double3(1.6,0,0),double3(1.4,0,0));
		CHECK_RESULT(0.1,-1,0,0)

		result = sp_centre.lineXsurface(double3(0.6,0,0),double3(-0.6,0,0));
		CHECK_RESULT(0.1,1,0,0)

		result = sp_centre.lineXsurface(double3(-0.6,0,0),double3(0.6,0,0));
		CHECK_RESULT(0.1,-1,0,0)
#undef CHECK_RESULT
	}

};


#endif /* GEOMETRY_H_ */
