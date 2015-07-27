/*
 * geometry.h
 *
 *  Created on: 23 Jan 2015
 *      Author: mrobins
 */

#ifndef GEOMETRY_TEST_H_
#define GEOMETRY_TEST_H_

#include <cxxtest/TestSuite.h>

#include <random>
typedef std::mt19937 generator_type;
generator_type generator;

#define LOG_LEVEL 1
#include "Aboria.h"

using namespace Aboria;


class GeometryTest : public CxxTest::TestSuite {
public:
	void test_diffusion_reflection(void) {
		const double tol = 1e-9;
		Sphere sp_centre(Vect3d(0,0,0),0.5,false);
		Vect3d p1(0,0,0);
		Vect3d p2;
		std::uniform_real_distribution<double> uni(0,0.1);
		for (int i = 0; i < 100; ++i) {
			p2 = p1 + Vect3d(uni(generator),uni(generator),uni(generator));
			reflect_once(p1,p2,sp_centre);
			TS_ASSERT_EQUALS(sp_centre.is_in(p2),false);
			p1 = p2;

		}
	}
	void testReflectOnce(void) {
		const double tol = GEOMETRY_TOLERANCE;

		Sphere sp_centre(Vect3d(0,0,0),0.5,false);
		Vect3d p1(0.4,0,0);
		Vect3d p2(0.6,0,0);
		TS_ASSERT_EQUALS(reflect_once(p1,p2,sp_centre),true);
		TS_ASSERT_DELTA(p2[0],0.4,tol);
		TS_ASSERT_DELTA(p2[1],0,tol);
		TS_ASSERT_DELTA(p2[2],0,tol);

		Vect3d p3(0.6,0,0);
		Vect3d p4(0.4,0,0);
		TS_ASSERT_EQUALS(reflect_once(p3,p4,sp_centre),false);

		Vect3d p5(0.6,0,0);
		Vect3d p6(0.4,100,0);
		TS_ASSERT_EQUALS(reflect_once(p5,p6,sp_centre),false);

		Vect3d p7(0,0.4,0);
		Vect3d p8(0,0.7,0);
		TS_ASSERT_EQUALS(reflect_once(p7,p8,sp_centre),true);
		TS_ASSERT_DELTA(p8[0],0,tol);
		TS_ASSERT_DELTA(p8[1],0.3,tol);
		TS_ASSERT_DELTA(p8[2],0,tol);
	}

    void testReflectOnceSymbolic(void) {
		const double tol = GEOMETRY_TOLERANCE;

        ABORIA_VARIABLE(radius,double,"radius")
    	typedef Particles<radius> ParticlesType;
       	ParticlesType particles;
       	auto position_ = get_vector<position>(particles);
       	auto radius_ = get_vector<radius>(particles);
       	particles.push_back(Vect3d(0,0,0));


		Sphere sp_centre(Vect3d(0,0,0),0.5,false);
		Vect3d p1(0.4,0,0);
		Vect3d p2(0.6,0,0);
        position_ =  reflect_(p1,p2,sp_centre);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0.4,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);

		GeometrySymbolic<Sphere> sphere_;		

        position_ = Vect3d(0,-0.4,0);
        position_ = Vect3d(0,-0.3,0) | sphere_(position_,0.5,false);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],-0.3,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);

        position_ = Vect3d(0,-0.4,0);
        position_ = Vect3d(0,-0.2,0) | sphere_(position_,0.3,true);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);

        position_ = Vect3d(0,-0.4,0);
        radius_ = 0.5;
        position_ = Vect3d(0,-0.3,0) | sphere_(position_,radius_,false);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],-0.3,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);

        position_ = Vect3d(0,-0.4,0);
        radius_ = 0.3;
        position_ = Vect3d(0,-0.2,0) | sphere_(position_,radius_,true);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);


       	ParticlesType sphere_particles;
       	auto sphere_position_ = get_vector<position>(sphere_particles);
       	auto sphere_radius_ = get_vector<radius>(sphere_particles);
        
       	sphere_particles.push_back(Vect3d(0,0,0));
        sphere_radius_ = 0.3;

    	sphere_particles.init_neighbour_search(Vect3d(-4,-4,-4),Vect3d(4,4,4),0.3,Vect3b(false,false,false));

        position_ = Vect3d(0,0.4,0);

		GeometriesSymbolic<Sphere> spheres_;		
       	Label<1> b;

        position_ = Vect3d(0,-0.2,0) | spheres_(b=sphere_particles,sphere_radius_);
		TS_ASSERT_DELTA(get<position>(particles[0])[0],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[1],0,tol);
		TS_ASSERT_DELTA(get<position>(particles[0])[2],0,tol);


	}

	void testSphereIsIn(void) {
		Sphere sp_centre(Vect3d(0,0,0),0.5,true);
		Sphere sp_off_centre(Vect3d(1.0,0,0),0.5,false);
		TS_ASSERT_EQUALS(sp_centre.is_in(Vect3d(0,0,0)),true);
		TS_ASSERT_EQUALS(sp_centre.is_in(Vect3d(0.4,0,0)),true);
		TS_ASSERT_EQUALS(sp_centre.is_in(Vect3d(0.6,0,0)),false);
		TS_ASSERT_EQUALS(sp_off_centre.is_in(Vect3d(1.0,0,0)),false);
		TS_ASSERT_EQUALS(sp_off_centre.is_in(Vect3d(1.4,0,0)),false);
		TS_ASSERT_EQUALS(sp_off_centre.is_in(Vect3d(1.6,0,0)),true);
	}
	void testSphereLineXSurface(void) {
		Sphere sp_centre(Vect3d(0,0,0),0.5,true);
		Sphere sp_off_centre(Vect3d(1.0,0,0),0.5,false);
		const double tol = GEOMETRY_TOLERANCE;
		TS_ASSERT_EQUALS(sp_centre.lineXsurface(Vect3d(0,0,0),Vect3d(0,0,0.1)).first,-1)

#define CHECK_RESULT(a,b,c,d) \
		TS_ASSERT_DELTA(result.first,a,tol); \
		TS_ASSERT_DELTA(result.second[0],b,tol); \
		TS_ASSERT_DELTA(result.second[1],c,tol); \
		TS_ASSERT_DELTA(result.second[2],d,tol); \

		std::pair<double,Vect3d> result = sp_centre.lineXsurface(Vect3d(0,0,0),Vect3d(0,0,1.0));
		CHECK_RESULT(0.5,0,0,1)

		result = sp_centre.lineXsurface(Vect3d(0.6,0,0),Vect3d(0.4,0,0));
		CHECK_RESULT(0.1,1,0,0)

		result = sp_centre.lineXsurface(Vect3d(0.6,0,0),Vect3d(-0.6,0,0));
		CHECK_RESULT(0.1,1,0,0)

		result = sp_off_centre.lineXsurface(Vect3d(1,0,0),Vect3d(1,0,1));
		CHECK_RESULT(0.5,0,0,-1)

		result = sp_off_centre.lineXsurface(Vect3d(1.6,0,0),Vect3d(1.4,0,0));
		CHECK_RESULT(0.1,-1,0,0)

		result = sp_centre.lineXsurface(Vect3d(0.6,0,0),Vect3d(-0.6,0,0));
		CHECK_RESULT(0.1,1,0,0)

		result = sp_centre.lineXsurface(Vect3d(-0.6,0,0),Vect3d(0.6,0,0));
		CHECK_RESULT(0.1,-1,0,0)
#undef CHECK_RESULT
	}

};


#endif /* GEOMETRY_H_ */
