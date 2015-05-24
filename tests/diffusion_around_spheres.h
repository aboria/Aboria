/*
 * diffusion_around_spheres.h
 *
 *  Created on: 6 Jul 2015
 *      Author: robinsonm
 */

#ifndef DIFFUSION_AROUND_SPHERES_H_
#define DIFFUSION_AROUND_SPHERES_H_


#include <cxxtest/TestSuite.h>

#include <random>
typedef std::mt19937 generator_type;
generator_type generator;

#define LOG_LEVEL 1
#include "Aboria.h"

using namespace Aboria;


class DiffusionAroundSpheres : public CxxTest::TestSuite {
public:
	void test_diffusion_around_spheres(void) {
		const double tol = GEOMETRY_TOLERANCE;

		ABORIA_VARIABLE(radius,double,"radius")

		Particles<radius> spheres;

		const double L = 10.0;
		const double D = 1.0;
		const double dt = 0.1;
		const double timesteps = 1000;

		spheres.push_back(Vect3d(0,0,0));
		spheres[0].set<radius>(1.0);
		spheres.push_back(Vect3d(5,0,0));
		spheres[1].set<radius>(2.0);
		spheres.push_back(Vect3d(0,-5,0));
		spheres[2].set<radius>(1.5);
		spheres.push_back(Vect3d(0,0,5));
		spheres[3].set<radius>(1.0);

    	spheres.init_neighbour_search(Vect3d(-L,-L,-L),Vect3d(L,L,L),4,Vect3b(true,true,true));

		Particles<> points;
		std::uniform_real_distribution<double> uni(-L,L);
		for (int i = 0; i < 1000; ++i) {
			points.push_back(Vect3d(uni(generator),uni(generator),uni(generator)));
		}

		auto spheres_position = get_vector<position>(spheres);
		auto spheres_radius = get_vector<radius>(spheres);

		auto points_position = get_vector<position>(points);
		auto points_alive = get_vector<alive>(points);


		Label<0> a;
		Label<1> b;
		Dx dx;
		Normal N;
		GeometriesSymbolic<Sphere> spheres_;		
		VectorSymbolic<double> vector;		

		/*
		 * Kill any points within spheres
		 */
		points_alive = first_(b=spheres, norm_(dx) < spheres_radius[b],false,true);

		/*
		 * Check no points within spheres
		 */
		for(auto i: points) {
			TS_ASSERT_RELATION(std::greater<double>, norm(get<position>(i) - Vect3d(0,0,0)), 1.0);
			TS_ASSERT_RELATION(std::greater<double>, norm(get<position>(i) - Vect3d(5,0,0)), 2.0);
			TS_ASSERT_RELATION(std::greater<double>, norm(get<position>(i) - Vect3d(0,-5,0)), 1.5);
			TS_ASSERT_RELATION(std::greater<double>, norm(get<position>(i) - Vect3d(0,0,5)), 1.0);
		}

		/*
		 * Diffusion step for points and reflect off spheres
		 */
		for (int i = 0; i < timesteps; ++i) {
			points_position += std::sqrt(2*D*dt)*vector(N,N,N) | spheres_(b=spheres,spheres_radius[b]);
		}

		/*
		 * Check still no points within spheres
		 */
		for(auto i: points) {
			TS_ASSERT_RELATION(std::greater<double>, norm(get<position>(i) - Vect3d(0,0,0)), 1.0);
			TS_ASSERT_RELATION(std::greater<double>, norm(get<position>(i) - Vect3d(5,0,0)), 2.0);
			TS_ASSERT_RELATION(std::greater<double>, norm(get<position>(i) - Vect3d(0,-5,0)), 1.5);
			TS_ASSERT_RELATION(std::greater<double>, norm(get<position>(i) - Vect3d(0,0,5)), 1.0);
		}
	}
};

#endif /* DIFFUSION_AROUND_SPHERES_H_ */
