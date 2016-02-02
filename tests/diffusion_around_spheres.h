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
		const double timesteps = 100;

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

        Symbol<position> p;
        Symbol<id> id_;
        Symbol<alive> alive_;
        Symbol<radius> r;
        Label<0,Particles<radius> > a_s(spheres);
        Label<1,Particles<radius> > b_s(spheres);
        Label<0,Particles<> > a_p(points);
        Label<1,Particles<> > b_p(points);

		Dx dx;
		Normal N;
		GeometriesSymbolic<Sphere> spheres_;		
		VectorSymbolic<double> vector;		
        Accumulate<std::bit_or<bool> > any;

		/*
		 * Kill any points within spheres
		 */
		alive_[a_p] = !any(b_s, norm(dx) < r[b_s],true);

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
			p[a_p] += std::sqrt(2*D*dt)*vector(N,N,N) | spheres_(b_s,r[b_s]);
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
