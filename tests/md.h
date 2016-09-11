/*
 * md.h
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */
#ifndef MD_TEST_H_
#define MD_TEST_H_

#include <cxxtest/TestSuite.h>

//[md
#include <random>

#include "Aboria.h"
using namespace Aboria;

#include <boost/math/constants/constants.hpp>
#include <math.h>

//<-
class MDTest : public CxxTest::TestSuite {
public:

    typedef std::mt19937 generator_type;
    generator_type generator;

    void test_md(void) {
//->
//=int main() {
        const double PI = boost::math::constants::pi<double>();

        /*
         * Create a 2d particle container with one additional variable
         * "velocity", represented by a 2d double vector
         */
        ABORIA_VARIABLE(velocity,double2,"velocity")
        typedef Particles<std::tuple<velocity>,2> container_type;
        typedef container_type::position position;
        container_type particles;

        /*
         * set parameters for the MD simulation
         */
        const int timesteps = 3000;
        const int nout = 200;
        const int timesteps_per_out = timesteps/nout;
        const double L = 31.0/1000.0;
        const int N = 30;
        const double diameter = 0.0022;
        const double k = 1.0e01;
        const double dens = 1160.0;
        //const double mass = PI*std::pow(0.5*diameter,2)*dens;
        const double mass = (1.0/6.0)*PI*std::pow(0.5*diameter,3)*dens;
        const double reduced_mass = 0.5*mass;
        const double dt = (1.0/50.0)*PI/sqrt(k/reduced_mass);
        const double v0 = L/(timesteps*dt);

        /*
         * initiate neighbour search on a periodic 2d domain of side length L
         */
        particles.init_neighbour_search(double2(0,0),double2(L,L),diameter,bool2(true,true));

        /*
         * create N particles, ensuring that they do not overlap, according 
         * to the set diameter. set the initial velocity in a random direction
         */
        std::uniform_real_distribution<double> uni(0,1);
        for (int i = 0; i < N; ++i) {
            bool free_position = false;

            /*
             * create new particle
             */
            container_type::value_type p;

            /*
             * set a random direction, and initialise velocity
             */
            const double theta = uni(generator)*2*PI;
            get<velocity>(p) = v0*double2(cos(theta),sin(theta));

            /*
             * randomly choose positions within the domain until one is 
             * found with no other particles within a range equal to diameter
             */
            while (free_position == false) {
                get<position>(p) = double2(uni(generator)*L,uni(generator)*L);
                free_position = true;

                /*
                 * loop over all neighbouring particles within a square with
                 * side length "diameter" (see init_neighbour_search call above)
                 */
                for (auto tpl: particles.get_neighbours(get<position>(p))) {

                    /*
                     * tpl variable is a tuple containing:
                     *  (0) -> neighbouring particle value_type
                     *  (1) -> relative position of neighbouring particle
                     *         from query point
                     */
                    const double2& dx = std::get<1>(tpl);
                    const container_type::value_type& j = std::get<0>(tpl);
                    if (dx.norm() < diameter) {
                        free_position = false;
                        break;
                    }
                }
            }
            particles.push_back(p);
        }

        /*
         * create symbols and labels in order to use the expression API
         */
        Symbol<position> p;
        Symbol<velocity> v;
        Symbol<id> id_;
        Label<0,container_type> a(particles);
        Label<1,container_type> b(particles);

        /*
         * dx is a symbol representing the difference in positions of 
         * particle a and b.
         */
        auto dx = create_dx(a,b);

        /*
         * sum is a symbolic function that sums a sequence of 2d vectors
         */
        Accumulate<std::plus<double2> > sum;
        
        /*
         * perform MD timestepping
         */
        for (int io = 0; io < nout; ++io) {

            /*
             * on every i/o step write particle container to a vtk
             * unstructured grid file
             */
            std::cout << "." << std::flush;
#ifdef HAVE_VTK
            vtkWriteGrid("particles",io,particles.get_grid(true));
#endif
            for (int i = 0; i < timesteps_per_out; i++) {
                /*
                 * leap frog integrator
                 */
                v[a] += dt * ( 
                        // spring force between particles
                        sum(b, id_[a]!=id_[b] && norm(dx)<diameter, 
                              -k*(diameter/norm(dx)-1)*dx)
                        /mass
                        );
                p[a] += dt*v[a];
            }
        }
        std::cout << std::endl;
    }
//]

};

#endif /* MD_TEST_H_ */

