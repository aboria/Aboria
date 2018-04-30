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
/*`
This example creates $N$ particles within a two-dimensional square domain,
with periodic boundary conditions.

There is a linear spring force
$\mathbf{f}\_{ij}$ between particles $i$ and $j$ with a
rest separation of $r$ (constant for all particles), and a cutoff at $r$. That
is, if
$\mathbf{r}\_i$ is the position of particle $i$ and
$\mathbf{dx}\_{ij}=\mathbf{r}\_j-\mathbf{r}\_j$, then

$$
\mathbf{f}\_{ij} = \begin{cases}
            \frac{r-|\mathbf{dx}\_{ij}|}{|\mathbf{dx}\_{ij}|}\mathbf{dx}\_{ij},
& \text{for }
              |\mathbf{dx}\_{ij}|<r \\\\
            0 & \text{otherwise}.
            \end{cases}
$$


We wish to use a leap frog integrator to evolve positions $\mathbf{r}\_i$ using
velocities $\mathbf{v}\_i$ and accelerations $\mathbf{a}\_i = \sum_j
\mathbf{f}\_{ij}$. This gives the following update equations for
each timestep $n$

\begin{align*}
\mathbf{v}^{n+1}\_i &= \mathbf{v}^n_i + \frac{dt}{m_i} \sum_j \mathbf{f}^n_{ij}
\\\\
\mathbf{r}^{n+1}_i &= \mathbf{r}^n_i + dt\, \mathbf{v}^{n+1}_i.
\end{align*}

We implement this in Aboria using the code given below. Firstly we create the
particle set data structure and add particles, ensuring that we have an initial
condition where all the spring forces are $\mathbf{f}\_{ij}=0$. Then we start
the timestep loop, using our update equations given above.
*/

#include <random>

#include "Aboria.h"
using namespace Aboria;

#include <boost/math/constants/constants.hpp>
#include <math.h>

//<-
class MDLevel1Test : public CxxTest::TestSuite {
public:
  ABORIA_VARIABLE(velocity, vdouble2, "velocity")
  ABORIA_VARIABLE(verlet, std::vector<int>, "verlet-list")

  template <template <typename, typename> class Vector,
            template <typename> class SearchMethod>
  void helper_md(void) {
    //->
    //=int main() {
    const double PI = boost::math::constants::pi<double>();

    /*
     * Create a 2d particle container with one additional variable
     * "velocity", represented by a 2d double vector
     */
    //<-
    typedef Particles<std::tuple<velocity, verlet>, 2, Vector, SearchMethod>
        container_type;
    //->
    //=        typedef Particles<std::tuple<velocity>,2> container_type;

    typedef typename container_type::position position;
    typedef typename container_type::query_type query_type;
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::raw_reference raw_reference;
    typedef typename container_type::raw_const_reference raw_const_reference;

    /*
     * set parameters for the MD simulation
     */
    const int timesteps = 300;
    const int nout = 200;
    const int timesteps_per_out = 10;
    const double L = 31.0 / 1000.0;
    const int N = 2000;
    const double diameter = 0.0022;
    const double diameter2 = std::pow(diameter, 2);
    const double k = 1.0e01;
    const double dens = 1160.0;
    const double mass = (1.0 / 6.0) * PI * std::pow(0.5 * diameter, 3) * dens;
    const double reduced_mass = 0.5 * mass;
    const double dt = (1.0 / 50.0) * PI / sqrt(k / reduced_mass);
    const double max_velocity = 30.0;
    const double max_velocity2 = std::pow(30.0, 2);
    const double max_step = max_velocity * dt;
    const int n_reuse_verlet_list = 5;
    const double buffer = 2 * n_reuse_verlet_list * max_step;
    std::cout << "setting search radius/L = " << (diameter + buffer) / L
              << std::endl;

    /*
     * construct N particles
     */
    container_type particles(N);

    /*
     * randomly set position for N particles
     */
    detail::for_each(particles.begin(), particles.end(),
                     [] CUDA_HOST_DEVICE(raw_reference i) {
#if defined(__CUDACC__)
                       thrust::uniform_real_distribution<double> uni(0, 1);
#else
                       std::uniform_real_distribution<double> uni(0, 1);
#endif
                       generator_type &gen = get<generator>(i);
                       get<position>(i) = vdouble2(uni(gen), uni(gen));
                       get<velocity>(i) = vdouble2(0, 0);
                     });

    /*
     * initiate neighbour search on a periodic 2d domain of side length L
     * set average number of particles per cell to 1
     */
    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(L, L),
                                    vbool2(true, true));

    /*
     * perform MD timestepping
     */
    for (int i = 0; i < timesteps / n_reuse_verlet_list; ++i) {
      // build verlet list
      for (auto i : particles) {
        std::vector<int> &verlet_list = get<verlet>(i);
        verlet_list.clear();
        for (auto j = euclidean_search(particles.get_query(), get<position>(i),
                                       diameter + buffer);
             j != false; ++j) {
          // hack to get index: look at id pointers
          const int j_index = &get<id>(*j) - &get<id>(particles[0]);
          verlet_list.push_back(j_index);
        }
      }
      // timesteps
      for (int vi = 0; vi < n_reuse_verlet_list; ++vi) {
        for (auto i : particles) {
          vdouble2 force_sum(0, 0);
          for (int j_index : get<verlet>(i)) {
            const_reference j = particles[j_index];
            const vdouble2 dx = particles.correct_dx_for_periodicity(
                get<position>(j) - get<position>(i));
            const double r2 = dx.squaredNorm();
            if (r2 > 0 && r2 < diameter2) {
              force_sum -= k * (diameter / std::sqrt(r2) - 1) * dx;
            }
          }
          get<velocity>(i) += dt * force_sum / mass;
          const double v2 = get<velocity>(i).squaredNorm();
          if (v2 > max_velocity2) {
            get<velocity>(i) *= max_velocity / std::sqrt(v2);
          }
        }
        for (auto i : particles) {
          get<position>(i) += dt * get<velocity>(i);
        }
        std::cout << "." << std::flush;
#ifdef HAVE_VTK
        vtkWriteGrid("particles", io, particles.get_grid(true));
#endif
      }
      particles.update_positions();
    }
  }
  std::cout << std::endl;
}
//]

    void test_std_vector_CellList(void) {
  helper_md<std::vector, CellList>();
}

void test_std_vector_CellListOrdered(void) {
  helper_md<std::vector, CellListOrdered>();
}

void test_std_vector_HyperOctree(void) { helper_md<std::vector, HyperOctree>(); }
}
;

#endif /* MD_TEST_H_ */
