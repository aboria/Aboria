/*
 * md.h
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */
#ifndef MD_TEST_H_
#define MD_TEST_H_

#include <cxxtest/TestSuite.h>

#include <boost/math/constants/constants.hpp>
#include <math.h>
#include <random>

#include "Aboria.h"
using namespace Aboria;

class MDTest : public CxxTest::TestSuite {
public:
  typedef std::mt19937 generator_type;
  generator_type generator;

  template <template <typename> class SearchMethod>
  void helper_md_iterator(void) {
    //[md_iterator
    /*`
    Simple Molecular Dynamics of $N$ particles interacting via a linear spring
    potential.

    This example creates $N$ particles within a two-dimensional square domain,
    with periodic boundary conditions.

    There is a linear spring force
    $\mathbf{f}\_{ij}$ between particles $i$ and $j$ with a
    rest separation of $r$ (constant for all particles), and a cutoff at $r$.
    That is, if
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


    We wish to use a leap frog integrator to evolve positions $\mathbf{r}\_i$
    using velocities $\mathbf{v}\_i$ and accelerations $\mathbf{a}\_i = \sum_j
    \mathbf{f}\_{ij}$. This gives the following update equations for
    each timestep $n$

    \begin{align*}
    \mathbf{v}^{n+1}\_i &= \mathbf{v}^n_i + \frac{dt}{m_i} \sum_j
    \mathbf{f}^n_{ij}
    \\\\
    \mathbf{r}^{n+1}_i &= \mathbf{r}^n_i + dt\, \mathbf{v}^{n+1}_i.
    \end{align*}

    We implement this in Aboria using the code given below. Firstly we create
    the particle set data structure and add particles, ensuring that we have an
    initial condition where all the spring forces are $\mathbf{f}\_{ij}=0$. Then
    we start the timestep loop, using our update equations given above.
    */
    //<-
    // clang-format off
    //->
//=#include <boost/math/constants/constants.hpp>
//=#include <math.h>
//=#include <random>
//=
//=#include "Aboria.h"
//=using namespace Aboria;
//=
//=int main() {
    //<-
    // clang-format on
    //->
    const double PI = boost::math::constants::pi<double>();

    // Create a 2d particle container with one additional variable
    // "velocity", represented by a 2d double vector
    ABORIA_VARIABLE(velocity, vdouble2, "velocity")
    //<-
    typedef Particles<std::tuple<velocity>, 2, std::vector, SearchMethod>
        container_type;
    //->
    //=        typedef Particles<std::tuple<velocity>,2> container_type;

    typedef typename container_type::position position;
    container_type particles;

    // set parameters for the MD simulation
    const int timesteps = 3000;
    const int nout = 200;
    const int timesteps_per_out = timesteps / nout;
    const double L = 31.0 / 1000.0;
    const int N = 30;
    const double diameter = 0.0022;
    const double k = 1.0e01;
    const double dens = 1160.0;
    const double mass = PI * std::pow(0.5 * diameter, 2) * dens;
    const double reduced_mass = 0.5 * mass;
    const double dt = (1.0 / 50.0) * PI / std::sqrt(k / reduced_mass);
    const double v0 = L / (timesteps * dt);

    // initiate neighbour search on a periodic 2d domain of side length L
    // set average number of particles per cell to 1
    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(L, L),
                                    vbool2(true, true));

    // create N particles, ensuring that they do not overlap, according
    // to the set diameter. set the initial velocity in a random direction
    std::uniform_real_distribution<double> uni(0, 1);
    for (int i = 0; i < N; ++i) {
      bool free_position = false;

      // create new particle
      typename container_type::value_type p;

      // set a random direction, and initialise velocity
      const double theta = uni(generator) * 2 * PI;
      get<velocity>(p) = v0 * vdouble2(cos(theta), sin(theta));

      // randomly choose positions within the domain until one is
      // found with no other particles within a range equal to diameter
      while (free_position == false) {
        get<position>(p) = vdouble2(uni(generator) * L, uni(generator) * L);
        free_position = true;

        // loop over all neighbouring particles within a euclidean distance
        // of size "diameter"
        for (auto j = euclidean_search(particles.get_query(), get<position>(p),
                                       diameter);
             j != false; ++j) {
          //<-
          (void)j;
          //->
          free_position = false;
          break;
        }
      }
      particles.push_back(p);
    }

    // perform MD timestepping
    for (int io = 0; io < nout; ++io) {

      // on every i/o step write particle container to a vtk
      // unstructured grid file
      std::cout << "." << std::flush;
#ifdef HAVE_VTK
      vtkWriteGrid("particles", io, particles.get_grid(true));
#endif
      for (int t = 0; t < timesteps_per_out; t++) {

        // loop through particles and calculate velocity update
        for (auto i : particles) {
          auto sum = vdouble2::Constant(0);

          // use a range search with radius `diameter` to
          // accelerate velocity update
          for (auto j = euclidean_search(particles.get_query(),
                                         get<position>(i), diameter);
               j != false; ++j) {
            const double r = j.dx().norm();
            if (r != 0) {
              sum += -k * (diameter / r - 1.0) * j.dx();
            }
          }
          get<velocity>(i) += dt * sum / mass;
        }

        // now loop through particles and use velocity
        // to update the particle positions
        for (auto i : particles) {
          get<position>(i) += dt * get<velocity>(i);
        }

        // since we have changed the particle positions, we need
        // to update the spatial data structures
        particles.update_positions();
      }
    }
    std::cout << std::endl;
  }
  //]

  template <template <typename> class SearchMethod>
  void helper_md_symbolic(void) {
    //[md_symbolic
    /*`
    Simple Molecular Dynamics of $N$ particles interacting via a linear spring
    potential.

    This example creates $N$ particles within a two-dimensional square domain,
    with periodic boundary conditions.

    There is a linear spring force
    $\mathbf{f}\_{ij}$ between particles $i$ and $j$ with a
    rest separation of $r$ (constant for all particles), and a cutoff at $r$.
    That is, if
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


    We wish to use a leap frog integrator to evolve positions $\mathbf{r}\_i$
    using velocities $\mathbf{v}\_i$ and accelerations $\mathbf{a}\_i = \sum_j
    \mathbf{f}\_{ij}$. This gives the following update equations for
    each timestep $n$

    \begin{align*}
    \mathbf{v}^{n+1}\_i &= \mathbf{v}^n_i + \frac{dt}{m_i} \sum_j
    \mathbf{f}^n_{ij}
    \\\\
    \mathbf{r}^{n+1}_i &= \mathbf{r}^n_i + dt\, \mathbf{v}^{n+1}_i.
    \end{align*}

    We implement this in Aboria using the code given below. Firstly we create
    the particle set data structure and add particles, ensuring that we have an
    initial condition where all the spring forces are $\mathbf{f}\_{ij}=0$. Then
    we start the timestep loop, using our update equations given above.
    */
    //<-
    // clang-format off
    //->
//=#include <boost/math/constants/constants.hpp>
//=#include <math.h>
//=#include <random>
//=
//=#include "Aboria.h"
//=using namespace Aboria;
//=
//=int main() {
    //<-
    // clang-format on
    //->

    const double PI = boost::math::constants::pi<double>();

    // Create a 2d particle container with one additional variable
    // "velocity", represented by a 2d double vector
    ABORIA_VARIABLE(velocity, vdouble2, "velocity")
    //<-
    typedef Particles<std::tuple<velocity>, 2, std::vector, SearchMethod>
        container_type;
    //->
    //=        typedef Particles<std::tuple<velocity>,2> container_type;

    typedef typename container_type::position position;

    // set parameters for the MD simulation
    const int timesteps = 3000;
    const int nout = 200;
    const int timesteps_per_out = timesteps / nout;
    const double L = 31.0 / 1000.0;
    const int N = 30;
    const double diameter = 0.0022;
    const double k = 1.0e01;
    const double dens = 1160.0;
    const double mass = PI * std::pow(0.5 * diameter, 2) * dens;
    const double reduced_mass = 0.5 * mass;
    const double dt = (1.0 / 50.0) * PI / std::sqrt(k / reduced_mass);

    // create N particles
    container_type particles(N);

    // create symbols and labels in order to use the expression API
    Symbol<position> p;
    Symbol<velocity> v;
    Symbol<id> id_;
    Uniform uniform;
    VectorSymbolic<double, 2> vector;
    Label<0, container_type> a(particles);
    Label<1, container_type> b(particles);

    // dx is a symbol representing the difference in positions of
    // particle a and b.
    auto dx = create_dx(a, b);

    // sum is a symbolic function that sums a sequence of 2d vectors
    AccumulateWithinDistance<std::plus<vdouble2>> sum(diameter);

    // initialise particles to a uniform distribution
    p[a] = L * vector(uniform[a], uniform[a]);

    // zero initial velocity
    v[a] = vector(0, 0);

    // initiate neighbour search on a periodic 2d domain of side length
    // L set average number of particles per cell to 1
    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(L, L),
                                    vbool2(true, true));

    // perform MD timestepping
    for (int io = 0; io < nout; ++io) {

      // on every i/o step write particle container to a vtk
      // unstructured grid file
      std::cout << "." << std::flush;
#ifdef HAVE_VTK
      vtkWriteGrid("particles", io, particles.get_grid(true));
#endif
      for (int i = 0; i < timesteps_per_out; i++) {
        // leap frog integrator
        v[a] += dt *
                // spring force between particles
                sum(b, if_else(id_[a] != id_[b],
                               -k * (diameter / norm(dx) - 1.0), 0.0) *
                           dx) /
                mass;

        // update positions
        p[a] += dt * v[a];
      }
    }
    std::cout << std::endl;
  }
  //]

  void test_CellListOrdered() {
    helper_md_iterator<CellListOrdered>();
    helper_md_symbolic<CellListOrdered>();
  }

  void test_CellList() {
    helper_md_iterator<CellList>();
    helper_md_symbolic<CellList>();
  }

  void test_Kdtree() {
    helper_md_iterator<Kdtree>();
    helper_md_symbolic<Kdtree>();
  }

  void test_KdtreeNanoflann() {
    helper_md_iterator<KdtreeNanoflann>();
    helper_md_symbolic<KdtreeNanoflann>();
  }

  void test_HyperOctree() {
    helper_md_iterator<HyperOctree>();
    helper_md_symbolic<HyperOctree>();
  }
};

#endif /* MD_TEST_H_ */
