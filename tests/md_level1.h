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

  /*
template <typename Query>
struct set_position_lambda {
  typedef typename Query::traits_type::raw_reference reference;
  typedef typename Query::traits_type::position position;

  CUDA_HOST_DEVICE
  void operator()(reference i) const {
#if defined(__CUDACC__)
      thrust::uniform_real_distribution<double> uni(0,1);
#else
      std::uniform_real_distribution<double> uni(0,1);
#endif
      generator_type& gen = get<generator>(i);
      get<position>(i) = vdouble2(uni(gen),uni(gen));
  }
};

template <typename Query>
struct timestep_lambda {
  typedef typename Query::traits_type::raw_reference reference;
  typedef typename Query::traits_type::raw_const_reference const_reference;
  typedef typename Query::traits_type::position position;

  Query query;
  double diameter;
  double k;
  double dt;
  double mass;

  timestep_lambda(const Query &query, double diameter, double k, double dt,
double mass): query(query),diameter(diameter),k(k),dt(dt),mass(mass) {}

  CUDA_HOST_DEVICE
  void operator()(reference i) const {
      vdouble2 force_sum(0,0);
      std::cout << "id = "<<get<id>(i)<<" old position = "<<get<position>(i)<<
std::endl; for (auto tpl: euclidean_search(query,get<position>(i),diameter)) {
          const vdouble2& dx = detail::get_impl<1>(tpl);
          const_reference j = detail::get_impl<0>(tpl);
          const double r = dx.norm();
          std::cout << "id = "<<get<id>(i)<<" dx = "<<dx<< " r =
"<<r<<std::endl; if (r != 0) { force_sum -= k*(diameter/r-1)*dx;
          }
      }
      get<velocity>(i) += dt*force_sum/mass;
      get<position>(i) += dt*get<velocity>(i);
      std::cout << "id = "<<get<id>(i)<<" new position = "<<get<position>(i)<<
std::endl;
  }
};
*/

#if defined(HAVE_THRUST)
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
    typedef Particles<std::tuple<velocity>, 2, Vector, SearchMethod>
        container_type;
    //->
    //=        typedef Particles<std::tuple<velocity>,2> container_type;

    typedef typename container_type::position position;
    typedef typename container_type::query_type query_type;
    typedef typename container_type::raw_reference raw_reference;

    /*
     * set parameters for the MD simulation
     */
    const int timesteps = 300;
    const int nout = 200;
    const int timesteps_per_out = timesteps / nout;
    const double L = 31.0 / 1000.0;
    const int N = 10000;
    const double diameter = 0.0022;
    const double k = 1.0e01;
    const double dens = 1160.0;
    const double mass = (1.0 / 6.0) * PI * std::pow(0.5 * diameter, 3) * dens;
    const double reduced_mass = 0.5 * mass;
    const double dt = (1.0 / 50.0) * PI / sqrt(k / reduced_mass);

    /*
     * construct N particles
     */
    std::cout << "create particles" << std::endl;
    container_type particles(N);

    /*
     * randomly set position for N particles
     */
    std::cout << "randomly set positions" << std::endl;
    thrust::tabulate(get<position>(particles).begin(),
                     get<position>(particles).end(),
                     [=] CUDA_HOST_DEVICE(const int i) {
                       thrust::default_random_engine gen;
                       thrust::uniform_real_distribution<float> uni(0, L);
                       gen.discard(i);
                       return vdouble2(uni(gen), uni(gen));
                     });

    std::cout << "init vel" << std::endl;
    // initialise velocity to zero
    thrust::fill(get<velocity>(particles).begin(),
                 get<velocity>(particles).end(), vdouble2(0, 0));

    /*
     * initiate neighbour search on a periodic 2d domain of side length L
     * set average number of particles per cell to 1
     */
    std::cout << "init search" << std::endl;
    particles.init_neighbour_search(vdouble2(0, 0), vdouble2(L, L),
                                    vbool2(true, true));

    /*
     * perform MD timestepping
     */
    std::cout << "time loop" << std::endl;
    for (int io = 0; io < nout; ++io) {

      /*
       * on every i/o step write particle container to a vtk
       * unstructured grid file
       */
      std::cout << "." << std::flush;
#ifdef HAVE_VTK
      vtkWriteGrid("particles", io, particles.get_grid(true));
#endif
      for (int i = 0; i < timesteps_per_out; i++) {
        const query_type query = particles.get_query();
        Aboria::detail::for_each(
            particles.begin(), particles.end(),
            [=] CUDA_HOST_DEVICE(raw_reference i) {
              vdouble2 force_sum(0, 0);
              for (auto j = euclidean_search(query, get<position>(i), diameter);
                   j != false; ++j) {
                const double r = j.dx().norm();
                /*
                printf("found pair at (%f,%f) and (%f,%f) with dx = (%f,%f)\n"
                                                               ,get<position>(i)[0],
                                                               get<position>(i)[1],
                                                               get<position>(*j)[0],
                                                               get<position>(*j)[1],
                                                               j.dx()[0],
                                                               j.dx()[1]);
                                                               */
                if (r != 0) {
                  force_sum -= k * (diameter / r - 1) * j.dx();
                }
              }
              get<velocity>(i) += dt * force_sum / mass;
            });

        Aboria::detail::for_each(particles.begin(), particles.end(),
                                 [=] CUDA_HOST_DEVICE(raw_reference i) {
                                   get<position>(i) += dt * get<velocity>(i);
                                 });

        particles.update_positions();
      }
    }
    std::cout << std::endl;
  }
  //]

#endif

  void test_std_vector_CellList(void) {
    // nvcc bug cant compile with std::tuple
    //
    // workaround - define ABORIA_THRUST_USE_THRUST_TUPLE
    //
    // https://devtalk.nvidia.com/default/topic/1028112/cuda-setup-and-installation/nvcc-bug-related-to-gcc-6-lt-tuple-gt-header-/
    // helper_md<std::vector, CellList>();
  }

  void test_std_vector_CellListOrdered(void) {
    // helper_md<std::vector, CellListOrdered>();
  }

  void test_std_vector_HyperOctree(void) {
    // helper_md<std::vector, HyperOctree>();
  }

  // void test_thrust_vector_CellList(void) {
  //#if //defined(HAVE_THRUST)
  //    helper_md<thrust::device_vector,CellList>();
  //#end//if
  //}

  void test_thrust_vector_CellListOrdered(void) {
#if defined(HAVE_THRUST)
    helper_md<thrust::device_vector, CellListOrdered>();
#endif
  }

  void test_thrust_vector_HyperOctree(void) {
#if defined(HAVE_THRUST)
    helper_md<thrust::device_vector, HyperOctree>();
#endif
  }
};

#endif /* MD_TEST_H_ */
