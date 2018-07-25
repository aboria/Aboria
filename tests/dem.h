/*
 * dem_symbolic.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */
#ifndef DEM_TEST_H_
#define DEM_TEST_H_

#include <cxxtest/TestSuite.h>

//[dem
/*`
We wish to describe a set of particles in 3D cuboid, which is periodic in the
Cartesian $x$, and $y$ directions. There is a linear spring force
$\mathbf{f}\_{ij}$ plus dissipation term between particles $i$ and $j$ with a
rest separation at $s_i + s_j$ and a cutoff at $s_i + s_j$. That is, if
$\mathbf{r}\_i$ is the position of particle $i$ and
$\mathbf{dx}\_{ij}=\mathbf{r}\_j-\mathbf{r}\_j$, then

$$
\mathbf{f}\_{ij} = \begin{cases}
            \frac{s_i+s_j-|\mathbf{dx}\_{ij}|}{|\mathbf{dx}\_{ij}|}\mathbf{dx}\_{ij}
            + \gamma(\mathbf{v}\_j-\mathbf{v}\_i), & \text{for }
              |\mathbf{dx}\_{ij}|<s_i+s_j \\\\
            0 & \text{otherwise}.
            \end{cases}
$$

We wish to use a leap frog integrator to evolve positions $\mathbf{r}\_i$ using
velocities $\mathbf{v}\_i$ and accelerations $\mathbf{a}\_i = \frac{1}{m_i}
\sum_j \mathbf{f}\_{ij} - \mathbf{g}$. This gives the following update equations
for each timestep $n$

\begin{align*}
\mathbf{v}^{n+1}\_i &= \mathbf{v}^n_i + \frac{dt}{m_i} \sum_j \mathbf{f}^n_{ij}
\\\\
\mathbf{r}^{n+1}_i &= \mathbf{r}^n_i + dt\, \mathbf{v}^{n+1}_i.
\end{align*}

[$images/dem/domain.svg [width 100%]  [align center]]

This figure above shows the simulation domain. As well as periodic in $x$, and
$y$, we wish to add a soft surface at $z=0$ that also interacts with the
particles using the same linear spring force. The domain is left open at $z=h$.
We wish to initialise $N$ particles within this domain, and ensure that they are
non-interacting at $t=0$.
 */
#include <random>

#include "Aboria.h"
using namespace Aboria;

#include <boost/math/constants/constants.hpp>

//<-
class DEMTest : public CxxTest::TestSuite {
public:
  void test_dem(void) {
    //->
    //=int main() {
    const double PI = boost::math::constants::pi<double>();

    /*
     * create variable types
     */
    ABORIA_VARIABLE(radius, double, "radius")
    ABORIA_VARIABLE(mass, double, "mass")
    ABORIA_VARIABLE(velocity, vdouble3, "velocity")
    ABORIA_VARIABLE(acceleration, vdouble3, "acceleration")

    /*
     * create particle container
     */
    typedef Particles<std::tuple<radius, velocity, acceleration, mass>>
        dem_type;
    typedef position_d<3> position;
    dem_type dem;

    /*
     * simulation parameters
     */
    const int timesteps = 3000;
    const int nout = 200;
    const int timesteps_per_out = timesteps / nout;
    const double L = 31.0 / 1000.0;
    const int N = 30;

    /*
     * dem parameters
     */
    const double dem_diameter = 0.0022;
    const double dem_gamma = 0.0004;
    const double dem_k = 1.0e01;
    const double dem_dens = 1160.0;
    const double dem_mass_min =
        (1.0 / 6.0) * PI * std::pow(0.5 * dem_diameter, 3) * dem_dens;
    const double dem_mass_max =
        (1.0 / 6.0) * PI * std::pow(dem_diameter, 3) * dem_dens;
    const double dem_min_reduced_mass =
        dem_mass_min * dem_mass_max / (dem_mass_min + dem_mass_max);
    const double dt = (1.0 / 50.0) * PI /
                      sqrt(dem_k / dem_min_reduced_mass -
                           std::pow(0.5 * dem_gamma / dem_min_reduced_mass, 2));

    /*
     * initialise neighbour search with 3d cuboid domain,
     * periodic in x and y, set cell size so that there is
     * an average of 2 particles within each cell
     */
    dem.init_neighbour_search(vdouble3(0, 0, -dem_diameter),
                              vdouble3(L / 3, L / 3, L + dem_diameter),
                              vbool3(true, true, false), 2);

    /*
     * create N random, non-overlaping particles with
     * random radius between 0.5*dem_diameter and dem_diameter
     */
    std::uniform_real_distribution<double> uni(0, 1);
    std::default_random_engine generator;
    for (int i = 0; i < N; ++i) {
      bool free_position = false;
      dem_type::value_type p;

      get<radius>(p) = 0.25 * (uni(generator) + 1) * dem_diameter;
      while (free_position == false) {
        get<position>(p) =
            vdouble3(uni(generator) * L / 3, uni(generator) * L / 3,
                     uni(generator) * (L - dem_diameter) + dem_diameter / 2);
        free_position = true;
        /*
         * loop over all neighbouring particles within "dem_diameter" distance
         */
        for (auto j = euclidean_search(dem.get_query(), get<position>(p),
                                       dem_diameter);
             j != false; ++j) {
          if (j.dx().norm() < get<radius>(*j) + get<radius>(p)) {
            free_position = false;
            break;
          }
        }
      }
      dem.push_back(p);
    }

    /*
     * create symbols and labels in order to use the expression API
     */
    Symbol<position> p;
    Symbol<radius> r;
    Symbol<mass> m;
    Symbol<velocity> v;
    Symbol<acceleration> dvdt;
    Symbol<id> id_;
    Label<0, dem_type> a(dem);
    Label<1, dem_type> b(dem);

    /*
     * dx is a symbol representing the difference in positions of
     * particle a and b.
     */
    auto dx = create_dx(a, b);

    /*
     * sum is a symbolic function that sums a sequence of 3d vectors
     * over neighbouring particles within "dem_diameter"
     */
    AccumulateWithinDistance<std::plus<vdouble3>> sum(dem_diameter);

    /*
     * vector creates a new double vector of dimension 3
     */
    VectorSymbolic<double, 3> vector;

    /*
     * perform MD timestepping
     */
    v[a] = vdouble3::Zero();
    m[a] = (1.0 / 6.0) * PI * 8 * r[a] * r[a] * r[a] * dem_dens;
    for (int io = 0; io < nout; ++io) {
      std::cout << "." << std::flush;

#ifdef HAVE_VTK
      vtkWriteGrid("dem", io, dem.get_grid(true), {{"time", io / 2.0}});
#endif
      for (int i = 0; i < timesteps_per_out; i++) {
        p[a] += v[a] * dt;

        dvdt[a] =
            ( // spring force between dem particles
                sum(b, if_else(norm(dx) < r[a] + r[b] && id_[a] != id_[b],
                               -dem_k * ((r[a] + r[b]) / norm(dx) - 1) * dx +
                                   dem_gamma * (v[b] - v[a]),
                               vdouble3::Zero()))

                // spring force between particles and bottom wall
                + if_else(r[a] - p[a][2] > 0, dem_k * (r[a] - p[a][2]), 0.0) *
                      vdouble3(0, 0, 1)

                // gravity force
                + vdouble3(0, 0, -9.81) * m[a]

                ) /
            m[a];

        v[a] += dvdt[a] * dt;
      }
    }
    std::cout << std::endl;
  }
  //]
};

#endif /* DEM_TEST_H_ */
