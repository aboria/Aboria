/*
 * bd_symbolic.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */
#ifndef BD_TEST_H_
#define BD_TEST_H_

#include <cxxtest/TestSuite.h>

//[bd
/*`

This example creates two particle sets, the first representing a set of point
particles undergoing brownian motion, the second representing a set of spheres
with a given radius, different for each sphere.

The point particles diffuse around the three-dimensional, periodic domain, and
reflect off the spheres whenever they encounter them.

Let indicies $i$ and $j$ refer to a pair of point particles, and $a$ and $b$
refer to a pair of spheres. Let $\mathbf{p}\_i$ be the position of particle $i$,
and let $\mathbf{dx}\_{ij}$ refer to the shortest difference between the
positions of particles $i$ and $j$. Let $r_b$ be the radius of sphere $b$.
Then the update equations used to evolve the system are given by

$$
\mathbf{p}\_i = \mathbf{p}\_i + \sqrt{2\\, D\\, dt}\\, \mathbf{N}
        + \sum_b \begin{cases}
                    -2\left(\frac{r_b}{||\mathbf{dx}\_{ib}||} - 1\right)
\mathbf{dx}\_{ib}, & \text{for } ||\mathbf{dx}\_{ib}||<r_b \\\ 0 &
\text{otherwise}. \end{cases}
$$

where $D$ is the diffusion constant, $dt$ is the timestep and $\mathbf{N}$
is a three-dimensional vector containing random samples from a normal
distribution
*/
#include "Aboria.h"
#include <random>
using namespace Aboria;

//<-
class BDTest : public CxxTest::TestSuite {
public:
  template <template <typename> class SearchMethod> void helper_bd(void) {
    //->
    //=int main() {
    ABORIA_VARIABLE(radius, double, "radius")

    //<-
    typedef Particles<std::tuple<radius>, 3, std::vector, SearchMethod>
        spheres_type;
    typedef Particles<std::tuple<>, 3, std::vector, SearchMethod> points_type;
    //->
    //=        typedef Particles<std::tuple<radius>,3,std::vector> spheres_type;
    //=        typedef Particles<std::tuple<>,3,std::vector> points_type;

    typedef position_d<3> position;
    spheres_type spheres;

    const double L = 10.0;
    const double D = 1.0;
    const double dt = 0.01;
    const double timesteps = 500;

    spheres.push_back(vdouble3(0, 0, 0));
    get<radius>(spheres[0]) = 1.0;
    spheres.push_back(vdouble3(5, 0, 0));
    get<radius>(spheres[1]) = 2.0;
    spheres.push_back(vdouble3(0, -5, 0));
    get<radius>(spheres[2]) = 1.5;
    spheres.push_back(vdouble3(0, 0, 5));
    get<radius>(spheres[3]) = 1.0;

    points_type points;
    std::uniform_real_distribution<double> uni(-L + L / 5, L - L / 5);
    std::default_random_engine generator;
    for (int i = 0; i < 1000; ++i) {
      points.push_back(
          vdouble3(uni(generator), uni(generator), uni(generator)));
    }

    points.init_neighbour_search(vdouble3(-L, -L, -L), vdouble3(L, L, L),
                                 vbool3(true, true, true));
    spheres.init_neighbour_search(vdouble3(-L, -L, -L), vdouble3(L, L, L),
                                  vbool3(false, false, false));

    Symbol<position> p;
    Symbol<radius> r;
    Symbol<alive> alive_;
    Label<0, spheres_type> a(spheres);
    Label<1, spheres_type> b(spheres);
    Label<0, points_type> i(points);
    Label<1, points_type> j(points);
    auto dx = create_dx(i, b);
    Normal N;
    VectorSymbolic<double, 3> vector;
    AccumulateWithinDistance<std::bit_or<bool>> any(4);
    AccumulateWithinDistance<std::plus<vdouble3>> sum(4);

    int count_before = 0;
    for (auto point : points) {
      if ((get<position>(point) - get<position>(spheres[0])).norm() <
          get<radius>(spheres[0])) {
        count_before++;
      }
    }

    /*
     * Kill any points within spheres
     */
    alive_[i] = !any(b, if_else(norm(dx) < r[b], true, false));

    int count_after = 0;
    for (auto point : points) {
      if ((get<position>(point) - get<position>(spheres[0])).norm() <
          get<radius>(spheres[0])) {
        count_after++;
      }
    }

    std::cout << " found " << count_before << " before and " << count_after
              << " after" << std::endl;

    /*
     * Diffusion step for points and reflect off spheres
     */
    for (int ts = 1; ts < timesteps; ++ts) {
      if (ts % 10 == 0) {

#ifdef HAVE_VTK
        vtkWriteGrid("bd", ts / 10, points.get_grid(true));
#endif
        std::cout << "." << std::flush;
      }
      p[i] += std::sqrt(2 * D * dt) * vector(N[i], N[i], N[i]);
      p[i] += sum(b, if_else(norm(dx) < r[b] && norm(dx) != 0,
                             -2 * (r[b] / norm(dx) - 1), 0) *
                         dx);
    }
    std::cout << std::endl;
  }
  //]

  void test_CellListOrdered() {
    std::cout << "CellListOrdered" << std::endl;
    helper_bd<CellListOrdered>();
  }

  void test_CellList() {
    std::cout << "CellList" << std::endl;
    helper_bd<CellList>();
  }

  void test_Kdtree() {
    std::cout << "Kdtree" << std::endl;
    helper_bd<Kdtree>();
  }

  void test_HyperOctree() {
    std::cout << "HyperOctree" << std::endl;
    helper_bd<HyperOctree>();
  }
};

#endif /* BD_TEST_H_ */
