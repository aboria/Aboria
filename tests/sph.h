/*
 * sph.cpp
 *
 *  Created on: 30 Jan 2014
 *      Author: mrobins
 */
#ifndef SPH_TEST_H_
#define SPH_TEST_H_

#include <cxxtest/TestSuite.h>

//[sph
/*`

Simulation of a water column in hydrostatic equilibrium. SPH discretises
the Navier-Stokes equations using radial interpolation kernels defined
over a given particle set. See the following papers for more details than
can be described here:

[@ http://onlinelibrary.wiley.com/doi/10.1002/fld.2677/abstract
M. Robinson, J. J. Monaghan, Direct numerical simulation of decaying
two-dimensional turbulence in a no-slip square box using smoothed par-
ticle hydrodynamics, International Journal for Numerical Methods in
Fluids 70 (1) (2012) 37–55]

[@http://www.sciencedirect.com/science/article/pii/S0301932213001882 M.
Robinson, M. Ramaioli, S. Luding, Fluid-particle flow simulations us- ing
two-way-coupled mesoscale SPH-DEM and validation, International Journal of
Multiphase Flow 59 (2014) 121–134.]

SPH is based on the idea of kernel interpolation. A fluid variable
$A(\mathbf{r})$ (such as velocity or density) is interpolated using a kernel
$W$, which depends on the smoothing length variable $h$.

\begin{equation}\label{Eq:integralInterpolant}
A(\mathbf{r}) = \int A(\mathbf{r'})W(\mathbf{r}-\mathbf{r'},h)d\mathbf{r'}.
\end{equation}

where $m_b$ and $\rho_b$ are the mass and density of particle $b$.

Here we have used the fifth-order Wendland kernel for $W$, which has the form

\begin{eqnarray}
W(q) = \frac{\beta}{h^d}
\begin{cases}
(2-q)^4(1+2q) & \text{for } 0 \leq q \le 2, \\
0             & \text{for } q > 2.
\end{cases}
\end{eqnarray}

where $q = ||\mathbf{r}-\mathbf{r'}||/h$.

The rate of change of density, or continuity equation, is given by

\begin{equation} \label{Eq:changeInDensity}
\frac{D\rho_a}{Dt} = \frac{1}{\Omega_a}\sum_b m_b \mathbf{v}\_{ab} \cdot
\nabla_a W_{ab}, \end{equation}

where $\mathbf{v}_{ab}=\mathbf{v}_a-\mathbf{v}_b$. The correction term
$\Omega_a$ due to variable $h$ is given by

\begin{equation}
\Omega_a = 1 - \frac{\partial h_a}{\partial \rho_a} \sum_b m_b \frac{\partial
W_{ab}(h_a)}{\partial h_a}. \end{equation}

The equation of state used is the standard quasi-compressible form

\begin{equation}
P = B \left ( \left ( \frac{\rho}{\rho_0} \right )^\gamma - 1 \right ),
\end{equation}

where $\gamma = 7$ is a typical value and $\rho_0$ is a reference density that
is normally set to the density of the fluid.

The SPH momentum equation is given as below. Viscosity is included by adding a
viscous term $\Pi$

\begin{equation}\label{Eq:sphJustPressureForce}
\frac{d\mathbf{v}_a}{dt} = -\sum_b m_b \left [ \left ( \frac{P_a}{\Omega_a
\rho_a^2} + \Pi_{ab} \right ) \nabla_a W_{ab}(h_a) + \left ( \frac{P_b}{\Omega_b
\rho_b^2} + \Pi_{ab} \right ) \nabla_a W_{ab}(h_b) \right ], \end{equation}

The SPH literature contains many different forms for $\Pi$. We have used the
term

\begin{equation}\label{Eq:monaghansViscousTerm}
\Pi_{ab} = - \alpha \frac{v_{sig} (\mathbf{v}\_{ab} \cdot \mathbf{r}\_{ab} )}{2
\overline{\rho}\_{ab} |\mathbf{r}\_{ab}|}, \end{equation}

where $v_{sig} = 2(c_s + |\mathbf{v}\_{ab} \cdot \mathbf{r}\_{ab}| /
|\mathbf{r}_{ab}| )$ is a signal velocity that represents the speed at which
information propagates between the particles.

The particle's position and velocity were integrated using the
Leapfrog second order method, which is also reversible in time in the absence of
viscosity. To preserve the reversibility of the simulation, $d\rho/dt$ was
calculated using the particle's position and velocity at the end of the
timestep, rather than the middle as is commonly done. The full integration
scheme is given by

\begin{align}
\mathbf{r}^{\frac{1}{2}} &= \mathbf{r}^{0} + \frac{\delta t}{2} \mathbf{v}^{0},
\\\\
\mathbf{v}^{\frac{1}{2}} &= \mathbf{v}^{0} + \frac{\delta t}{2}
F(\mathbf{r}^{-\frac{1}{2}},\mathbf{v}^{-\frac{1}{2}},\rho^{-\frac{1}{2}}), \\\\
\rho^{\frac{1}{2}} &= \rho^{0} + \frac{\delta t}{2}
D(\mathbf{r}^0,\mathbf{v}^0), \label{Eq:timestepDensity1} \\\\
\mathbf{v}^{1} &= \mathbf{v}^{0} + \delta t
F(\mathbf{r}^{\frac{1}{2}},\mathbf{v}^{\frac{1}{2}},\rho^{\frac{1}{2}}), \\\\
\mathbf{r}^{1} &= \mathbf{r}^{\frac{1}{2}} + \frac{\delta t}{2} \mathbf{v}^{1},
\\\\ \rho^{1} &= \rho^{\frac{1}{2}} + \frac{\delta t}{2}
D(\mathbf{r}^1,\mathbf{v}^1), \label{Eq:timestepDensity2} \end{align}

where $\mathbf{r}^0$, $\mathbf{r}^{1/2}$ and $\mathbf{r}^1$ is $\mathbf{r}$ at
the start, mid-point and end of the timestep respectively. The functions $F$ and
$D$ are respectivly the momentum and continuity equations given above. The
timestep $\delta t$ is bounded by the standard Courant condition

\begin{equation}
\delta t_1 \le \min_a \left ( 0.8 \frac{h_a}{v_{sig}} \right ),
\end{equation}

where the minimum is taken over all the particles.

 */
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>

const double PI = boost::math::constants::pi<double>();
const double WCON_WENDLAND = 21.0 / (256.0 * PI);
const unsigned int NDIM = 3;

/*
 * Note that we are using standard C++ function objects to implement
 * the kernel W and its gradient F. We can wrap these as lazy Aboria
 * functions that can be used within the symbolic Level 3 API.
 */

struct F_fun {
  typedef double result_type;
  double operator()(const double r, const double h) const {
    if (r == 0)
      return 0;
    const double q = r / h;
    if (q <= 2.0) {
      return (1 / std::pow(h, NDIM + 2)) * WCON_WENDLAND *
             (-4 * std::pow(2 - q, 3) * (1 + 2 * q) + 2 * std::pow(2 - q, 4)) /
             q;
    } else {
      return 0.0;
    }
  }
};

struct W_fun {
  typedef double result_type;
  double operator()(const double r, const double h) const {
    const double q = r / h;
    if (q <= 2.0) {
      return (1 / std::pow(h, NDIM)) * WCON_WENDLAND * std::pow(2.0 - q, 4) *
             (1.0 + 2.0 * q);
    } else {
      return 0.0;
    }
  }
};

ABORIA_BINARY_FUNCTION(F, F_fun, SymbolicDomain);
ABORIA_BINARY_FUNCTION(W, W_fun, SymbolicDomain);

//<-
class SPHTest : public CxxTest::TestSuite {
public:
  template <template <typename> class SearchMethod> void helper_sph(void) {
    //->
    //=int main() {
    ABORIA_VARIABLE(kernel_radius, double, "kernel radius");
    ABORIA_VARIABLE(velocity, vdouble3, "velocity");
    ABORIA_VARIABLE(velocity_tmp, vdouble3, "temp velocity");
    ABORIA_VARIABLE(varh_omega, double, "varh omega");
    ABORIA_VARIABLE(density, double, "density");
    ABORIA_VARIABLE(total_force, vdouble3, "total force");
    ABORIA_VARIABLE(is_fixed, uint8_t, "fixed boundary");
    ABORIA_VARIABLE(pressure_div_density2, double, "pressure div density2");

    //<-
    typedef Particles<
        std::tuple<kernel_radius, velocity, velocity_tmp, varh_omega, density,
                   total_force, is_fixed, pressure_div_density2>,
        3, std::vector, SearchMethod>
        sph_type;
    //->
    //=typedef Particles<
    //=    std::tuple<kernel_radius, velocity, velocity_tmp, varh_omega,
    //=               density, total_force, is_fixed,
    //=               pressure_div_density2>, 3>
    //=    sph_type;

    typedef position_d<3> position;
    sph_type sph;

    Symbol<position> p;
    Symbol<id> id_;
    Symbol<velocity> v;
    Symbol<velocity_tmp> v0;
    Symbol<density> rho;
    Symbol<total_force> dvdt;
    Symbol<is_fixed> fixed;
    Symbol<varh_omega> omega;
    Symbol<kernel_radius> h;
    Symbol<pressure_div_density2> pdr2;
    Label<0, sph_type> a(sph);
    Label<1, sph_type> b(sph);

    const int timesteps = 200;
    const int nout = 10;
    const int timesteps_per_out = timesteps / nout;
    const double L = 31.0 / 1000.0;
    const int nx = 5;

    /*
     * sph parameters
     */
    const double hfac = 1.5;
    const double visc = 8.9e-07;
    const double refd = 1000.0;
    const double dens = 1000.0;
    const double gamma = 7;
    const double VMAX = 2.0 * sqrt(2 * 9.81 * L);
    const double CSFAC = 10.0;
    const double spsound = CSFAC * VMAX;
    const double prb = std::pow(refd / dens, gamma - 1.0) *
                       std::pow(spsound, 2) * refd / gamma;
    const double psep = L / nx;
    double dt = std::min(0.25 * hfac * psep / spsound,
                         0.125 * std::pow(hfac * psep, 2) / visc);
    const double mass = dens * std::pow(psep, NDIM);

    std::cout << "h = " << hfac * psep << " vmax = " << VMAX << std::endl;

    const double time_damping = dt * 500;
    double t = 0;

    const vdouble3 low(0, 0, -3.0 * psep);
    const vdouble3 high(L, L, L);
    const vbool3 periodic(true, true, false);

    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < nx; j++) {
        for (int k = 0; k < nx + 3; k++) {
          typename sph_type::value_type p;
          get<position>(p) = low + vdouble3((i + 0.5) * psep, (j + 0.5) * psep,
                                            (k + 0.5) * psep);
          get<kernel_radius>(p) = hfac * psep;
          get<velocity>(p) = vdouble3(0, 0, 0);
          get<velocity_tmp>(p) = vdouble3(0, 0, 0);
          get<varh_omega>(p) = 1.0;
          get<density>(p) = dens;
          get<total_force>(p) = vdouble3(0, 0, 0);
          get<is_fixed>(p) = get<position>(p)[2] < 0;
          sph.push_back(p);
        }
      }
    }

    std::cout << "starting...." << std::endl;
    sph.init_neighbour_search(low, high, periodic);

    auto dx = create_dx(a, b);
    AccumulateWithinDistance<std::plus<vdouble3>> sum_vect(2 * hfac * psep);
    AccumulateWithinDistance<std::plus<double>> sum(2 * hfac * psep);
    Accumulate<Aboria::max<double>> max;
    max.set_init(0);
    Accumulate<Aboria::min<double>> min;
    min.set_init(L);
    VectorSymbolic<double, 3> vector;

    for (int i = 0; i < nout; ++i) {
#ifdef HAVE_VTK
      vtkWriteGrid("sph", i, sph.get_grid(true));
#endif
      for (int k = 0; k < timesteps_per_out; ++k) {
        /*
         * 0 -> 1/2 step
         */
        v[a] += if_else(fixed[a] == false, dt / 2, 0) * dvdt[a];
        if (t < time_damping)
          v[a] *= 0.98;
        p[a] += dt / 2 * v[a];

        /*
         * Calculate omega
         */
        omega[a] =
            1.0 - (mass / (rho[a] * NDIM)) *
                      sum(b, if_else(norm(dx) < 2 * h[a],
                                     pow(norm(dx), 2) * F(norm(dx), h[a]) +
                                         NDIM * W(norm(dx), h[a]),
                                     0));

        /*
         * 1/2 -> 1 step
         */

        /*
         * calculate change in density and kernel radius
         */
        rho[a] += dt * (mass / omega[a]) *
                  sum(b, if_else(norm(dx) < 2 * h[a],
                                 dot(v[b] - v[a], dx * F(norm(dx), h[a])), 0));
        h[a] = pow(mass / rho[a], 1.0 / NDIM);

        /*
         * reset neighbour search for new kernel radius
         */
        const double maxh = eval(max(a, h[a]));
        const double minh = eval(min(a, h[a]));
        sum_vect.set_max_distance(2 * maxh);
        sum.set_max_distance(2 * maxh);

        /*
         * advance velocity, position and calculate pdr2
         */
        v0[a] = v[a];
        v[a] += if_else(fixed[a] == false, dt / 2.0, 0.0) * dvdt[a];
        pdr2[a] = prb * (pow(rho[a] / refd, gamma) - 1.0) / pow(rho[a], 2);
        p[a] += dt / 2 * v0[a];

        /*
         * acceleration due to gravity
         */
        dvdt[a] = vdouble3(0, 0, -9.81);

        /*
         * acceleration on SPH calculation
         */
        dvdt[a] +=
            mass *
            sum_vect(b,
                     if_else(norm(dx) < 2 * h[a],

                             // pressure force
                             ((1.0 / omega[a]) * pdr2[a] * F(norm(dx), h[a]) +
                              (1.0 / omega[b]) * pdr2[b] * F(norm(dx), h[b])) *
                                     dx

                                 // viscous force
                                 + (v[a] - v[b]) * visc * (rho[a] + rho[b]) /
                                       (rho[a] * rho[b]) * F(norm(dx), h[a])

                                 ,
                             vector(0, 0, 0)));

        /*
         * 1/2 -> 1 step for velocity
         */
        v[a] = if_else(fixed[a] == false, v0[a] + dt / 2 * dvdt[a],
                       vector(0, 0, 0));

        t += dt;

        /*
         * sph timestep condition
         */
        dt = std::min(0.25 * minh / spsound, 0.125 * std::pow(minh, 2) / visc);
      }
      std::cout << "iteration " << i << std::endl;
    }
  }
  //]

  void test_CellListOrdered() { helper_sph<CellListOrdered>(); }

  void test_CellList() { helper_sph<CellList>(); }

  void test_Kdtree() { helper_sph<Kdtree>(); }
};

#endif /* SPH_TEST_H_ */
