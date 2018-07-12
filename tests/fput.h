#ifndef FPUT_TEST_H_
#define FPUT_TEST_H_

#include <cxxtest/TestSuite.h>

//[fput
/*`
The first-ever numerical experiment was performed in 1955 by Fermi, Pasta, Ulam,
and Tsingou. They simulated a lattice of particles, with nearest neighbour
particles connected by springs. The motion of the particles over time is chaotic
and quasi-periodic.

This experiment was the first to demonstrate the importance of computer
simulation in the analysis of non-linear systems, and demonstrates the paradox
that many apparently chaotic systems exhibit periodic behaviour.

We reproduce this experiment using Aboria below. Of course, being a 1D lattice
code, Aboria is overkill for this type of simulation, but given that the first
numerical experiment was a particle simulation I had to include it in the
examples.

The equations of motion for each particle in terms of its displacement from a
rest configuration are

$$
\frac{d^2 u_j}{dt^2} = \frac{c^2}{h^2} (u_{j+1}+u_{j-1}-2u_{j}) + \alpha
[(u_{j+1}-u_{j})^2 - (u_{j}-u_{j-1})^2]
$$

If we calculate the energy in the lowest mode

$$
E_1 = 0.5(\dot{A_1}^2 + w_1^2 A_1^2)
$$

where $w_1^2 = 4 \sin^2(\pi/2N)$, and $A_1 = \sqrt{2/(N+1)}\sum_{n=1}^N u_n
\sin(n\pi/(N+1))$ and set up an initial condition with all the energy in this
lowest mode, then we see that while the energy initially diffuses away to the
higher modes, after a long enough time the energy returns to the lowest modes
until the system is in its original state.

If you compile the code below with Cairo enabled, then it will output svg files
showing the image below. The circles are the particles, coloured red by their
displacement $u$, and the line shows the amount of energy in the lowest mode
versus time.

[$images/fput/fput999.svg [width 100%]  [align center]]

 */
#include "Aboria.h"
#include <random>
using namespace Aboria;

#include <boost/math/constants/constants.hpp>

//=int main() {
//<-
class FPUTTest : public CxxTest::TestSuite {
public:
  void test_fput(void) {
    //->
    const double PI = boost::math::constants::pi<double>();

    // setup types
    ABORIA_VARIABLE(velocity, vdouble1, "velocity")
    ABORIA_VARIABLE(acceleration, vdouble1, "acceleration")
    ABORIA_VARIABLE(acceleration0, vdouble1, "old acceleration")
    typedef Particles<std::tuple<velocity, acceleration, acceleration0>, 1>
        particles_t;
    typedef particles_t::position position;

    // simulation parameters
    const double c = 1.0;
    const double alpha = 0.25;
    const double h = 1.0;
    const double scale = std::pow(c, 2) / std::pow(h, 2);
    const int N = 32;

    const double w1 = 2.0 * std::sin(PI / (2.0 * N));
    const double Tf = 160 * 2 * PI / w1;
    const long timesteps = 1000000;
    const int nout = 1000;
    const int timesteps_per_out = timesteps / nout;
    const double dt = static_cast<double>(Tf) / timesteps;

    // create particles
    particles_t particles(N);

    // set intial variables, position here is simply a displacement
    // for each particle away from its resting position
    for (size_t i = 0; i < N; ++i) {
      get<position>(particles)[i][0] = std::sin((i + 1) * PI / (N + 1));
      get<acceleration>(particles)[i][0] = 0;
      get<velocity>(particles)[i][0] = 0;
    }

    // time loop
    std::vector<double> energy(nout, 0);
    for (size_t out = 0; out < nout; ++out) {
      for (size_t t = 0; t < timesteps_per_out; ++t) {
        // advance position: x1 = x0 + dt*v + 0.5*dt^2*a
        for (size_t i = 0; i < N; ++i) {
          get<position>(particles)[i] +=
              dt * get<velocity>(particles)[i] +
              0.5 * std::pow(dt, 2) * get<acceleration>(particles)[i];
        }

        // calculate acceleration: a1 = f(x1)
        for (size_t i = 0; i < N; ++i) {
          // get displacements, taking into account boundary conditions
          // of x0 = xN = 0
          const double &x = get<position>(particles)[i][0];
          const double xp1 =
              i != N - 1 ? get<position>(particles)[i + 1][0] : 0;
          const double xm1 = i != 0 ? get<position>(particles)[i - 1][0] : 0;

          get<acceleration0>(particles)[i] = get<acceleration>(particles)[i];
          get<acceleration>(particles)[i][0] =
              scale * (xp1 + xm1 - 2 * x) +
              alpha * (std::pow(xp1 - x, 2) - std::pow(x - xm1, 2));
        }

        // advance velocity: v1 = v0 + 0.5*dt*(a0 + a1)
        for (size_t i = 0; i < N; ++i) {
          get<velocity>(particles)[i] += 0.5 * dt *
                                         (get<acceleration0>(particles)[i] +
                                          get<acceleration>(particles)[i]);
        }
      }
      // calculate energy in the lowest mode
      // (w1*A1)^2 = potential energy
      // (Adot1)^2 = kinetic energy
      double A1 = 0;
      double Adot1 = 0;
      for (size_t i = 0; i < N; ++i) {
        const double &u = get<position>(particles)[i][0];
        const double &v = get<velocity>(particles)[i][0];
        A1 += std::sqrt(2.0 / (N + 1)) * u * std::sin((i + 1) * PI / (N + 1));
        Adot1 +=
            std::sqrt(2.0 / (N + 1)) * v * std::sin((i + 1) * PI / (N + 1));
      }
      energy[out] = 0.5 * (std::pow(Adot1, 2) + std::pow(w1 * A1, 2));

      // visualise system if compiled with cairo
#ifdef HAVE_CAIRO
      // create surface
      const int image_size = 512;
      const double ratio = 1.0 / 5.0;
      cairo_surface_t *surface = cairo_svg_surface_create(
          ("fput" + std::to_string(out) + ".svg").c_str(), image_size,
          image_size * ratio);
      cairo_svg_surface_restrict_to_version(surface, CAIRO_SVG_VERSION_1_2);
      cairo_t *cr = cairo_create(surface);

      // set source
      cairo_scale(cr, image_size, image_size);
      cairo_set_source_rgba(cr, 0, 0, 0, 1.0);
      const double lw = 0.01;
      cairo_set_line_width(cr, lw);

      // draw particles
      for (size_t i = 0; i < N; ++i) {
        const double &u = get<position>(particles)[i][0];
        const double pos = 0.1 * u + (i + 1) * dx;
        cairo_set_source_rgba(cr, std::abs(u), 0, 0, 1.0);
        cairo_arc(cr, pos, 0.5 * ratio, lw, 0, 2 * PI);
        cairo_fill(cr);
      }

      // draw energy
      cairo_set_source_rgba(cr, 1, 0, 0, 0.8);
      cairo_move_to(cr, 0, energy[0] * ratio);
      for (size_t i = 1; i < out; ++i) {
        cairo_line_to(cr, static_cast<double>(i) / nout,
                      (1.0 - 10 * energy[i]) * ratio);
      }
      cairo_stroke(cr);

      // clean up
      cairo_destroy(cr);
      cairo_surface_destroy(surface);
#endif // HAVE_CAIRO
    }
  }
  //]
};

#endif /* FPUT_TEST_H_ */
