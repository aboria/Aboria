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

#ifndef RBF_INTERPOLATION_TEST_H_
#define RBF_INTERPOLATION_TEST_H_

#include <sys/resource.h>
#include <sys/time.h>

#include <chrono>
#include <cxxtest/TestSuite.h>
#include <fstream>
typedef std::chrono::system_clock Clock;

//[rbf_interpolation_global
/*`

Here we interpolate a function $f(x,y)$ over the two-dimensional unit cube from
$(0,0)$ to $(1,1)$, using the multiquadric Radial Basis Function (RBF).

The function to be interpolated is

\begin{equation}
f(x,y) = \exp(-9(x-\frac{1}{2})^2 - 9(y-\frac{1}{4})^2),
\end{equation}

and the multiquadric RBF $K(r)$ is given by

\begin{equation}
K(r) = \sqrt{r^2 + c^2}.
\end{equation}

We create a set of basis functions around a set of $N$ particles with positions
$\mathbf{x}_i$. In this case the radial coordinate around each point is $r =
||\mathbf{x}_i - \mathbf{x}||$. The function $f$ is approximated using the sum
of these basis functions, with the addition of a constant factor $\beta$

\begin{equation}
\overline{f}(\mathbf{x}) = \beta + \sum_i \alpha_j
K(||\mathbf{x}_i-\mathbf{x}||). \end{equation}

We exactly interpolate the function at $\mathbf{x}_i$ (i.e. set
$\overline{f}(\mathbf{x}_i)=f(\mathbf{x}_i)$), leading to

\begin{equation}
f(\mathbf{x}_i) = \beta + \sum_j \alpha_j K(||\mathbf{x}_j-\mathbf{x}_i||).
\end{equation}

Note that the sum $j$ is over the same particle set.

We also need one additional constraint to find $\beta$

\begin{equation}
0 = \sum_j \alpha_j.
\end{equation}

We rewrite the two previous equations as a linear algebra expression

\begin{equation}
\mathbf{\Phi} = \begin{bmatrix} \mathbf{G}&\mathbf{P} \\\\ \mathbf{P}^T &
\mathbf{0} \end{bmatrix} \begin{bmatrix} \mathbf{\alpha} \\\\ \mathbf{\beta}
\end{bmatrix} = \mathbf{W} \mathbf{\gamma}, \end{equation}

where $\mathbf{\Phi} = [f(\mathbf{x}_1),f(\mathbf{x}_2),...,f(\mathbf{x}_N),0]$,
$\mathbf{G}$ is an $N \times N$ matrix with elements $G_{ij} =
K(||\mathbf{x}_j-\mathbf{x}_i||)$ and $\mathbf{P}$ is a $N \times 1$ vector of
ones.

The basis function coefficients $\mathbf{\gamma}$ are found by solving this
equation. To do this, we use the iterative solvers found in the Eigen package
and Aboria's ability to wrap C++ function objects as Eigen matricies.

*/
#include "Aboria.h"
#include <random>

using namespace Aboria;

//<-
class RbfInterpolationTest : public CxxTest::TestSuite {
public:
  template <template <typename> class SearchMethod> void helper_global(void) {
#ifdef HAVE_EIGEN
    std::cout << "---------------------\n"
              << "Running global test....\n"
              << "---------------------" << std::endl;
    //->
    //=int main() {
    auto funct = [](const double x, const double y) {
      return std::exp(-9 * std::pow(x - 0.5, 2) - 9 * std::pow(y - 0.25, 2));
      // return x;
    };

    ABORIA_VARIABLE(alpha, double, "alpha value")
    ABORIA_VARIABLE(interpolated, double, "interpolated value")
    ABORIA_VARIABLE(constant2, double, "c2 value")

    typedef Particles<std::tuple<alpha, constant2, interpolated>, 2,
                      std::vector, SearchMethod>
        ParticlesType;
    typedef position_d<2> position;
    typedef typename ParticlesType::const_reference const_particle_reference;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
    ParticlesType knots;
    ParticlesType augment;
    ParticlesType test;

    const double c = 0.5;

    const int N = 1000;

    const int max_iter = 100;
    const int restart = 101;
    typename ParticlesType::value_type p;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
      get<position>(p) =
          vdouble2(distribution(generator), distribution(generator));
      get<constant2>(p) = std::pow(c, 2);
      knots.push_back(p);

      get<position>(p) =
          vdouble2(distribution(generator), distribution(generator));
      get<constant2>(p) = std::pow(c, 2);
      test.push_back(p);
    }

    // knots.init_neighbour_search(min,max,periodic);

    augment.push_back(p);

    auto kernel = [](const_particle_reference a, const_particle_reference b) {
      return std::sqrt((get<position>(b) - get<position>(a)).squaredNorm() +
                       get<constant2>(b));
    };

    auto one = [](const_particle_reference a, const_particle_reference b) {
      return 1.0;
    };

    auto G = create_dense_operator(knots, knots, kernel);
    auto P = create_dense_operator(knots, augment, one);
    auto Pt = create_dense_operator(augment, knots, one);
    auto Zero = create_zero_operator(augment, augment);

    auto W = create_block_operator<2, 2>(G, P, Pt, Zero);

    auto G_test = create_dense_operator(test, knots, kernel);
    auto W_test = create_block_operator<2, 2>(G_test, P, Pt, Zero);

    vector_type phi(N + 1), gamma(N + 1);
    for (size_t i = 0; i < knots.size(); ++i) {
      const double x = get<position>(knots[i])[0];
      const double y = get<position>(knots[i])[1];
      phi[i] = funct(x, y);
    }
    phi[knots.size()] = 0;

    matrix_type W_matrix(N + 1, N + 1);
    W.assemble(W_matrix);

    gamma = W_matrix.ldlt().solve(phi);

    Eigen::GMRES<matrix_type> gmres;
    gmres.setMaxIterations(max_iter);
    gmres.set_restart(restart);
    gmres.compute(W_matrix);
    gamma = gmres.solve(phi);
    std::cout << "GMRES:       #iterations: " << gmres.iterations()
              << ", estimated error: " << gmres.error() << std::endl;

    phi = W * gamma;
    double rms_error = 0;
    double scale = 0;
    for (size_t i = 0; i < knots.size(); ++i) {
      const double x = get<position>(knots[i])[0];
      const double y = get<position>(knots[i])[1];
      const double truth = funct(x, y);
      const double eval_value = phi[i];
      rms_error += std::pow(eval_value - truth, 2);
      scale += std::pow(truth, 2);
      // TS_ASSERT_DELTA(eval_value,truth,2e-3);
    }

    std::cout << "rms_error for global support, at centers  = "
              << std::sqrt(rms_error / scale) << std::endl;
    TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-6);

    phi = W_test * gamma;
    rms_error = 0;
    scale = 0;
    for (size_t i = 0; i < test.size(); ++i) {
      const double x = get<position>(test[i])[0];
      const double y = get<position>(test[i])[1];
      const double truth = funct(x, y);
      const double eval_value = phi[i];
      rms_error += std::pow(eval_value - truth, 2);
      scale += std::pow(truth, 2);
      // TS_ASSERT_DELTA(eval_value,truth,2e-3);
    }
    std::cout << "rms_error for global support, away from centers  = "
              << std::sqrt(rms_error / scale) << std::endl;
    TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-5);

//=}
//]
#endif // HAVE_EIGEN
  }

  template <template <typename> class SearchMethod> void helper_compact(void) {
#ifdef HAVE_EIGEN
    std::cout << "---------------------\n"
              << "Running compact test....\n"
              << "---------------------" << std::endl;

    //[rbf_interpolation_compact
    //=#include "Aboria.h"
    //=#include <random>
    //=int main() {
    auto funct = [](const double x, const double y) {
      return std::exp(-9 * std::pow(x - 0.5, 2) - 9 * std::pow(y - 0.25, 2));
      // return x;
    };

    ABORIA_VARIABLE(alpha, double, "alpha value")
    ABORIA_VARIABLE(interpolated, double, "interpolated value")
    ABORIA_VARIABLE(constant2, double, "c2 value")

    typedef Particles<std::tuple<alpha, constant2, interpolated>, 2,
                      std::vector, SearchMethod>
        ParticlesType;
    typedef position_d<2> position;
    typedef typename ParticlesType::const_reference const_particle_reference;
    typedef typename position::value_type const &const_position_reference;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    typedef Eigen::SparseMatrix<double> matrix_type;
    ParticlesType knots, augment;
    ParticlesType test;

    const double hfac = 4.0;
    vdouble2 min = vdouble2::Constant(0);
    vdouble2 max = vdouble2::Constant(1);
    vbool2 periodic = vbool2::Constant(false);

    const int N = 1000;

    const int max_iter = 100;
    const int restart = 101;
    const double delta = std::pow(double(N), -1.0 / 2.0);
    const double h = hfac * delta;
    std::cout << "using h = " << h << std::endl;
    const double RASM_size = 2 * h;
    const int RASM_n = N * std::pow(RASM_size, 2) / (max - min).prod();
    const double RASM_buffer = 0.9 * RASM_size;
    std::cout << "RASM_size = " << RASM_size << " RASM_n = " << RASM_n
              << " RASM_buffer = " << RASM_buffer << std::endl;

    typename ParticlesType::value_type p;

    std::default_random_engine generator(123);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
      get<position>(p) =
          vdouble2(distribution(generator), distribution(generator));
      knots.push_back(p);

      get<position>(p) =
          vdouble2(distribution(generator), distribution(generator));
      test.push_back(p);
    }
    augment.push_back(p);

    knots.init_neighbour_search(min, max, periodic);

    auto kernel = [h](const_position_reference dx, const_particle_reference a,
                      const_particle_reference b) {
      return std::pow(2.0 - dx.norm() / h, 4) * (1.0 + 2.0 * dx.norm() / h);
    };

    auto one = [](const_particle_reference a, const_particle_reference b) {
      return 1.0;
    };

    auto G = create_sparse_operator(knots, knots, 2 * h, kernel);
    auto P = create_dense_operator(knots, augment, one);
    auto Pt = create_dense_operator(augment, knots, one);
    auto Zero = create_zero_operator(augment, augment);

    auto W = create_block_operator<2, 2>(G, P, Pt, Zero);

    auto G_test = create_sparse_operator(test, knots, 2 * h, kernel);
    auto W_test = create_block_operator<2, 2>(G_test, P, Pt, Zero);

    vector_type phi(N + 1), gamma(N + 1);
    for (size_t i = 0; i < knots.size(); ++i) {
      const double x = get<position>(knots[i])[0];
      const double y = get<position>(knots[i])[1];
      phi[i] = funct(x, y);
    }
    phi[knots.size()] = 0;

    matrix_type W_matrix(N + 1, N + 1);
    W.assemble(W_matrix);

    Eigen::ConjugateGradient<matrix_type, Eigen::Lower | Eigen::Upper,
                             Eigen::DiagonalPreconditioner<double>>
        cg_test;
    cg_test.setMaxIterations(max_iter);
    cg_test.compute(W_matrix);
    gamma = cg_test.solve(phi);
    std::cout << "CG:          #iterations: " << cg_test.iterations()
              << ", estimated error: " << cg_test.error() << std::endl;

    Eigen::ConjugateGradient<
        matrix_type, Eigen::Lower | Eigen::Upper,
        SchwartzPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
        cg;
    cg.setMaxIterations(max_iter);
    cg.preconditioner().analyzePattern(W);
    cg.compute(W_matrix);
    gamma = cg.solve(phi);
    std::cout << "CG-RASM:     #iterations: " << cg.iterations()
              << ", estimated error: " << cg.error() << std::endl;

    Eigen::GMRES<matrix_type,
                 SchwartzPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
        gmres;
    gmres.setMaxIterations(max_iter);
    gmres.preconditioner().analyzePattern(W);
    gmres.set_restart(restart);
    gmres.compute(W_matrix);
    gamma = gmres.solve(phi);
    std::cout << "GMRES-RASM:  #iterations: " << gmres.iterations()
              << ", estimated error: " << gmres.error() << std::endl;

    Eigen::BiCGSTAB<matrix_type,
                    SchwartzPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
        bicg;
    bicg.setMaxIterations(max_iter);
    bicg.preconditioner().analyzePattern(W);
    bicg.compute(W_matrix);
    gamma = bicg.solve(phi);
    std::cout << "BiCGSTAB-RASM:#iterations: " << bicg.iterations()
              << ", estimated error: " << bicg.error() << std::endl;

    double rms_error = 0;
    double scale = 0;
    phi = W * gamma;
    for (size_t i = 0; i < knots.size(); ++i) {
      const double x = get<position>(knots[i])[0];
      const double y = get<position>(knots[i])[1];
      const double truth = funct(x, y);
      const double eval_value = phi[i];
      rms_error += std::pow(eval_value - truth, 2);
      scale += std::pow(truth, 2);
      // TS_ASSERT_DELTA(eval_value,truth,2e-3);
    }
    std::cout << "rms_error for compact support, at centers  = "
              << std::sqrt(rms_error / scale) << std::endl;
    TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-4);

    rms_error = 0;
    scale = 0;
    phi = W_test * gamma;
    for (size_t i = 0; i < test.size(); ++i) {
      const double x = get<position>(test[i])[0];
      const double y = get<position>(test[i])[1];
      const double truth = funct(x, y);
      const double eval_value = phi[i];
      rms_error += std::pow(eval_value - truth, 2);
      scale += std::pow(truth, 2);
      // TS_ASSERT_DELTA(eval_value,truth,2e-3);
    }

    std::cout << "rms_error for compact support, away from centers  = "
              << std::sqrt(rms_error / scale) << std::endl;
    TS_ASSERT_LESS_THAN(std::sqrt(rms_error / scale), 1e-2);

//=}
//]
#endif // HAVE_EIGEN
  }

  template <template <typename> class SearchMethod> void helper_h2(void) {
#ifdef HAVE_H2LIB
    std::cout << "---------------------\n"
              << "Running h2 test....\n"
              << "---------------------" << std::endl;
    //[rbf_interpolation_h2
    //=#include "Aboria.h"
    //=#include <random>
    //=int main() {
    auto funct = [](const double x, const double y) {
      return std::exp(-9 * std::pow(x - 0.5, 2) - 9 * std::pow(y - 0.25, 2));
      // return x;
    };

    ABORIA_VARIABLE(alpha, double, "alpha value")
    ABORIA_VARIABLE(interpolated, double, "interpolated value")

    typedef Particles<std::tuple<alpha, interpolated>, 2, std::vector,
                      SearchMethod>
        ParticlesType;
    typedef position_d<2> position;
    typedef typename ParticlesType::const_reference const_particle_reference;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
    ParticlesType knots;
    ParticlesType augment;
    ParticlesType test;

    const double c = 1.0 / 0.50;
    vdouble2 min = vdouble2::Constant(0);
    vdouble2 max = vdouble2::Constant(1);
    vbool2 periodic = vbool2::Constant(false);

    const int N = 2000;

    // const double RASM_size = 0.3 / c;
    // const int RASM_n = N * std::pow(RASM_size, 2) / (max - min).prod();
    // const double RASM_buffer = 0.9 * RASM_size;
    // std::cout << "RASM_size = "<<RASM_size<<" RASM_n = "<<RASM_n<<"
    // RASM_buffer = "<<RASM_buffer<<std::endl;

    // const int nx = 3;
    const int max_iter = 200;
    // const int restart = 101;
    // const double delta = 1.0 / nx;
    typename ParticlesType::value_type p;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
      get<position>(p) =
          vdouble2(distribution(generator), distribution(generator));
      knots.push_back(p);

      get<position>(p) =
          vdouble2(distribution(generator), distribution(generator));
      test.push_back(p);
    }

    const size_t order = 6;
    knots.init_neighbour_search(min, max, periodic, std::pow(order, 2));
    test.init_neighbour_search(min, max, periodic, std::pow(order, 2));

    augment.push_back(p);

    const double mscale = std::sqrt(3.0) * c;
    auto kernel = [&](const auto &a, const auto &b) {
      const double r = (b - a).norm();
      return (1.0 + mscale * r) * std::exp(-r * mscale);
    };

    auto self_kernel = [&](const_particle_reference a,
                           const_particle_reference b) {
      double result = kernel(get<position>(a), get<position>(b));
      if (get<id>(a) == get<id>(b)) {
        result += 1e-5;
      }
      return result;
    };

    auto G = create_h2_operator(knots, knots, order, kernel, self_kernel, 1.0);
    G.get_first_kernel().compress(1e-10);
    matrix_type G_matrix(N, N);
    G.assemble(G_matrix);
    auto Gtest = create_fmm_operator<order>(test, knots, kernel, self_kernel);
    // knots.init_neighbour_search(min, max, periodic, 35);

    vector_type phi(N), gamma(N);
    for (size_t i = 0; i < knots.size(); ++i) {
      const double x = get<position>(knots[i])[0];
      const double y = get<position>(knots[i])[1];
      phi[i] = funct(x, y);
    }
    // phi[knots.size()] = 0;

    Eigen::BiCGSTAB<decltype(G), Eigen::DiagonalPreconditioner<double>> dgmres;
    dgmres.setMaxIterations(max_iter);
    dgmres.compute(G);
    gamma = dgmres.solve(phi);
    std::cout << "BiCGSTAB:  #iterations: " << dgmres.iterations()
              << ", estimated error: " << dgmres.error()
              << " h2 error = " << (G * gamma - phi).norm()
              << " true error = " << (G_matrix * gamma - phi).norm()
              << std::endl;

    Eigen::BiCGSTAB<decltype(G),
                    SchwartzPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
        bicg_rasm2;
    bicg_rasm2.setMaxIterations(max_iter);
    bicg_rasm2.compute(G);
    gamma = bicg_rasm2.solve(phi);
    std::cout << "BiCGSTAB-RASM-local:#iterations: " << bicg_rasm2.iterations()
              << ", estimated error: " << bicg_rasm2.error()
              << " h2 error = " << (G * gamma - phi).norm()
              << " true error = " << (G_matrix * gamma - phi).norm()
              << std::endl;

    vector_type phi_test = G * gamma;
    double rms_error = 0;
    double scale = 0;
    for (size_t i = 0; i < knots.size(); ++i) {
      const double x = get<position>(knots[i])[0];
      const double y = get<position>(knots[i])[1];
      const double truth = funct(x, y);
      const double eval_value = phi_test[i];
      rms_error += std::pow(eval_value - truth, 2);
      scale += std::pow(truth, 2);
      // TS_ASSERT_DELTA(eval_value,truth,2e-3);
    }

    std::cout << "rms_error for global support, at centers  = "
              << std::sqrt(rms_error / scale) << std::endl;

    rms_error = 0;
    scale = 0;
    phi_test = Gtest * gamma;
    for (size_t i = 0; i < test.size(); ++i) {
      const vdouble2 p = get<position>(test)[i];
      const double eval_value = phi_test[i];
      const double truth = funct(p[0], p[1]);
      rms_error += std::pow(eval_value - truth, 2);
      scale += std::pow(truth, 2);
      // TS_ASSERT_DELTA(eval_value,truth,2e-3);
    }
    std::cout << "rms_error for global support, away from centers  = "
              << std::sqrt(rms_error / scale) << std::endl;

//=}
//]
#endif // HAVE_H2LIB
  }

  template <unsigned int D, template <typename> class SearchMethod>
  void helper_param_sweep(const double sigma, const int N, const size_t order,
                          std::ofstream &out_it, std::ofstream &out_err,
                          std::ofstream &out_setup, std::ofstream &out_solve,
                          std::ofstream &out_h2_error, std::ofstream &out_mem,
                          int which) {
#ifdef HAVE_H2LIB

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    long mem = usage.ru_maxrss;

    const int width = 11;
    char *argv[] = {"program name", NULL};
    int argc = sizeof(argv) / sizeof(char *) - 1;
    char **argv2 = &argv[0];
    init_h2lib(&argc, &argv2);

    out_it << std::setw(width) << N << " " << std::setw(width) << sigma << " "
           << std::setw(width) << D << " " << std::setw(width) << order;
    out_err << std::setw(width) << N << " " << std::setw(width) << sigma << " "
            << std::setw(width) << D << " " << std::setw(width) << order;
    out_setup << std::setw(width) << N << " " << std::setw(width) << sigma
              << " " << std::setw(width) << D << " " << std::setw(width)
              << order;
    out_solve << std::setw(width) << N << " " << std::setw(width) << sigma
              << " " << std::setw(width) << D << " " << std::setw(width)
              << order;
    out_h2_error << std::setw(width) << N << " " << std::setw(width) << sigma
                 << " " << std::setw(width) << D << " " << std::setw(width)
                 << order;
    out_mem << std::setw(width) << N << " " << std::setw(width) << sigma << " "
            << std::setw(width) << D << " " << std::setw(width) << order;

    auto funct = [](const auto &x) {
      // rosenbrock functions, need even D
      double ret = 0;
      for (size_t i = 0; i < D / 2; ++i) {
        ret += 100 * std::pow(std::pow(x[2 * i], 2) - x[2 * i + 1], 2) +
               std::pow(x[2 * i] - 1, 2);
      }
      return ret;
    };

    ABORIA_VARIABLE(alpha, double, "alpha value")
    ABORIA_VARIABLE(interpolated, double, "interpolated value")

    typedef Particles<std::tuple<alpha, interpolated>, D, std::vector,
                      SearchMethod>
        ParticlesType;
    typedef position_d<D> position;
    typedef typename ParticlesType::const_reference const_particle_reference;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_type;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_type;
    ParticlesType knots[3];
    knots[which].resize(N);

    const double c = 1.0 / sigma;
    auto min = Vector<double, D>::Constant(0);
    auto max = Vector<double, D>::Constant(1);
    auto periodic = Vector<bool, D>::Constant(false);

    const int max_iter = 1000;

    std::default_random_engine generator(0);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (int i = 0; i < N; ++i) {
      for (size_t d = 0; d < D; ++d) {
        get<position>(knots[which])[i][d] = distribution(generator);
      }
    }

    // const size_t n_subdomain = static_cast<size_t>(std::pow(0.12,D)*N);
    // const size_t order = static_cast<size_t>(std::pow(n_subdomain,1.0/D));
    const size_t n_subdomain = std::pow(order, D);
    knots[which].init_neighbour_search(min, max, periodic, n_subdomain);

    const double jitter = 1e-5;

    const double gscale = std::pow(c, 2);
    auto gaussian_kernel = [&](const auto &a, const auto &b) {
      return std::exp(-(b - a).squaredNorm() * gscale);
    };

    const double mscale = std::sqrt(3.0) * c;
    auto matern_kernel = [&](const auto &a, const auto &b) {
      const double r = (b - a).norm();
      return (1.0 + mscale * r) * std::exp(-r * mscale);
    };

    auto wendland_kernel = [&](const auto &a, const auto &b) {
      const double r = (b - a).norm();
      if (r < 2 * sigma) {
        return std::pow(2.0 - r * c, 4) * (1.0 + 2.0 * r * c);
      } else {
        return 0.0;
      }
    };

    auto matern_self_kernel = [&](const_particle_reference a,
                                  const_particle_reference b) {
      double ret = matern_kernel(get<position>(a), get<position>(b));
      if (get<id>(a) == get<id>(b)) {
        ret += jitter;
      }
      return ret;
    };

    auto gaussian_self_kernel = [&](const_particle_reference a,
                                    const_particle_reference b) {
      double ret = gaussian_kernel(get<position>(a), get<position>(b));
      if (get<id>(a) == get<id>(b)) {
        ret += jitter;
      }
      return ret;
    };

    auto wendland_self_kernel = [&](const_particle_reference a,
                                    const_particle_reference b) {
      double ret = wendland_kernel(get<position>(a), get<position>(b));
      if (get<id>(a) == get<id>(b)) {
        ret += jitter;
      }
      return ret;
    };

    const double eta = D > 3 ? -1 : 1;
    auto Ggaussian = create_h2_operator(
        knots[0], knots[0], order, gaussian_kernel, gaussian_self_kernel, eta);
    Ggaussian.get_first_kernel().compress(1e-10);

    auto Gmatern = create_h2_operator(knots[1], knots[1], order, matern_kernel,
                                      matern_self_kernel, eta);
    Gmatern.get_first_kernel().compress(1e-10);

    auto Gwendland = create_h2_operator(
        knots[2], knots[2], order, wendland_kernel, wendland_self_kernel, eta);
    Gwendland.get_first_kernel().compress(1e-10);

    // const int nleaf = std::pow(sigma,2)*N;
    // knots.init_neighbour_search(min, max, periodic, nleaf);

    vector_type phi(N), phi_matrix(N), gamma(N);
    gamma = vector_type::Random(N);
    if (N < 5000) {
      matrix_type Gmatrix(N, N);
      if (which == 0) {
        Ggaussian.assemble(Gmatrix);
        phi = Ggaussian * gamma;
      } else if (which == 1) {
        Gmatern.assemble(Gmatrix);
        phi = Gmatern * gamma;
      } else {
        Gwendland.assemble(Gmatrix);
        phi = Gwendland * gamma;
      }
      phi_matrix = Gmatrix * gamma;
      out_h2_error << " " << std::setw(width)
                   << (phi - phi_matrix).norm() / phi_matrix.norm()
                   << std::endl;
    } else {
      out_h2_error << std::endl;
    }

    for (size_t i = 0; i < knots[which].size(); ++i) {
      const auto x = get<position>(knots[which][i]);
      phi[i] = funct(x);
    }

    if (which == 0) {
      // Eigen::BiCGSTAB<decltype(Ggaussian),
      // Eigen::DiagonalPreconditioner<double>>
      Eigen::BiCGSTAB<decltype(Ggaussian), Eigen::IdentityPreconditioner> bicg;
      bicg.setMaxIterations(max_iter);
      auto t0 = Clock::now();
      bicg.compute(Ggaussian);
      auto t1 = Clock::now();
      gamma = bicg.solve(phi);
      auto t2 = Clock::now();
      out_it << " " << std::setw(width) << bicg.iterations();
      out_err << " " << std::setw(width) << bicg.error();
      out_setup
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
      out_solve
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();

      getrusage(RUSAGE_SELF, &usage);
      out_mem << " " << std::setw(width) << (usage.ru_maxrss - mem) / 1000;
    } else if (which == 1) {
      // Eigen::BiCGSTAB<decltype(Gmatern),
      // Eigen::DiagonalPreconditioner<double>>
      Eigen::BiCGSTAB<decltype(Gmatern), Eigen::IdentityPreconditioner> bicg;
      bicg.setMaxIterations(max_iter);
      auto t0 = Clock::now();
      bicg.compute(Gmatern);
      auto t1 = Clock::now();
      gamma = bicg.solve(phi);
      auto t2 = Clock::now();
      out_it << " " << std::setw(width) << bicg.iterations();
      out_err << " " << std::setw(width) << bicg.error();
      out_setup
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
      out_solve
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
      getrusage(RUSAGE_SELF, &usage);
      out_mem << " " << std::setw(width) << (usage.ru_maxrss - mem) / 1000;
    } else {
      // Eigen::BiCGSTAB<decltype(Gwendland),
      // Eigen::DiagonalPreconditioner<double>>
      Eigen::BiCGSTAB<decltype(Gwendland), Eigen::IdentityPreconditioner> bicg;
      bicg.setMaxIterations(max_iter);
      auto t0 = Clock::now();
      bicg.compute(Gwendland);
      auto t1 = Clock::now();
      gamma = bicg.solve(phi);
      auto t2 = Clock::now();
      out_it << " " << std::setw(width) << bicg.iterations();
      out_err << " " << std::setw(width) << bicg.error();
      out_setup
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
      out_solve
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
      getrusage(RUSAGE_SELF, &usage);
      out_mem << " " << std::setw(width) << (usage.ru_maxrss - mem) / 1000;
    }

    if (D < 3) {
      if (which == 0) {
        Eigen::BiCGSTAB<decltype(Ggaussian),
                        SchwartzPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
            bicg_rasm2;
        bicg_rasm2.setMaxIterations(max_iter);
        // bicg_rasm2.preconditioner().set_neighbourhood_buffer_size(sigma);
        // bicg_rasm2.preconditioner().set_coarse_grid_n(20);
        auto t0 = Clock::now();
        bicg_rasm2.compute(Ggaussian);
        auto t1 = Clock::now();
        gamma = bicg_rasm2.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_rasm2.iterations();
        out_err << " " << std::setw(width) << bicg_rasm2.error();
        out_setup << " " << std::setw(width)
                  << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0)
                         .count();
        out_solve << " " << std::setw(width)
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1)
                         .count();
        getrusage(RUSAGE_SELF, &usage);
        out_mem << " " << std::setw(width) << (usage.ru_maxrss - mem) / 1000;
      } else if (which == 1) {
        Eigen::BiCGSTAB<decltype(Gmatern),
                        SchwartzPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
            bicg_rasm2;
        bicg_rasm2.setMaxIterations(max_iter);
        // bicg_rasm2.preconditioner().set_neighbourhood_buffer_size(sigma);
        // bicg_rasm2.preconditioner().set_coarse_grid_n(20);
        auto t0 = Clock::now();
        bicg_rasm2.compute(Gmatern);
        auto t1 = Clock::now();
        gamma = bicg_rasm2.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_rasm2.iterations();
        out_err << " " << std::setw(width) << bicg_rasm2.error();
        out_setup << " " << std::setw(width)
                  << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0)
                         .count();
        out_solve << " " << std::setw(width)
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1)
                         .count();
        getrusage(RUSAGE_SELF, &usage);
        out_mem << " " << std::setw(width) << (usage.ru_maxrss - mem) / 1000;
      } else {
        Eigen::BiCGSTAB<decltype(Gwendland),
                        SchwartzPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
            bicg_rasm2;
        bicg_rasm2.setMaxIterations(max_iter);
        // bicg_rasm2.preconditioner().set_neighbourhood_buffer_size(sigma);
        // bicg_rasm2.preconditioner().set_coarse_grid_n(20);
        auto t0 = Clock::now();
        bicg_rasm2.compute(Gwendland);
        auto t1 = Clock::now();
        gamma = bicg_rasm2.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_rasm2.iterations();
        out_err << " " << std::setw(width) << bicg_rasm2.error();
        out_setup << " " << std::setw(width)
                  << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0)
                         .count();
        out_solve << " " << std::setw(width)
                  << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1)
                         .count();
        getrusage(RUSAGE_SELF, &usage);
        out_mem << " " << std::setw(width) << (usage.ru_maxrss - mem) / 1000;
      }
    } else {
      out_it << " " << std::setw(width) << " ";
      out_err << " " << std::setw(width) << " ";
      out_setup << " " << std::setw(width) << " ";
      out_solve << " " << std::setw(width) << " ";
      out_mem << " " << std::setw(width) << " ";
    }

    if (which == 0) {
      Eigen::BiCGSTAB<decltype(Ggaussian), SchwartzSamplingPreconditioner<
                                               Eigen::LLT<Eigen::MatrixXd>>>
          bicg_rasm;
      bicg_rasm.setMaxIterations(max_iter);
      bicg_rasm.preconditioner().set_number_of_random_particles(300);
      bicg_rasm.preconditioner().set_sigma(1.0 / c);
      bicg_rasm.preconditioner().set_rejection_sampling_scale(1.0);
      auto t0 = Clock::now();
      bicg_rasm.compute(Ggaussian);
      auto t1 = Clock::now();
      gamma = bicg_rasm.solve(phi);
      auto t2 = Clock::now();
      out_it << " " << std::setw(width) << bicg_rasm.iterations();
      out_err << " " << std::setw(width) << bicg_rasm.error();
      out_setup
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
      out_solve
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    } else if (which == 1) {
      Eigen::BiCGSTAB<decltype(Gmatern), SchwartzSamplingPreconditioner<
                                             Eigen::LLT<Eigen::MatrixXd>>>
          bicg_rasm;
      bicg_rasm.setMaxIterations(max_iter);
      bicg_rasm.preconditioner().set_number_of_random_particles(300);
      bicg_rasm.preconditioner().set_sigma(1.0 / (2 * c));
      bicg_rasm.preconditioner().set_rejection_sampling_scale(1.0);
      auto t0 = Clock::now();
      bicg_rasm.compute(Gmatern);
      auto t1 = Clock::now();
      gamma = bicg_rasm.solve(phi);
      auto t2 = Clock::now();
      out_it << " " << std::setw(width) << bicg_rasm.iterations();
      out_err << " " << std::setw(width) << bicg_rasm.error();
      out_setup
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
      out_solve
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    } else {
      Eigen::BiCGSTAB<decltype(Gwendland), SchwartzSamplingPreconditioner<
                                               Eigen::LLT<Eigen::MatrixXd>>>
          bicg_rasm;
      bicg_rasm.setMaxIterations(max_iter);
      bicg_rasm.preconditioner().set_number_of_random_particles(300);
      bicg_rasm.preconditioner().set_sigma(1.0 / c);
      bicg_rasm.preconditioner().set_rejection_sampling_scale(1.0);
      auto t0 = Clock::now();
      bicg_rasm.compute(Gwendland);
      auto t1 = Clock::now();
      gamma = bicg_rasm.solve(phi);
      auto t2 = Clock::now();
      out_it << " " << std::setw(width) << bicg_rasm.iterations();
      out_err << " " << std::setw(width) << bicg_rasm.error();
      out_setup
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
      out_solve
          << " " << std::setw(width)
          << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    }

    /*
    if (which == 0) {
        Eigen::BiCGSTAB<decltype(Ggaussian),
                        NystromPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
            bicg_nystrom;
        bicg_nystrom.setMaxIterations(max_iter);
        bicg_nystrom.preconditioner().set_number_of_random_particles(500);
        bicg_nystrom.preconditioner().set_lambda(1e-4);
        auto t0 = Clock::now();
        bicg_nystrom.compute(Ggaussian);
        auto t1 = Clock::now();
        gamma = bicg_nystrom.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_nystrom.iterations();
        out_err << " " << std::setw(width) << bicg_nystrom.error();
        out_setup << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t1
    - t0).count(); out_solve << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t2
    - t1).count(); } else if (which == 1) { Eigen::BiCGSTAB<decltype(Gmatern),
                        NystromPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
            bicg_nystrom;
        bicg_nystrom.setMaxIterations(max_iter);
        bicg_nystrom.preconditioner().set_number_of_random_particles(500);
        bicg_nystrom.preconditioner().set_lambda(1e-4);
        auto t0 = Clock::now();
        bicg_nystrom.compute(Gmatern);
        auto t1 = Clock::now();
        gamma = bicg_nystrom.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_nystrom.iterations();
        out_err << " " << std::setw(width) << bicg_nystrom.error();
        out_setup << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t1
    - t0).count(); out_solve << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t2
    - t1).count(); } else { Eigen::BiCGSTAB<decltype(Gwendland),
                        NystromPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
            bicg_nystrom;
        bicg_nystrom.setMaxIterations(max_iter);
        bicg_nystrom.preconditioner().set_number_of_random_particles(500);
        bicg_nystrom.preconditioner().set_lambda(1e-4);
        auto t0 = Clock::now();
        bicg_nystrom.compute(Gwendland);
        auto t1 = Clock::now();
        gamma = bicg_nystrom.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_nystrom.iterations();
        out_err << " " << std::setw(width) << bicg_nystrom.error();
        out_setup << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t1
    - t0).count(); out_solve << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t2
    - t1).count();
    }

    if (which == 0) {
        Eigen::BiCGSTAB<decltype(Ggaussian),
                        ChebyshevPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
            bicg_cheb;
        bicg_cheb.setMaxIterations(max_iter);
        bicg_cheb.preconditioner().set_order(30);
        auto t0 = Clock::now();
        bicg_cheb.compute(Ggaussian);
        auto t1 = Clock::now();
        gamma = bicg_cheb.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_cheb.iterations();
        out_err << " " << std::setw(width) << bicg_cheb.error();
        out_setup << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t1
    - t0).count(); out_solve << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t2
    - t1).count(); } else if (which == 1) { Eigen::BiCGSTAB<decltype(Gmatern),
                        ChebyshevPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
            bicg_cheb;
        bicg_cheb.setMaxIterations(max_iter);
        bicg_cheb.preconditioner().set_order(22);
        auto t0 = Clock::now();
        bicg_cheb.compute(Gmatern);
        auto t1 = Clock::now();
        gamma = bicg_cheb.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_cheb.iterations();
        out_err << " " << std::setw(width) << bicg_cheb.error();
        out_setup << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t1
    - t0).count(); out_solve << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t2
    - t1).count(); } else { Eigen::BiCGSTAB<decltype(Gwendland),
                        ChebyshevPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
            bicg_cheb;
        bicg_cheb.setMaxIterations(max_iter);
        bicg_cheb.preconditioner().set_order(22);
        auto t0 = Clock::now();
        bicg_cheb.compute(Gwendland);
        auto t1 = Clock::now();
        gamma = bicg_cheb.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_cheb.iterations();
        out_err << " " << std::setw(width) << bicg_cheb.error();
        out_setup << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t1
    - t0).count(); out_solve << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t2
    - t1).count();
    }

    if (which == 0) {
        Eigen::BiCGSTAB<decltype(Ggaussian),
                        ReducedOrderPreconditioner<H2LibCholeskyDecomposition>>
            bicg_ro;
        bicg_ro.setMaxIterations(max_iter);
        bicg_ro.preconditioner().set_tolerance(1e-6);
        auto t0 = Clock::now();
        bicg_ro.compute(Ggaussian);
        auto t1 = Clock::now();
        gamma = bicg_ro.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_ro.iterations();
        out_err << " " << std::setw(width) << bicg_ro.error();
        out_setup << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t1
    - t0).count(); out_solve << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t2
    - t1).count(); } else if (which == 1) { Eigen::BiCGSTAB<decltype(Gmatern),
                        ReducedOrderPreconditioner<H2LibCholeskyDecomposition>>
            bicg_ro;
        bicg_ro.setMaxIterations(max_iter);
        bicg_ro.preconditioner().set_tolerance(1e-6);
        auto t0 = Clock::now();
        bicg_ro.compute(Gmatern);
        auto t1 = Clock::now();
        gamma = bicg_ro.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_ro.iterations();
        out_err << " " << std::setw(width) << bicg_ro.error();
        out_setup << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t1
    - t0).count(); out_solve << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t2
    - t1).count(); } else { Eigen::BiCGSTAB<decltype(Gwendland),
                        ReducedOrderPreconditioner<H2LibCholeskyDecomposition>>
            bicg_ro;
        bicg_ro.setMaxIterations(max_iter);
        bicg_ro.preconditioner().set_tolerance(1e-6);
        auto t0 = Clock::now();
        bicg_ro.compute(Gwendland);
        auto t1 = Clock::now();
        gamma = bicg_ro.solve(phi);
        auto t2 = Clock::now();
        out_it << " " << std::setw(width) << bicg_ro.iterations();
        out_err << " " << std::setw(width) << bicg_ro.error();
        out_setup << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t1
    - t0).count(); out_solve << " " << std::setw(width)
                         << std::chrono::duration_cast<std::chrono::seconds>(t2
    - t1).count();
    }
    */

    /*
        Eigen::BiCGSTAB<decltype(Ggaussian),
                        ReducedOrderPreconditioner<H2LibCholeskyDecomposition>>
            bicg_pre;
        bicg_pre.preconditioner().set_tolerance(1e-2);
        bicg_pre.setMaxIterations(max_iter);
        bicg_pre.compute(Ggaussian);
        gamma = bicg_pre.solve(phi);
        out_it << " " << bicg_pre.iterations();
        out_err << " " << bicg_pre.error();
        */
    out_it << std::endl;
    out_err << std::endl;
    out_setup << std::endl;
    out_solve << std::endl;
    out_mem << std::endl;

#endif // HAVE_H2LIB
  }

  void test_param_sweep(void) {
    std::cout << "-------------------------------------------\n"
              << "Running precon param sweep             ....\n"
              << "------------------------------------------" << std::endl;
    std::ofstream out_it;
    out_it.open("iterations_gaussian.txt", std::ios::out);
    std::ofstream out_error;
    out_error.open("error_gaussian.txt", std::ios::out);
    std::ofstream out_setup;
    out_setup.open("setup_gaussian.txt", std::ios::out);
    std::ofstream out_solve;
    out_solve.open("solve_gaussian.txt", std::ios::out);
    std::ofstream out_h2_error;
    out_h2_error.open("h2_error_gaussian.txt", std::ios::out);
    std::ofstream out_mem;
    out_mem.open("mem_gaussian.txt", std::ios::out);

    const int width = 12;
    auto header = [](auto &out) {
      out << std::setw(width) << "N " << std::setw(width) << "sigma "
          << std::setw(width) << "D " << std::setw(width) << "order "
          << std::setw(width) << "diag " << std::setw(width) << "srtz "
          << std::setw(width)
          << "srtsamp "
          //<< std::setw(width) << "nystrom "
          //<< std::setw(width) << "cheby "
          //<< std::setw(width) << "reduce "
          << std::endl;
    };

    header(out_it);
    header(out_error);
    header(out_setup);
    header(out_solve);
    header(out_h2_error);
    header(out_mem);
    for (int N = 1000; N < -300000; N *= 2) {
      for (double sigma = 0.001; sigma < -0.1; sigma += 0.1) {
        for (size_t order = 2; order < 8; order += 2) {
          helper_param_sweep<2, KdtreeNanoflann>(
              sigma, N, order, out_it, out_error, out_setup, out_solve,
              out_h2_error, out_mem, 0);
          if (order < 6)
            helper_param_sweep<4, KdtreeNanoflann>(
                sigma, N, order, out_it, out_error, out_setup, out_solve,
                out_h2_error, out_mem, 0);
          // if (order < 4)
          //    helper_param_sweep<6,KdtreeNanoflann>(sigma, N, order, out_it,
          //    out_error,out_setup,out_solve,out_h2_error,out_mem,2);
        }
      }
    }
    out_it.close();
    out_error.close();
    out_setup.close();
    out_solve.close();
    out_h2_error.close();
    out_mem.close();

    out_it.open("iterations_matern.txt", std::ios::out);
    out_error.open("error_matern.txt", std::ios::out);
    out_setup.open("setup_matern.txt", std::ios::out);
    out_solve.open("solve_matern.txt", std::ios::out);
    out_h2_error.open("h2_error_matern.txt", std::ios::out);
    out_mem.open("mem_matern.txt", std::ios::out);

    header(out_it);
    header(out_error);
    header(out_setup);
    header(out_solve);
    header(out_h2_error);
    header(out_mem);
    for (int N = 1000; N < 300000; N *= 2) {
      for (double sigma = 0.1; sigma < 2.0; sigma += 0.2) {
        for (size_t order = 2; order < 8; order += 2) {
          helper_param_sweep<2, KdtreeNanoflann>(
              sigma, N, order, out_it, out_error, out_setup, out_solve,
              out_h2_error, out_mem, 1);
          if (order < 6)
            helper_param_sweep<4, KdtreeNanoflann>(
                sigma, N, order, out_it, out_error, out_setup, out_solve,
                out_h2_error, out_mem, 1);
          // if (order < 4)
          //    helper_param_sweep<6,KdtreeNanoflann>(sigma, N, order, out_it,
          //    out_error,out_setup,out_solve,out_h2_error,out_mem,1);
        }
      }
    }
    out_it.close();
    out_error.close();
    out_setup.close();
    out_solve.close();
    out_h2_error.close();
    out_mem.close();

    out_it.open("iterations_wendland.txt", std::ios::out);
    out_error.open("error_wendland.txt", std::ios::out);
    out_setup.open("setup_wendland.txt", std::ios::out);
    out_solve.open("solve_wendland.txt", std::ios::out);
    out_h2_error.open("h2_error_wendland.txt", std::ios::out);
    out_mem.open("mem_wendland.txt", std::ios::out);

    header(out_it);
    header(out_error);
    header(out_setup);
    header(out_solve);
    header(out_h2_error);
    header(out_mem);
    for (int N = 1000; N < 300000; N *= 2) {
      for (double sigma = 0.1; sigma < 2.0; sigma += 0.2) {
        for (size_t order = 2; order < 8; order += 2) {
          helper_param_sweep<2, KdtreeNanoflann>(
              sigma, N, order, out_it, out_error, out_setup, out_solve,
              out_h2_error, out_mem, 2);
          if (order < 6)
            helper_param_sweep<4, KdtreeNanoflann>(
                sigma, N, order, out_it, out_error, out_setup, out_solve,
                out_h2_error, out_mem, 2);
          // if (order < 4)
          // helper_param_sweep<6,KdtreeNanoflann>(sigma, N, order, out_it,
          // out_error,out_setup,out_solve,out_h2_error,out_mem,2);
        }
      }
    }
    out_it.close();
    out_error.close();
    out_setup.close();
    out_solve.close();
    out_h2_error.close();
    out_mem.close();
  }

  /*
  void n_sweep(void) {
    std::cout << "-------------------------------------------\n"
              << "Running precon n sweep             ....\n"
              << "------------------------------------------" << std::endl;
    std::ofstream out_it;
    out_it.open("iterations_gaussian_n.txt", std::ios::out);
    std::ofstream out_error;
    out_error.open("error_gaussian_n.txt", std::ios::out);
    std::ofstream out_setup;
    out_setup.open("setup_gaussian_n.txt", std::ios::out);
    std::ofstream out_solve;
    out_solve.open("solve_gaussian_n.txt", std::ios::out);

    const int width = 12;
    auto header = [](auto& out) {
        out<< std::setw(width) << "N "
           << std::setw(width) << "sigma "
           << std::setw(width) << "diag "
           << std::setw(width) << "srtsamp "
           << std::setw(width) << "srtz "
           << std::setw(width) << "nystrom "
           << std::setw(width) << "cheby "
           << std::setw(width) << "reduce "
           << std::endl;
    };

    header(out_it);
    header(out_error);
    header(out_setup);
    header(out_solve);
    const double sigma = 0.1;
    for (int N = 1000; N < -1e6; N *= 2) {
      helper_param_sweep<KdtreeNanoflann>(sigma, N, out_it,
  out_error,out_setup,out_solve,0);
    }
    out_it.close();
    out_error.close();
    out_setup.close();
    out_solve.close();

    out_it.open("iterations_matern_n.txt", std::ios::out);
    out_error.open("error_matern_n.txt", std::ios::out);
    out_setup.open("setup_matern_n.txt", std::ios::out);
    out_solve.open("solve_matern_n.txt", std::ios::out);

    header(out_it);
    header(out_error);
    header(out_setup);
    header(out_solve);
    for (int N = 1000; N < 1e6; N *= 2) {
      helper_param_sweep<KdtreeNanoflann>(sigma, N, out_it,
  out_error,out_setup,out_solve,1);
    }
    out_it.close();
    out_error.close();
    out_setup.close();
    out_solve.close();

    out_it.open("iterations_wendland_n.txt", std::ios::out);
    out_error.open("error_wendland_n.txt", std::ios::out);
    out_setup.open("setup_wendland_n.txt", std::ios::out);
    out_solve.open("solve_wendland_n.txt", std::ios::out);

    header(out_it);
    header(out_error);
    header(out_setup);
    header(out_solve);
    for (int N = 1000; N < 1e6; N *= 2) {
      helper_param_sweep<KdtreeNanoflann>(sigma, N, out_it,
  out_error,out_setup,out_solve,2);
    }
    out_it.close();
    out_error.close();
    out_setup.close();
    out_solve.close();

  }
  */

  void test_CellListOrdered() {
    std::cout << "-------------------------------------------\n"
              << "Running tests on CellListOrdered....\n"
              << "------------------------------------------" << std::endl;
    helper_global<CellListOrdered>();
    helper_compact<CellListOrdered>();
  }

  void test_CellList() {
    std::cout << "-------------------------------------------\n"
              << "Running tests on CellList....\n"
              << "------------------------------------------" << std::endl;
    helper_global<CellList>();
    helper_compact<CellList>();
  }

  void test_kdtree() {
    std::cout << "-------------------------------------------\n"
              << "Running tests on kdtree....\n"
              << "------------------------------------------" << std::endl;
    helper_h2<Kdtree>();
    helper_compact<Kdtree>();
  }

  void test_HyperOctree() {
    std::cout << "-------------------------------------------\n"
              << "Running tests on HyperOctree....\n"
              << "------------------------------------------" << std::endl;
    helper_h2<HyperOctree>();
    helper_compact<Kdtree>();
  }
};

#endif /* RBF_INTERPOLATION_TEST_H_ */
       //->
