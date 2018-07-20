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

#ifndef RBF_PARAM_SWEEP_H_
#define RBF_PARAM_SWEEP_H_

//#include "hilbert/hilbert.h"
#include <sys/resource.h>
#include <sys/time.h>

#include <chrono>
#include <cxxtest/TestSuite.h>
#include <fstream>
typedef std::chrono::system_clock Clock;

#ifdef HAVE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "Aboria.h"
using namespace Aboria;

class RbfParamSweepTest : public CxxTest::TestSuite {
public:
  ABORIA_VARIABLE(function, double, "function_eval");

  struct matern_kernel {
    double m_scale;
    double m_sigma;
    static constexpr const char *m_name = "matern";
    void set_sigma(const double sigma) {
      m_sigma = sigma;
      m_scale = std::sqrt(3.0) / sigma;
    }
    template <unsigned int D>
    CUDA_HOST_DEVICE double operator()(const Vector<double, D> &a,
                                       const Vector<double, D> &b) const {
      const double r = (b - a).norm();
      return (1.0 + m_scale * r) * std::exp(-r * m_scale);
    };
  };

  struct gaussian_kernel {
    double m_scale;
    double m_sigma;
    static constexpr const char *m_name = "gaussian";
    void set_sigma(const double sigma) {
      m_sigma = sigma;
      m_scale = 1.0 / std::pow(sigma, 2);
    }
    template <unsigned int D>
    CUDA_HOST_DEVICE double operator()(const Vector<double, D> &a,
                                       const Vector<double, D> &b) const {
      return std::exp(-(b - a).squaredNorm() * m_scale);
    };
  };

  struct exponential_kernel {
    double m_scale;
    double m_sigma;
    static constexpr const char *m_name = "exp";
    void set_sigma(const double sigma) {
      m_sigma = sigma;
      m_scale = 1.0 / sigma;
    }
    template <unsigned int D>
    CUDA_HOST_DEVICE double operator()(const Vector<double, D> &a,
                                       const Vector<double, D> &b) const {
      return std::exp(-(b - a).norm() * m_scale);
    };
  };

  struct rational_quadratic_kernel {
    double m_scale;
    double m_sigma;
    static constexpr const char *m_name = "rational";
    void set_sigma(const double sigma) {
      m_sigma = sigma;
      m_scale = std::pow(sigma, 2);
    }
    template <unsigned int D>
    CUDA_HOST_DEVICE double operator()(const Vector<double, D> &a,
                                       const Vector<double, D> &b) const {
      const double r2 = (b - a).squaredNorm();
      return 1.0 - r2 / (r2 + m_scale);
    };
  };

  struct inverse_multiquadric_kernel {
    double m_scale;
    double m_sigma;
    static constexpr const char *m_name = "inv_mq";
    void set_sigma(const double sigma) {
      m_sigma = sigma;
      m_scale = std::pow(sigma, 2);
    }
    template <unsigned int D>
    CUDA_HOST_DEVICE double operator()(const Vector<double, D> &a,
                                       const Vector<double, D> &b) const {
      const double r2 = (b - a).squaredNorm();
      return 1.0 / std::sqrt(r2 + m_scale);
    };
  };

  struct output_files {
    std::ofstream out_op_setup_time;
    std::ofstream out_op_apply_time;
    std::ofstream out_op_apply_error;
    std::ofstream out_op_memory;
    std::ofstream out_ds_setup_time;

    const int width = 12;
    static const int Nops = 2;
    std::ofstream out_solve_setup_time[Nops];
    std::ofstream out_solve_solve_time[Nops];
    std::ofstream out_solve_solve_memory[Nops];
    std::ofstream out_solve_iterations[Nops];
    std::ofstream out_solve_error[Nops];
    std::ofstream out_solve_test_error[Nops];

    output_files(const std::string &name) {

      auto solve_header = [width = width + 1](auto &out) {
        out << std::setw(width) << "N " << std::setw(width) << "sigma "
            << std::setw(width) << "D " << std::setw(width) << "order "
            << std::setw(width) << "Nsubdomain " << std::setw(width) << "chol "
            << std::setw(width) << "diag " << std::setw(width) << "srtz "
            << std::setw(width)
            << "nystrom "
            //<< std::setw(width) << "nys_srtz "
            << std::endl;
      };

      auto op_header = [width = width + 1](auto &out) {
        out << std::setw(width) << "N " << std::setw(width) << "sigma "
            << std::setw(width) << "D " << std::setw(width) << "order "
            << std::setw(width) << "Nsubdomain " << std::setw(width)
            << "matrix " << std::setw(width) << "dense " << std::setw(width)
            << "fmm " << std::setw(width) << "h2 " << std::endl;
      };

      auto ds_header = [width = width + 1](auto &out) {
        out << std::setw(width) << "N " << std::setw(width) << "sigma "
            << std::setw(width) << "D " << std::setw(width) << "order "
            << std::setw(width) << "Nsubdomain " << std::setw(width) << "kdree "
            << std::setw(width) << std::endl;
      };

      out_op_setup_time.open(name + "_op_setup_time.txt", std::ios::out);
      op_header(out_op_setup_time);
      out_op_apply_time.open(name + "_op_apply_time.txt", std::ios::out);
      op_header(out_op_apply_time);
      out_op_apply_error.open(name + "_op_apply_error.txt", std::ios::out);
      op_header(out_op_apply_error);
      out_op_memory.open(name + "_op_memory.txt", std::ios::out);
      op_header(out_op_memory);

      out_ds_setup_time.open(name + "_ds_setup_time.txt", std::ios::out);
      ds_header(out_ds_setup_time);

      std::string op_names[2] = {"_matrix_", "_h2_"};
      for (int i = 0; i < Nops; ++i) {
        out_solve_setup_time[i].open(
            name + op_names[i] + "solve_setup_time.txt", std::ios::out);
        solve_header(out_solve_setup_time[i]);
        out_solve_solve_time[i].open(
            name + op_names[i] + "solve_solve_time.txt", std::ios::out);
        solve_header(out_solve_solve_time[i]);
        out_solve_solve_memory[i].open(
            name + op_names[i] + "solve_solve_memory.txt", std::ios::out);
        solve_header(out_solve_solve_memory[i]);

        out_solve_iterations[i].open(
            name + op_names[i] + "solve_iterations.txt", std::ios::out);
        solve_header(out_solve_iterations[i]);
        out_solve_error[i].open(name + op_names[i] + "solve_error.txt",
                                std::ios::out);
        solve_header(out_solve_error[i]);
        out_solve_test_error[i].open(
            name + op_names[i] + "solve_test_error.txt", std::ios::out);
        solve_header(out_solve_test_error[i]);
      }
    }

    void new_line_end() {
      out_op_setup_time << std::endl;
      out_op_apply_time << std::endl;
      out_op_apply_error << std::endl;
      out_op_memory << std::endl;
      out_ds_setup_time << std::endl;
      for (int i = 0; i < Nops; ++i) {
        out_solve_setup_time[i] << std::endl;
        out_solve_solve_time[i] << std::endl;
        out_solve_solve_memory[i] << std::endl;
        out_solve_iterations[i] << std::endl;
        out_solve_error[i] << std::endl;
        out_solve_test_error[i] << std::endl;
      }
    }

    void new_line_start(const size_t N, const double sigma, const size_t D,
                        const size_t order, const size_t Nsubdomain) {
      auto start = [&](auto &out) {
        out << std::setw(width) << N << " " << std::setw(width) << sigma << " "
            << std::setw(width) << D << " " << std::setw(width) << order << " "
            << std::setw(width) << Nsubdomain;
      };

      start(out_op_setup_time);
      start(out_op_apply_time);
      start(out_op_apply_error);
      start(out_op_memory);
      start(out_ds_setup_time);
      for (int i = 0; i < Nops; ++i) {
        start(out_solve_setup_time[i]);
        start(out_solve_solve_time[i]);
        start(out_solve_solve_memory[i]);
        start(out_solve_iterations[i]);
        start(out_solve_error[i]);
        start(out_solve_test_error[i]);
      }
    }

    void close() {
      out_op_setup_time.close();
      out_op_apply_time.close();
      out_op_apply_error.close();
      out_op_memory.close();
      out_ds_setup_time.close();
      for (int i = 0; i < Nops; ++i) {
        out_solve_setup_time[i].close();
        out_solve_solve_time[i].close();
        out_solve_solve_memory[i].close();
        out_solve_iterations[i].close();
        out_solve_error[i].close();
        out_solve_test_error[i].close();
      }
    }
  };

  template <unsigned int D> auto rosenbrock(size_t N, size_t Ntest) {
    typedef Particles<std::tuple<function>, D, std::vector, KdtreeNanoflann>
        Particles_t;
    typedef position_d<D> position;

    Particles_t knots(N + Ntest);

    auto funct = [](const auto &x) {
      double ret = 0;
      // rosenbrock functions all D
      for (size_t i = 0; i < D - 1; ++i) {
        ret += 100 * std::pow(x[i + 1] - std::pow(x[i], 2), 2) +
               std::pow(1 - x[i], 2);
      }
      return ret;
    };

    auto funct1d = [](const auto &x) { return std::pow(1.0 - x[0], 2); };

    std::default_random_engine generator(0);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (size_t i = 0; i < N + Ntest; ++i) {
      for (size_t d = 0; d < D; ++d) {
        get<position>(knots)[i][d] = distribution(generator);
      }
    }
    if (D == 1) {
      for (size_t i = 0; i < N + Ntest; ++i) {
        get<function>(knots)[i] = funct1d(get<position>(knots)[i]);
      }
    } else {
      for (size_t i = 0; i < N + Ntest; ++i) {
        get<function>(knots)[i] = funct(get<position>(knots)[i]);
      }
    }
    return knots;
  }

  template <typename Op, typename TrueOp, typename TestOp>
  void helper_operator_matrix(Op &G, const TrueOp &Gtrue,
                              const Eigen::VectorXd &phi, TestOp &Gtest,
                              const Eigen::VectorXd &phi_test,
                              output_files &out, bool do_solve) {

    if (do_solve) {
      std::cout << "SOLVING CHOL" << std::endl;
      auto t0 = Clock::now();
      auto solver = G.llt();
      auto t1 = Clock::now();
      Eigen::VectorXd gamma = solver.solve(phi);
      auto t2 = Clock::now();
      const double solve_error = (G * gamma - phi).norm() / phi.norm();

      for (int i = 0; i < 2; ++i) {
        out.out_solve_iterations[i] << " " << std::setw(out.width) << 1;
        out.out_solve_error[i] << " " << std::setw(out.width) << solve_error;
        out.out_solve_setup_time[i]
            << " " << std::setw(out.width)
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                   .count();
        out.out_solve_solve_time[i]
            << " " << std::setw(out.width)
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
                   .count();
        out.out_solve_solve_memory[i]
            << " " << std::setw(out.width)
            << G.rows() * G.cols() * sizeof(double) / 1e9;

        Eigen::VectorXd phi_proposed = Gtest * gamma;

        out.out_solve_test_error[i]
            << " " << std::setw(out.width)
            << (phi_proposed - phi_test).norm() / phi_test.norm();
      }
    } else {
      for (int i = 0; i < 2; ++i) {
        out.out_solve_iterations[i] << " " << std::setw(out.width) << -1;
        out.out_solve_error[i] << " " << std::setw(out.width) << -1;
        out.out_solve_setup_time[i] << " " << std::setw(out.width) << -1;
        out.out_solve_solve_time[i] << " " << std::setw(out.width) << -1;
        out.out_solve_solve_memory[i] << " " << std::setw(out.width) << -1;
        out.out_solve_test_error[i] << " " << std::setw(out.width) << -1;
      }
    }
  }

  template <typename Op, typename TrueOp, typename TestOp>
  void helper_operator(const Op &G, const TrueOp &Gtrue,
                       const Eigen::VectorXd &phi, const TestOp &Gtest,
                       const Eigen::VectorXd &phi_test, const int max_iter,
                       const size_t Nbuffer, output_files &out, int do_solve) {

    const size_t N = G.cols();
    Eigen::VectorXd gamma = Eigen::VectorXd::Random(N);
    auto t0 = Clock::now();
    Eigen::VectorXd phi_random = G * gamma;
    auto t1 = Clock::now();

    std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    Eigen::VectorXd phi_random_true = Gtrue * gamma;

    out.out_op_apply_error << " " << std::setw(out.width)
                           << (phi_random - phi_random_true).norm() /
                                  phi_random_true.norm();
    out.out_op_apply_time
        << " " << std::setw(out.width)
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
               .count();

    if (do_solve) {

      std::cout << "SOLVING IDENTITY" << std::endl;
      Eigen::BiCGSTAB<Op, Eigen::IdentityPreconditioner> bicg;

      bicg.setMaxIterations(max_iter);
      // bicg.setTolerance(X);
      t0 = Clock::now();
      bicg.compute(G);
      t1 = Clock::now();
      Eigen::VectorXd gamma = bicg.solve(phi);
      auto t2 = Clock::now();
      out.out_solve_iterations[do_solve - 1] << " " << std::setw(out.width)
                                             << bicg.iterations();
      out.out_solve_error[do_solve - 1] << " " << std::setw(out.width)
                                        << bicg.error();
      out.out_solve_setup_time[do_solve - 1]
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                 .count();
      out.out_solve_solve_time[do_solve - 1]
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
                 .count();
      out.out_solve_solve_memory[do_solve - 1]
          << " " << std::setw(out.width) << G.rows() * sizeof(double) / 1e9;

      Eigen::VectorXd phi_proposed = Gtest * gamma;

      out.out_solve_test_error[do_solve - 1]
          << " " << std::setw(out.width)
          << (phi_proposed - phi_test).norm() / phi_test.norm();

      std::cout << "SOLVING Schwartz" << std::endl;
      Eigen::BiCGSTAB<Op, SchwartzPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
          bicg2;
      bicg2.setMaxIterations(max_iter);
      bicg2.preconditioner().set_max_buffer_n(Nbuffer);
      // bicg2.preconditioner().set_decimate_factor(2);
      t0 = Clock::now();
      bicg2.compute(G);

#ifdef HAVE_GPERFTOOLS
      ProfilerStart("schwartz_solve");
#endif

      t1 = Clock::now();
      gamma = bicg2.solve(phi);
      t2 = Clock::now();

#ifdef HAVE_GPERFTOOLS
      ProfilerStop();
#endif

      out.out_solve_iterations[do_solve - 1] << " " << std::setw(out.width)
                                             << bicg2.iterations();
      out.out_solve_error[do_solve - 1] << " " << std::setw(out.width)
                                        << bicg2.error();
      out.out_solve_setup_time[do_solve - 1]
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                 .count();
      out.out_solve_solve_time[do_solve - 1]
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
                 .count();
      auto &knots = G.get_first_kernel().get_row_elements();
      const size_t ndomains = knots.size() / knots.get_max_bucket_size();
      const size_t domain_size = 0.75 * knots.get_max_bucket_size() + Nbuffer;
      out.out_solve_solve_memory[do_solve - 1]
          << " " << std::setw(out.width)
          << ndomains * std::pow(domain_size, 2) * sizeof(double) / 1e9;

      phi_proposed = Gtest * gamma;

      out.out_solve_test_error[do_solve - 1]
          << " " << std::setw(out.width)
          << (phi_proposed - phi_test).norm() / phi_test.norm();

      std::cout << "SOLVING Nystrom" << std::endl;
      Eigen::BiCGSTAB<Op, NystromPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
          bicg3;
      bicg3.setMaxIterations(max_iter);
      const size_t Ninducing =
          0.5 * (-2.0 * knots.size() +
                 std::sqrt(4.0 * std::pow(knots.size(), 2) +
                           4.0 * ndomains * std::pow(domain_size, 2)));
      bicg3.preconditioner().set_number_of_random_particles(Ninducing);
      bicg3.preconditioner().set_lambda(1e-5);
      t0 = Clock::now();
      bicg3.compute(G);
      t1 = Clock::now();
      gamma = bicg3.solve(phi);
      t2 = Clock::now();

      out.out_solve_iterations[do_solve - 1] << " " << std::setw(out.width)
                                             << bicg3.iterations();
      out.out_solve_error[do_solve - 1] << " " << std::setw(out.width)
                                        << bicg3.error();
      out.out_solve_setup_time[do_solve - 1]
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                 .count();
      out.out_solve_solve_time[do_solve - 1]
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
                 .count();
      out.out_solve_solve_memory[do_solve - 1]
          << " " << std::setw(out.width)
          << (Ninducing * G.rows() + std::pow(Ninducing, 2)) * sizeof(double) /
                 1e9;

      phi_proposed = Gtest * gamma;

      out.out_solve_test_error[do_solve - 1]
          << " " << std::setw(out.width)
          << (phi_proposed - phi_test).norm() / phi_test.norm();

      /*
      std::cout << "SOLVING Nystrom Swartz" << std::endl;
      Eigen::BiCGSTAB<Op,
                      NystromSwartzPreconditioner<Eigen::LLT<Eigen::MatrixXd>>>
          bicg4;
      bicg4.setMaxIterations(max_iter);
      bicg4.preconditioner().nystrom().set_number_of_random_particles(
          Ninducing);
      bicg4.preconditioner().nystrom().set_lambda(1e-5);
      bicg4.preconditioner().swartz().set_max_buffer_n(Nbuffer);
      t0 = Clock::now();
      bicg4.compute(G);
      t1 = Clock::now();
      gamma = bicg4.solve(phi);
      t2 = Clock::now();

      out.out_solve_iterations[do_solve - 1] << " " << std::setw(out.width)
                                             << bicg4.iterations();
      out.out_solve_error[do_solve - 1] << " " << std::setw(out.width)
                                        << bicg4.error();
      out.out_solve_setup_time[do_solve - 1]
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                 .count();
      out.out_solve_solve_time[do_solve - 1]
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1)
                 .count();
      out.out_solve_solve_memory[do_solve - 1]
          << " " << std::setw(out.width)
          << (Ninducing * G.rows() + std::pow(Ninducing, 2)) * sizeof(double) /
                     1e9 +
                 ndomains * std::pow(domain_size, 2) * sizeof(double) / 1e9;

      phi_proposed = Gtest * gamma;

      out.out_solve_test_error[do_solve - 1]
          << " " << std::setw(out.width)
          << (phi_proposed - phi_test).norm() / phi_test.norm();
          */
    }
  }

  template <size_t Order, typename Particles_t, typename Kernel>
  void helper_param_sweep(const Particles_t &particles, const size_t Ntest,
                          const Kernel &kernel, const double jitter,
                          const size_t Nsubdomain, output_files &out) {
#ifdef HAVE_EIGEN
    typedef typename Particles_t::raw_const_reference raw_const_reference;

    char name[50] = "program name";
    char *argv[] = {name, NULL};
    int argc = sizeof(argv) / sizeof(char *) - 1;
    char **argv2 = &argv[0];
    init_h2lib(&argc, &argv2);

    out.new_line_start(particles.size() - Ntest, kernel.m_sigma,
                       Particles_t::dimension, Order, Nsubdomain);

    const int max_iter = 1000;
    const unsigned int D = Particles_t::dimension;
    const size_t n_subdomain = Nsubdomain;
    const int Nbuffer = 4 * n_subdomain;
    typedef position_d<D> position;

    Particles_t knots(particles.size() - Ntest);
    Particles_t test(Ntest);

    bbox<D> knots_box;
    for (size_t i = 0; i < knots.size(); ++i) {
      knots[i] = particles[i];
      knots_box = knots_box + bbox<D>(get<position>(knots)[i]);
    }

    bbox<D> test_box;
    for (size_t i = 0; i < test.size(); ++i) {
      test[i] = particles[knots.size() + i];
      test_box = test_box + bbox<D>(get<position>(test)[i]);
    }

    auto t0 = Clock::now();
    knots.init_neighbour_search(knots_box.bmin, knots_box.bmax,
                                Vector<bool, D>::Constant(false), n_subdomain);
    auto t1 = Clock::now();
    out.out_ds_setup_time
        << " " << std::setw(out.width)
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
               .count();

    test.init_neighbour_search(test_box.bmin, test_box.bmax,
                               Vector<bool, D>::Constant(false), n_subdomain);
    std::cout << "FINISHED INIT NEIGHBOUR" << std::endl;

    auto self_kernel = [=] CUDA_HOST_DEVICE(raw_const_reference a,
                                            raw_const_reference b) {
      double ret = kernel(get<position>(a), get<position>(b));
      if (get<id>(a) == get<id>(b)) {
        ret += jitter;
      }
      return ret;
    };

    typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> map_type;
    map_type phi(get<function>(knots).data(), knots.size());
    map_type phi_test(get<function>(test).data(), test.size());

    t0 = Clock::now();
    auto Gmatrix = create_matrix_operator(knots, knots, self_kernel);
    t1 = Clock::now();
    out.out_op_setup_time
        << " " << std::setw(out.width)
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
               .count();
    out.out_op_memory << " " << std::setw(out.width)
                      << std::pow(knots.size(), 2) * sizeof(double) / 1e9;
    auto Gmatrix_test = create_matrix_operator(test, knots, self_kernel);

    std::cout << "APPLYING MATRIX OPERATOR" << std::endl;
    helper_operator_matrix(Gmatrix.get_first_kernel().get_matrix(), Gmatrix,
                           phi, Gmatrix_test, phi_test, out, true);
    helper_operator(Gmatrix, Gmatrix, phi, Gmatrix_test, phi_test, max_iter,
                    Nbuffer, out, 1);

    std::cout << "CREATING DENSE OPERATOR" << std::endl;
    t0 = Clock::now();
    auto Gdense = create_dense_operator(knots, knots, self_kernel);
    t1 = Clock::now();
    out.out_op_setup_time
        << " " << std::setw(out.width)
        << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
               .count();
    out.out_op_memory << " " << std::setw(out.width) << 0;

    auto Gdense_test = create_dense_operator(test, knots, self_kernel);

    std::cout << "APPLYING DENSE OPERATOR" << std::endl;
    helper_operator(Gdense, Gmatrix, phi, Gdense_test, phi_test, max_iter,
                    Nbuffer, out, 0);

    if (D < 1) {
      const size_t n_subdomain = std::pow(Order, D);
      const int Nbuffer = 4 * n_subdomain;

      knots.init_neighbour_search(knots_box.bmin, knots_box.bmax,
                                  Vector<bool, D>::Constant(false),
                                  n_subdomain);
      test.init_neighbour_search(test_box.bmin, test_box.bmax,
                                 Vector<bool, D>::Constant(false), n_subdomain);
      std::cout << "FINISHED FMM and H2 INIT NEIGHBOUR" << std::endl;

      // need to recreate Gmatrix as order has changed
      auto Gmatrix = create_matrix_operator(knots, knots, self_kernel);

      // need to remap phi and phi_test as the order has probably changed
      map_type phi(get<function>(knots).data(), knots.size());
      map_type phi_test(get<function>(test).data(), test.size());

      std::cout << "CREATING FMM OPERATOR" << std::endl;
      t0 = Clock::now();
      auto G_FMM =
          create_fmm_operator<Order>(knots, knots, kernel, self_kernel);
      t1 = Clock::now();
      out.out_op_setup_time
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                 .count();
      out.out_op_memory << " " << std::setw(out.width) << 0;

      auto G_FMM_test =
          create_fmm_operator<Order>(test, knots, kernel, self_kernel);

      std::cout << "APPLYING FMM OPERATOR" << std::endl;
      helper_operator(G_FMM, Gmatrix, phi, G_FMM_test, phi_test, max_iter,
                      Nbuffer, out, 0);

      std::cout << "CREATING H2 OPERATOR" << std::endl;
      const double eta = 2.0;
      const double beta = 1.0; // std::max(2.0 / D, 1.0);
      t0 = Clock::now();
      auto G_H2 = create_h2_operator(knots, knots, Order, kernel, self_kernel,
                                     eta, beta);
      t1 = Clock::now();
      out.out_op_setup_time
          << " " << std::setw(out.width)
          << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                 .count();

      size_t size = G_H2.get_first_kernel().get_h2_matrix().get_size();
      /*
      auto near_size =
          GmaternH2.get_first_kernel().get_h2_matrix().get_near_size();
          out.out_op_memory << " " << std::setw(out.width) << 0;
          */
      out.out_op_memory << " " << std::setw(out.width) << size / 1e9;

      auto G_H2_test = create_h2_operator(test, knots, Order, kernel,
                                          self_kernel, eta, beta);

      helper_operator(G_H2, Gmatrix, phi, G_H2_test, phi_test, max_iter,
                      Nbuffer, out, 2);
    }

    out.new_line_end();

#endif // HAVE_EIGEN
  }

  template <typename Kernel>
  void helper_param_sweep_per_kernel(const int minN) {
    Kernel kernel;

    std::cout << "-------------------------------------------\n"
              << "Running precon param sweep with kernel = " << kernel.m_name
              << "....\n"
              << "------------------------------------------" << std::endl;

    output_files out(kernel.m_name);

    const size_t Ntest = 1000;
    const double jitter = 1e-5;

    for (int N = 8000; N < 30000; N *= 2) {
      for (double sigma = 0.9; sigma < 2.0; sigma += 0.4) {
        kernel.set_sigma(sigma);
        for (size_t n_subdomain = 50; n_subdomain < 400; n_subdomain += 100) {
          /*
          helper_param_sweep<2>(rosenbrock<14>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
          helper_param_sweep<2>(rosenbrock<10>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
          helper_param_sweep<2>(rosenbrock<8>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
          helper_param_sweep<3>(rosenbrock<5>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
          helper_param_sweep<4>(rosenbrock<4>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
                                */
          helper_param_sweep<6>(rosenbrock<3>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
          /*
          helper_param_sweep<5>(rosenbrock<3>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
          helper_param_sweep<4>(rosenbrock<3>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
          helper_param_sweep<10>(rosenbrock<2>(N, Ntest), Ntest, kernel, jitter,
                                 n_subdomain, out);
          helper_param_sweep<8>(rosenbrock<2>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
          helper_param_sweep<6>(rosenbrock<2>(N, Ntest), Ntest, kernel, jitter,
                                n_subdomain, out);
          helper_param_sweep<10>(rosenbrock<1>(N, Ntest), Ntest, kernel, jitter,
                                 n_subdomain, out);
                                 */
        }
      }
    }
  }

  void test_gaussian(void) {
    helper_param_sweep_per_kernel<gaussian_kernel>(16000);
  }
  void test_matern(void) {
    helper_param_sweep_per_kernel<matern_kernel>(32000);
  }
  void test_exponential(void) {
    helper_param_sweep_per_kernel<exponential_kernel>(32000);
  }
  void test_rational_quadratic(void) {
    helper_param_sweep_per_kernel<rational_quadratic_kernel>(1000);
  }
  void test_inverse_multiquadric(void) {
    helper_param_sweep_per_kernel<inverse_multiquadric_kernel>(1000);
  }
};

#endif /* RBF_PARAM_SWEEP_H_ */
