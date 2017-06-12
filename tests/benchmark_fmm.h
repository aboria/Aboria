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


#ifndef BENCHMARK_FMM_H_
#define BENCHMARK_FMM_H_

#include <cxxtest/TestSuite.h>
#ifdef HAVE_EIGEN
#include <Eigen/Core>
#endif
#include "Aboria.h"
#include <chrono>
typedef std::chrono::system_clock Clock;
#include <fstream>      // std::ofstream
#include <thread>
#ifdef HAVE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif


#include <boost/math/constants/constants.hpp>
const double PI = boost::math::constants::pi<double>();

using namespace Aboria;

class BenchmarkFMM: public CxxTest::TestSuite {
public:
    template<unsigned int N, template <typename> class SearchMethod
            ,unsigned int D, typename Kernel>
    double multiquadric_fmm(
            const std::vector<Vector<double,D>>& position_vect, 
            const std::vector<double>& source_vect, 
            const std::vector<double>& target_vect, 
            const Kernel& kernel,
            const size_t repeats) {
        typedef Vector<double,D> double_d;
        typedef Vector<bool,D> bool_d;
        const size_t n = source_vect.size();
        std::cout << "multiquadric_fmm: D = "<<D<<" N = "<<N<<" n = "<<n<<std::endl;
        ABORIA_VARIABLE(source,double,"source")
        ABORIA_VARIABLE(target,double,"target")
    	typedef Particles<std::tuple<source,target>,D,std::vector,SearchMethod> particles_type;
        typedef typename particles_type::position position;
       	particles_type particles(n);

        double_d min(-0.1);
        double_d max(1.1);
        bool_d periodic(false);
        
        for (int i=0; i<n; ++i) {
            get<position>(particles)[i] = position_vect[i];
            get<source>(particles)[i] = source_vect[i];
            get<target>(particles)[i] = 0;
        }

        particles.init_neighbour_search(min,max,periodic);

        auto fmm = make_fmm_query(particles.get_query(),
                make_black_box_expansion<D,N>(kernel));

        fmm.gemv(get<target>(particles),get<source>(particles));

        auto t0 = Clock::now();
#ifdef HAVE_GPERFTOOLS
        ProfilerStart("multiquadric_fmm");
#endif
        for (int ii=0; ii<repeats; ++ii) {
            fmm.gemv(get<target>(particles),get<source>(particles));
        }
#ifdef HAVE_GPERFTOOLS
        ProfilerStop();
#endif
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;

        const double L2 = std::inner_product(
                std::begin(get<target>(particles)), std::end(get<target>(particles)),
                std::begin(target_vect), 
                0.0,
                [](const double t1, const double t2) { return t1 + t2; },
                [](const double t1, const double t2) { return (t1-t2)*(t1-t2); }
                );
        TS_ASSERT_LESS_THAN(L2,1e-2);

        return dt.count()/repeats;
    }

    template<unsigned int D, typename Kernel>
    double multiquadric_vector(
            std::vector<double>& target_vect, 
            const std::vector<Vector<double,D>>& position_vect, 
            const std::vector<double>& source_vect, 
            const Kernel& kernel,
            const size_t repeats) {

        typedef Vector<double,D> double_d;
        const size_t n = position_vect.size();
        target_vect.resize(N);
        std::fill(std::begin(target_vect),std::end(target_vect),0.0);
        std::cout << "multiquadric_vector: D = "<<D<<" N = "<<n<<" repeats = "<<repeats<<std::endl;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                target_vect[i] += kernel(position_vect[j]-position_vect[i],
                                         position_vect[i],
                                         position_vect[j])
                                    *source_vect[j];
            }
        }
        
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    target_vect[i] += kernel(position_vect[j]-position_vect[i],
                                             position_vect[i],
                                             position_vect[j])
                                        *source_vect[j];
                }
            }
        }
        auto t1 = Clock::now();

        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;
        return dt.count()/repeats;
    }

    template <unsigned int D>
    void helper_multiquadric() {
        std::ofstream file;
        const size_t base_repeats = 5e6;
        file.open("benchmark_fmm_" + std::to_string(D) + ".csv");
        file <<"#"<< std::setw(14) << "N" 
             << std::setw(15) << "vector" 
             << std::setw(15) << "fmm-kdtree-N-2" 
             << std::setw(15) << "fmm-octtree-N-2" 
             << std::setw(15) << "fmm-kdtree-N-3" 
             << std::setw(15) << "fmm-octtree-N-3" 
             << std::endl;
        for (double i = 100; i < 10000; i *= 1.1) {
            const size_t N = i;
            // randomly generate a bunch of positions over a range 
            const double pos_min = 0;
            const double pos_max = 1;
            std::uniform_real_distribution<double> U(pos_min,pos_max);
            generator_type generator(time(NULL));
            auto gen = std::bind(U, generator);
            typedef Vector<double,D> double_d;

            std::vector<double_d> positions(N);
            std::vector<double> source(N);
            std::vector<double> target(N);

            // generate a source vector using a smooth cosine
            auto source_fn = [&](const double_d &p) {
                double ret=1.0;
                const double scale = 2.0*detail::PI/(pos_max-pos_min); 
                for (int i=0; i<D; i++) {
                    ret *= cos((p[i]-pos_min)*scale);
                }
                return ret/N;
            };

            for (int i=0; i<N; i++) {
                double_d& p = positions[i];
                for (int d=0; d<D; ++d) {
                    p[d] = gen();
                }
                source[i] = source_fn(p);
            }

            // multiquadric kernel
            const double c = 0.1;
            auto kernel = [&c](const double_d &dx, const double_d &pa, const double_d &pb) {
                return std::sqrt(dx.squaredNorm() + c); 
            };

            //const size_t repeats = base_repeats/std::pow(N,2) + 1;
            const size_t repeats = 1.0;
            file << std::setw(15) << N;
#ifdef HAVE_OPENMP
            omp_set_num_threads(1);
#ifdef HAVE_EIGEN
            Eigen::setNbThreads(1);
#endif
#endif
            file << std::setw(15) << multiquadric_vector(target,positions,source,kernel,repeats);
            file << std::setw(15) << 
                multiquadric_fmm<2,nanoflann_adaptor>(positions,source,
                                                        target,kernel,repeats);
            file << std::setw(15) << 
                multiquadric_fmm<2,octtree>(positions,source,
                                                        target,kernel,repeats);
            file << std::setw(15) << 
                multiquadric_fmm<3,nanoflann_adaptor>(positions,source,
                                                        target,kernel,repeats);
            file << std::setw(15) << 
                multiquadric_fmm<3,octtree>(positions,source,
                                                        target,kernel,repeats);


            file << std::endl;
        }
        file.close();
    }

    void test_multiquadric() {
        helper_multiquadric<2>();
    }


};

#endif 
