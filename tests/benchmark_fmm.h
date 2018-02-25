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
    double multiquadric_h2(
            const std::vector<Vector<double,D>>& position_vect, 
            const std::vector<double>& source_vect, 
            const std::vector<double>& target_vect, 
            const Kernel& kernel,
            const size_t repeats,
            const size_t nbucket) {
        typedef Vector<double,D> double_d;
        typedef Vector<bool,D> bool_d;
        const size_t n = source_vect.size();
        std::cout << "multiquadric_h2: D = "<<D<<" N = "<<N<<" n = "<<n<<" nbucket = "<<nbucket<<" repeats = "<<repeats<<std::endl;
        ABORIA_VARIABLE(source,double,"source")
        ABORIA_VARIABLE(target,double,"target")
    	typedef Particles<std::tuple<source,target>,D,std::vector,SearchMethod> particles_type;
        typedef typename particles_type::position position;
       	particles_type particles(n);

        double_d min = double_d::Constant(-0.1);
        double_d max = double_d::Constant(1.1);
        bool_d periodic = bool_d::Constant(false);
        
        for (int i=0; i<n; ++i) {
            get<position>(particles)[i] = position_vect[i];
            get<source>(particles)[i] = source_vect[i];
            get<target>(particles)[i] = 0;
        }

        particles.init_neighbour_search(min,max,periodic,nbucket);

        auto kernel2 = [&](typename particles_type::const_reference a,
                           typename particles_type::const_reference b) {
            return kernel(get<position>(a),get<position>(b));
        };

        auto h2 = create_h2_operator(particles,particles,N,
                                     kernel,kernel2);

        typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> map_type; 
        map_type targetv(get<target>(particles).data(),particles.size());
        map_type sourcev(get<source>(particles).data(),particles.size());
        targetv += h2*sourcev;

        auto t0 = Clock::now();
#ifdef HAVE_GPERFTOOLS
        typedef typename particles_type::traits_type traits_type;
        if (std::is_same<SearchMethod<traits_type>,
                Kdtree<traits_type>>::value) {
            ProfilerStart(("multiquadric_h2_kdtree"+std::to_string(N)+"_"+std::to_string(D)).c_str());
        } else {
            ProfilerStart(("multiquadric_h2_octtree"+std::to_string(N)+"_"+std::to_string(D)).c_str());
        }
#endif
        for (int ii=0; ii<repeats; ++ii) {
            targetv += h2*sourcev;
        }
#ifdef HAVE_GPERFTOOLS
        ProfilerStop();
#endif
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;

        /*
        double L2 = 0;
        for (int i=0; i<n; ++i) {
            const double e = get<target>(particles)[i]-target_vect[get<id>(particles)[i]];
            L2 += std::pow(e,2);
        }
            
        TS_ASSERT_LESS_THAN(L2,1e-2);
        */

        return dt.count()/repeats;
    }

    template<unsigned int N, template <typename> class SearchMethod
            ,unsigned int D, typename Kernel>
    double multiquadric_fmm(
            const std::vector<Vector<double,D>>& position_vect, 
            const std::vector<double>& source_vect, 
            const std::vector<double>& target_vect, 
            const Kernel& kernel,
            const size_t repeats,
            const size_t nbucket) {
        typedef Vector<double,D> double_d;
        typedef Vector<bool,D> bool_d;
        const size_t n = source_vect.size();
        std::cout << "multiquadric_fmm: D = "<<D<<" N = "<<N<<" n = "<<n<<" nbucket = "<<nbucket<<" repeats = "<<repeats<<std::endl;
        ABORIA_VARIABLE(source,double,"source")
        ABORIA_VARIABLE(target,double,"target")
    	typedef Particles<std::tuple<source,target>,D,std::vector,SearchMethod> particles_type;
        typedef typename particles_type::position position;
       	particles_type particles(n);

        double_d min = double_d::Constant(-0.1);
        double_d max = double_d::Constant(1.1);
        bool_d periodic = bool_d::Constant(false);
        
        for (int i=0; i<n; ++i) {
            get<position>(particles)[i] = position_vect[i];
            get<source>(particles)[i] = source_vect[i];
            get<target>(particles)[i] = 0;
        }

        particles.init_neighbour_search(min,max,periodic,nbucket);
        auto kernel2 = [&](typename particles_type::const_reference a,
                           typename particles_type::const_reference b) {
            return kernel(get<position>(a),get<position>(b));
        };
        auto fmm = create_fmm_operator<N>(particles,particles,
                                    kernel,kernel2);

        typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> map_type; 
        map_type targetv(get<target>(particles).data(),particles.size());
        map_type sourcev(get<source>(particles).data(),particles.size());
        targetv += fmm*sourcev;

        auto t0 = Clock::now();
#ifdef HAVE_GPERFTOOLS
        typedef typename particles_type::traits_type traits_type;
        if (std::is_same<SearchMethod<traits_type>,
                Kdtree<traits_type>>::value) {
            ProfilerStart(("multiquadric_fmm_kdtree"+std::to_string(N)+"_"+std::to_string(D)).c_str());
        } else {
            ProfilerStart(("multiquadric_fmm_octtree"+std::to_string(N)+"_"+std::to_string(D)).c_str());
        }
#endif
        for (int ii=0; ii<repeats; ++ii) {
            targetv += fmm*sourcev;
        }
#ifdef HAVE_GPERFTOOLS
        ProfilerStop();
#endif
        auto t1 = Clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cout << "time = "<<dt.count()/repeats<<std::endl;

        /*
        double L2 = 0;
        for (int i=0; i<n; ++i) {
            const double e = get<target>(particles)[i]-target_vect[get<id>(particles)[i]];
            L2 += std::pow(e,2);
        }
            
        TS_ASSERT_LESS_THAN(L2,1e-2);
        */

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

        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                target_vect[i] += kernel(position_vect[i],
                                         position_vect[j])
                                    *source_vect[j];
            }
        }
        
        auto t0 = Clock::now();
        for (int r=0; r<repeats; ++r) {
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    target_vect[i] += kernel(position_vect[i],
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
        const size_t nbucket_min = 10;
        const size_t nbucket_max = 11;
        const size_t nbucket_incr = 10;
        const size_t base_repeatsN2 = 1e7;
        const size_t base_repeatsN = 2e5;
        file.open("benchmark_fmm_" + std::to_string(D) + ".csv");
        file <<"#"<< std::setw(25) << "N" 
             << std::setw(25) << "vector";
        for (size_t i = nbucket_min; i < nbucket_max; i += nbucket_incr) {
             file << std::setw(25) << "fmm-kdtree-N-2-nb-"+std::to_string(i) 
                  << std::setw(25) << "fmm-octtree-N-2-nb-"+std::to_string(i) 
                  << std::setw(25) << "h2-kdtree-N-2-nb-"+std::to_string(i) 
                  << std::setw(25) << "h2-octtree-N-2-nb-"+std::to_string(i) 
                  << std::setw(25) << "fmm-kdtree-N-3-nb-"+std::to_string(i) 
                  << std::setw(25) << "fmm-octtree-N-3-nb-"+std::to_string(i) 
                  << std::setw(25) << "h2-kdtree-N-3-nb-"+std::to_string(i) 
                  << std::setw(25) << "h2-octtree-N-3-nb-"+std::to_string(i) 
                  << std::setw(25) << "fmm-kdtree-N-4-nb-"+std::to_string(i) 
                  << std::setw(25) << "fmm-octtree-N-4-nb-"+std::to_string(i)
                  << std::setw(25) << "h2-kdtree-N-4-nb-"+std::to_string(i) 
                  << std::setw(25) << "h2-octtree-N-4-nb-"+std::to_string(i);
        }
        file << std::endl;
        for (double i = 1000; i < 1000000; i *= 1.1) {
        //{
        //    double i = 10000;
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
                for (size_t i = 0; i < D;  ++i) {
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
            auto kernel = [&c](const double_d &pa, const double_d &pb) {
                return std::sqrt((pb-pa).squaredNorm() + c); 
            };

            const size_t repeatsN2 = base_repeatsN2/std::pow(N,2) + 1;
            const size_t repeatsN = base_repeatsN/N + 1;
            file << std::setw(15) << N;
#ifdef HAVE_OPENMP
            omp_set_num_threads(1);
#ifdef HAVE_EIGEN
            Eigen::setNbThreads(1);
#endif
#endif
            file << std::setw(25) << multiquadric_vector(target,positions,source,kernel,repeatsN2);

            for (size_t i = nbucket_min; i < nbucket_max; i += nbucket_incr) {
            //{
            //    size_t i = nbucket_min;
                file << std::setw(25) << 
                    multiquadric_fmm<2,Kdtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_fmm<2,octtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_h2<2,Kdtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_h2<2,octtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_fmm<3,Kdtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_fmm<3,octtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_h2<3,Kdtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_h2<3,octtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_fmm<4,Kdtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_fmm<4,octtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_h2<4,Kdtree>(positions,source,
                            target,kernel,repeatsN,i);
                file << std::setw(25) << 
                    multiquadric_h2<4,octtree>(positions,source,
                            target,kernel,repeatsN,i);

            }

            file << std::endl;
        }
        file.close();
    }

    void test_multiquadric() {
        //helper_multiquadric<1>();
        helper_multiquadric<2>();
        //helper_multiquadric<3>();
        //helper_multiquadric<4>();
    }



};

#endif 
