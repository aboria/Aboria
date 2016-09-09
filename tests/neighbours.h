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


#ifndef NEIGHBOURS_H_
#define NEIGHBOURS_H_

#include <cxxtest/TestSuite.h>

#define LOG_LEVEL 3
#include "Aboria.h"

using namespace Aboria;

class NeighboursTest : public CxxTest::TestSuite {
public:

    template<template <typename,typename> class Vector,template <typename> class SearchMethod>
    void helper_single_particle(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<std::tuple<scalar>,3,Vector,SearchMethod> Test_type;
        typedef position_d<3> position;
    	Test_type test;
    	double3 min(-1,-1,-1);
    	double3 max(1,1,1);
    	double3 periodic(true,true,true);
    	double diameter = 0.1;
    	test.init_neighbour_search(min,max,diameter,periodic);
    	typename Test_type::value_type p;

        get<position>(p) = double3(0,0,0);
    	test.push_back(p);

    	int count = 0;
    	for (auto tpl: test.get_neighbours(double3(diameter/2,diameter/2,0))) {
    		count++;
    	}
    	TS_ASSERT_EQUALS(count,1);

    	auto tpl = test.get_neighbours(double3(diameter/2,diameter/2,0));
    	TS_ASSERT_EQUALS(std::distance(tpl.begin(),tpl.end()),1);

    	tpl = test.get_neighbours(double3(2*diameter,0,0));
    	TS_ASSERT_EQUALS(std::distance(tpl.begin(),tpl.end()),0);
    }

    template<template <typename,typename> class Vector,template <typename> class SearchMethod>
    void helper_two_particles(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<std::tuple<scalar>,3,Vector,SearchMethod> Test_type;
        typedef position_d<3> position;
    	Test_type test;
    	double3 min(-1,-1,-1);
    	double3 max(1,1,1);
    	double3 periodic(true,true,true);
    	double diameter = 0.1;
    	test.init_neighbour_search(min,max,diameter,periodic);
    	typename Test_type::value_type p;

        get<position>(p) = double3(0,0,0);
    	test.push_back(p);

        get<position>(p) = double3(diameter/2,0,0);
    	test.push_back(p);

    	auto tpl = test.get_neighbours(double3(1.1*diameter,0,0));
    	TS_ASSERT_EQUALS(std::distance(tpl.begin(),tpl.end()),1);
    	const typename Test_type::value_type &pfound = tuple_ns::get<0>(*tpl.begin());
    	TS_ASSERT_EQUALS(get<id>(pfound),get<id>(test[1]));

    	tpl = test.get_neighbours(double3(0.9*diameter,0,0));
    	TS_ASSERT_EQUALS(std::distance(tpl.begin(),tpl.end()),2);

    	tpl = test.get_neighbours(double3(1.6*diameter,0,0));
    	TS_ASSERT_EQUALS(std::distance(tpl.begin(),tpl.end()),0);

    	tpl = test.get_neighbours(double3(0.25*diameter,0.99*diameter,0));
    	TS_ASSERT_EQUALS(std::distance(tpl.begin(),tpl.end()),2);

    	tpl = test.get_neighbours(double3(0.25*diameter,1.01*diameter,0));
    	TS_ASSERT_EQUALS(std::distance(tpl.begin(),tpl.end()),0);
    }

    template <typename Search>
    struct has_n_neighbours {
        unsigned int n;
        Search search;

        CUDA_HOST_DEVICE 
        has_n_neighbours(const Search& search, unsigned int n):search(search),n(n) {}

        CUDA_HOST_DEVICE 
        void operator()(typename Search::reference i) {
            auto tpl = search.get_neighbours(get<typename Search::position>(i));
            TS_ASSERT_EQUALS(tpl.end()-tpl.begin(),n);
        }
    };

    template<unsigned int D, 
             template <typename,typename> class VectorType,
             template <typename> class SearchMethod>
    void helper_d(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	typedef Particles<std::tuple<scalar>,D,VectorType,SearchMethod> Test_type;
        typedef position_d<D> position;
        typedef Vector<double,D> double_d;
        typedef Vector<bool,D> bool_d;
        typedef Vector<unsigned int,D> uint_d;
    	Test_type test;
    	double_d min(-1);
    	double_d max(1);
    	bool_d periodic(true);
    	typename Test_type::value_type p;
        uint_d index(0);
        unsigned int n = 10;
        double dx = 2.0/n;

    	double diameter = dx*1.00001;
        unsigned int expect_n = std::pow(3,D);

        bool finished = false;
        while (finished != true) {
            double_d pos = index*dx+min+dx/2;
            get<position>(p) = pos;
            test.push_back(p);
            index[0]++;
            for (int i=0; i<D; i++) {
                if (index[i] >= n) {
                    if (i==D-1) {
                        finished = true;
                    } else {
                        index[i+1]++;
                        index[i] = 0;
                    }
                } else {
                    break;
                }
            }
        }

    	test.init_neighbour_search(min,max,diameter,periodic);
        for_each(test.begin(),test.end(),
                has_n_neighbours<typename Test_type::query_type>(test.get_query(),expect_n));
    }

    void test_std_vector_bucket_search_serial(void) {
        helper_single_particle<std::vector,bucket_search_serial>();
        helper_two_particles<std::vector,bucket_search_serial>();
        helper_d<1,std::vector,bucket_search_serial>();
        helper_d<2,std::vector,bucket_search_serial>();
        helper_d<3,std::vector,bucket_search_serial>();
        helper_d<4,std::vector,bucket_search_serial>();
    }

    void test_std_vector_bucket_search_parallel(void) {
        helper_single_particle<std::vector,bucket_search_parallel>();
        helper_two_particles<std::vector,bucket_search_parallel>();
        helper_d<1,std::vector,bucket_search_parallel>();
        helper_d<2,std::vector,bucket_search_parallel>();
        helper_d<3,std::vector,bucket_search_parallel>();
        helper_d<4,std::vector,bucket_search_parallel>();
    }

    void test_thrust_vector_bucket_search_serial(void) {
#if defined(__CUDACC__)
        helper_single_particle<thrust::device_vector,bucket_search_serial>();
        helper_two_particles<thrust::device_vector,bucket_search_serial>();
        helper_d<1,thrust::device_vector,bucket_search_serial>();
        helper_d<2,thrust::device_vector,bucket_search_serial>();
        helper_d<3,thrust::device_vector,bucket_search_serial>();
        helper_d<4,thrust::device_vector,bucket_search_serial>();
#endif
    }

    void test_thrust_vector_bucket_search_parallel(void) {
#if defined(__CUDACC__)
        helper_single_particle<thrust::device_vector,bucket_search_parallel>();
        helper_two_particles<thrust::device_vector,bucket_search_parallel>();
        helper_d<1,thrust::device_vector,bucket_search_parallel>();
        helper_d<2,thrust::device_vector,bucket_search_parallel>();
        helper_d<3,thrust::device_vector,bucket_search_parallel>();
        helper_d<4,thrust::device_vector,bucket_search_parallel>();
#endif
    }

};



#endif /* NEIGHBOURS_H_ */
