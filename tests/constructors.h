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


#ifndef CONSTRUCTORS_H_
#define CONSTRUCTORS_H_

#include <cxxtest/TestSuite.h>

#include "Level1.h"

using namespace Aboria;


class ConstructorsTest : public CxxTest::TestSuite {
public:
    template<template <typename,typename> class V>
    void helper_OneDouble(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
        typedef std::tuple<scalar> variables_type;
    	Particles<variables_type,3,V> test;
    }

    template<template <typename,typename> class V>
    void helper_OneVect3d(void) {
        ABORIA_VARIABLE(vector,vdouble3,"vector")
        typedef std::tuple<vector> variables_type;
    	Particles<variables_type,3,V> test;
    }

    template<template <typename,typename> class V>
    void helper_NoData(void) {
    	Particles<> test;
    }

    template<template <typename,typename> class V>
    void helper_Dimension(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
        typedef std::tuple<scalar> variable_type;
    	Particles<variable_type,8,V> test8d;
    	Particles<variable_type,7,V> test7d;
    	Particles<variable_type,6,V> test6d;
    	Particles<variable_type,5,V> test5d;
    	Particles<variable_type,4,V> test4d;
    	Particles<variable_type,3,V> test3d;
    	Particles<variable_type,2,V> test2d;
    	Particles<variable_type,1,V> test1d;
    }

    template<template <typename,typename> class V>
    void helper_Multiple(void) {
        ABORIA_VARIABLE(vector,vdouble3,"vector")
        ABORIA_VARIABLE(var1,double,"var1")
        ABORIA_VARIABLE(var2,int,"var2")
        ABORIA_VARIABLE(var3,unsigned int,"var3")
        ABORIA_VARIABLE(var4,bool,"var4")
        typedef std::tuple<vector,var1,var2,var3,var4> variables_type;
    	Particles<variables_type,3,V> test;
    }

    void test_std_vector(void) {
        helper_OneDouble<std::vector>();
        helper_OneVect3d<std::vector>();
        helper_NoData<std::vector>();
        helper_Dimension<std::vector>();
        helper_Multiple<std::vector>();
    }

    void test_thrust_vector(void) {
#if defined(__CUDACC__)
        helper_OneDouble<thrust::device_vector>();
        helper_OneVect3d<thrust::device_vector>();
        helper_NoData<thrust::device_vector>();
        helper_Dimension<thrust::device_vector>();
        helper_Multiple<thrust::device_vector>();
#endif
    }


};


#endif /* CONSTRUCTORS_H_ */
