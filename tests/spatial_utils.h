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


#ifndef SPATIAL_UTILS_H_
#define SPATIAL_UTILS_H_

#include <cxxtest/TestSuite.h>

#include "Aboria.h"

using namespace Aboria;


class SpatialUtilsTest : public CxxTest::TestSuite {
public:
    void test_bucket_indicies(void) {
        typedef Vector<unsigned int,3> vect;
        detail::bucket_index<3> bi(vect(4,7,2));
        unsigned int index = bi.collapse_index_vector(vect(1,2,1));
        // index = 1 + 2*2 + 1*2*7 = 19 
    	TS_ASSERT_EQUALS(index,19);
        vect vindex = bi.reassemble_index_vector(index);
    	TS_ASSERT_EQUALS(vindex[0],1);
    	TS_ASSERT_EQUALS(vindex[1],2);
    	TS_ASSERT_EQUALS(vindex[2],1);

        typedef Vector<unsigned int,4> vect2;
        detail::bucket_index<4> bi2(vect2(4,7,1,6));
        index = bi2.collapse_index_vector(vect2(1,2,0,4));
        // index = 4 + 0*6 + 2*6*1 + 1*6*1*7 = 58
    	TS_ASSERT_EQUALS(index,58);
        vect2 vindex2 = bi2.reassemble_index_vector(index);
    	TS_ASSERT_EQUALS(vindex2[0],1);
    	TS_ASSERT_EQUALS(vindex2[1],2);
    	TS_ASSERT_EQUALS(vindex2[2],0);
    	TS_ASSERT_EQUALS(vindex2[3],4);

    }


};


#endif /* CONSTRUCTORS_H_ */
