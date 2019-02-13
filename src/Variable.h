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


#ifndef VARIABLE_H_ 
#define VARIABLE_H_ 

#include <boost/preprocessor/cat.hpp>
#include "Vector.h"
#include "Random.h"
#include <vector>

namespace Aboria {

/// \brief variables are attached to particles and have a name and container type
///
/// \param NAME a type with a \c const char* member variable containing the name of the variable 
/// \param T a type used to contain the data of the variable, e.g. \c int, double.
template<typename T, typename NAME>
struct Variable {
    const char *name = NAME().name;
    typedef T value_type;
};

/// \brief a macro to conveniently define variable types
/// \param NAME the name of the generated type 
/// \param DATA_TYPE the type used to contain the data of the variable, e.g. \c int, double
/// \param NAME_STRING a string used to name or describe the variable, e.g. "scalar", "velocity"
#define ABORIA_VARIABLE(NAME,DATA_TYPE,NAME_STRING)      \
    struct BOOST_PP_CAT(NAME,_description) {                            \
    	const char* name = NAME_STRING; \
    };                                                   \
    typedef Aboria::Variable<DATA_TYPE,BOOST_PP_CAT(NAME,_description)> NAME;   \

#define ABORIA_VARIABLE_VECTOR(NAME,DATA_TYPE,NAME_STRING)      \
    struct BOOST_PP_CAT(NAME,_description) {                            \
    	const char* name = NAME_STRING; \
    };                                                   \
    template <unsigned int N> \
    using NAME = Aboria::Variable<Vector<DATA_TYPE,N>,BOOST_PP_CAT(NAME,_description)>;   \


ABORIA_VARIABLE_VECTOR(position_d,double,"position")
ABORIA_VARIABLE_VECTOR(particles_d,size_t,"particles_id")
ABORIA_VARIABLE(alive,uint8_t,"is_alive")
ABORIA_VARIABLE(id,size_t,"id")
ABORIA_VARIABLE(generator,generator_type,"random_generator_seed")

}
#endif /* VARIABLE_H_ */
