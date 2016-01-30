/*
 * Variable.h
 *
 *  Created on: 23 Feb 2015
 *      Author: robinsonm
 */

#ifndef VARIABLE_H_ 
#define VARIABLE_H_ 

#include <boost/preprocessor/cat.hpp>

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
    typedef Variable<DATA_TYPE,BOOST_PP_CAT(NAME,_description)> NAME;   \

ABORIA_VARIABLE(position,Vect3d,"position")
ABORIA_VARIABLE(alive,bool,"is alive")
ABORIA_VARIABLE(id,size_t,"id")

}
#endif /* VARIABLE_H_ */
