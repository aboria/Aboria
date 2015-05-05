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

template<typename T, typename NAME>
struct Variable {
    const char *name = NAME::name;
    typedef T value_type;
};

#define ABORIA_VARIABLE(NAME,DATA_TYPE,NAME_STRING)      \
    struct BOOST_PP_CAT(NAME,_description) {                            \
    	const char *name = NAME_STRING; \
    };                                                   \
    typedef Variable<DATA_TYPE,BOOST_PP_CAT(NAME,_description)> NAME;   \

ABORIA_VARIABLE(position,Vect3d,"position")
ABORIA_VARIABLE(alive,bool,"is alive")
ABORIA_VARIABLE(id,size_t,"id")

}
#endif /* VARIABLE_H_ */
