/*
 * variables.h
 * 
 * Copyright 2014 Martin Robinson
 *
 * This file is part of Aboria.
 *
 * Aboria is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Aboria is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Aboria.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 27 Nov 2014
 *      Author: robinsonm
 */

#ifndef VARIABLES_H_ 
#define VARIABLES_H_ 

#include <cxxtest/TestSuite.h>

#define LOG_LEVEL 1
#include "Aboria.h"

using namespace Aboria;


class VariablesTest: public CxxTest::TestSuite {
public:
    void testScalar(void) {

        ABORIA_VARIABLE(a,double,"a")
        ABORIA_VARIABLE(b,float,"b")
        ABORIA_VARIABLE(c,int,"c")

        Particles<a,b,c>::value_type p;
        set<a>(p,1.2);
        set<b>(p,1.3);
        set<c>(p,1);

        TS_ASSERT_DELTA(get<a>(p),1.2,std::numeric_limits<double>::epsilon());
        TS_ASSERT_DELTA(get<b>(p),1.3,std::numeric_limits<float>::epsilon());
        TS_ASSERT_EQUALS(get<c>(p),1);

        get<a>(p) = 2.2;
        get<b>(p) = 2.3;
        get<c>(p) = 2;

        TS_ASSERT_DELTA(get<a>(p),2.2,std::numeric_limits<double>::epsilon());
        TS_ASSERT_DELTA(get<b>(p),2.3,std::numeric_limits<float>::epsilon());
        TS_ASSERT_EQUALS(get<c>(p),2);

    }

    void testVector(void) {

        ABORIA_VARIABLE(a,Vect3d,"a")
        ABORIA_VARIABLE(b,Vect3i,"b")

        Particles<a,b>::value_type p;

        set<a>(p,Vect3d(1.1,1.2,1.3));
        set<b>(p,Vect3d(1,2,3));

        std::cout << "get<a>(p) = " << get<a>(p) <<std::endl;
        TS_ASSERT_EQUALS(get<a>(p),Vect3d(1.1,1.2,1.3));
        TS_ASSERT_EQUALS(get<b>(p),Vect3d(1,2,3));

        get<a>(p) = Vect3d(2.1,2.2,2.3);
        get<b>(p) = Vect3d(2,3,4);

        TS_ASSERT_EQUALS(get<a>(p),Vect3d(2.1,2.2,2.3));
        TS_ASSERT_EQUALS(get<b>(p),Vect3d(2,3,4));
    }

};


#endif /* CONSTRUCTORS_H_ */
