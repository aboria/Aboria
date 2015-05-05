/*
 * constructors.h
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

#ifndef CONSTRUCTORS_H_
#define CONSTRUCTORS_H_

#include <cxxtest/TestSuite.h>

#define LOG_LEVEL 1
#include "Aboria.h"

using namespace Aboria;


class ConstructorsTest : public CxxTest::TestSuite {
public:
    void testOneDouble(void) {
        ABORIA_VARIABLE(scalar,double,"scalar")
    	Particles<scalar> test;
    }

    void testOneVect3d(void) {
        ABORIA_VARIABLE(vector,Vect3d,"vector")
    	Particles<vector> test;
    }

    void testNoData(void) {
    	Particles<> test;
    }

    void testMultiple(void) {
        ABORIA_VARIABLE(vector,Vect3d,"vector")
        ABORIA_VARIABLE(var1,double,"var1")
        ABORIA_VARIABLE(var2,int,"var2")
        ABORIA_VARIABLE(var3,unsigned int,"var3")
        ABORIA_VARIABLE(var4,bool,"var4")
    	Particles<vector,var1,var2,var3,var4> test;
    }

};


#endif /* CONSTRUCTORS_H_ */
