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
    	Particles<std::tuple<double> > test;
    }

    void testOneVect3d(void) {
    	Particles<std::tuple<Vect3d> > test;
    }

    void testNoData(void) {
    	Particles<> test;
    }

    void testMultiple(void) {
    	Particles<std::tuple<Vect3d,double,int,unsigned int, bool> > test;
    }

};


#endif /* CONSTRUCTORS_H_ */
