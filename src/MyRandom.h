/*
 * MyRandom.h
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of RD_3D.
 *
 * RD_3D is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RD_3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with RD_3D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 12 Oct 2012
 *      Author: robinsonm
 */

#ifndef MYRANDOM_H_
#define MYRANDOM_H_

#include <boost/random.hpp>

namespace Aboria {
//typedef boost::minstd_rand base_generator_type;
typedef boost::mt19937  base_generator_type;
extern base_generator_type generator;
}


#endif /* MYRANDOM_H_ */
