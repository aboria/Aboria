/* 
 * Tyche.cpp
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of Tyche.
 *
 * Tyche is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Tyche is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Tyche.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Mar 28, 2013
 *      Author: mrobins
 */

#include <time.h>
#include "Tyche.h"

namespace Tyche {

base_generator_type generator;

void init(int argc, char *argv[]) {
	generator.seed(time(NULL));
}

void random_seed(unsigned int seed) {
  generator.seed(seed);
}

}
