/*
 * RD_3D.h
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
 *  Created on: 30 Oct 2012
 *      Author: robinsonm
 */

#ifndef TYCHE_H_
#define TYCHE_H_

#include "Species.h"
#include "Diffusion.h"
#include "NextSubvolumeMethod.h"
#include "MyRandom.h"
#include "Boundary.h"
#include "Geometry.h"
#include "Reaction.h"
#include "ReactionEquation.h"
#include "Run.h"
#include "Io.h"
#include "Vector.h"
#include "Union.h"
#include "Output.h"
#include "Control.h"
#include "Visualisation.h"

namespace Tyche {

void init(int argc, char *argv[]);
void random_seed(unsigned int);

}



#endif /* TYCHE_H_ */
