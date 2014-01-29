/*
 * vector.h   
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
 *  Created on: 18 Oct 2012
 *      Author: robinsonm
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include "Eigen/Dense"

namespace Tyche {
typedef Eigen::Vector3d Vect3d;
typedef Eigen::Vector3i Vect3i;
typedef Eigen::Matrix<int, 6, 1> Vect6i;
typedef Eigen::Matrix<bool, 3, 1> Vect3b;
typedef Eigen::Vector2d Vect2d;
typedef Eigen::Array2i Array2i;
typedef Eigen::Array3Xf Array3Xf;

template<typename T,int N>
std::ostream& operator<< (std::ostream& out, const Eigen::Matrix<T,N,1>& v) {
	out << "(";
	for (int i = 0; i < N; ++i) {
		out << v[i];
		if (i != N-1) out << ",";
	}
	return out << ")";
}
}
#endif /* VECTOR_H_ */
