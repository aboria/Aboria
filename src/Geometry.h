/*
 * Geometry.h
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
 *  Created on: 19 Oct 2012
 *      Author: robinsonm
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "Vector.h"
#include "Log.h"


namespace Aboria {

const double GEOMETRY_TOLERANCE = 1.0/100000.0;

class Geometry {
public:
	virtual bool is_in(const Vect3d &point) const = 0;
	virtual std::pair<double,Vect3d> lineXsurface(const Vect3d &p1, const Vect3d &p2) const = 0;
	virtual const Vect3d shortest_vector_to_boundary(const Vect3d &point) const = 0;
	friend std::ostream& operator<<( std::ostream& out, const Geometry& b ) {
		b.print(out);
		return out;
	}
	virtual void print(std::ostream& out) const = 0;
};


class Sphere : public Geometry {
public:
	Sphere(const Vect3d& position,
			const double radius,
			const bool in) :
				position(position),radius(radius),radius_sq(radius*radius),in(in) {
	}
	static std::auto_ptr<Sphere> New(const Vect3d& position,
			const double radius,
			const bool in) {
		return std::auto_ptr<Sphere>(new Sphere(position,radius,in));
	}

	bool is_in(const Vect3d& point) const {
		bool inside;
		const double radial_dist_sq = radial_distance_to_boundary_sq(point);
		if (in) {
			inside = (radial_dist_sq < pow(radius-GEOMETRY_TOLERANCE,2));
		} else {
			inside = (radial_dist_sq < pow(radius+GEOMETRY_TOLERANCE,2));
		}
		return inside == in;
	}

	std::pair<double,Vect3d> lineXsurface(const Vect3d &p1, const Vect3d &p2) {
		std::pair<double,Vect3d> ret;
		Vect3d dx = (p2-p1);
		dx.normalize();
		const Vect3d dr = p1-position;
		const double dr_norm2 = dr.squaredNorm();
		const double dx_dot_dr = dx.dot(dr);

		const double in_sqrt = pow(dx_dot_dr,2) - dr_norm2 + radius_sq;
		if (in_sqrt < 0) {
			//no intersect
			ret.first = -1;
		} else if (in_sqrt == 0) {
			//one intersect
			ret.first = -dx_dot_dr;
			ret.second =  p1 + ret.first*dx - position;
			ret.second.normalize();
		} else {
			//two intersect, take the first one
			ret.first = -dx_dot_dr - sqrt(in_sqrt);
			ret.second =  p1 + ret.first*dx - position;
			ret.second.normalize();
		}
		return ret;
	}

	const Vect3d shortest_vector_to_boundary(const Vect3d& point) const {
		const double radial_dist = sqrt(radial_distance_to_boundary_sq(point));
		return (radius/radial_dist-1.)*point;
	}

	double distance_to_boundary(const Vect3d& point) const {
		double dist = shortest_vector_to_boundary(point).norm();
		// Boundaries expect the distance to be negative if outside of
		// geometry
		if (is_in(point)) {
			return dist;
		} else {
			return -dist;
		}
	}

	virtual void print(std::ostream& out) const {
		out << "Sphere with radius " << radius << " at position " << position;
	}
private:
	friend std::ostream& operator<< (std::ostream& out, const Sphere& p);
	Vect3d position;
	double radius,radius_sq;
	bool in;

	inline double radial_distance_to_boundary_sq(const Vect3d& point) const {
		return (position-point).squaredNorm();
	}
};


}

#endif /* GEOMETRY_H_ */
