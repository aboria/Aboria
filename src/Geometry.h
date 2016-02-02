/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Aboria.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "Vector.h"
#include "Log.h"
#include "Ptr.h"


namespace Aboria {

const double GEOMETRY_TOLERANCE = 1.0/100000.0;

class GeometryBase {
public:
	virtual bool is_in(const Vect3d &point) const = 0;
	virtual std::pair<double,Vect3d> lineXsurface(const Vect3d &p1, const Vect3d &p2) const = 0;
	virtual const Vect3d shortest_vector_to_boundary(const Vect3d &point) const = 0;
	friend std::ostream& operator<<( std::ostream& out, const GeometryBase& b ) {
		b.print(out);
		return out;
	}
	virtual void print(std::ostream& out) const = 0;
};

template<typename T>
bool reflect_once(const Vect3d& p1, Vect3d &p2, const T& geometry) {
	std::pair<double,Vect3d> res = geometry.lineXsurface(p1,p2);
	if (res.first >= 0) {
		/*
		 * if line going through surface with normal then don't reflect
		 */
		const Vect3d dx = p2 - p1;
		const double dxCosTheta = res.second.dot(dx);
		if (dxCosTheta > 0) {
			return false;
		}

		/*
		 * calculate new position after reflection
		 */
		const double ratio = res.first/dx.norm();
		const Vect3d v_reflect = dx - 2*dxCosTheta*res.second;
		p2 = p1 + (dx-v_reflect)*ratio + v_reflect;
		return true;
	} else {
		return false;
	}
}

template <typename T>
Vect3d operator| (const Vect3d& vector, const T& geometry) {
    Vect3d result = vector;
    reflect_once(Vect3d(0,0,0),result,geometry);
    return result;
}

template <typename T>
Vect3d operator| (const T& geometry, const Vect3d& vector) {
    Vect3d result = vector;
    reflect_once(Vect3d(0,0,0),result,geometry);
    return result;
}

class Sphere : public GeometryBase {
public:
	Sphere(const Vect3d& position,
			const double radius,
			const bool in) :
				position(position),radius(radius),radius_sq(radius*radius),in(in) {
	}
	static ptr<Sphere> New(const Vect3d& position,
			const double radius,
			const bool in) {
		return ptr<Sphere>(new Sphere(position,radius,in));
	}

	bool is_in(const Vect3d& point) const {
		bool inside;
		const double radial_dist_sq = radial_distance_to_boundary_sq(point);
		if (in) {
			inside = (radial_dist_sq < std::pow(radius-GEOMETRY_TOLERANCE,2));
		} else {
			inside = (radial_dist_sq < std::pow(radius+GEOMETRY_TOLERANCE,2));
		}
		return inside == in;
	}

	std::pair<double,Vect3d> lineXsurface(const Vect3d &p1, const Vect3d &p2) const {
		std::pair<double,Vect3d> ret;
		Vect3d dx = (p2-p1);
		const double maxd = dx.norm();
		dx /= maxd;
		const Vect3d dr = p1-position;
		const double dr_norm2 = dr.squaredNorm();
		const double dx_dot_dr = dx.dot(dr);
		const double in_sqrt = std::pow(dx_dot_dr,2) - dr_norm2 + radius_sq;
		if (in_sqrt < 0) {
			//no intersect
			ret.first = -1;
		} else if (in_sqrt == 0) {
			//one intersect
			ret.first = -dx_dot_dr;
		} else {
			//two intersect
			ret.first = -dx_dot_dr - std::sqrt(in_sqrt);
			if (ret.first < 0) {
				ret.first = -dx_dot_dr + std::sqrt(in_sqrt);

			}
		}
		if ((ret.first > 0) && (ret.first <= maxd)) {
			ret.second =  p1 + ret.first*dx - position;
			if (in) {
				ret.second /= ret.second.norm();
			} else {
				ret.second /= -ret.second.norm();
			}
		} else {
			ret.first = -1;
		}
		return ret;
	}

	const Vect3d shortest_vector_to_boundary(const Vect3d& point) const {
		const double radial_dist = std::sqrt(radial_distance_to_boundary_sq(point));
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
	//friend std::ostream& operator<< (std::ostream& out, const Sphere& p);
	Vect3d position;
	double radius,radius_sq;
	bool in;

	inline double radial_distance_to_boundary_sq(const Vect3d& point) const {
		return (position-point).squaredNorm();
	}
};




}

#endif /* GEOMETRY_H_ */
