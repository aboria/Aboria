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


namespace Aboria {

const double GEOMETRY_TOLERANCE = 1.0/100000.0;

class Geometry {
public:
  virtual bool is_in(const Vect3d &point) const = 0;
  virtual bool lineXsurface(const Vect3d &p1, const Vect3d &p2, double *intersect_point=NULL, Vect3d *intersect_normal=NULL) const = 0;
  virtual const Vect3d shortest_vector_to_boundary(const Vect3d &point) const = 0;
  friend std::ostream& operator<<( std::ostream& out, const Geometry& b ) {
	  b.print(out);
	  return out;
  }
  virtual void print(std::ostream& out) const = 0;
};

class NullGeometry : public Geometry {
        bool is_in(const Vect3d &point) const {
	  return false;
	}

	bool lineXsurface(const Vect3d &p1, const Vect3d &p2, double *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
	  if ((intersect_point!=NULL) || (intersect_normal != NULL)) {
	    ASSERT(intersect_point==NULL, "intersect point not supported for NullGeometry");
	    ASSERT(intersect_normal==NULL, "intersect normal not supported for NullGeometry");
	  }
	  return false;
	}
	
	virtual const Vect3d shortest_vector_to_boundary(const Vect3d &point) const {
	  return std::nan("")*Vect3d::Ones();
	}

	virtual void print(std::ostream& out) const {
		out << "Null Geometry";
	}
};


template<unsigned int DIM>
class AxisAlignedRectangle;

static const int dim_map[][2] = {{1,2}, {0,2}, {0,1}};

template<unsigned int DIM>
class AxisAlignedPlane: public Geometry {
public:
	typedef AxisAlignedRectangle<DIM> SurfaceElementType;
	static const int dim = DIM;

	AxisAlignedPlane(const double coord, const int normal):
		coord(coord),
		normal(normal)
	{}
	AxisAlignedPlane(const double coord):
			coord(coord),
			normal(1.0)
		{}
	AxisAlignedPlane(const AxisAlignedPlane& arg):
			coord(arg.coord),
			normal(arg.normal)
	{}

	static std::auto_ptr<AxisAlignedPlane<DIM> > New(const double coord, const int normal) {
		return std::auto_ptr<AxisAlignedPlane<DIM> >(new AxisAlignedPlane<DIM>(coord,normal));
	}
	bool lineXsurface(const Vect3d &p1, const Vect3d &p2, double *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
		if (((p2[DIM]>=coord)&&(p1[DIM]<coord))||((p2[DIM]<coord)&&(p1[DIM]>=coord))) {
			if (intersect_point != NULL) {
				(*intersect_point) = (coord-p1[DIM])/(p2[DIM]-p1[DIM]);
//				(*intersect_point)[dim_map[DIM][0]] = 0.5*(p1[dim_map[DIM][0]] + p2[dim_map[DIM][0]]);
//				(*intersect_point)[dim_map[DIM][1]] = 0.5*(p1[dim_map[DIM][1]] + p2[dim_map[DIM][1]]);
			}
			if (intersect_normal != NULL) {
				(*intersect_normal)[DIM] = 1.0;
				(*intersect_normal)[dim_map[DIM][0]] = 0.0;
				(*intersect_normal)[dim_map[DIM][1]] = 0.0;
			}
			return true;
		} else {
			return false;
		}
	}


	bool is_in(const Vect3d& point) const {
		return normal*(point[DIM]-coord) < GEOMETRY_TOLERANCE;
	}

	const Vect3d shortest_vector_to_boundary(const Vect3d& r) const {
	   Vect3d shortest = Vect3d::Zero();
		shortest[DIM] = coord-r[DIM];
		return shortest;
	}
	inline double distance_to_boundary(const Vect3d& r) const {
		return normal*(r[DIM]-coord);
	}


	virtual void print(std::ostream& out) const {
		const std::string options[3] = {"x","y","z"};
		out << options[DIM] << " = " << get_coord() << " with normal " << get_normal();
	}

protected:
	double coord;
	int normal;

};


typedef AxisAlignedPlane<0> xplane;
typedef AxisAlignedPlane<1> yplane;
typedef AxisAlignedPlane<2> zplane;



template<unsigned int DIM>
class AxisAlignedRectangle: public AxisAlignedPlane<DIM> {
public:
	AxisAlignedRectangle(const Vect3d _low, const Vect3d _high, const int _normal):
	   AxisAlignedPlane<DIM>(_low[DIM], _normal),
	   low(_low),high(_high),
	   normal_vector(Vect3d::Zero()) {
	   high[DIM] = low[DIM];
	   normal_vector[DIM] = this->normal;
	}
	AxisAlignedRectangle(const AxisAlignedRectangle<DIM>& arg):
		AxisAlignedPlane<DIM>(arg),
		low(arg.low),
		high(arg.high),
		normal_vector(arg.normal_vector)
						{}

	static std::auto_ptr<AxisAlignedRectangle<DIM> > New(const Vect3d _low, const Vect3d _high, const int _normal) {
		return std::auto_ptr<AxisAlignedRectangle<DIM> >(new AxisAlignedRectangle<DIM>(_low,_high,_normal));
	}

	AxisAlignedRectangle<DIM>& operator=(const AxisAlignedRectangle<DIM>& arg) {
		AxisAlignedPlane<DIM>::operator=(arg);
		low = arg.low;
		high = arg.high;
		normal_vector = arg.normal_vector;
		return *this;
	}

	bool lineXsurface(const Vect3d& p1, const Vect3d& p2, double *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
		if (((p2[DIM]>=low[DIM])&&(p1[DIM]<low[DIM]))||((p2[DIM]<low[DIM])&&(p1[DIM]>=low[DIM]))) {
			const double intersect0 = 0.5*(p1[dim_map[DIM][0]] + p2[dim_map[DIM][0]]);
			if ((intersect0 >= low[dim_map[DIM][0]]) && (intersect0 < high[dim_map[DIM][0]])) {
				const double intersect1 = 0.5*(p1[dim_map[DIM][1]] + p2[dim_map[DIM][1]]);
				if ((intersect1 >= low[dim_map[DIM][1]]) && (intersect1 < high[dim_map[DIM][1]])) {
					if (intersect_point != NULL) {
						(*intersect_point) = (low[DIM]-p1[DIM])/(p2[DIM]-p1[DIM]);
					}
					if (intersect_normal != NULL) {
						(*intersect_normal)[DIM] = 1.0;
						(*intersect_normal)[dim_map[DIM][0]] = 0.0;
						(*intersect_normal)[dim_map[DIM][1]] = 0.0;
					}
					return true;
				}
			}
		}
		return false;
	}

	bool is_in(const Vect3d& point) const {
		return (normal_vector[DIM]*(point[DIM]-low[DIM]) < GEOMETRY_TOLERANCE) &&
				(point[dim_map[DIM][0]] > low[dim_map[DIM][0]]+GEOMETRY_TOLERANCE) &&
				(point[dim_map[DIM][1]] > low[dim_map[DIM][1]]+GEOMETRY_TOLERANCE) &&
				(point[dim_map[DIM][0]] < high[dim_map[DIM][0]]-GEOMETRY_TOLERANCE) &&
				(point[dim_map[DIM][1]] < high[dim_map[DIM][1]]-GEOMETRY_TOLERANCE)
				;
	}

	void get_random_point_and_normal(Vect3d& p, Vect3d& n) {
	   p = get_random_point();
	   n = normal_vector;
	}
	Vect3d get_random_point() {
		Vect3d ret;
		ret[DIM] = this->coord;
		boost::variate_generator<base_generator_type&, boost::uniform_real<> >
					uni1(generator,boost::uniform_real<>(low[dim_map[DIM][0]],high[dim_map[DIM][0]]));
		ret[dim_map[DIM][0]] = uni1();
		boost::variate_generator<base_generator_type&, boost::uniform_real<> >
					uni2(generator,boost::uniform_real<>(low[dim_map[DIM][1]],high[dim_map[DIM][1]]));
		ret[dim_map[DIM][1]] = uni2();
		return ret;
	}

	void get_random_point_and_normal_triangle(Vect3d& p, Vect3d& n) {
		p = get_random_point_triangle();
		n = normal_vector;
	}
	Vect3d get_random_point_triangle() {
		Vect3d ret;
		ret[DIM] = this->coord;
		boost::variate_generator<base_generator_type&, boost::triangle_distribution<> >
							tri1(generator,boost::triangle_distribution<>(low[dim_map[DIM][0]],
									0.5*(low[dim_map[DIM][0]]+high[dim_map[DIM][0]]),
									high[dim_map[DIM][0]]));
		ret[dim_map[DIM][0]] = tri1();
		boost::variate_generator<base_generator_type&, boost::triangle_distribution<> >
								tri2(generator,boost::triangle_distribution<>(low[dim_map[DIM][1]],
										0.5*(low[dim_map[DIM][1]]+high[dim_map[DIM][1]]),
										high[dim_map[DIM][1]]));
		ret[dim_map[DIM][1]] = tri2();
		return ret;
	}

	const Vect3d& get_low() const {return low;}
	const Vect3d& get_high() const {return high;}

	virtual void print(std::ostream& out) const {
		out << "Rectangle aligned with dimension " << DIM << ". Lower point in other dimensions is "<<get_low()<<". Upper point in other dimensions is "<<get_high()<<".";
	}

private:
	Vect3d low,high,normal_vector;
};

typedef AxisAlignedRectangle<0> xrect;
typedef AxisAlignedRectangle<1> yrect;
typedef AxisAlignedRectangle<2> zrect;


class Rectangle {
public:
	Rectangle(const Vect3d& lower_corner,
			const Vect3d& upper_left_corner,
			const Vect3d& lower_right_corner):
				low(lower_corner),high(upper_left_corner) {
		l = upper_left_corner-lower_corner;
		r = lower_right_corner-lower_corner;
		normal = l.cross(r);
		normal.normalize();
	}

	bool lineXsurface(const Vect3d& p1, const Vect3d& p2, double *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
		for (int d = 0; d < 3; ++d) {
			if (((p1[d]<low[d]) && (p2[d]<low[d])) || ((p1[d]>=high[d]) && (p2[d]>=high[d]))) {
				return false;
			}
		}
		const double denominator = (p2-p1).dot(normal);
		if (denominator==0) return false;
		const double numerator = (low-p1).dot(normal);

		Vect3d vintersect_point = (numerator/denominator) * (p2-p1) + p1;
		for (int d = 0; d < 3; ++d) {
			if (((vintersect_point)[d]>=high[d]) && ((vintersect_point)[d]<low[d])) {
				return false;
			}
		}
		if (intersect_point != NULL) {
			*intersect_point = numerator/denominator;
		}
		if (intersect_normal != NULL) {
			*intersect_normal = normal;
		}
		return true;
	}


	void get_random_point_and_normal(Vect3d& p, Vect3d& n) {
		p = get_random_point();
		n = normal;
	}
	Vect3d get_random_point() {
		boost::variate_generator<base_generator_type&, boost::uniform_real<> >
			uni(generator,boost::uniform_real<>(0,1));
		return low + l*uni() + r*uni();
	}

	void get_random_point_and_normal_triangle(Vect3d& p, Vect3d& n) {
		p = get_random_point_triangle();
		n = normal;
	}
	Vect3d get_random_point_triangle() {
		boost::variate_generator<base_generator_type&, boost::triangle_distribution<> >
					tri(generator,boost::triangle_distribution<>(0,0.5,1));
		return low + l*tri() + r*tri();
	}
	const Vect3d& get_low() const {return low;}
	const Vect3d& get_l() const {return l;}
	const Vect3d& get_r() const {return r;}
	const Vect3d& get_normal() const {return normal;}

	virtual void print(std::ostream& out) const {
		 out << "Rectangle with equation x = "<<get_low()<<" + s*"<<get_l()<<" + t*"<<get_r()<<" and normal = "<<get_normal();
	}
private:
	Vect3d low,high;
	Vect3d l,r;
	Vect3d normal;
};


class Box : public Geometry {
public:
	Box(const Vect3d& lower_corner,
			const Vect3d& upper_corner,
			const bool in):low(lower_corner),high(upper_corner),in(in) {
	}
	static std::auto_ptr<Box> New(const Vect3d& lower_corner,
			const Vect3d& upper_corner, const bool in) {
		return std::auto_ptr<Box>(new Box(lower_corner,upper_corner,in));
	}
	bool is_in(const Vect3d& point) const {
		bool inside;
		if (in) {
			inside = ((point.array() > low.array()+GEOMETRY_TOLERANCE).all() && (point.array() < high.array()-GEOMETRY_TOLERANCE).all());
		} else {
			inside = ((point.array() > low.array()-GEOMETRY_TOLERANCE).all() && (point.array() < high.array()+GEOMETRY_TOLERANCE).all());
		}
		//std::cout << "testing point "<<point<<". result = "<< (inside==in) << std::endl;
		return inside == in;
	}
	bool lineXsurface(const Vect3d& p1, const Vect3d& p2, double *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
		ASSERT(intersect_point==NULL, "intersect point not supported for Box geometry");
		ASSERT(intersect_normal==NULL, "intersect normal not supported for Box geometry");
		const bool is_in1 = is_in(p1);
		const bool is_in2 = is_in(p2);
		return is_in1 != is_in2;
	}

	const Vect3d shortest_vector_to_boundary(const Vect3d& point) const {
		Vect3d box_point;
		if ((point.array() >= low.array()).all() && (point.array() <= high.array()).all()) {
			double min_value = 100000*(high[0]-low[0]);
			int min_d;
			bool is_low;
			for (int d = 0; d < 3; ++d) {
				if (point[d]-low[d] < min_value) {
					min_value = point[d]-low[d];
					is_low = true;
					min_d = d;
				}
				if (high[d]-point[d] < min_value) {
					min_value = high[d]-point[d];
					is_low = false;
					min_d = d;
				}
			}
			box_point = point;
			if (is_low) {
				box_point[min_d] = low[min_d];
			} else {
				box_point[min_d] = high[min_d];
			}
		} else {
			for (int d = 0; d < 3; ++d) {
				if (point[d]<=low[d]) {
					box_point[d] = low[d];
				} else if (point[d]>=high[d]) {
					box_point[d] = high[d];
				} else {
					box_point[d] = point[d];
				}
			}
		}
		return box_point - point;
	}

	double distance_to_boundary(const Vect3d& point) const {
		return shortest_vector_to_boundary(point).norm();
	}

	virtual void print(std::ostream& out) const {
		out << "Box";
	}
private:
	Vect3d low,high;
    bool in;
};


class MultipleBoxes : public Geometry {
public:
	MultipleBoxes(const bool in):in(in) {
	}
	static std::auto_ptr<MultipleBoxes> New( const bool in) {
		return std::auto_ptr<MultipleBoxes>(new MultipleBoxes(in));
	}
	void add_box(const Vect3d& lower_corner, const Vect3d& upper_corner) {
		boxes.push_back(Box(lower_corner,upper_corner,in));
	}
	bool is_in(const Vect3d& point) const {
		bool ret = true;
		const int n = boxes.size();
		if (n==0) return !in;
		for (int i = 0; i < n; ++i) {
			ret &= boxes[i].is_in(point);
		}
		return ret;
	}
	bool lineXsurface(const Vect3d& p1, const Vect3d& p2, double *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
		ASSERT(intersect_point==NULL, "intersect point not supported for BoxWithHoles geometry");
		ASSERT(intersect_normal==NULL, "intersect normal not supported for BoxWithHoles geometry");
		const bool is_in1 = is_in(p1);
		const bool is_in2 = is_in(p2);
		return is_in1 != is_in2;
	}

	const Vect3d shortest_vector_to_boundary(const Vect3d& point) const {
	  ERROR("NOT IMPLEMENTED!");
	  return Vect3d::Zero();
	}
	virtual void print(std::ostream& out) const {
		out << "Multiple Boxes";
	}
private:
	std::vector<Box> boxes;
    bool in;
};

  
template<unsigned int DIM>
class AxisAlignedCylinder : public Geometry {
public:
  AxisAlignedCylinder(const Vect3d& base,
		      const double radius,
		      const bool in) : 
    base(base),radius(radius),radius_sq(radius*radius),in(in) {
  }
  static std::auto_ptr<AxisAlignedCylinder> New(const Vect3d& base,
				     const double radius,
				     const bool in) {
    return std::auto_ptr<AxisAlignedCylinder>(new AxisAlignedCylinder(base,radius,in));
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
  bool lineXsurface(const Vect3d& p1, const Vect3d& p2, double *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
    ASSERT(intersect_point==NULL, "intersect point not supported for Cylinder geometry");
    ASSERT(intersect_normal==NULL, "intersect normal not supported for Cylinder geometry");
    const bool is_in1 = is_in(p1);
    const bool is_in2 = is_in(p2);
    return is_in1 != is_in2;
  }
  
  const Vect3d shortest_vector_to_boundary(const Vect3d& point) const {
    const double radial_dist = sqrt(radial_distance_to_boundary_sq(point));

    Vect3d cyl_point = point;
    cyl_point *= radius/radial_dist;
    cyl_point[DIM] = point[DIM];
    return cyl_point-point;
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
	  out << "Cylinder aligned with dimension " << DIM;
  }
private:
  Vect3d base;
  double radius,radius_sq;
  bool in;

  inline double radial_distance_to_boundary_sq(const Vect3d& point) const {
    double radial_dist_sq = 0;
    for (int i=0; i<3; i++) {
      if (i!=DIM) {
	double d = point[i]-base[i];
	radial_dist_sq += d*d;
      }
    }
    return radial_dist_sq;
  }
};
 
typedef AxisAlignedCylinder<0> xcylinder;
typedef AxisAlignedCylinder<1> ycylinder;
typedef AxisAlignedCylinder<2> zcylinder;

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
  bool lineXsurface(const Vect3d& p1, const Vect3d& p2, double *intersect_point=NULL, Vect3d *intersect_normal=NULL) const {
    ASSERT(intersect_point==NULL, "intersect point not supported for Cylinder geometry");
    ASSERT(intersect_normal==NULL, "intersect normal not supported for Cylinder geometry");
    const bool is_in1 = is_in(p1);
    const bool is_in2 = is_in(p2);
    return is_in1 != is_in2;
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
