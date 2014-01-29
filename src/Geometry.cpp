/*
 * Geometry.cpp
 *
 *  Created on: 19 Oct 2012
 *      Author: robinsonm
 */

#include "Geometry.h"
#include "Constants.h"

namespace Tyche {


std::ostream& operator<< (std::ostream& out, const NullGeometry& b) {
	return out << "Null Geometry";
}

std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<0>& p) {
	return out << "x = " << p.get_coord() << " with normal " << p.get_normal();
}

std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<1>& p) {
	return out << "y = " << p.get_coord() << " with normal " << p.get_normal();
}

std::ostream& operator<< (std::ostream& out, const AxisAlignedPlane<2>& p) {
	return out << "z = " << p.get_coord() << " with normal " << p.get_normal();
}

std::ostream& operator<< (std::ostream& out, const Rectangle& p) {
	return out << "Rectangle with equation x = "<<p.get_low()<<" + s*"<<p.get_l()<<" + t*"<<p.get_r()<<" and normal = "<<p.get_normal();
}

std::ostream& operator<< (std::ostream& out, const Box& p) {
	return out << "Box";
}

std::ostream& operator<< (std::ostream& out, const MultipleBoxes& p) {
	return out << "Multiple Boxes";
}

std::ostream& operator<< (std::ostream& out, const Sphere& p) {
	return out << "Sphere";

}

std::ostream& operator<< (std::ostream& out, const Shell& p) {
	return out << "Shell";
}

}
