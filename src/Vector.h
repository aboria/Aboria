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

#include <math.h>
#include <iostream>


namespace Aboria {


template<typename T,int N>
class Vector {
public:
	Vector() {}
	Vector(T arg1) {
		mem[0] = arg1;
	}
	Vector(T arg1,T arg2) {
		mem[0] = arg1;
		mem[1] = arg2;
	}
	Vector(T arg1,T arg2,T arg3) {
		mem[0] = arg1;
		mem[1] = arg2;
		mem[2] = arg3;
	}
	template<typename T2>
	Vector(const Vector<T2,N> &arg) {
		for (int i = 0; i < N; ++i) {
			mem[i] = arg[i];
		}
	}
	template<typename T2>
	Vector<T,N> &operator =(Vector<T2,N> &arg) {
		for (int i = 0; i < N; ++i) {
			mem[i] = arg[i];
		}
		return *this;
	}
	template<typename T2>
	Vector<T,N> &operator =(T2 *arg) {
		for (int i = 0; i < N; ++i) {
			mem[i] = arg[i];
		}
		return *this;
	}
	const T &operator[](unsigned int n) const {
		return mem[n];
	}
	T &operator[](unsigned int n) {
		return mem[n];
	}
	template<typename T2>
	double dot(const Vector<T2,N> &arg) const {
		double ret = 0;
		for (int i = 0; i < N; ++i) {
			ret += arg[i]*mem[i];
		}
		return ret;
	}
	double squaredNorm() const {
		double ret = 0;
		for (int i = 0; i < N; ++i) {
			ret += pow(mem[i],2);
		}
		return ret;
	}
		
	double norm() const {
		return sqrt(squaredNorm());
	}

	void normalize() {
		double n = norm();
		for (int i = 0; i < N; ++i) {
			mem[i] /= n;
		}
	} 
		
	bool all() const {
		bool ret = true;
		for (int i = 0; i < N; ++i) {
			ret &= mem[i];
		}
		return ret;
	}
	bool any() const {
		bool ret = false;
		for (int i = 0; i < N; ++i) {
			ret |= mem[i];
		}
		return ret;
	}
	T minCoeff() const {
		T min = mem[0];
		for (int i = 1; i < N; ++i) {
			if (mem[i]<min) {
				min = mem[i];
			}
		}
		return min;
	}
	T maxCoeff() const {
		T max = mem[0];
		for (int i = 1; i < N; ++i) {
			if (mem[i]>max) {
				max = mem[i];
			}
		}
		return max;
	}
	T prod() const {
		T ret = 1;
		for (int i = 0; i < N; ++i) {
			ret *= mem[i];
		}
		return ret;
	}
	T *data() {
		return mem;
	}
private:
	T mem[N];
};

#define OPERATOR(the_op) \
		template<typename T1,typename T2,int N> \
				Vector<double,N> operator the_op(const Vector<T1,N> &arg1, const Vector<T2,N> &arg2) { \
					Vector<double,N> ret; \
					for (int i = 0; i < N; ++i) { \
						ret[i] = arg1[i] the_op arg2[i]; \
					} \
					return ret; \
				} \
				template<typename T1,typename T2,int N> \
								Vector<double,N> operator the_op(const Vector<T1,N> &arg1, const T2 &arg2) { \
									Vector<double,N> ret; \
									for (int i = 0; i < N; ++i) { \
										ret[i] = arg1[i] the_op arg2; \
									} \
									return ret; \
								} \
								template<typename T1,typename T2,int N> \
												Vector<double,N> operator the_op(const T1 &arg1, const Vector<T2,N> &arg2) { \
													Vector<double,N> ret; \
													for (int i = 0; i < N; ++i) { \
														ret[i] = arg1 the_op arg2[i]; \
													} \
													return ret; \
												} \
 \
		template<int,int,int N> \
		Vector<int,N> operator the_op(const Vector<int,N> &arg1, const Vector<int,N> &arg2) { \
			Vector<int,N> ret; \
			for (int i = 0; i < N; ++i) { \
				ret[i] = arg1[i] the_op arg2[i]; \
			} \
			return ret; \
		} \
		template<int,int,int N> \
				Vector<int,N> operator the_op(const int &arg1, const Vector<int,N> &arg2) { \
					Vector<int,N> ret; \
					for (int i = 0; i < N; ++i) { \
						ret[i] = arg1 the_op arg2[i]; \
					} \
					return ret; \
				} \
				template<int,int,int N> \
						Vector<int,N> operator the_op(const Vector<int,N> &arg1, const int &arg2) { \
							Vector<int,N> ret; \
							for (int i = 0; i < N; ++i) { \
								ret[i] = arg1[i] the_op arg2; \
							} \
							return ret; \
						} \


OPERATOR(+)
OPERATOR(-)
OPERATOR(/)
OPERATOR(*)

#define COMPARISON(the_op) \
		template<typename T1,typename T2,int N> \
				Vector<bool,N> operator the_op(const Vector<T1,N> &arg1, const Vector<T2,N> &arg2) { \
					Vector<bool,N> ret; \
					for (int i = 0; i < N; ++i) { \
						ret[i] = arg1[i] the_op arg2[i]; \
					} \
					return ret; \
				} \
				template<typename T1,typename T2,int N> \
								Vector<bool,N> operator the_op(const Vector<T1,N> &arg1, const T2 &arg2) { \
									Vector<bool,N> ret; \
									for (int i = 0; i < N; ++i) { \
										ret[i] = arg1[i] the_op arg2; \
									} \
									return ret; \
								} \
								template<typename T1,typename T2,int N> \
																Vector<bool,N> operator the_op(const T1 &arg1, const T2 &arg2) { \
																	Vector<bool,N> ret; \
																	for (int i = 0; i < N; ++i) { \
																		ret[i] = arg1 the_op arg2; \
																	} \
																	return ret; \
																} \

COMPARISON(>)
COMPARISON(<)
COMPARISON(<=)
COMPARISON(>=)

#define COMPOUND_ASSIGN(the_op) \
		template<typename T1,typename T2,int N> \
				Vector<double,N> &operator the_op(Vector<T1,N> &arg1, const Vector<T2,N> &arg2) { \
					for (int i = 0; i < N; ++i) { \
						arg1[i] the_op arg2[i]; \
					} \
					return arg1; \
				} \
				template<typename T1,typename T2,int N> \
								Vector<double,N> &operator the_op(Vector<T1,N> &arg1, const T2 &arg2) { \
									for (int i = 0; i < N; ++i) { \
										arg1[i] the_op arg2; \
									} \
									return arg1; \
								} \
 \
		template<int,int,int N> \
		Vector<int,N> &operator the_op(Vector<int,N> &arg1, const Vector<int,N> &arg2) { \
			for (int i = 0; i < N; ++i) { \
				arg1[i] the_op arg2[i]; \
			} \
			return arg1; \
		} \
				template<int,int,int N> \
						Vector<int,N> &operator the_op(Vector<int,N> &arg1, const int &arg2) { \
							for (int i = 0; i < N; ++i) { \
								arg1[i] the_op arg2; \
							} \
							return arg1; \
						} \


COMPOUND_ASSIGN(+=)
COMPOUND_ASSIGN(-=)
COMPOUND_ASSIGN(*=)
COMPOUND_ASSIGN(/=)

template<typename T>
Vector<T,3> cross(const Vector<T,3> &arg1,const Vector<T,3> &arg2) {
	Vector<T,3> ret;
	ret[0] = arg1[1]*arg2[2] - arg1[2]*arg2[1];
	ret[1] = -arg1[0]*arg2[2] + arg1[2]*arg2[0];
	ret[2] = arg1[0]*arg2[1] - arg1[1]*arg2[0];
	return ret;
}

template<typename T,int N>
std::ostream& operator<< (std::ostream& out, const Vector<T,N>& v) {
	out << "(";
	for (int i = 0; i < N; ++i) {
		out << v[i];
		if (i != N-1) out << ",";
	}
	return out << ")";
}

typedef Vector<double,3> Vect3d;
typedef Vector<int,3> Vect3i;
typedef Vector<int,6> Vect6i;
typedef Vector<bool,3> Vect3b;
typedef Vector<double,2> Vect2d;
typedef Vector<int,2> Array2i;


}
#endif /* VECTOR_H_ */
