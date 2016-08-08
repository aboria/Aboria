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


#ifndef VECTOR_H_
#define VECTOR_H_

//#include <math.h>
#include <cmath>

#include <iostream>


namespace Aboria {


template<typename T,unsigned int N>
class Vector {
public:
	typedef T value_type;
	const static int size = N;

	Vector() {}
	Vector(T arg1) {
        for (unsigned int i=0; i<N; i++) {
		    mem[i] = arg1;
        }
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
	double inner_product(const Vector<T2,N> &arg) const {
		double ret = 0;
		for (int i = 0; i < N; ++i) {
			ret += arg[i]*mem[i];
		}
		return ret;
	}

    template <typename T2>
    Vector<T2,N> cast() {
        Vector<T2,N> ret;
        for (int i = 0; i < N; ++i) {
            ret[i] = static_cast<T2>(mem[i]);
		}
		return ret;
	}


	template<typename T2>
	double dot(const Vector<T2,N> &arg) const {
		return inner_product(arg);
	}
	double squaredNorm() const {
		double ret = 0;
		for (int i = 0; i < N; ++i) {
			ret += std::pow(mem[i],2);
		}
		return ret;
	}

		
	double norm() const {
		return std::sqrt(squaredNorm());
	}

	template<typename EXP_T>
	Vector<T,N> pow(const EXP_T exponent) {
		Vector<T,N> n = *this;
		for (int i = 0; i < N; ++i) {
			n[i] = std::pow(n[i],exponent);
		}
		return n;
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

template<typename T, unsigned int N, typename EXP_T>
Vector<T,N> pow(Vector<T,N> arg, EXP_T exponent) {
	return arg.pow(exponent);
}

/*
template<typename T1,typename T2,unsigned int N> 
bool operator ==(const Vector<T1,N> &arg1, const Vector<T2,N> &arg2) { 
    bool ret = true; 
    for (int i = 0; i < N; ++i) { 
        ret &= arg1[i] == arg2[i]; 
    } 
    return ret; 
} 
*/


#define UNARY_OPERATOR(the_op) \
		template<typename T,unsigned int N> \
				Vector<double,N> operator the_op(const Vector<T,N> &arg1) { \
					Vector<double,N> ret; \
					for (int i = 0; i < N; ++i) { \
						ret[i] = the_op arg1[i]; \
					} \
					return ret; \
				} \

UNARY_OPERATOR(-)

#define OPERATOR(the_op) \
		template<typename T1,typename T2,unsigned int N> \
				Vector<double,N> operator the_op(const Vector<T1,N> &arg1, const Vector<T2,N> &arg2) { \
					Vector<double,N> ret; \
					for (int i = 0; i < N; ++i) { \
						ret[i] = arg1[i] the_op arg2[i]; \
					} \
					return ret; \
				} \
				template<typename T1,typename T2,unsigned int N> \
								Vector<double,N> operator the_op(const Vector<T1,N> &arg1, const T2 &arg2) { \
									Vector<double,N> ret; \
									for (int i = 0; i < N; ++i) { \
										ret[i] = arg1[i] the_op arg2; \
									} \
									return ret; \
								} \
								template<typename T1,typename T2,unsigned int N> \
												Vector<double,N> operator the_op(const T1 &arg1, const Vector<T2,N> &arg2) { \
													Vector<double,N> ret; \
													for (int i = 0; i < N; ++i) { \
														ret[i] = arg1 the_op arg2[i]; \
													} \
													return ret; \
												} \
 \
		template<int,int,unsigned int N> \
		Vector<int,N> operator the_op(const Vector<int,N> &arg1, const Vector<int,N> &arg2) { \
			Vector<int,N> ret; \
			for (int i = 0; i < N; ++i) { \
				ret[i] = arg1[i] the_op arg2[i]; \
			} \
			return ret; \
		} \
		template<int,int,unsigned int N> \
				Vector<int,N> operator the_op(const int &arg1, const Vector<int,N> &arg2) { \
					Vector<int,N> ret; \
					for (int i = 0; i < N; ++i) { \
						ret[i] = arg1 the_op arg2[i]; \
					} \
					return ret; \
				} \
				template<int,int,unsigned int N> \
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
		template<typename T1,typename T2,unsigned int N> \
				Vector<bool,N> operator the_op(const Vector<T1,N> &arg1, const Vector<T2,N> &arg2) { \
					Vector<bool,N> ret; \
					for (int i = 0; i < N; ++i) { \
						ret[i] = arg1[i] the_op arg2[i]; \
					} \
					return ret; \
				} \
				template<typename T1,typename T2,unsigned int N> \
								Vector<bool,N> operator the_op(const Vector<T1,N> &arg1, const T2 &arg2) { \
									Vector<bool,N> ret; \
									for (int i = 0; i < N; ++i) { \
										ret[i] = arg1[i] the_op arg2; \
									} \
									return ret; \
								} \
								template<typename T1,typename T2,unsigned int N> \
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
COMPARISON(==)

#define COMPOUND_ASSIGN(the_op) \
		template<typename T1,typename T2,unsigned int N> \
				Vector<double,N> &operator the_op(Vector<T1,N> &arg1, const Vector<T2,N> &arg2) { \
					for (int i = 0; i < N; ++i) { \
						arg1[i] the_op arg2[i]; \
					} \
					return arg1; \
				} \
				template<typename T1,typename T2,unsigned int N> \
								Vector<double,N> &operator the_op(Vector<T1,N> &arg1, const T2 &arg2) { \
									for (int i = 0; i < N; ++i) { \
										arg1[i] the_op arg2; \
									} \
									return arg1; \
								} \
 \
		template<int,int,unsigned int N> \
		Vector<int,N> &operator the_op(Vector<int,N> &arg1, const Vector<int,N> &arg2) { \
			for (int i = 0; i < N; ++i) { \
				arg1[i] the_op arg2[i]; \
			} \
			return arg1; \
		} \
				template<int,int,unsigned int N> \
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

#define UFUNC(the_op) \
	template<typename T,unsigned int N> \
	Vector<T,N> the_op(const Vector<T,N> &arg1) { \
		Vector<T,N> ret; \
	    for (int i = 0; i < N; ++i) { \
		    ret[i] = std::the_op(arg1[i]); \
        }  \
		return ret; \
	} \

UFUNC(floor)
UFUNC(ceil)
UFUNC(round)

template<typename T, int I>
double norm(const Vector<T,I> &arg1) {
	return arg1.norm();
}


template<typename T, int I>
double squaredNorm(const Vector<T,I> &arg1) {
	return arg1.squaredNorm();
}



template<typename T1, typename T2, int I>
double dot(const Vector<T1,I> &arg1, const Vector<T2,I> &arg2) {
	return arg1.inner_product(arg2);
}

template<typename T>
Vector<T,3> cross(const Vector<T,3> &arg1,const Vector<T,3> &arg2) {
	Vector<T,3> ret;
	ret[0] = arg1[1]*arg2[2] - arg1[2]*arg2[1];
	ret[1] = -arg1[0]*arg2[2] + arg1[2]*arg2[0];
	ret[2] = arg1[0]*arg2[1] - arg1[1]*arg2[0];
	return ret;
}



template<typename T,unsigned int N>
std::ostream& operator<< (std::ostream& out, const Vector<T,N>& v) {
	out << "(";
	for (int i = 0; i < N; ++i) {
		out << v[i];
		if (i != N-1) out << ",";
	}
	return out << ")";
}

template<typename T,unsigned int N>
std::istream& operator>> (std::istream& out, Vector<T,N>& v) {
    out.get();
	for (int i = 0; i < N; ++i) {
		out >> v[i];
		out.get();
	}
    return out;
}

typedef Vector<double,1> aboria_double1;
typedef Vector<double,2> aboria_double2;
typedef Vector<double,3> aboria_double3;
typedef Vector<double,4> aboria_double4;
typedef Vector<double,5> aboria_double5;
typedef Vector<double,6> aboria_double6;
typedef Vector<double,7> aboria_double7;

typedef Vector<int,1> aboria_int1;
typedef Vector<int,2> aboria_int2;
typedef Vector<int,3> aboria_int3;
typedef Vector<int,4> aboria_int4;
typedef Vector<int,5> aboria_int5;
typedef Vector<int,6> aboria_int6;
typedef Vector<int,7> aboria_int7;

typedef Vector<bool,1> aboria_bool1;
typedef Vector<bool,2> aboria_bool2;
typedef Vector<bool,3> aboria_bool3;
typedef Vector<bool,4> aboria_bool4;
typedef Vector<bool,5> aboria_bool5;
typedef Vector<bool,6> aboria_bool6;
typedef Vector<bool,7> aboria_bool7;

#define double1 aboria_double1
#define double2 aboria_double2
#define double3 aboria_double3
#define double4 aboria_double4
#define double5 aboria_double5
#define double6 aboria_double6
#define double7 aboria_double7

#define int1 aboria_int1
#define int2 aboria_int2
#define int3 aboria_int3
#define int4 aboria_int4
#define int5 aboria_int5
#define int6 aboria_int6
#define int7 aboria_int7

#define bool1 aboria_bool1
#define bool2 aboria_bool2
#define bool3 aboria_bool3
#define bool4 aboria_bool4
#define bool5 aboria_bool5
#define bool6 aboria_bool6
#define bool7 aboria_bool7

}
#endif /* VECTOR_H_ */
