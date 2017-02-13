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


#ifndef CHEBYSHEV_H_
#define CHEBYSHEV_H_

#include <boost/math/constants/constants.hpp>

namespace Aboria {
namespace detail {

// constexpr sin and cos routines from (Boost Licence):
// https://github.com/pkeir/ctfft/blob/master/const_math.hpp
// Based on the triple-angle formula: sin 3x = 3 sin x - 4 sin ^3 x

const double PI = boost::math::constants::pi<double>();
const double PI_2 = boost::math::constants::pi<double>()/2;

constexpr double tol = 0.001;

constexpr double cube(const double x) { return x*x*x; }

constexpr
double sin_helper(const double x) {
  return x < tol ? x : 3*(sin_helper(x/3.0)) - 4*cube(sin_helper(x/3.0));
}

constexpr
double sin(const double x) {
  return sin_helper(x < 0 ? -x+PI : x);
}

//sinh 3x = 3 sinh x + 4 sinh ^3 x
constexpr
double sinh_helper(const double x) {
  return x < tol ? x : 3*(sinh_helper(x/3.0)) + 4*cube(sinh_helper(x/3.0));
}

//sinh 3x = 3 sinh x + 4 sinh ^3 x
constexpr
double sinh(const double x) {
  return x < 0 ? -sinh_helper(-x) : sinh_helper(x);
}

constexpr double cos (const double x) { return sin(PI_2 - x); }

// chebychev polynomial of order k, evaluated at the i-th root of a chebychev
// polynomial of order n
template <typename T>
constexpr T chebyshev_at_node(unsigned int i, unsigned int k, unsigned int n) {
    return 
}



template <typename T>
T chebyshev_polynomial(const T &x, unsigned int n) {
    if (n==0) {
        return 1.0;
    } else if (n==1) {
        return x;
    } else {
        T Tn_2 = 1.0;
        T Tn_1 = x;
        T Tn;
        for (int i=2; i<n; ++i) {
            Tn = 2.0*x*Tn_1 - Tn_2;
            Tn_2 = Tn_1;
            Tn_1 = Tn;
        }
    }
}

// evaluate 1/n + 2/n * sum_{k=1}^{n-1} T_k(y_i) T_k(x), where y_i is the i-th
// root of the the chebyshev polynomial of order n
template <typename T>
T chebyshev_Sn(const T &x, unsigned int i, unsigned int n) {
    // Clenshaw algorithm: \alpha = 2x, \beta = -1, T0=1, T1=x
    //                     a_0 = 1/n, a_k = 2/n * cos(k*(2i-1)/(2n) * pi)
    T bk_1 = 0;
    T bk_2 = 0;
    for (unsigned int k=n-1; k>=1; --k) {
        T bk = cos(k*(2*i-1)*PI/2) + 2*x*bk_1 - bk_2; 
        bk_2 = bk_1;
        bk_1 = bk;
    }
    // one more step with 2*a_0, then s(x) = 0.5*(b0 - b2)
    return (1.0/n + x*bk_1 - bk_2);
}

template <typename T, unsigned int N>
T chebyshev_Rn(const Vector<T,N> &x, const Vector<unsigned int i,N> &i, unsigned int n) {
    T Rn = 1.0;
    for (unsigned int d = 0; d < N; ++d) {
        Rn *= chebyshev_Sn(x[d],i[d],n);
    }
}

template <typename T, typename Traits>
struct Chebyshev_Rn {
    Vector Sn;
    unsigned int n;
    Chebyshev_Rn(Vector unsigned int n):n(n),Sn(n) {
        for (int i=0; i<n; ++i) {
            Sn[i] = chebyshev_Sn();
    }

}
 


 
    
}
}


