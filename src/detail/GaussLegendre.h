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

#ifndef GAUSS_LEGENDRE_H_
#define GAUSS_LEGENDRE_H_

#include <array>

namespace Aboria {

namespace detail {

template <size_t N> struct GaussLegendre {};

template <> struct GaussLegendre<2> {

  static constexpr std::array<double, 2> nodes = {
      {-0.5773502691896257645091488, 0.5773502691896257645091488}};
  static constexpr std::array<double, 2> weights = {{1.0, 1.0}};
};

template <> struct GaussLegendre<3> {
  static constexpr std::array<double, 3> nodes = {
      {-0.7745966692414833770358531, 0.0, 0.7745966692414833770358531}};

  static constexpr std::array<double, 3> weights = {
      {0.5555555555555555555555556, 0.8888888888888888888888889,
       0.5555555555555555555555556}};
};

template <> struct GaussLegendre<4> {
  static constexpr std::array<double, 4> nodes = {
      {-0.8611363115940525752239465, -0.3399810435848562648026658,
       0.3399810435848562648026658, 0.8611363115940525752239465}};

  static constexpr std::array<double, 4> weights = {
      {0.3478548451374538573730639, 0.6521451548625461426269361,
       0.6521451548625461426269361, 0.3478548451374538573730639}};
};

template <> struct GaussLegendre<5> {
  static constexpr std::array<double, 5> nodes = {
      {-0.9061798459386639927976269, -0.5384693101056830910363144, 0.0,
       0.5384693101056830910363144, 0.9061798459386639927976269}};

  static constexpr std::array<double, 5> weights = {
      {0.2369268850561890875142640, 0.4786286704993664680412915,
       0.5688888888888888888888889, 0.4786286704993664680412915,
       0.2369268850561890875142640}};
};

template <> struct GaussLegendre<6> {
  static constexpr std::array<double, 6> nodes = {
      {-0.9324695142031520278123016, -0.6612093864662645136613996,
       -0.2386191860831969086305017, 0.2386191860831969086305017,
       0.6612093864662645136613996, 0.9324695142031520278123016}};

  static constexpr std::array<double, 6> weights = {
      {0.1713244923791703450402961, 0.3607615730481386075698335,
       0.4679139345726910473898703, 0.4679139345726910473898703,
       0.3607615730481386075698335, 0.1713244923791703450402961}};
};

template <> struct GaussLegendre<7> {
  static constexpr std::array<double, 7> nodes = {
      {-0.9491079123427585245261897, -0.7415311855993944398638648,
       -0.4058451513773971669066064, 0.0, 0.4058451513773971669066064,
       0.7415311855993944398638648, 0.9491079123427585245261897}};

  static constexpr std::array<double, 7> weights = {
      {0.1294849661688696932706114, 0.2797053914892766679014678,
       0.3818300505051189449503698, 0.4179591836734693877551020,
       0.3818300505051189449503698, 0.2797053914892766679014678,
       0.1294849661688696932706114}};
};

template <> struct GaussLegendre<8> {
  static constexpr std::array<double, 8> nodes = {
      {-0.9602898564975362316835609, -0.7966664774136267395915539,
       -0.5255324099163289858177390, -0.1834346424956498049394761,
       0.1834346424956498049394761, 0.5255324099163289858177390,
       0.7966664774136267395915539, 0.9602898564975362316835609}};

  static constexpr std::array<double, 8> weights = {
      {0.1012285362903762591525314, 0.2223810344533744705443560,
       0.3137066458778872873379622, 0.3626837833783619829651504,
       0.3626837833783619829651504, 0.3137066458778872873379622,
       0.2223810344533744705443560, 0.1012285362903762591525314}};
};

} // namespace detail

} // namespace Aboria

#endif
