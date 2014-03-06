/*
 * sph_common.h
 * 
 * Copyright 2014 Martin Robinson
 *
 * This file is part of Aboria.
 *
 * Aboria is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Aboria is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Aboria.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 7 Feb 2014
 *      Author: robinsonm
 */

#ifndef SPH_COMMON_H_
#define SPH_COMMON_H_

#define _3D_
const int NDIM=3;

#ifdef _1D_
const double WCON = 1;
const double WCON_HANN = 0.25;
const double WCON_PRICE2 = 6.0/(2.0*(PRICE2_A*pow(PRICE2_ALPHA,6)+PRICE2_B*pow(PRICE2_BETA,6)+pow(2.0,6)));
#endif
#ifdef _2D_
const double WCON = 15.0/(7.0*PI);
const double WCON_QUINTIC = 7.0/(478.0*PI);
const double WCON_HANN = 1.0/(4.0*PI-16.0/PI);
const double WCON_PRICE2 = 21.0/(PI*(PRICE2_A*pow(PRICE2_ALPHA,7)+PRICE2_B*pow(PRICE2_BETA,7)+pow(2.0,7)));
const double WCON_MY_CUBIC = 10.0/(PI*(32.0+MY_CUBIC_A*pow(MY_CUBIC_ALPHA,5)));
const double WCON_WENDLAND = 7.0/(64.0*PI);
#endif
#ifdef _3D_
const double WCON = 3.0/(2.0*PI);
const double WCON_WENDLAND = 21.0/(256.0*PI);
#endif


#define WENDLAND

inline double F(const double q,const double h) {
#ifdef WENDLAND
         if (q<=2.0) {
            return (1/pow(h,NDIM+2))*WCON_WENDLAND*(-4*pow(2-q,3)*(1+2*q) + 2*pow(2-q,4))/q;
         } else {
            return 0.0;
         }
#else
#ifdef MY_CUBIC
         if (q<=MY_CUBIC_ALPHA) {
            return -(1/pow(h,NDIM+2))*WCON_MY_CUBIC*3*(pow(2-q,2) + MY_CUBIC_A*pow(MY_CUBIC_ALPHA-q,2))/q;
         } else if (q<=2.0) {
            return -(1/pow(h,NDIM+2))*WCON_MY_CUBIC*3*pow(2-q,2)/q;
         } else {
            return 0.0;
         }
#else
#ifdef PRICE2
         if (q<=PRICE2_BETA) {
            return -(1/pow(h,NDIM+2))*WCON_PRICE2*5*(pow(2-q,4) + PRICE2_A*pow(PRICE2_ALPHA-q,4) + PRICE2_B*pow(PRICE2_BETA-q,4))/q;
         } else if (q<=PRICE2_ALPHA) {
            return -(1/pow(h,NDIM+2))*WCON_PRICE2*5*(pow(2-q,4) + PRICE2_A*pow(PRICE2_ALPHA-q,4))/q;
         } else if (q<=2.0) {
            return -(1/pow(h,NDIM+2))*WCON_PRICE2*5*pow(2-q,4)/q;
         } else {
            return 0.0;
         }
#else
#ifdef QUINTIC
         if (q<=1.0) {
            return -(1/pow(h,NDIM+2))*WCON_QUINTIC*5*(pow(3-q,4) - 6*pow(2-q,4) + 15*pow(1-q,4))/q;
         } else if (q<=2.0) {
            return -(1/pow(h,NDIM+2))*WCON_QUINTIC*5*(pow(3-q,4) - 6*pow(2-q,4))/q;
         } else if (q<=3.0) {
            return -(1/pow(h,NDIM+2))*WCON_QUINTIC*5*pow(3-q,4)/q;
         } else {
            return 0.0;
         }
#else
#ifdef HANN
         if (q<=2.0) {
            return (1/pow(h,NDIM+2))*WCON_HANN*0.5*PI*sin(0.5*PI*(q-2.0))/q;
         } else {
            return 0.0;
         }
#else
         if (q<=1.0) {
             return (1/pow(h,NDIM+2))*WCON*(-2.0+ 1.5*q);
         }
         else if (q<=2.0) {
            return -(1/pow(h,NDIM+2))*3.0*(WCON/6.0)*pow(2.0-q,2)/q;
         }
         else {
            return 0.0;
         }
#endif
#endif
#endif
#endif
#endif
      }


      inline double K(const double q,const double h) {
#ifdef WENDLAND
         if (q<=2.0) {
            return WCON_WENDLAND*pow(2.0-q,4)*(1.0+2.0*q);
         } else {
            return 0.0;
         }
#else
#ifdef MY_CUBIC
         if (q<=MY_CUBIC_ALPHA) {
            return WCON_MY_CUBIC*(pow(2-q,3) + MY_CUBIC_A*pow(MY_CUBIC_ALPHA-q,3));
         } else if (q<=2.0) {
            return WCON_MY_CUBIC*pow(2-q,3);
         } else {
            return 0.0;
         }
#else
#ifdef PRICE2
         if (q<=PRICE2_BETA) {
            return WCON_PRICE2*(pow(2-q,5) + PRICE2_A*pow(PRICE2_ALPHA-q,5) + PRICE2_B*pow(PRICE2_BETA-q,5));
         } else if (q<=PRICE2_ALPHA) {
            return WCON_PRICE2*(pow(2-q,5) + PRICE2_A*pow(PRICE2_ALPHA-q,5));
         } else if (q<=2.0) {
            return WCON_PRICE2*pow(2-q,5);
         } else {
            return 0.0;
         }
#else
#ifdef QUINTIC
         if (q<=1.0) {
            return WCON_QUINTIC*(pow(3-q,5) - 6*pow(2-q,5) + 15*pow(1-q,5));
         } else if (q<=2.0) {
            return WCON_QUINTIC*(pow(3-q,5) - 6*pow(2-q,5));
         } else if (q<=3.0) {
            return WCON_QUINTIC*pow(3-q,5);
         } else {
            return 0.0;
         }
#else
#ifdef HANN
         if (q<=2.0) {
            return WCON_HANN*(1.0-cos(0.5*PI*(q-2.0)));
         } else {
            return 0.0;
         }
#else
         if (q<=1.0) {
            double q2 = pow(q,2);
            double q3 = q*q2;
            return WCON*(2.0/3.0 - q2 + 0.5*q3);
         }
         else if (q<=2.0) {
            return (WCON/6.0)*pow(2.0-q,3);
         }
         else {
            return 0.0;
         }
#endif
#endif
#endif
#endif
#endif
      }

      inline double W(const double q,const double h) {
         return (1/pow(h,NDIM))*K(q,h);
      }


#endif /* SPH_COMMON_H_ */
