/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2025 Wolfgang Betz
 *
 * Fesslix is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fesslix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Fesslix.  If not, see <http://www.gnu.org/licenses/>. 
 */

#pragma once

#include "flxmath_interp.h"

#include <cmath>


typedef tdouble (*flx_math_fun_ptr)(const tdouble, void*);
    

/**
* @brief throw away all numbers before the dot
*/
inline const tdouble frac(const tdouble x) {  // digits after the dot
  if ( x >= ZERO ) return x - std::floor(x);
  return std::ceil(x) - x;
}
inline const tdouble arccot(const tdouble x) {
  return PI * ONE/2 - std::atan(x);
}
inline const tdouble arctanh(const tdouble x) {
  return (ONE/2) * std::log((ONE + x) / (ONE - x));
}

/**
* @brief returns the larger value of x and y
*/
inline const tdouble max(const tdouble x, const tdouble y) { return (y>x)?y:x; }
/**
* @brief returns the smaller value of x and y
*/
inline const tdouble min(const tdouble x, const tdouble y) { return (y<x)?y:x; }
/**
* @brief returns ONE if x>0, ZERO if x<0 and ZERO otherwise
*/
inline const tdouble sign(const tdouble x) { return (x>ZERO)?ONE:((x<ZERO)?-ONE:ZERO); }
/**
* @brief returns x with the sign of y
*/
inline const tdouble sign(const tdouble x, const tdouble y) {
  return (y >= ZERO ? (x >= ZERO ? x : -x) : (x >= ZERO ? -x : x));
}
/**
* @brief returns true if x and y have the same sign and false otherwise
*/
inline const bool sign_same(const tdouble x, const tdouble y) {
  return (x*y)>=ZERO;
}

template <typename T>
T pow_int_g(T x, tuint n)
{
  T v = ONE;
  do {
    if(n & 1) v *= x;        // for odd n
    n >>= 1;
    x *= x;
  } while (n);
  return v;
}

template <typename T>
inline T pow_int(T x, const tuint n)
{
  switch(n) {
    case 0:
      return 1;
    case 1:
      return x;
    case 2:
      return x*x;
    case 3:
      return x*x*x;
    case 4:
    {
      T v = x*x;
      return v*v;
    }
    case 5:
    {
      T v = x*x;
      return v*v*x;
    }
    case 6:
    {
      T v = x*x;
      return v*v*v;
    }
    case 7:
    {
      T v = x*x*x;
      return v*v*x;
    }
    case 8:
    {
      T v = x*x;
      T v2 = v*v;
      return v2*v2;
    }
    case 9:
    {
      T v = x*x*x;
      return v*v*v;
    }
    case 10:
    {
      T v = x*x;
      T v2 = v*v;
      return v2*v2*v;
    }
    default:
    {
      return pow_int_g(x,n);
    }
  };
}

template <typename T>
inline T pow2(T x)
{
  return x*x;
}

template <typename T>
inline T pow_b2(T n)
{
  return 1<<n;
}

/**
* @brief calculate the Logarithm of the Gamma Function
* @brief returns the value ln(GAMMA(xx)] for xx > 0
*/
FLXLIB_EXPORT const tdouble GammaLn(const tdouble &xx);
/**
* @brief solves GammaLn(alpha+beta)-GammaLn(beta)
*/
FLXLIB_EXPORT const tdouble GammaLn_diff(const tdouble alpha, const tdouble beta);

const tdouble BetaFunLn(const tdouble &z, const tdouble &w);

inline const tdouble BetaFun(const tdouble &z, const tdouble &w) {
  return std::exp(BetaFunLn(z,w));
}

/**
* @brief error function
*/
inline const tdouble flxerf(const tdouble x) { return erf(x); }
/**
* @brief inverse error function
*/
FLXLIB_EXPORT const tdouble flxerf_inv(const tdouble p);

/**
* @brief gamma function
*/
FLXLIB_EXPORT const tdouble flxgamma(const tdouble x);
/**
* @brief incomplete upper gamma function
*/
FLXLIB_EXPORT const tdouble flxgamma(const tdouble a, const tdouble z);
/**
* @brief incomplete lower gamma function
*/
FLXLIB_EXPORT const tdouble flxgamma_l(const tdouble a, const tdouble z);
/**
* @brief regularized incomplete upper gamma function
*/
FLXLIB_EXPORT const tdouble flxgamma_ru(const tdouble a, const tdouble z);
/**
* @brief regularized incomplete lower gamma function
*/
FLXLIB_EXPORT const tdouble flxgamma_rl(const tdouble a, const tdouble z);
/**
* @brief inverse regularized incomplete upper gamma function
*/
FLXLIB_EXPORT const tdouble flxgamma_ru_inv(const tdouble a, const tdouble q);
/**
* @brief inverse regularized incomplete lower gamma function
*/
FLXLIB_EXPORT const tdouble flxgamma_rl_inv(const tdouble a, const tdouble p);
/**
* @brief regularized incomplete Beta function
*/
FLXLIB_EXPORT const tdouble iBeta_reg(const tdouble alpha, const tdouble beta, const tdouble x);
/**
* @brief complement of regularized incomplete Beta function
*/
FLXLIB_EXPORT const tdouble iBeta_regc(const tdouble alpha, const tdouble beta, const tdouble x);
/**
* @brief inverse of the regularized incomplete Beta function
*/
FLXLIB_EXPORT const tdouble iBeta_reg_inv(const tdouble alpha, const tdouble beta, const tdouble p);
/**
* @brief inverse of the complement of the regularized incomplete Beta function
*/
FLXLIB_EXPORT const tdouble iBetac_reg_inv(const tdouble alpha, const tdouble beta, const tdouble p);
/**
* @brief Digamma function
*/
FLXLIB_EXPORT const tdouble flxdigamma(const tdouble x);


/**
* @brief Factorial Function
*/
FLXLIB_EXPORT const tdouble Factorial(const int n);
FLXLIB_EXPORT const tdouble DoubleFactorial(const int n);
FLXLIB_EXPORT const tdouble FactorialLn(const int n);

/**
* @brief calculate the log of the Binomial Coefficient (n, m)
*/
inline const tdouble BinomCoeff_ln(const int n, const int m) {
  return GammaLn(n+1)-GammaLn(m+1)-GammaLn(n-m+1);
}

/**
* @brief calculate the Binomial Coefficient (n, m)
*/
inline const tdouble BinomCoeff(const int n, const int m)
{
  //return std::floor(ONE/2+std::exp(FactorialLn(n)-FactorialLn(m)-FactorialLn(n-m)));
  return exp(BinomCoeff_ln(n,m));
}

/**
* @brief rounds number to integer
*/
inline const tdouble round_flx(const tdouble r) {
  return (r >= ZERO) ? std::floor(r + ONE/2) : std::ceil(r - ONE/2);
}

/**
* @brief rounds number to 10^-n basis
*/
inline const tdouble round_flx(tdouble r, const tuint n) {
  const tdouble p = std::pow(tdouble(10),n);
  r *= p;
  return ((r >= ZERO) ? std::floor(r + ONE/2) : std::ceil(r - ONE/2))/p;
}

/**
* @brief rounds number to floating point basis
*/
inline const tdouble round_flx_fb(tdouble r, const tuint n) {
  int p10 = std::log10(std::fabs(r));
  if (fabs(r)<ONE) --p10;
  const tdouble p = pow(10,int(n)-p10);
  r *= p;
  return ((r >= ZERO) ? std::floor(r + ONE/2) : std::ceil(r - ONE/2))/p;
}


// ------------------------------------ optimization routines ----------------------------------

/**
* finds the minimum of a 1D function
* @param x_lb lower bound of search interval
* @param x_ub upper bound of search interval
* @param x_min an initial guess for the minimum
* @param fun a pointer to the function to minimize
* @param dp data-pointer that is passed to fun
* @param use_brent true: brent's method, false: golden section search
* @param use_initial true: x_min is a first guess for the optimum
* @param NITER maximum number of iterations
* @param NEX maximum number of iterations when trying to bracket the minimum
* @param eps1 absolute error
* @param eps2 relative error
* @param osp a pointer to a output-stream (can be NULL)
* @returns parameter x_min: the location of the minimum
*/
FLXLIB_EXPORT void flx_optim(tdouble x_lb, tdouble x_ub, tdouble &x_min, flx_math_fun_ptr fun, void* dp, const bool use_brent, const bool use_initial, const tuint NITER=1000, const tuint NEX=100, const tdouble eps1=1e-6, const tdouble eps2=1e-6, ostreamp osp=NULL);


// ------------------------------------ root search routines ----------------------------------

/**
* @brief finds a root of a function using the bisection algorithm
* @return the result of the root-search
* @param fun a pointer to the function to analyze
* @param dp data-pointer that is passed to fun
* @param theconst the root-search variable (const-variable)
* @param start start of the search-interval
* @param end end of the search-interval
* @param dx epsilon for interval size
* @param dy epsilon for function-value
*/
FLXLIB_EXPORT const tdouble flx_RootSearch_Bisec(flx_math_fun_ptr fun, void* dp, tdouble start, tdouble end, const tdouble dx=1e-6, const tdouble dy=1e-8, ostreamp streamp=NULL);
FLXLIB_EXPORT const tdouble flx_RootSearch_RegulaFalsi(flx_math_fun_ptr fun, void* dp, tdouble start, tdouble end, const tdouble dx=1e-6, const tdouble dy=1e-8, ostreamp streamp=NULL);



