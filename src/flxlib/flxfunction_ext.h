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

#ifndef fesslix_flxfunction_ext_H
#define fesslix_flxfunction_ext_H

#include "flxfunction.h"


struct flx_fun_void_data {
  tdouble* theconst;
  FunBase* funR;
  flx_fun_void_data() : theconst(NULL), funR(NULL) {}
};

inline tdouble flx_fun_void(const tdouble x, void *params) {
  flx_fun_void_data *p = (flx_fun_void_data *)params;
  *(p->theconst) = x;
  return p->funR->calc();
}


/**
* @brief creates a FlxPoint out of FlxFunctions
*/
class FLXLIB_EXPORT FlxFunPoint {
  public:
    /**
    * @brief frame of reference (coordinate system)
    */
    enum FoR { 
      cartesian, cylindrical, spherical
    };
    static FoR get_FoR(char strlFoR,bool errSerious=true);
    static char get_FoR(FoR elFoR);
  private:
    FoR FrORef;
    FlxFunction *d1, *d2, *d3;
  public:
    FlxFunPoint(const FlxFunPoint& ref);
    FlxFunPoint(ReadStream *reader, FlxFunctionReader *funReader,bool errSerious=true);
    FlxFunPoint(const FoR FrORef, FlxFunction *d1, FlxFunction  *d2, FlxFunction *d3);
    ~FlxFunPoint();
    FlxFunPoint& operator=(const FlxFunPoint& ref);
    /**
    * @brief returns a pointer to a point
    * note: memory has to be deallocated
    */
    const flxPoint calc() const;
    /**
    * @brief frame of reference (coordinate system)
    * @param dir x, y or z
    */
    const bool is_direction_zero(const char dir);
    
    friend std::ostream& operator<<(std::ostream& os, const FlxFunPoint& val);
};

std::ostream& operator<<(std::ostream& os, const FlxFunPoint& val);


/**
* @brief integrates a function numerically using Gaussian quadrature
* @return the result of the integration
* @param funI the function to integrate
* @param theconst the integration variable (const-variable)
* @param start start of the integral
* @param end end of the integral
* @param GPN the number of Gauss-points to use
*/
FLXLIB_EXPORT const tdouble FlxFun_GaussIntegration(FunBase* funI, tdouble* theconst, const tdouble& start, const tdouble& end, const tuint& GPN, GaussIntegration& GI);

/**
* @brief estimates the number of Gauss-points required in order to integrate funI with an error less than err
* @return number of Gauss-points 
*   returns 0 if funI is constant w.r.t. theconst
*   returns GPmax+1 if the iteration did not converge within the specified bounds
* @param funI the function to integrate
* @param theconst the integration variable (const-variable)
* @param start start of the integral
* @param end end of the integral
* @param GPmax the maximum number of Gauss-points to check the integration-error
* @param GPtest the number of Gauss-points to use as reference solution
* @param err the error (relative error) required to stop the iteration
*/
FLXLIB_EXPORT const tuint FlxFun_EstimateGaussPoint(FunBase* funI, tdouble* theconst, const tdouble& start, const tdouble& end, const tuint& GPmax, const tuint& GPtest, const tdouble& err, GaussIntegration& GI);

/**
* @brief finds a root of a function using the bisection algorithm
* @return the result of the root-search
* @param funR the function to analyze
* @param theconst the root-search variable (const-variable)
* @param start start of the search-interval
* @param end end of the search-interval
* @param dx epsilon for interval size
* @param dy epsilon for function-value
*/
FLXLIB_EXPORT const tdouble FlxFun_RootSearch_Bisec(FunBase* funR, tdouble* theconst, tdouble start, tdouble end, const tdouble dx, const tdouble dy, ostreamp streamp);
FLXLIB_EXPORT const tdouble FlxFun_RootSearch_RegulaFalsi(FunBase* funR, tdouble* theconst, tdouble start, tdouble end, const tdouble dx, const tdouble dy, ostreamp streamp);


#endif // fesslix_flxfunction_ext_H
