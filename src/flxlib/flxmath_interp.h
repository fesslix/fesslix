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

#include "flxglobal.h"


class FLXLIB_EXPORT flx_interp {
  private:
    size_t Nreserve;
    size_t Nsmpl;
    tdouble* dptr;
    
    const tdouble get_x(const size_t i) const { return dptr[2*i]; }
    const tdouble get_fx(const size_t i) const { return dptr[2*i+1]; }
    
    const bool find_3p(const tdouble f,const size_t pos, tdouble& xn1, tdouble& xn2) const;
    const tdouble interpolate_3p(const tdouble x,const size_t pos) const;
  public:
    flx_interp(size_t Nreserve);
    ~flx_interp();
    
    /**
    * @returns the index of the first element larger or equal than x
    */
    const size_t find_larger_eq(const tdouble x) const;
    /**
    * @brief Appends x and fx=f(x) to the list
    * @returns true, if the value could be appended
    */
    const bool append(const tdouble x, const tdouble fx);
    /**
    * @brief interpolates the function at position x
    */
    const tdouble interpolate(const tdouble x) const;
    
    const tdouble find_1st_x_before_xs_smaller_than_f(const tdouble xs, const tdouble f, const bool smaller=true) const;
        
    void reset() { Nsmpl=0; }
};

