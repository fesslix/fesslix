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

#include "flxglobaldef.h"

#include <cstddef>
#include <cmath>

/**
* @brief a 'double' class that tries to minimize potential numerical summation errors by quantifying this errors
* reference: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
*/
class pdouble {
  private: 
    tdouble v;   // the running sum
    tdouble c;   // the accumulated error; the best-guess for the actual number is 'v+c' !!! (NOT 'v-c', because the sign of c is changed compared to the Wikipedia version)
  public:
    pdouble() : v(ZERO), c(ZERO) {}
    pdouble(const tdouble& t) : v(t), c(ZERO) {}
    pdouble(const pdouble& p) : v(p.v), c(p.c) {}
    pdouble(const tdouble& v, const tdouble &c) : v(v), c(c) {}
    pdouble& operator=(const tdouble& p) {
      v=p;
      c=ZERO;
      return *this;
    }
    pdouble& operator=(const pdouble& p) {
      if (this!=&p) {
        v=p.v;
        c=p.c;
      }
      return *this;
    }
    tdouble cast2double() const { return v+c; }
  
    pdouble& operator+=(const tdouble& right);
    pdouble& operator+=(const pdouble& right) { 
      if (this!=&right) {
        this->operator+=(right.c);
        this->operator+=(right.v);
      } else {
        v *= 2;
        c *= 2;
      }
      return *this;
    }
    pdouble& operator-=(const tdouble& right) {
      return operator+=(-right);
    }
    pdouble& operator-=(const pdouble& right) { 
      if (this!=&right) {
        this->operator-=(right.c);
        this->operator-=(right.v);
      } else {
        v = ZERO;
        c = ZERO;
      }
      return *this;
    }
    pdouble& operator*=(const tdouble& right) {
      v*=right;
      c*=right;
      return *this;
    }
    pdouble& operator*=(const pdouble& right) {
      const tdouble v1 = v;
      const tdouble c1 = c;
      const tdouble v2 = right.v;
      const tdouble c2 = right.c;
      v = c1*c2;
      c = ZERO;
      operator+=(v1*c2);
      operator+=(v2*c1);
      operator+=(v1*v2);
      return *this;
    }
    pdouble& operator/=(const tdouble& right) {
      v/=right;
      c/=right;
      return *this;
    }
    pdouble& operator/=(const pdouble& right) {
      if (this!=&right) {
        v/=right.v;
        c/=right.v;
      } else {
        v = ONE;
        c = ZERO;
      }
      return *this;
    }
    const pdouble operator*(const tdouble& right) const {
      pdouble result(this->v*right,this->c*right);
      return result;
    }
    const pdouble operator*(const pdouble& right) const {
      pdouble result(*this);
      result *= right;
      return result;
    }
    const pdouble operator/(const tdouble& right) {
      pdouble result(this->v/right,this->c/right);
      return result;
    }
    const pdouble operator/(const pdouble& right) {
      return operator/(right.v);
    }
    const pdouble operator+(const pdouble& right) const {
      pdouble result(*this);
      result += right;
      return result;
    }
    const pdouble operator-(const pdouble& right) const {
      pdouble result(*this);
      result -= right;
      return result;
    }

    const tdouble get_v() const { return v; }
    const tdouble get_c() const { return c; }
};


/**
* @brief a 'double' class that tries to minimize potential numerical summation errors by performing the summation in bins
*/
class qdouble {
  private:
    /**
    * @brief if the bin is full, add to this value
    */
    pdouble final_sum;
    /**
    * @brief temporary sum (the bin)
    */
    pdouble t_sum;
    /**
    * @brief number of points per bins
    */
    const size_t Np;
    /**
    * @brief current number of points in the bin
    */
    size_t cNp;
    /**
    * @brief total number of recorded points
    */
    size_t c_total;
    
    inline void increase_it() {
      ++cNp;
      ++c_total;
      if (cNp>=Np) {
        final_sum+=t_sum;
        t_sum = ZERO;
        cNp = 0;
      }
    }
  public:
    /**
    * @param Np ppb=true: number of points per bin; otherwise: estimated total number of points
    */
    qdouble(const size_t NpV,const bool ppb) : final_sum(ZERO), t_sum(ZERO), Np(ppb?NpV:size_t(sqrt(tdouble(NpV)))), cNp(0), c_total(0) {}
    
    const tdouble cast2double() const { return final_sum.cast2double()+t_sum.cast2double(); }
    const pdouble cast2pdouble() const { pdouble t(final_sum);t+=t_sum;return t; }
    const tdouble get_average() const { return cast2double()/c_total; }
    
    qdouble& operator+=(const tdouble& right) {
      t_sum+=right;
      increase_it();
      return *this;
    }
    
    qdouble& operator+=(const pdouble& right) {
      t_sum+=right;
      increase_it();
      return *this;
    }
};

/**
* @brief a 'double' class that estimates the mean and the variance in a summation 
*        it also records minimum and maximum values
*/
class vdouble {
  private: 
    pdouble s_p_ref;      // basis to center summations
    pdouble s_p1;          // sum of numbers
    pdouble s_p2;          // sum of squared numbers centered around s_p2
    tdouble min;           // smallest recorded value 
    tdouble max;           // largest recorded value
    /**
    * @brief number of total points contained in summation so far
    */
    size_t Np;
    bool corrected;        // if s_p2_ref equals s_p1/Np
    bool nan_warning;
    
    void correct_p_ref();
  public:
    vdouble(const bool nan_warning=false) : min(ZERO), max(ZERO), Np(0), corrected(false), nan_warning(nan_warning) {}

    vdouble& operator=(const vdouble& p) {
      if (this!=&p) {
        s_p_ref=p.s_p_ref;
        s_p1=p.s_p1;
        s_p2=p.s_p2;
        min = p.min;
        max = p.max;
        Np = p.Np;
        corrected = p.corrected;
      }
      return *this;
    }
    
    tdouble get_mean() { 
      correct_p_ref();
      return s_p_ref.cast2double(); 
    }
    const pdouble& get_mean_p() {
      correct_p_ref();
      return s_p_ref;
    }
    tdouble get_variance() {
      correct_p_ref();
      return s_p2.cast2double()/(Np-1);
    };
    const pdouble get_variance_p() {
      correct_p_ref();
      return s_p2/(Np-1);
    };
    tdouble get_min() const { return min; }
    tdouble get_max() const { return max; }
    tdouble get_sum() const { return (s_p_ref*Np + s_p1).cast2double(); }
    pdouble get_sum_p() const { return (s_p_ref*Np + s_p1); }
    pdouble get_sum_of_squares_p() { return ( s_p_ref*( (s_p_ref*Np + s_p1)*2 - s_p_ref*Np)  + s_p2);  }
    tdouble get_sum_of_squares() { return get_sum_of_squares_p().cast2double();  }
    size_t get_size() const { return Np; }
    tdouble get_mean_sample(const tdouble y);
    tdouble get_var_sample(const tdouble y);
    
    vdouble& operator+=(const tdouble& right);
    
    void clear();
};



