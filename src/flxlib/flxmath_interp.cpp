/* Fesslix - Stochastic Analysis
 * Copyright (C) 2010-2026 Wolfgang Betz
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


#include "flxmath.h"
#include <string.h>


const tdouble flx_interpolate_linear(const tdouble x_val, const tdouble* x_ptr, const tdouble* f_ptr, const size_t N)
{
  size_t Nel = N;
  size_t cpos = 0;
  // check beginning and ending of matrix
    if (x_val<=x_ptr[0]) {
      return f_ptr[0];
    }
    if (x_val>=x_ptr[Nel-1]) {
      return f_ptr[Nel-1];
    }
  // find the lower bound of the interpolation interval
    while (Nel>1) {
      size_t npos = cpos + Nel/2;
      if (x_ptr[npos]>x_val) {
        Nel = npos - cpos;
      } else {
        Nel = Nel - (npos - cpos);
        cpos = npos;
      }
    };
  // interpolate
    return f_ptr[cpos] + (f_ptr[cpos+1]-f_ptr[cpos])*((x_val-x_ptr[cpos])/(x_ptr[cpos+1]-x_ptr[cpos]));
}

const size_t flx_interpolate_find_larger_eq(const tdouble x_val, const tdouble* x_ptr, const size_t N)
{
  if (N==0) return 0;
  if (x_val<x_ptr[0]) return 0;
  if (x_val>x_ptr[N-1]) return N;
  size_t start = 0;
  size_t length = N;
  while (length>1) {
    const size_t check = start + (length/2);
    if (x_val<x_ptr[check]) {
      // start = start;
      length = check-start;
    } else {
      length = (start+length)-check;
      start = check;
    }
  };
  return start+1;
}


flx_interp::flx_interp(size_t Nreserve)
: Nreserve(Nreserve), Nsmpl(0), dptr(new tdouble[2*Nreserve])
{

}

flx_interp::~flx_interp()
{
  delete [] dptr;
}

const size_t flx_interp::find_larger_eq(const tdouble x) const
{
  if (Nsmpl==0) return 0;
  if (x<get_x(0)) return 0;
  if (x>get_x(Nsmpl-1)) return Nsmpl;
  size_t start = 0;
  size_t length = Nsmpl;
  while (length>1) {
    const size_t check = start + (length/2);
    if (x<get_x(check)) {
      // start = start;
      length = check-start;
    } else {
      length = (start+length)-check;
      start = check;
    }
  };
  return start+1; 
}

const bool flx_interp::append(const tdouble x, const tdouble fx)
{
  if (Nsmpl>=Nreserve) return false;
  const size_t pos = find_larger_eq(x);
  if (pos<Nsmpl) {
    if (fabs(x-get_x(pos))<1e-6) {
      if (fabs(fx-get_fx(pos))>1e-6) throw FlxException("flx_interp::append","Same value x with different values for fx.");
      return true;
    } else {
      memmove(dptr+2*pos+2,dptr+2*pos,(Nsmpl-pos)*2*sizeof(tdouble));
      dptr[2*pos] = x;
      dptr[2*pos+1] = fx;
    }
  } else {
    dptr[2*pos] = x;
    dptr[2*pos+1] = fx;
  }
  ++Nsmpl;
  return true;
}

const bool flx_interp::find_3p(const tdouble f, const size_t pos, tdouble& xn1, tdouble& xn2) const
{
  if (pos==0) {                        // linear interpolation
    const tdouble x0 = dptr[0];
    const tdouble y0 = dptr[1];
    const tdouble x1 = dptr[2];
    const tdouble y1 = dptr[3];
    xn1 = (f-y0)/(y1-y0)*(x1-x0)+x0;
    if (xn1>=x1) throw FlxException("flx_interp::find_3p_01");
    return false;
  } else if (pos+1==Nsmpl) {        // linear interpolation
    const tdouble x0 = dptr[2*Nsmpl-4];
    const tdouble y0 = dptr[2*Nsmpl-3];
    const tdouble x1 = dptr[2*Nsmpl-2];
    const tdouble y1 = dptr[2*Nsmpl-1];
    xn1 = (f-y0)/(y1-y0)*(x1-x0)+x0;
    if (xn1<x1) throw FlxException("flx_interp::find_3p_02");
    return false;
  } else {                        // quadratic interpolation
    const tdouble x0 = dptr[2*(pos-1)];
    const tdouble y0 = dptr[2*(pos-1)+1]-f;
    const tdouble x1 = dptr[2*pos];
    const tdouble y1 = dptr[2*pos+1]-f;
    const tdouble x2 = dptr[2*pos+2];
    const tdouble y2 = dptr[2*pos+3]-f;
    const tdouble c = y1;
    const tdouble a = (y2-c-(y0-c)/(x0-x1)*(x2-x1))/((x2-x1)*(x2-x0));
    if (fabs(a)<GlobalVar.TOL()) {
      xn1 = y0*(x1-x0)/(y0-y1)+x0;
      return false;
    }
    const tdouble b = (y0-c)/(x0-x1)-a*(x0-x1);
    const tdouble st = pow2(b)-4*a*c;
    if (st<ZERO) throw FlxException("flx_interp::find_3p_03","No root found");
    xn1 = (-b-sqrt(st))/(2*a)+x1;
    xn2 = (-b+sqrt(st))/(2*a)+x1;
    if (xn1>x2||xn1<x0) {
      if (xn2>x2||xn2<x0) {
        throw FlxException("flx_interp::find_3p_04","No root found");
      } else {
        xn1 = xn2;
        return false;
      }
    }
    if (xn2>x2||xn2<x0) {
      return false;
    }
    return true;
  }
}

const tdouble flx_interp::interpolate_3p(const tdouble x, const size_t pos) const
{
  if (pos==0) {                        // linear interpolation
    const tdouble x0 = dptr[0];
    const tdouble y0 = dptr[1];
    const tdouble x1 = dptr[2];
    const tdouble y1 = dptr[3];
    return (x-x1)/(x0-x1)*y0 + (x-x0)/(x1-x0)*y1;
  } else if (pos+1==Nsmpl) {        // linear interpolation
    const tdouble x0 = dptr[2*Nsmpl-4];
    const tdouble y0 = dptr[2*Nsmpl-3];
    const tdouble x1 = dptr[2*Nsmpl-2];
    const tdouble y1 = dptr[2*Nsmpl-1];
    return (x-x1)/(x0-x1)*y0 + (x-x0)/(x1-x0)*y1;
  } else {                        // quadratic interpolation
    const tdouble x0 = dptr[2*(pos-1)];
    const tdouble y0 = dptr[2*(pos-1)+1];
    const tdouble x1 = dptr[2*pos];
    const tdouble y1 = dptr[2*pos+1];
    const tdouble x2 = dptr[2*pos+2];
    const tdouble y2 = dptr[2*pos+3];
    return (x-x1)*(x-x2)/((x0-x1)*(x0-x2))*y0 
          +(x-x0)*(x-x2)/((x1-x0)*(x1-x2))*y1
          +(x-x0)*(x-x1)/((x2-x0)*(x2-x1))*y2;
  }
}

const tdouble flx_interp::interpolate(const tdouble x) const
{
  if (Nsmpl<2) throw FlxException("flx_interp::interpolate_1","Not enough points in the set to interpolate.");
  const size_t pos = find_larger_eq(x);
  if (x==get_x(pos)) {
    return get_fx(pos);
  }
  if (pos==0) {
    return interpolate_3p(x,pos);
  } else if (pos>=Nsmpl) {
    return interpolate_3p(x,Nsmpl-1);
  } else {
    const tdouble w = (x-get_x(pos-1))/(get_x(pos)-get_x(pos-1));
    #if FLX_DEBUG
      if (w>ONE || w<ZERO) throw FlxException_Crude("flx_interp::interpolate_2");
    #endif
    return (ONE-w)*interpolate_3p(x,pos-1)+w*interpolate_3p(x,pos);
  }
}

const tdouble flx_interp::find_1st_x_before_xs_smaller_than_f(const tdouble xs, const tdouble f, const bool smaller) const
{
  if (Nsmpl<2) throw FlxException("flx_interp::find_1st_x_after_xs_smaller_than_f_01","Not enough points in the set.");
  const tdouble fs = interpolate(xs);
  if (smaller?(fs<=f):(fs>=f)) return xs;
  size_t pos = find_larger_eq(xs);
  if (pos==0) throw FlxException_Crude("flx_interp::find_1st_x_after_xs_smaller_than_f_02");
  --pos;
  while (smaller?(get_fx(pos)>f):(get_fx(pos)<f)) {
    if (pos==0) {
      return xs;
    }
    --pos;
  }
  tdouble x1=ZERO,x2=ZERO;
  if (pos==Nsmpl) {
    find_3p(f,pos,x1,x2);
    return x1;
  } else {
    bool b = find_3p(f,pos,x1,x2);
    // ensure that upper root is within bounds
      tdouble xn1;
      if (b) {
        if (x2>get_x(pos+1)) xn1 = x1;
        else xn1 = x2;
      } else {
        xn1 = x1;
      }
    b = find_3p(f,pos+1,x1,x2);
    // ensure that upper root is within bounds
      tdouble xn2;
      if (b) {
        if (x2>get_x(pos+1)) xn2 = x1;
        else xn2 = x2;
      } else {
        xn2 = x1;
      }
    // make sure that both roots are within bounds
      if (xn1>get_x(pos+1) || xn2>get_x(pos+1)) throw FlxException_Crude("flx_interp::find_1st_x_after_xs_smaller_than_f_03");
    const tdouble xn = (xn1+xn2)/2;
    const tdouble w = (xn-get_x(pos))/(get_x(pos+1)-get_x(pos));
    #if FLX_DEBUG
      if (w>ONE || w<ZERO) throw FlxException_Crude("flx_interp::find_1st_x_after_xs_smaller_than_f_04");
    #endif
    return (ONE-w)*xn1+w*xn2;
  }
}





