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

#define FLXLIB_CPP

#include "pdouble.h"
#include "flxexception.h"
#include "flxglobal.h"
#include "flxmath_rnd.h"


pdouble& pdouble::operator+=(const tdouble& right) {
    if (std::isfinite(v+right)) {
        if (fabs(v)>=fabs(right)) {
            tdouble y = right + c;
            tdouble t = v + y;
            c = (v-t)+y;
            v = t;
        } else {
            tdouble y = v + c;
            tdouble t = right + y;
            c = (right-t)+y;
            v = t;
        }
    } else {
        v += right;
        c = ZERO;
    }
    return *this;
}

void vdouble::correct_p_ref() {
    if (!corrected) {
        if (std::isinf(s_p1.cast2double())) {
            s_p_ref = s_p1;
            s_p2 = std::numeric_limits<double>::infinity();
        } else {
            // determin the corrected mean value
                const pdouble b(s_p1/Np + s_p_ref);
            // correct the centered mean sum
                const pdouble s_diff = (s_p_ref-b);
                s_p1 += (s_p_ref-b)*Np;   // ideally, this should be zero at this point
            // correct the sum of squares
                s_p2 += s_diff*( s_p1*2 - s_diff*Np );
                s_p_ref = b;
                if (nan_warning) {
                    if (std::isnan(s_p_ref.cast2double())) {
                        throw FlxException("vdouble::correct_p_ref_01", "Mean is nan.");
                    }
                }
        }
        corrected = true;
        #if FLX_DEBUG
            if (s_p2.cast2double() < ZERO) {
                throw FlxException_Crude("vdouble::correct_p_ref_02");
            }
        #endif
    }
}

tdouble vdouble::get_mean_sample(const tdouble y)
{
    // get properties
        const tdouble m = this->get_mean();
        const tdouble s = sqrt(this->get_variance());
        const tdouble N = this->get_size();
    // get sample of Student's t distribution
        tdouble t;
        if (y<=ZERO) {
            t = rv_InvCDF_Studentst(N-ONE,rv_Phi(y));
        } else {
            t = -rv_InvCDF_Studentst(N-ONE,rv_Phi(-y));
        }
    // generate final sample
        return m + t*s/sqrt(N);
}

tdouble vdouble::get_var_sample(const tdouble y)
{
    // assuming that sample variance is correlated with sample mean (which is actually not the case)
        // const tdouble var = this->get_variance() * pow2(mu/this->get_mean());
    // assuming that the sample variance and the sample mean are independent (which only holds if the underlying sample follow a Normal distribution)
        const tdouble var = this->get_variance();
    // assuming an underlying Normal distribution (which is only the case for large N)
        const tdouble delta = sqrt(ONE*2)/(this->get_size()-1);
        return var *(ONE+delta*y);
}

vdouble& vdouble::operator+=(const tdouble& right) {
    ++Np;
    corrected = false;
    if (nan_warning) {
        if (std::isnan((right))) {
            throw FlxException("vdouble::operator+=_01", "Right-hand side is nan.");
        }
    }
    // update min & max
    if (Np==1) {
        min = right;
        max = right;
    } else {
        if (right<min) {
        min = right;
        }
        if (right>max) {
        max = right;
        }
    }
    // update the centered sum
    if (!(std::isinf(s_p1.cast2double()))) {
        pdouble p(s_p_ref);
        p -= right;
        s_p1 -= p;
        // update the sum of squares
        if (std::isinf(right)) {
            s_p2 = std::numeric_limits<double>::infinity();
        } else {
            p *= p;
            #if FLX_DEBUG
                if (p.cast2double()<ZERO) {
                    throw FlxException_Crude("vdouble::operator+=_02");
                }
            #endif
            if (std::isfinite(p.cast2double())) {
                s_p2 += p;
            } else {
                s_p2 = std::numeric_limits<double>::infinity();
            }
        }
    }
    // correct the basis of the sum of squares
    if (Np==1 || Np==10 || Np==100 || Np==1000 || Np==10000) {
        correct_p_ref();
    }
    return *this;
}

void vdouble::clear() {
      s_p_ref = ZERO;
      s_p1 = ZERO;
      s_p2 = ZERO;
      min = ZERO;
      max = ZERO;
      Np = 0;
      corrected = false;
    }
