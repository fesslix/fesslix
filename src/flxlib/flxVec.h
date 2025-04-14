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

#include "flxmath.h"
#include <string.h>


class flxpVec;

class FLXLIB_EXPORT flxVec {
  private:
    const tuint N;
    tdouble* vptr;
    const bool is_ptr;
    
  public:
    flxVec(const tuint N);
    flxVec(const flxVec& rhs);
    flxVec(const flxpVec& rhs);
    flxVec() = delete;
    flxVec(flxVec&& rhs);
    /**
    * @param copy_mem: true: copy memory, false: memory is only reference (and not managed)
    * @param manage_mem: true: memory is managed; otherwise, not (only relevant if copy_mem is false!)
    */
    flxVec(tdouble* ptr, const tuint& NV, const bool copy_mem=false, const bool manage_mem=false);
    flxVec(const tdouble* ptr, const tuint& NV, const bool copy_mem=false);
    ~flxVec();
    
    const tuint get_N() const { return N; }
    tdouble* get_tmp_vptr() { return vptr; }
    const tdouble* get_tmp_vptr_const() const { return vptr; }
    
    void set_zero() { for (tuint i=0;i<N;++i) vptr[i] = ZERO; }
    void set_nan();
    
    flxVec& operator=(const flxVec& rhs) {
      if (vptr!=rhs.vptr && N>0) {
        #if FLX_DEBUG
          if (N!=rhs.N) {
            throw FlxException_Crude("flxVec::operator=_1");
          }
        #endif
        memcpy(vptr,rhs.vptr,N*sizeof(tdouble));      
      }
      return *this;
    }
    
    flxVec& operator=(const flxpVec& rhs);
    
    flxVec& operator=(const tdouble& rhs) {
      for (tuint i=0;i<N;++i) vptr[i] = rhs;
      return *this;
    }
    
    flxVec& operator*=(const tdouble& rhs) {
      for (tuint i=0;i<N;++i) vptr[i]*=rhs;
      return *this;
    }
    flxVec& operator/=(const tdouble& rhs) {
      for (tuint i=0;i<N;++i) vptr[i]/=rhs;
      return *this;
    }
    
    flxVec& operator+=(const flxVec& rhs) {
      const tdouble* tp = rhs.get_tmp_vptr_const();
      for (tuint i=0;i<N;++i) vptr[i]+=tp[i];
      return *this;
    }
    
    flxVec& operator-=(const flxVec& rhs) {
      const tdouble* tp = rhs.get_tmp_vptr_const();
      for (tuint i=0;i<N;++i) vptr[i]-=tp[i];
      return *this;
    }
    
    flxVec& operator+=(const tdouble& rhs) {
      for (tuint i=0;i<N;++i) vptr[i]+=rhs;
      return *this;
    }
    
    flxVec& operator-=(const tdouble& rhs) {
      for (tuint i=0;i<N;++i) vptr[i]-=rhs;
      return *this;
    }
    
    /**
    * @brief *this += (vector*scalar)
    */
    flxVec& add(const flxVec& v, const tdouble s) {
      #if FLX_DEBUG
        if (N!=v.N) throw FlxException_Crude("flxVec::add");
      #endif
      const tdouble* tp = v.get_tmp_vptr_const();
      for (tuint i=0;i<N;++i) vptr[i]+=tp[i]*s;
      return *this;
    }
    
    flxVec& mult_coeff(const flxVec& rhs1, const flxVec& rhs2);
    
    const tdouble operator*(const flxVec& rhs) const {
      #if FLX_DEBUG
        if (N!=rhs.N) {
          throw FlxException_Crude("flxVec::operator*");
        }
      #endif
      #if FLX_KAHAN_DOT
        pdouble tm=ZERO;
      #else
        tdouble tm=ZERO;
      #endif
      const tdouble* v1p = vptr;
      const tdouble* v2p =rhs.get_tmp_vptr_const();
      for (tuint i=0;i<N;++i) {
        tm+=v1p[i]*v2p[i];
      }
      #if FLX_KAHAN_DOT
        return tm.cast2double();
      #else
        return tm;
      #endif
    }
    
    tdouble& operator[](const tuint idx) { 
      #if FLX_DEBUG
        if (idx>=N) { 
          throw FlxException_Crude("flxVec::operator[]_1"); 
        }
      #endif
      return vptr[idx];
    }
    
    const tdouble& operator[](const tuint idx) const { 
      #if FLX_DEBUG
        if (idx>=N) {
          throw FlxException_Crude("flxVec::operator[]_2");
        }
      #endif
      return vptr[idx];
    }
    
    const bool operator==(const flxVec& rhs) const;
    const bool operator!=(const flxVec& rhs) const { return !operator==(rhs); }
    
    const tdouble get_NormMax() const {
      tdouble tm=fabs(vptr[0]);
      for (tuint i=1;i<N;++i) {
        if (tm < fabs(vptr[i])) {
          tm=fabs(vptr[i]);
        }
      }
      return tm;
    }
    
    /**
    * @brief computes the distance between this vector and rhs
    * @returns distance squared
    */
    const tdouble comp_dist_NOroot(const flxVec& rhs) const {
      #if FLX_KAHAN_2NORM
        pdouble tm=ZERO;
      #else
        tdouble tm=ZERO;
      #endif
      const tdouble* rhsp = rhs.get_tmp_vptr_const();
      for (tuint i=0;i<N;++i) {
        tm+=pow2(vptr[i]-rhsp[i]);
      }
      #if FLX_KAHAN_2NORM
        return tm.cast2double();
      #else
        return tm;
      #endif
    }
    
    const tdouble get_Norm2_NOroot() const {
      #if FLX_KAHAN_2NORM
        pdouble tm=ZERO;
      #else
        tdouble tm=ZERO;
      #endif
      for (tuint i=0;i<N;++i) {
        tm+=vptr[i]*vptr[i];
      }
      #if FLX_KAHAN_2NORM
        return tm.cast2double();
      #else
        return tm;
      #endif
    }

    /**
    * @brief sorts the components of this vector
    */
    void sort();
    
    const tdouble get_Norm2() const { return sqrt(get_Norm2_NOroot()); }
    /**
    * @brief computes the squared distance of the difference vector (*this-rhs)
    */
    const tdouble get_Norm2_NOroot_of_distance(const flxVec& rhs) const;
    const tdouble get_max() const;
    const tdouble get_min() const;
    const size_t get_maxID() const;
    const size_t get_minID() const;
    const tdouble get_sum() const;
    const tdouble get_absMean() const;
    const tdouble get_Mean() const;
    const tdouble get_Var(const tdouble& mean) const;
    const tdouble get_sd(const tdouble& mean) const;
    flxVec& normalize();

    const size_t count_nan() const;
    void copy_vals_without_nan(flxVec& rhs) const;
    
    void check_TOL() const {
      const tdouble GT = get_NormMax() * GlobalVar.TOL();
      for (tuint i=0;i<N;++i) {
        if (fabs(vptr[i]) <= GT) vptr[i] = ZERO;
      }
    }
    
    void swap(flxVec& sV);
    /**
    * @brief pics N elements with distance slicegap from the vector ptr
    */
    flxVec& slice(const tdouble* ptr, const tuint stride);
    void cast2tuint(iVec& R, const bool No0);
    /**
    * @brief sets the vector equal to 'rhs' -> check dimensions of both vectors first!
    */
    void assign_save(const flxVec& rhs);
    const bool check_size(const tuint NV, const bool throwErr=true) const;

    struct Iterator {   // based on https://www.internalpointers.com/post/writing-custom-iterators-modern-cpp
      using iterator_category  = std::forward_iterator_tag;
      using difference_type     = std::ptrdiff_t;
      using value_type         = tdouble;
      using pointer            = tdouble*;
      using reference          = tdouble&;

      Iterator(pointer ptr) : m_ptr(ptr) {}

      reference operator*() const { return *m_ptr; }
      pointer operator->() { return m_ptr; }

      // Prefix increment
      Iterator& operator++() { m_ptr++; return *this; }

      // Postfix increment
      Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }

      friend bool operator== (const Iterator& a, const Iterator& b) { return a.m_ptr == b.m_ptr; };
      friend bool operator!= (const Iterator& a, const Iterator& b) { return a.m_ptr != b.m_ptr; };

      private:
        pointer m_ptr;
    };

    Iterator begin() { return Iterator(&vptr[0]); }
    Iterator end()   { return Iterator(&vptr[N]); }

};

typedef flxVec* flxVec_ptr;


class FLXLIB_EXPORT flxpVec {
  private:
    const tuint N;
    pdouble* vptr;
    const bool is_ptr;
    
  public:
    flxpVec(const tuint N);
    flxpVec(const flxpVec& rhs);
    flxpVec(const flxVec& rhs);
    flxpVec(pdouble* ptr, const tuint N);
    flxpVec(const pdouble* ptr, const tuint N);
    ~flxpVec();
    
    const tuint get_N() const { return N; }
    pdouble* get_tmp_vptr() const { return vptr; }
    const pdouble* get_tmp_vptr_const() const { return vptr; }
    
    void set_zero() { for (tuint i=0;i<N;++i) vptr[i] = ZERO; }
    
    flxpVec& operator=(const flxpVec& rhs);
    
    flxpVec& operator+=(const flxVec& rhs);
    flxpVec& operator+=(const flxpVec& rhs);
    flxpVec& operator-=(const flxVec& rhs);
    flxpVec& operator-=(const flxpVec& rhs);
    
    /**
    * @brief *this += (vector*scalar)
    */
    flxpVec& add(const flxVec& v, const tdouble s);
    flxpVec& add(const flxpVec& v, const tdouble s);
    
    flxpVec& operator*=(const tdouble& rhs);
    flxpVec& operator*=(const pdouble& rhs);
    
    const pdouble operator*(const flxpVec& rhs) const;
    
    pdouble& operator[](const tuint idx) { 
      #if FLX_DEBUG
        if (idx>=N) throw FlxException_Crude("flxpVec::operator[]");
      #endif
      return vptr[idx];
    }
    
    const pdouble& operator[](const tuint idx) const { 
      #if FLX_DEBUG
        if (idx>=N) throw FlxException_Crude("flxpVec::operator[]_2");
      #endif
      return vptr[idx];
    }

    const tdouble get_NormMax() const;
    
    const tdouble get_Norm2_NOroot() const;
    
    const tdouble get_Norm2() const { return sqrt(get_Norm2_NOroot()); }
    
    void check_TOL() const;
    
};


// ---------------------------------------------------------------------------------------------------------------

/**
* @brief finds the smallest index in smpl_ref_list which is larger than e
* @param vp a sorted tdouble-vector
* @param N length of vp
* @param p the percentile to return
* 
*/
FLXLIB_EXPORT const tdouble flx_percentile(const tdouble* const vp, const tuint N, const tdouble p, const bool inverse=false);

/**
* @brief finds the smallest index in smpl_ref_list which is larger than e
* @param vp a sorted tdouble-vector
* @param N length of vp
* @param e the value to search for
* @returns N if g>=largest entry
*/
FLXLIB_EXPORT const tuint flx_find_pos(const tdouble* const vp, const tuint N, const tdouble e);

/**
* @brief finds the smallest index in smpl_ref_list which is larger or equal than e
* @param vp a sorted tdouble-vector
* @param N length of vp
* @param e the value to search for
* @returns N if g>=largest entry
*/
FLXLIB_EXPORT const tuint flx_find_pos2(const tdouble* const vp, const tuint N, const tdouble e);

FLXLIB_EXPORT const tdouble calc_distance(const tdouble* const v1, const tdouble* const v2, const tuint N);


// -------------------------- VECTOR - tVec ----------------------------------------------------------------------

inline const tdouble Norm2(const tdouble* v, const tuint N) {
  #if FLX_KAHAN_2NORM
    pdouble tm=0.0;
  #else
    tdouble tm=0.0;
  #endif
  for (size_t i=0;i<N;++i) {
    tm+=v[i]*v[i];
  }
  #if FLX_KAHAN_2NORM
    return std::sqrt(tm.cast2double());
  #else
    return std::sqrt(tm);
  #endif
}



FLXLIB_EXPORT void flxVec_simple_plot(std::ostream& os, const flxVec& V, const bool checkTOL, const int prec, const int fixW, const bool brackets=false);
FLXLIB_EXPORT void flxVec_totalPrec_plot(std::ostream& os, const flxVec& V);
FLXLIB_EXPORT std::ostream& operator<<(std::ostream& os, const flxVec& V);
FLXLIB_EXPORT std::ostream& operator<<(std::ostream& os, const bVec& V);
FLXLIB_EXPORT std::ostream& operator<<(std::ostream&os, const std::vector< tuint >& V);

FLXLIB_EXPORT std::ostream& operator<<(std::ostream& os, const pdouble& p);



