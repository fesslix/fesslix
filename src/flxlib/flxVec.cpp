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

#define FLXLIB_CPP

#include "flxVec.h"



flxVec::flxVec(const tuint N)
: N(N), vptr((N>0)?(new tdouble[N]):nullptr), is_ptr(false)
{
  set_zero();
}

flxVec::flxVec(const flxVec& rhs)
: N(rhs.N), vptr((N>0)?(new tdouble[N]):nullptr), is_ptr(false)
{
  memcpy(vptr,rhs.vptr,N*sizeof(tdouble));     
}

flxVec::flxVec(flxVec&& rhs)
: N(rhs.N), vptr(rhs.vptr), is_ptr(rhs.is_ptr)
{
  rhs.vptr = nullptr;
}

flxVec::flxVec(const flxpVec& rhs)
: N(rhs.get_N()), vptr((N>0)?(new tdouble[N]):nullptr), is_ptr(false)
{
  const pdouble* vp = rhs.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) {
    vptr[i] = vp[i].cast2double();
  }
}

flxVec::flxVec(tdouble* ptr, const tuint& NV, const bool copy_mem, const bool manage_mem)
: N(NV), vptr(ptr), is_ptr(copy_mem?false:(!manage_mem))
{
  if (copy_mem) {
    vptr = (N>0)?(new tdouble[N]):nullptr;
    memcpy(vptr,ptr,N*sizeof(tdouble));    
  }
}

flxVec::flxVec(const tdouble* ptr, const tuint& NV, const bool copy_mem)
:N(NV), vptr(const_cast<tdouble*>(ptr)), is_ptr(!copy_mem)
{
  if (copy_mem) {
    vptr = (N>0)?(new tdouble[N]):nullptr;
    memcpy(vptr,ptr,N*sizeof(tdouble));    
  }
}

flxVec::~flxVec()
{
  if (!is_ptr) if (vptr) delete [] vptr;
}

flxVec& flxVec::operator=(const flxpVec& rhs)
{
  #if FLX_DEBUG
    if (N!=rhs.get_N()) throw FlxException_Crude("flxVec::operator=_2");
  #endif
  const pdouble* tp = rhs.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) {
    vptr[i] = tp[i].cast2double();
  }
  return *this;
}

void flxVec::set_nan()
{
  for (tuint i=0;i<N;++i) vptr[i] = NAN;
}

flxVec& flxVec::mult_coeff(const flxVec& rhs1, const flxVec& rhs2)
{
  #if FLX_DEBUG
    if (N!=rhs1.get_N() || N!=rhs2.get_N()) throw FlxException_Crude("flxVec::operator=_2");
  #endif
  const tdouble* p1 = rhs1.get_tmp_vptr_const();
  const tdouble* p2 = rhs2.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) {
    vptr[i] = p1[i]*p2[i];
  }
  return *this;
}

const bool flxVec::operator==(const flxVec& rhs) const
{
  #if FLX_DEBUG
    if (N!=rhs.get_N()) throw FlxException_Crude("flxVec::operator==");
  #endif
  for (tuint i=0;i<N;++i) {
    if ( fabs(vptr[i]-rhs.vptr[i])>GlobalVar.TOL() ) return false;
  }
  return true;
}

void flxVec::sort()
{
  std::sort(vptr,vptr+N);
}

const tdouble flxVec::get_Norm2_NOroot_of_distance(const flxVec& rhs) const
{
  #if FLX_DEBUG
    if (N!=rhs.get_N() ) throw FlxException_Crude("flxVec::get_Norm2_NOroot_of_distance");
  #endif
  #if FLX_KAHAN_2NORM
      pdouble tm=ZERO;
    #else
      tdouble tm=ZERO;
    #endif
    for (tuint i=0;i<N;++i) {
      tm+=pow2(vptr[i]-rhs[i]);
    }
    #if FLX_KAHAN_2NORM
      return tm.cast2double();
    #else
      return tm;
    #endif
}

const tdouble flxVec::get_max() const
{
  tdouble m = vptr[0];
  for (size_t i=1;i<N;++i) {
    if (vptr[i]>m) m = vptr[i];
  }
  return m;
}

const tdouble flxVec::get_min() const
{
  tdouble m = vptr[0];
  for (size_t i=1;i<N;++i) {
    if (vptr[i]<m) m = vptr[i];
  }
  return m;
}

const size_t flxVec::get_maxID() const
{
  tdouble m = vptr[0];
  size_t c = 0;
  for (size_t i=1;i<N;++i) {
    if (vptr[i]>m) {
      m = vptr[i];
      c = i;
    }
  }
  return c;
}

const size_t flxVec::get_minID() const
{
  tdouble m = vptr[0];
  size_t c = 0;
  for (size_t i=1;i<N;++i) {
    if (vptr[i]<m) {
      m = vptr[i];
      c = i;
    }
  }
  return c;
}

const tdouble flxVec::get_sum() const
{
//   qdouble q(N,false);   (very inefficient for TMCMC!!!)
  tdouble q=ZERO;
  for (size_t i=0;i<N;++i) {
    q += vptr[i];
  }
  return q;
//   return q.cast2double();
}

const tdouble flxVec::get_absMean() const
{
  tdouble m = fabs(vptr[0]);
  for (tuint i=1;i<N;++i) {
    m += fabs(vptr[i]);
  }
  m /= N;
  return m;
}

const tdouble flxVec::get_Mean() const
{
  qdouble q(N,false);
  for (size_t i=0;i<N;++i) {
    q += vptr[i];
  }
  tdouble p = q.cast2double();
  p /= N;
  return p;
}

const tdouble flxVec::get_Var(const tdouble& mean) const
{
  const tdouble TOL(1e-150);
  qdouble q(N,false);
  for (size_t i=0;i<N;++i) {
    q += pow2(vptr[i]-mean);
  }
  pdouble p = q.cast2pdouble();
  p /= (N-ONE);
  const tdouble var = p.cast2double();
  if (var<TOL && mean<TOL && mean > ZERO) {
    qdouble q2(N,false);
    for (size_t i=0;i<N;++i) {
      q2 += pow2((vptr[i]-mean)/mean);
    }
    pdouble p2 = q2.cast2pdouble();
    p2 /= (N-ONE);
    const tdouble var2 = p2.cast2double()*pow2(mean);
    return var2;
  } else {
    return var;
  }
}

const tdouble flxVec::get_sd(const tdouble& mean) const
{
  const tdouble TOL(1e-150);
  qdouble q(N,false);
  for (size_t i=0;i<N;++i) {
    q += pow2(vptr[i]-mean);
  }
  pdouble p = q.cast2pdouble();
  p /= (N-ONE);
  const tdouble sd = sqrt(p.cast2double());
  if (sd<TOL && mean<TOL && mean > ZERO) {
    qdouble q2(N,false);
    for (size_t i=0;i<N;++i) {
      q2 += pow2((vptr[i]-mean)/mean);
    }
    pdouble p2 = q2.cast2pdouble();
    p2 /= (N-ONE);
    const tdouble sd2 = sqrt(p2.cast2double())*mean;
    return sd2;
  } else {
    return sd;
  }
}

flxVec& flxVec::normalize()
{
  const tdouble l = get_Norm2();
  return operator/=(l);
}

const size_t flxVec::count_nan() const
{
  size_t c = 0;
  for (tuint i=0;i<N;++i) {
    if (std::isnan(vptr[i])) {
      ++c;
    }
  }
  return c;
}

void flxVec::copy_vals_without_nan(flxVec& rhs) const
{
  const tuint N_rhs = rhs.get_N();
  tdouble* p_rhs = rhs.get_tmp_vptr();
  size_t c = 0;
  for (tuint i=0;i<N;++i) {
    if (std::isnan(vptr[i])==false) {
      if (c>=N_rhs) throw FlxException_Crude("flxVec::copy_vals_without_nan_01");
      p_rhs[c++] = vptr[i];
    }
  }
  if (c!=N_rhs) throw FlxException_Crude("flxVec::copy_vals_without_nan_02");
}

void flxVec::swap(flxVec& sV)
{
  #if FLX_DEBUG
    if (get_N()!=sV.get_N()) throw FlxException_Crude("flxVec::swap");
  #endif
  tdouble* y1 = get_tmp_vptr();
  tdouble* y2 = sV.get_tmp_vptr();
  tdouble t;
  for (tuint i=0;i<N;++i) {
    t = y1[i];
    y1[i] = y2[i];
    y2[i] = t;
  }
}

flxVec& flxVec::slice(const tdouble* ptr, const tuint stride)
{
  const tdouble* ps = ptr;
  for (tuint i=0;i<N;++i) {
    vptr[i] = *ps;
    ps += stride;
  }
  return *this;
}

void flxVec::cast2tuint(iVec& R, const bool No0)
{
  #if FLX_DEBUG
    if (N != R.size()) {
      std::ostringstream ssV;
      ssV << "Vector sizes do not match.";
      throw FlxException("tVec_cast2tuintNo0_1", ssV.str() );
    }
  #endif
  tuint* Rp = &R[0];
  for (tuint i=0;i<N;++i) {
    if ( vptr[i] < 0.0) { 
      std::ostringstream ssV; ssV << "Number must not be negative ([" << i+1 << "]=" << vptr[i] << ").";
      FlxException("tVec_cast2tuintNo0_2", "Expected unsigned integer!", ssV.str() );     
    } else {
      Rp[i] = static_cast<tuint>(vptr[i]);
      if (No0 && Rp[i]==0) {
        std::ostringstream ssV; ssV << "Number must not be zero ([" << i+1 << "]).";
        FlxException("tVec_cast2tuintNo0_3", "Expected non-zero integer!", ssV.str() );     
      }
    }
  }
}

void flxVec::assign_save(const flxVec& rhs)
{
  rhs.check_size(get_N());
  operator=(rhs);
}

const bool flxVec::check_size(const tuint NV, const bool throwErr) const
{
  const bool res = (get_N()==NV);
  if (!res && throwErr) {
    std::ostringstream ssV;
    ssV << "The size of the vector (" << get_N() << ") does not match the required size (" << NV << ").";
    throw FlxException("flxVec::check_size", "Vector has wrong size", ssV.str() );
  }
  return res;
}




flxpVec::flxpVec(const tuint N)
: N(N), vptr(new pdouble[N]), is_ptr(false)
{
  set_zero();
}

flxpVec::flxpVec(const flxpVec& rhs)
: N(rhs.N), vptr(new pdouble[N]), is_ptr(false)
{
  memcpy((char*)vptr,(char*)rhs.vptr,N*sizeof(pdouble));     
}

flxpVec::flxpVec(const flxVec& rhs)
: N(rhs.get_N()), vptr(new pdouble[N]), is_ptr(false)
{
  const tdouble* vp = rhs.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) {
    vptr[i] = vp[i];
  }
}

flxpVec::flxpVec(pdouble* ptr, const tuint N)
: N(N), vptr(ptr), is_ptr(true)
{

}

flxpVec::flxpVec(const pdouble* ptr, const tuint N)
:N(N), vptr(const_cast<pdouble*>(ptr)), is_ptr(true)
{
  
}

flxpVec::~flxpVec()
{
  if (!is_ptr) delete [] vptr;
}

flxpVec& flxpVec::operator=(const flxpVec& rhs)
{
  if (this!=&rhs) {
    #if FLX_DEBUG
      if (N!=rhs.N) throw FlxException_Crude("flxpVec::operator=");
    #endif
    memcpy((char*)vptr,(char*)rhs.vptr,N*sizeof(pdouble));
  }
  return *this;
}

flxpVec& flxpVec::operator+=(const flxVec& rhs)
{
  #if FLX_DEBUG
    if (N!=rhs.get_N()) throw FlxException_Crude("flxpVec::operator+=_1");
  #endif
  const tdouble* tp = rhs.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) vptr[i]+=tp[i];
  return *this;
}
    
flxpVec& flxpVec::operator+=(const flxpVec& rhs)
{
  #if FLX_DEBUG
    if (N!=rhs.N) throw FlxException_Crude("flxpVec::operator+=_2");
  #endif
  const pdouble* tp = rhs.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) vptr[i]+=tp[i];
  return *this;
}

flxpVec& flxpVec::operator-=(const flxVec& rhs)
{
  #if FLX_DEBUG
    if (N!=rhs.get_N()) throw FlxException_Crude("flxpVec::operator-=_1");
  #endif
  const tdouble* tp = rhs.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) vptr[i]-=tp[i];
  return *this;
}

flxpVec& flxpVec::operator-=(const flxpVec& rhs)
{
  #if FLX_DEBUG
    if (N!=rhs.N) throw FlxException_Crude("flxpVec::operator-=_2");
  #endif
  const pdouble* tp = rhs.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) vptr[i]-=tp[i];
  return *this;
}

flxpVec& flxpVec::add(const flxVec& v, const tdouble s)
{
  #if FLX_DEBUG
    if (N!=v.get_N()) throw FlxException_Crude("flxpVec::add_1");
  #endif
  const tdouble* tp = v.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) vptr[i]+=tp[i]*s;
  return *this;
}

flxpVec& flxpVec::add(const flxpVec& v, const tdouble s)
{
  #if FLX_DEBUG
    if (N!=v.get_N()) throw FlxException_Crude("flxpVec::add_2");
  #endif
  const pdouble* tp = v.get_tmp_vptr_const();
  for (tuint i=0;i<N;++i) vptr[i]+=tp[i]*s;
  return *this;
}

flxpVec& flxpVec::operator*=(const tdouble& rhs)
{
  for (tuint i=0;i<N;++i) vptr[i]*=rhs;
  return *this;
}

flxpVec& flxpVec::operator*=(const pdouble& rhs)
{
  for (tuint i=0;i<N;++i) vptr[i]*=rhs;
  return *this;
}

const pdouble flxpVec::operator*(const flxpVec& rhs) const
{
  #if FLX_DEBUG
    if (N!=rhs.N) throw FlxException_Crude("flxpVec::operator*");
  #endif
  pdouble tm=0.0;
  const pdouble* v1p = vptr;
  const pdouble* v2p =rhs.get_tmp_vptr();
  for (tuint i=0;i<N;++i) {
    tm+=v1p[i]*v2p[i];
  }
  return tm;
}

const tdouble flxpVec::get_NormMax() const
{
  tdouble tm=fabs(vptr[0].cast2double());
  for (tuint i=1;i<N;++i) {
    if (tm < fabs(vptr[i].cast2double())) {
      tm=fabs(vptr[i].cast2double());
    }
  }
  return tm;
}

const tdouble flxpVec::get_Norm2_NOroot() const
{
  pdouble tm=0.0;
  for (tuint i=0;i<N;++i) {
    tm+=vptr[i]*vptr[i];
  }
  return tm.cast2double();
}

void flxpVec::check_TOL() const
{
  const tdouble GT = get_NormMax() * GlobalVar.TOL();
  for (tuint i=0;i<N;++i) {
    if (fabs(vptr[i].cast2double()) <= GT) vptr[i] = 0.0;
  }
}


const tdouble flx_percentile(const tdouble*const vp, const tuint N, const tdouble p, const bool inverse)
{
  // the method implemented is: https://en.wikipedia.org/wiki/Percentile#NIST_method
  // calculate the rank
    const tdouble n = p*(N+1);
  // split up the rank in an integer component and a decimal part
    tdouble intpart;
    tdouble fractpart = modf (n , &intpart);
    tuint ip2 = (tuint)intpart;
  // calculate the percentile
    if (inverse) {
      if (ip2>=N) ip2 = 0;
      else ip2 = N-ip2-1;


      if (ip2 == 0) {
        return vp[N-1];
      } else if (ip2>=N) {
        return vp[0];
      } else {
        return vp[ip2] + fractpart*(vp[ip2-1]-vp[ip2]);
      }
    } else {
      if (ip2 == 0) {
        return vp[0];
      } else if (ip2>=N) {
        return vp[N-1];
      } else {
        return vp[ip2-1] + fractpart*(vp[ip2]-vp[ip2-1]);
      }
    }
}

const tuint flx_find_pos(const tdouble*const vp, const tuint N, const tdouble e)
{
  if (N==0) return 0;
  if (e>=vp[N-1]) return N;
  tuint start = 0;
  tuint length = N;
  while (length>1) {
    const tuint check = start + ((length+1)/2)-1;
    if (e<vp[check]) {
      // start = start;
      length = check+1-start;
    } else {
      length = start+length-check-1;
      start = check+1;
    }
  };
  return start; 
}

const tuint flx_find_pos2(const tdouble*const vp, const tuint N, const tdouble e)
{
  if (N==0) return 0;
  if (e>vp[N-1]) return N;
  tuint start = 0;
  tuint length = N;
  while (length>1) {
    const tuint check = start + (length/2);
    if (e<vp[check]) {
      // start = start;
      length = check-start;
    } else {
      length = (start+length)-check;
      start = check;
    }
  };
  return start; 
}

const tdouble calc_distance(const tdouble*const v1, const tdouble*const v2, const tuint N)
{
  pdouble sum = 0;
  for (tuint i=0;i<N;++i) {
    sum += pow2(v1[i]-v2[i]);
  }
  return sqrt(sum.cast2double());
}





void flxVec_simple_plot(std::ostream& os, const flxVec& V, const bool checkTOL, const int prec, const int fixW, const bool brackets)
{
  if (brackets) os << "( ";
  for (size_t i = 0; i < V.get_N(); i++) {
    os << GlobalVar.Double2String(V[i],checkTOL,prec,fixW) << " ";
  }
  if (brackets) os << " )";
}

void flxVec_totalPrec_plot(std::ostream& os, const flxVec& V)
{
  for (size_t i=0; i<V.get_N(); ++i) {
    if (i!=0) os << ", ";
    os << GlobalVar.D2S_totalPrec(V[i]);
  }
}

std::ostream& operator<<(std::ostream& os, const flxVec& V)
{
  os << "( ";
  for (size_t i = 0; i < V.get_N(); i++) {
    if (i!=0) os << ", ";
    os << GlobalVar.Double2String( V[i] );
  }
  os << " )";
  return os;
}

std::ostream& operator<<(std::ostream& os, const bVec& V)
{
  os << "(";
  for (size_t i = 0; i < V.size(); i++) {
    if (i!=0) os << ", ";
    os << ((V[i]==true)?("1"):("0"));
  }
  os << " )";
  return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector< tuint >& V)
{
  os << "(";
  for (size_t i = 0; i < V.size(); i++) {
    if (i!=0) os << ",";
    os << V[i];
  }
  os << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const pdouble& p)
{
  os << "(" << GlobalVar.Double2String(p.get_v()) << "+" << GlobalVar.Double2String(p.get_c()) << ")";
  return os;
}
