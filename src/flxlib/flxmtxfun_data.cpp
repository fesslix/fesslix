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

#include "flxmtxfun_data.h"
#include "flxmtxfun.h"
#include "flxmtxfun_fun_calc.h"

using namespace std;

#if FLX_DEBUG
  int FlxConstMtxBox::Cinst = 0;
#endif

void FlxSMtx::check_1(const FlxSMtx& rhs) const
{
  if (rhs.get_ncols() != get_ncols() || rhs.get_nrows() != get_nrows() ) {
    std::ostringstream ssV;
    ssV << "Rows and columns of both matrices do not match.";
    throw FlxException("FlxSMtx::check_1", ssV.str() );
  }
}

void FlxSMtx::check_2(const FlxSMtx& rhs) const
{
  if ((rhs.mtx&&mtx==NULL)||(rhs.mtx==NULL&&mtx)) {
    std::ostringstream ssV;
    ssV << "Types of both matrices do not match.";
    throw FlxException("FlxSMtx::check_2", ssV.str() );
  }
}

FlxSMtx::FlxSMtx(FlxSMtx& rhs)
: nrows(rhs.nrows), ncols(rhs.ncols), dptr(rhs.dptr), mtx((rhs.mtx)?(rhs.mtx->copy()):NULL)
{
  
}

FlxSMtx::FlxSMtx(const tuint nrows, const tuint ncols, const tdouble val)
: nrows(nrows), ncols(ncols), dptr(nrows*ncols), mtx(NULL)
{
  dptr = val;
}

FlxSMtx::FlxSMtx(tdouble* dp, const tuint& nrows, const tuint& ncols)
: nrows(nrows), ncols(ncols), dptr(dp,nrows*ncols,false,true), mtx(NULL)
{

}

FlxSMtx::FlxSMtx(const flxVec& rhs)
: nrows(rhs.get_N()), ncols(1), dptr(rhs), mtx(NULL)
{
  
}

FlxSMtx::~FlxSMtx()
{
  if (mtx) delete mtx;
}

FlxSMtx& FlxSMtx::operator=(const FlxSMtx& rhs)
{
  check_1(rhs);
  check_2(rhs);
  if (mtx) {
    *mtx = *(rhs.mtx);
  } else {
    dptr = rhs.dptr;
  }
  return *this;
}

FlxSMtx& FlxSMtx::operator=(const tdouble& rhs)
{
  if (mtx) {
    throw FlxException_NotImplemented("FlxSMtx::operator=");
  } else {
    dptr = rhs;
  }
  return *this;
}

FlxSMtx& FlxSMtx::operator+=(const FlxSMtx& rhs)
{
  check_1(rhs);
  check_2(rhs);
  if (mtx) {
    throw FlxException_NotImplemented("FlxSMtx::operator+=");
  } else {
    dptr += rhs.dptr;
  }
  return *this;
}

FlxSMtx& FlxSMtx::operator*=(const tdouble& rhs)
{
  if (mtx) {
    *mtx *= rhs;
  } else {
    dptr *= rhs;
  }
  return *this;
}

const tdouble FlxSMtx::max() const
{
  if (mtx) {
    return mtx->max();
  } else {
    return dptr.get_max();
  }
}

const tdouble FlxSMtx::min() const
{
  if (mtx) {
    return mtx->min();
  } else {
    return dptr.get_min();
  }
}

const size_t FlxSMtx::maxID() const
{
  if (mtx) {
    return mtx->maxID();
  } else {
    return dptr.get_maxID();
  }
}

const size_t FlxSMtx::minID() const
{
  if (mtx) {
    return mtx->minID();
  } else {
    return dptr.get_minID();
  }
}

tdouble* FlxSMtx::get_internalPtr(const bool throwErr)
{
  if (mtx) {
    return mtx->get_internalPtr(throwErr);
  } else {
    return dptr.get_tmp_vptr();
  }
}

void FlxSMtx::insert(const tuint& i, const tuint& j, const tdouble& d)
{
  if (mtx) {
    mtx->add_value(i, j, -mtx->operator()(i,j)); 
    mtx->add_value(i, j, d);
  } else {
    dptr[i*ncols+j] = d;
  }
}

const tdouble FlxSMtx::operator()(const tuint& i, const tuint& j) const
{
  if (mtx) {
    return mtx->operator()(i,j);
  } else {
    return dptr[i*ncols+j];
  }
}

const tdouble FlxSMtx::operator()(const tuint& i) const
{
  if (mtx) {
    const tuint i_ = i/ncols;
    const tuint j_ = i%ncols;
    return mtx->operator()(i_,j_);
  } else {
    return dptr[i];
  }
}

void FlxSMtx::mult(const FlxSMtx& m1, const FlxSMtx& m2)
{
  // error checks
    if (m1.get_ncols()!=m2.get_nrows()) {
      std::ostringstream ssV;
      ssV << "Matrices can not be multiplied.";
      throw FlxException("FlxSMtx::mult_1", ssV.str() );
    }
    if (mtx) throw FlxException_Crude("FlxSMtx::mult_2");
    if (nrows!=m1.get_nrows() || ncols!=m2.get_ncols()) {
      throw FlxException_Crude("FlxSMtx::mult_3");
    }
  // 'create' matrices
    FlxMtx m_res(nrows,ncols,dptr.get_tmp_vptr());
    FlxMtx_base* mtx1 = m1.mtx;
    FlxMtx_base* mtx2 = m2.mtx;
    if (mtx1==NULL) {
      mtx1 = new FlxMtx(m1.get_nrows(),m1.get_ncols(),m1.dptr.get_tmp_vptr_const());
    }
    if (mtx2==NULL) {
      mtx2 = new FlxMtx(m2.get_nrows(),m2.get_ncols(),m2.dptr.get_tmp_vptr_const());
    }
  // perform the multiplication
    mtx1->MultMtx(*mtx2,m_res);
  // deallocate
    if (m1.mtx==NULL) delete mtx1;
    if (m2.mtx==NULL) delete mtx2;
}

ostream& operator<<(ostream& os, const FlxSMtx& M)
{
  tuint n = M.get_nrows();
  tuint m = M.get_ncols();
  for (tuint i = 0; i < n; ++i) {
    for (tuint j = 0; j < m; ++j) {
      os << " " << GlobalVar.Double2String( M(i,j) );
      if (j+1<m) os << ',';
    }
    if ( i+1 < n ) os << ';' << std::endl;
  }
  return os;
}

void SMtxBase_write_fullPrec(ostream& os, const FlxSMtx& M)
{
  os << "{";
  tuint n = M.get_nrows();
  tuint m = M.get_ncols();
  for (tuint i = 0; i < n; ++i) {
    for (tuint j = 0; j < m; ++j) {
      os << " " << GlobalVar.D2S_totalPrec( M(i,j) );
      if (j+1<m) os << ',';
    }
    if ( i+1 < n ) os << ';' << std::endl;
  }
  os << " }";
}


FlxConstMtxBox::FlxConstMtxBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxConstMtxBox' created ...";
      throw FlxException("FlxConstMtxBox::FlxConstMtxBox", ssV.str() );
    }
  #endif
  FunBaseFun_MtxConst::set_dataMtxConsts(this);
  FlxMtxBoxBase::set_ConstMtxBox(this);
}

FlxConstMtxBox::~FlxConstMtxBox()
{
  for (map<string, FlxSMtx*>::iterator pos = box.begin(); pos != box.end(); ++pos) {
    delete pos->second;
  }
}

void FlxConstMtxBox::insert(const string& name, FlxSMtx* value)
{
  pair<std::string, FlxSMtx*>Element(name, value);
  if ( ! box.insert(Element).second ) {  // ... and if the entry already exists ...
    // find it 
      map<std::string, FlxSMtx*>::iterator pos;
      pos = box.find(name);
      if ( !(pos != box.end()) ) {
        return; // that's not okay ... but should never happen
      }
    // delete the current value
      delete pos->second;
    // insert the new entry
      pos->second = value;    
  }
}

FlxSMtx* FlxConstMtxBox::get(const string& name, const bool err_if_unknown)
{
  map<std::string, FlxSMtx*>::iterator pos;
  pos = box.find(name);
  if ( pos != box.end() ) {
    return pos->second;
  } else {
    if (err_if_unknown) {
      std::ostringstream ssV;
      ssV << "Matrix-constant '" << name << "' has not yet been defined.";
      throw FlxException("FlxConstMtxBox::get_1", ssV.str() );
    } else {
      return NULL;
    }
  }
}

FlxSMtx* FlxConstMtxBox::get(const string& name, const tuint Nrows, const tuint Ncols, const bool forceSize)
{
  FlxSMtx* res = get(name,forceSize);
  if (res==NULL) {
    res = new FlxSMtx(Nrows,Ncols,ZERO);
    insert(name,res);
  } else if (res->get_nrows()!=Nrows || res->get_ncols()!=Ncols) {
    if (forceSize) {
      std::ostringstream ssV;
      ssV << "Dimension is " << res->get_nrows() << "x" << res->get_ncols() << " and not " << Nrows << "x" << Ncols;
      throw FlxException("FlxConstMtxBox::get_2", "Matrix-constant '" + name + "' has wrong dimension.", ssV.str() );
    } else {
      res = new FlxSMtx(Nrows,Ncols,ZERO);
      insert(name,res);
    }
  }
  return res;
}

tdouble* FlxConstMtxBox::get_Vec(const tuint N, const string& name, const bool forceSize)
{
  tuint NV = N;
  return get_Vec(name,NV,forceSize);
}

tdouble* FlxConstMtxBox::get_Vec(const string& name, tuint& N, const bool forceSize)
{
  #if FLX_DEBUG
    if (N==0 && forceSize) throw FlxException_Crude("FlxConstMtxBox::get_Vec_0");
  #endif
  if (N==0 || forceSize) {
    FlxSMtx* mp = get(name,true);
    const tuint Nr = mp->get_nrows();
    const tuint Nc = mp->get_ncols();
    if (Nr==1 || Nc==1) {        // qualifies as vector
      if (forceSize) {
        if ( (Nr==1 && Nc!=N) || (Nc==1 && Nr!=N) ) {
          std::ostringstream ssV;
          ssV << "Matrix-constant '" << name << "' has a size (" << Nr << "x" << Nc << ") "
          << "different from the one requested (vector of size " << N << ").";
          throw FlxException("FlxConstMtxBox::get_Vec_1",ssV.str());
        }
      } else {
        if (Nr==1) N = Nc;
        else N = Nr;
      }
      return mp->get_internalPtr(true);
    }
    // matrix!
      std::ostringstream ssV;
      ssV << "Matrix-constant '" << name << "' contains a matrix and not a vector.";
      throw FlxException("FlxConstMtxBox::get_Vec_2",ssV.str());
  }
  // the other case
  {
    FlxSMtx* mp = get(name,false);
    try {
      if (mp) {
        const tuint Nr = mp->get_nrows();
        const tuint Nc = mp->get_ncols();
        if (Nr==1 || Nc==1) {        // qualifies as vector
          if ( (Nr==1 && Nc==N) || (Nc==1 && Nr==N) ) {
            return mp->get_internalPtr(true);
          }
        }
      }
    } catch (FlxException &e) {
      // do nothing
    }
    // mp is not valid
    mp = new FlxSMtx(N,1,ZERO);
    insert(name,mp);
    return mp->get_internalPtr(true);
  }
}

const iVec FlxConstMtxBox::get_Vec_cast2tuint(const string& name, const bool No0)
{
  tuint N = 0;
  tdouble *dp = get_Vec(name,N);
  flxVec t(dp,N);
  iVec r(N);
  t.cast2tuint(r,No0);
  return r;
}

tdouble* FlxConstMtxBox::get_Mtx(const tuint Nr, const tuint Nc, const string& name, const bool forceSize)
{
  tuint NrV = Nr, NcV = Nc;
  return get_Mtx(name,NrV,NcV,forceSize);
}

tdouble* FlxConstMtxBox::get_Mtx(const string& name, tuint& Nr, tuint& Nc, const bool forceSize )
{
  #if FLX_DEBUG
    if (Nr*Nc==0 && Nr+Nc>0) throw FlxException_Crude("FlxConstMtxBox::get_Mtx_1");
    if (Nr+Nc==0 && forceSize) throw FlxException_Crude("FlxConstMtxBox::get_Mtx_2");
  #endif
  // in case the object must exist AND we need its size
  if (Nr+Nc==0) {
    FlxSMtx* mp = get(name,true);
    Nr = mp->get_nrows();
    Nc = mp->get_ncols();
    return mp->get_internalPtr(true);
  }
  if (forceSize) {
    FlxSMtx* mp = get(name,true);
    if (Nr==1 || Nc==1) {                // a vector is requested
      if (mp->get_nrows()==1 || mp->get_ncols()==1) {
        if (Nr*Nc == mp->get_nrows()*mp->get_ncols()) {
          return mp->get_internalPtr(true);
        }
      }
    } else {                                // a matrix is requested
      if (Nr==mp->get_nrows() && Nc==mp->get_ncols()) {
        return mp->get_internalPtr(true);
      }
    }
    std::ostringstream ssV;
    ssV << "Matrix-constant '" << name << "' has a size (" << mp->get_nrows() << "x" << mp->get_ncols() << ") "
    << "different from the one requested (" << Nr << "x" << Nc << ").";
    throw FlxException("FlxConstMtxBox::get_Mtx_3",ssV.str());
  }
  // the third case
  {
    FlxSMtx* mp = get(name,false);
    try {
      if (mp) {
        if (Nr==mp->get_nrows() && Nc==mp->get_ncols()) {
          return mp->get_internalPtr(true);
        }
      }
    } catch (FlxException &e) {
      // do nothing
    }
    // mp is not valid
    mp = new FlxSMtx(Nr,Nc,ZERO);
    insert(name,mp);
    return mp->get_internalPtr(true);
  }  
}

void FlxConstMtxBox::declareC(const string& name)
{
  if ( get(name) == NULL ) {
    insert(name, new FlxSMtx(1,1,ZERO));
  }
}

void FlxConstMtxBox::freeC(const string& name)
{
  delete get(name,true);
  box.erase(name);
}
