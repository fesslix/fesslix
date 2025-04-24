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

#include <cmath>

#include "flxrbrv.h"
#include "flxrbrv_rvs.h"
#include "flxdata.h"
#include "flxMtx_Eigen.h"
#include "flxobjrbrv.h"

tuint RBRV_set_base::static_ID = 0;


RBRV_entry::RBRV_entry(const std::string& name)
: 
  #if FLX_DEBUG
    valid(false),
  #endif
  value(ZERO),parent(NULL),name(name)
{
   
}

void RBRV_entry::eval_para()
{

}

const tdouble RBRV_entry::transform_x2y(const tdouble& x_val)
{
  std::ostringstream ssV;
  ssV << "This transformation is not available for this type of random variable (rv-name: " << name << ").";
  throw FlxException("RBRV_entry::transform_x2y", ssV.str() );
}

const tdouble RBRV_entry::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  std::ostringstream ssV;
  ssV << "This operation is not available for this type of random variable (rv-name: " << name << ").";
  throw FlxException("RBRV_entry::calc_pdf_x", ssV.str() );
}

const tdouble RBRV_entry::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  return log(calc_pdf_x(x_val,safeCalc));
}

const tdouble RBRV_entry::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  std::ostringstream ssV;
  ssV << "This operation is not available for this type of random variable (rv-name: " << name << ").";
  throw FlxException("RBRV_entry::calc_cdf_x", ssV.str() );
}

const tdouble RBRV_entry::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  return ONE-this->calc_cdf_x(x_val,safeCalc);
}

const tdouble RBRV_entry::calc_entropy()
{
  std::ostringstream ssV;
  ssV << "This operation is not available for this type of random variable (rv-name: " << name << ").";
  throw FlxException("RBRV_entry::calc_entropy", ssV.str() );
}

const tdouble RBRV_entry::get_median_current_config()
{
  std::ostringstream ssV;
  ssV << "This operation is not available for this type of random variable (rv-name: " << name << ").";
  throw FlxException("RBRV_entry::get_median_current_config", ssV.str() );
}

const tdouble RBRV_entry::get_mode_current_config()
{
  std::ostringstream ssV;
  ssV << "This operation is not available for this type of random variable (rv-name: " << name << ").";
  throw FlxException("RBRV_entry::get_mode_current_config", ssV.str() );
}

void RBRV_entry::set_x(const tdouble& xV)
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_entry::set_x");
    }
  #endif
  value = xV;
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_entry::set_parent(RBRV_set* parentV)
{
  if (parent) {
    throw FlxException_Crude("RBRV_entry::set_parent");
  }
  parent = parentV;
}

#if FLX_DEBUG
const tdouble& RBRV_entry::get_value() const
{
  if (valid) {
    return value;
  }
  {
    std::ostringstream ssV;
    ssV << "The RBRV '" << name << "' does not have a valid realization.";
    throw FlxException("RBRV_entry::get_value", ssV.str() );
  }
}
#endif

void RBRV_entry::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "* " << name << " [" << get_type() << "]";
  if (printID) {
    sout << " (" << counter++ << ")";
  }
  sout << std::endl;
}

py::dict RBRV_entry::info()
{
  py::dict res;
  res["type"] = get_type();
  res["name"] = name;
  return res;
}

RBRV_entry_fun::RBRV_entry_fun(const std::string& name, FlxFunction* fun)
 : RBRV_entry(name), fun(fun), fun_val(ZERO)
{
    eval_para();
}

void RBRV_entry_fun::eval_para()
{
  fun_val = fun->calc();
}

void RBRV_entry_fun::transform_y2x(const tdouble*const y_vec)
{
  #if FLX_DEBUG
    if (valid) throw FlxException_Crude("RBRV_entry_fun::transform_y2x");
  #endif
  value = fun_val;
  #if FLX_DEBUG
    valid = true;
  #endif
}

const bool RBRV_entry_fun::check_x(const tdouble xV)
{
  return true;
}


RBRV_entry_RV_base::RBRV_entry_RV_base(const std::string& name, const tuint iID)
: RBRV_entry(name), iID(iID),corr_rv(NULL), corr_valF(NULL), corr_val(ZERO), throwErrors(true)
{

}

RBRV_entry_RV_base::~RBRV_entry_RV_base()
{
  if (corr_valF) delete corr_valF;
}

void RBRV_entry_RV_base::init()
{
  this->eval_para();
  value = this->transform_y2x(ZERO);
  #if FLX_DEBUG
    valid = true;                // true, if value has been set
  #endif
}

void RBRV_entry_RV_base::transform_y2x(const tdouble*const y_vec)
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_entry_RV_base::transform_y2x");
    }
  #endif
  if (corr_rv) {
    eval_corr();
    const tdouble y1 = y_vec[corr_rv->iID];
    const tdouble y2 = y_vec[iID];
    const tdouble cond_mu = corr_val*y1;
    const tdouble cond_sd = sqrt(ONE-pow2(corr_val));
    const tdouble y = cond_sd*y2+cond_mu;
    value = transform_y2x(y);
  } else {
    const tdouble& y = y_vec[iID];
    value = transform_y2x(y);
  }
  #if FLX_DEBUG
    valid = true;
  #endif
}

struct flx_RBRV_entry_RV_data {
  tdouble p_i;                // probability of HPD interval
  RBRV_entry_RV_base* rbrv;
  tdouble om;                // smallest observed value
  flx_RBRV_entry_RV_data() : p_i(ZERO), rbrv(NULL), om(std::numeric_limits<double>::infinity()) {}
};

tdouble perfFun1D_RBRV_entry_RV_base(const tdouble x, void *params) {
  flx_RBRV_entry_RV_data *dat = (flx_RBRV_entry_RV_data *)params;
  // calculate lower bound
    tdouble lp = x;
    if (lp<=ZERO) return 1e90;
    if (lp >= ONE - dat->p_i) return 1e90;
    const tdouble lb = dat->rbrv->transform_y2x( rv_InvPhi( lp ) );
  // calculate upper bound
    tdouble up = lp + dat->p_i;
    if (up>ONE) throw FlxException_Crude("perfFun1D_RBRV_entry_RV_base");
    const tdouble ub = dat->rbrv->transform_y2x( rv_InvPhi( up ) );
  // remember minimum?
    if (ub-lb < dat->om) dat->om = ub-lb;
  return ub-lb;
}

const tdouble RBRV_entry_RV_base::get_HPD(const tdouble p)
{
  // initial settings
    const tdouble bd = 10.;        // number of intervals
    flx_RBRV_entry_RV_data dat;
    dat.p_i = p;
    dat.rbrv = this;
    tdouble x_min = (ONE-p)/2;
    const tdouble lmin = 1e-30;
    const tdouble lmax = (ONE-p)-1e-9;
  // perform optimization
    flx_optim(x_min-x_min/bd,x_min+x_min/bd,x_min,&perfFun1D_RBRV_entry_RV_base,&dat,true,true,1000,1000,1e-6,1e-6);        // ,&(GlobalVar.slogcout(1))
  if (perfFun1D_RBRV_entry_RV_base(lmin,&dat)<=dat.om) return ZERO;
  if (perfFun1D_RBRV_entry_RV_base(lmax,&dat)<=dat.om) return ONE-p;
  return x_min;
}

void RBRV_entry_RV_base::set_corr(RBRV_entry_RV_base* corr_rv_, FlxFunction* corrVal_, bool corrFixed_, const bool throwErrorsV)
{
  if (corr_rv_==NULL) throw FlxException_Crude("RBRV_entry_RV_base::set_corr_1");
  corr_rv = corr_rv_;
  throwErrors = throwErrorsV;
  corr_valF = new FunRBRV_calc_R_for_rhoPrime( corr_rv,this,new FlxFunction(*corrVal_), true );
  if (corrFixed_) {
    eval_corr();
    delete corr_valF;
    corr_valF = NULL;
  }
}

void RBRV_entry_RV_base::eval_corr()
{
  if (corr_valF) {
    FunRBRV_calc_R_for_rhoPrime* ftmp_ = dynamic_cast<FunRBRV_calc_R_for_rhoPrime*>(corr_valF);
    #if FLX_DEBUG
      if (ftmp_==NULL) throw FlxException_Crude("RBRV_entry_RV_base::eval_corr");
    #endif
    corr_val = ftmp_->calc_(throwErrors);
  }
}



RBRV_set_base::RBRV_set_base(const bool internal, const tuint sRV, const std::string& name, const bool noID)
: ID(noID?0:static_ID++),internal(internal),sRV(sRV),y_of_set(sRV),name(name)
{

}

const flxVec& RBRV_set_base::propose_y()
{
  RndCreator->gen_smp(y_of_set);
  return y_of_set;
}

void RBRV_set_base::get_y(tdouble*const y_vec)
{
  if (sRV>0) {
    flxVec y(y_vec,sRV);
    y = y_of_set;
  }
}

void RBRV_set_base::set_y(const tdouble*const y_vec)
{
  if (sRV>0) {
    const flxVec y(y_vec,sRV);
    y_of_set = y;
  }
}

void RBRV_set_base::transform_x2y()
{
  std::ostringstream ssV;
  ssV << "This transformation is not available for the set '" << name << "'.";
  throw FlxException("RBRV_set_base::transform_x2y", ssV.str() );
}

void RBRV_set_base::transform_y2w(const tdouble*const y_vec, tdouble*const w_vec)
{
  std::ostringstream ssV;
  ssV << "This transformation is not available for the set '" << name << "'.";
  throw FlxException("RBRV_set_base::transform_y2w", ssV.str() );
}

const FlxMtxSparsLTri* RBRV_set_base::calc_Jinv_1()
{
  throw FlxException_NotImplemented("RBRV_set_base::calc_Jinv_1");
}

void RBRV_set_base::calc_Jinv_2(tdouble* dxdw)
{
  throw FlxException_NotImplemented("RBRV_set_base::calc_Jinv_2");
}

const tdouble RBRV_set_base::get_pdf_x_eval_log()
{
  std::ostringstream ssV;
  ssV << "This operation is not available for the set '" << name << "'.";
  throw FlxException("RBRV_set_base::get_pdf_x_eval_log", ssV.str() );
}

void RBRV_set_base::add_covMTX(FlxMtxSym& cm)
{
  std::ostringstream ssV;
  ssV << "This operation is not available for the set '" << name << "'.";
  throw FlxException("add_covMTX::add_covMTX", ssV.str() );
}

const std::string RBRV_set_base::get_rv_name(const tuint index)
{
  if (index>=get_NOX()) throw FlxException_Crude("RBRV_set_base::get_rv_name");
  std::ostringstream ssV;
  ssV << name << "::" << index;
  return ssV.str();
}

const tdouble RBRV_set_base::get_pdf_y_eval_log() const
{
  tdouble s = ZERO;
  for (tuint i=0;i<sRV;++i) {
    s += rv_phi_log(y_of_set[i]);
  }
  return s;
}

RBRV_set_parents::RBRV_set_parents(const bool internal, const tuint sRV, const std::string& name, const tuint Nparents, RBRV_set_base**const parents, const bool noID)
: RBRV_set_base(internal,sRV,name,noID), Nparents(Nparents), parents(parents)
{

}

RBRV_set_parents::~RBRV_set_parents()
{
  if (parents) delete [] parents;
}

void RBRV_set_parents::find_dependent_sets(std::vector< RBRV_set_base* >& setvec)
{
  // make sure that it is not already in the list
    for (tuint i=0;i<setvec.size();++i) {
      if (setvec[i]==this) return;
    }
  // add the parents of this set first
    for (tuint i=0;i<Nparents;++i) {
      parents[i]->find_dependent_sets(setvec);
    }
  setvec.push_back(this);
}

void RBRV_set_parents::print_parents(std::ostream& sout)
{
  if (Nparents==0) return;
  sout << " (";
  for (tuint i=0;i<Nparents;++i) {
    if (i>0) sout << " ,";
    sout << parents[i]->name;
  }
  sout << ")";
}


RBRV_set::RBRV_set(const bool internal, const tuint sRV, const std::string& name, const bool noID, const tuint Nentries, RBRV_entry**const entries, 
                   const tuint Nparents, RBRV_set_base**const parents, const bool x2y_allowedV)
: RBRV_set_parents(internal, sRV, name,Nparents,parents,noID), Nentries(Nentries), entries(entries), x2y_allowed((Nparents==0)?(x2y_allowedV):false)
{
  for (tuint i=0;i<Nentries;++i) entries[i]->set_parent(this);
  if (x2y_allowed) {
    for (tuint i=0;i<Nentries;++i) {
      RBRV_entry_RV_base *rve = dynamic_cast<RBRV_entry_RV_base*>(entries[i]);
      if (rve==NULL) {
        x2y_allowed = false;
      } else {
        x2y_allowed = rve->allow_x2y();
      }
      if (!x2y_allowed) break;
    }
  }
}

RBRV_set::~RBRV_set()
{
  // delete entries
    for (tuint i=0;i<Nentries;++i) delete entries[i];
    delete [] entries;
}

void RBRV_set::set_is_valid(const bool is_valid)
{
  #if FLX_DEBUG
    for (tuint i=0;i<Nentries;++i) entries[i]->set_is_valid(is_valid);
  #endif
}

void RBRV_set::transform_y2x()
{
  const tdouble*const y_vec = y_of_set.get_tmp_vptr_const();
  for (tuint i=0;i<Nentries;++i) {
    entries[i]->eval_para();
    entries[i]->transform_y2x(y_vec);
  }
}

void RBRV_set::transform_x2y()
{
  if (!x2y_allowed) {
    std::ostringstream ssV;
    ssV << "The transformation from original-space to standard Normal space is not allowed for the set '" << name << "'.";
    throw FlxException("RBRV_set::transform_x2y_1", ssV.str() );
  }
  tdouble* y_vec = y_of_set.get_tmp_vptr();
  for (tuint i=0;i<Nentries;++i) {
    RBRV_entry_RV_base *rve = dynamic_cast<RBRV_entry_RV_base*>(entries[i]);
    if (rve) {
      y_vec[i] = rve->transform_x2y(rve->get_value());
    } else {
      throw FlxException_Crude("RBRV_set::transform_x2y_2");
    }
  }
}

void RBRV_set::transform_y2w(const tdouble*const y_vec, tdouble*const w_vec)
{
  if (!x2y_allowed) {
    std::ostringstream ssV;
    ssV << "The transformation from original-space to standard Normal space is not allowed for the set '" << name << "'.";
    throw FlxException("RBRV_set::transform_y2w_01", ssV.str() );
  }
  for (tuint i=0;i<Nentries;++i) {
    RBRV_entry_RV_base *rve = dynamic_cast<RBRV_entry_RV_base*>(entries[i]);
    if (rve==NULL) {
      throw FlxException("RBRV_set::transform_y2w_02","RBRV-type is not allowed for this operation.");
    } else {
      if (rve->allow_x2y()==false) {
        throw FlxException("RBRV_set::transform_y2w_03","RBRV-type is not allowed for this operation.");
      }
    }
    w_vec[i] = y_vec[i];
  }
}



void RBRV_set::set_x(const tdouble*const x_vec)
{
  for (tuint i=0;i<Nentries;++i) {
    entries[i]->eval_para();
    entries[i]->set_x(x_vec[i]);
  }
}

void RBRV_set::get_x(tdouble*const x_vec)
{
  for (tuint i=0;i<Nentries;++i) x_vec[i] = entries[i]->get_value();
}

void RBRV_set::get_mean(tdouble*const m_vec)
{
  for (tuint i=0;i<Nentries;++i) m_vec[i] = entries[i]->get_mean_current_config();
}

void RBRV_set::get_sd(tdouble*const s_vec)
{
  for (tuint i=0;i<Nentries;++i) s_vec[i] = entries[i]->get_sd_current_config();
}

const bool RBRV_set::check_xVec(const tdouble* xp)
{
  for (tuint i=0;i<Nentries;++i) {
    if ( entries[i]->check_x(xp[i])==false ) return false;
  }
  return true;
}

const FlxMtxSparsLTri* RBRV_set::calc_Jinv_1()
{
  return NULL;
}

void RBRV_set::calc_Jinv_2(tdouble* dxdw)
{
  // compute w-vector
    flxVec wv(dxdw,Nentries);
    wv = y_of_set;
  // calc dx/dw (diagonal matrix!)
    for (tuint i=0;i<Nentries;++i) {
      dxdw[i] = rv_phi(dxdw[i])/(entries[i]->calc_pdf_x(entries[i]->get_value()));
    }
}

void RBRV_set::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  for (tuint i=0;i<Nentries;++i) {
    entries[i]->eval_para();
    entries[i]->print(sout,prelim+"  ",counter,printID);
  }
}

const tdouble RBRV_set::get_pdf_x_eval_log()
{
  // transform x2y -> we will need it later
    transform_x2y();
  // calculate the determinant of the Jacobi matrix
    tdouble detJ = ZERO;        // log-transform!
    for (tuint i=0;i<Nentries;++i) {
      const tdouble xi = entries[i]->get_value();
      if (std::isinf(xi)) continue;
      detJ += entries[i]->calc_pdf_x_log(xi);
    }
  // transform y2x -> so we can check the error
    set_is_valid(false);
    transform_y2x();
  return detJ;
}

const std::string RBRV_set::get_rv_name(const tuint index)
{
  if (index>=get_NOX()) throw FlxException_Crude("RBRV_set_base::get_rv_name");
  return entries[index]->name;
}


RBRV_set_Nataf::RBRV_set_Nataf(const bool internal, const tuint sRV, const std::string& name, const bool noID, const tuint Nentries, RBRV_entry**const entries, FlxMtxSparsLTri* L)
: RBRV_set_base(internal, sRV, name, noID), Nentries(Nentries), w_of_set(Nentries), entries(entries), L(L)
{
  // for (tuint i=0;i<Nentries;++i) entries[i]->set_parent(this);
  for (tuint i=0;i<Nentries;++i) {
    RBRV_entry_RV_base *rve = dynamic_cast<RBRV_entry_RV_base*>(entries[i]);
    if (rve==NULL) {
      throw FlxException("RBRV_set_Nataf::RBRV_set_Nataf_1","RBRV-type is not allowed in a Nataf-set.");
    }
  }
  if (Nentries!=sRV) throw FlxException_Crude("RBRV_set_Nataf::RBRV_set_Nataf_2");
}

RBRV_set_Nataf::~RBRV_set_Nataf()
{
  // delete entries
    for (tuint i=0;i<Nentries;++i) delete entries[i];
    delete [] entries;
  if (L) delete L;
}

void RBRV_set_Nataf::set_is_valid(const bool is_valid)
{
  #if FLX_DEBUG
    for (tuint i=0;i<Nentries;++i) entries[i]->set_is_valid(is_valid);
  #endif
}

void RBRV_set_Nataf::transform_y2x()
{
  if (L) {
    L->MultMv(y_of_set,w_of_set);
  } else {
    w_of_set = y_of_set;
  }
  const tdouble*const w_vec = w_of_set.get_tmp_vptr_const();
  for (tuint i=0;i<Nentries;++i) {
    entries[i]->eval_para();
    entries[i]->transform_y2x(w_vec);
  }
}

void RBRV_set_Nataf::transform_x2y()
{
  tdouble* w_vec = w_of_set.get_tmp_vptr();
  for (tuint i=0;i<Nentries;++i) {
    RBRV_entry_RV_base *rve = dynamic_cast<RBRV_entry_RV_base*>(entries[i]);
    if (rve) {
      w_vec[i] = rve->transform_x2y(rve->get_value());
    } else {
      throw FlxException_Crude("RBRV_set_Nataf::transform_x2y");
    }
  }
  if (L) {
    L->MultInv(w_of_set,y_of_set);
  } else {
    y_of_set = w_of_set;
  }
}

void RBRV_set_Nataf::transform_y2w(const tdouble*const y_vec, tdouble*const w_vec)
{
  const flxVec yv(y_vec,get_NRV());
  flxVec wv(w_vec,get_NRV());
  if (L) {
    L->MultMv(yv,wv);
  } else {
    wv = yv;
  }
}

void RBRV_set_Nataf::set_x(const tdouble*const x_vec)
{
  for (tuint i=0;i<Nentries;++i) {
    entries[i]->eval_para();
    entries[i]->set_x(x_vec[i]);
  }
}

void RBRV_set_Nataf::get_x(tdouble*const x_vec)
{
  for (tuint i=0;i<Nentries;++i) x_vec[i] = entries[i]->get_value();
}

const bool RBRV_set_Nataf::check_xVec(const tdouble* xp)
{
  for (tuint i=0;i<Nentries;++i) {
    if ( entries[i]->check_x(xp[i])==false ) return false;
  }
  return true;
}

void RBRV_set_Nataf::get_mean(tdouble*const m_vec)
{
  for (tuint i=0;i<Nentries;++i) m_vec[i] = entries[i]->get_mean_current_config();
}

void RBRV_set_Nataf::get_sd(tdouble*const s_vec)
{
  for (tuint i=0;i<Nentries;++i) s_vec[i] = entries[i]->get_sd_current_config();
}

const FlxMtxSparsLTri* RBRV_set_Nataf::calc_Jinv_1()
{
  return L;
}

void RBRV_set_Nataf::calc_Jinv_2(tdouble* dxdw)
{
  // compute w-vector
    flxVec wv(dxdw,Nentries);
    if (L) {
      L->MultMv(y_of_set,wv);
    } else {
      wv = y_of_set;
    }
  // calc dx/dw (diagonal matrix!)
    for (tuint i=0;i<Nentries;++i) {
      dxdw[i] = rv_phi(dxdw[i])/(entries[i]->calc_pdf_x(entries[i]->get_value()));
    }
}

void RBRV_set_Nataf::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << " (Nataf-set)";
  if (printID) {
    sout << " (RV-ID: [" << counter << ";" << counter+get_NOX_only_this() << "[ )";
  }
  sout << std::endl;
  counter += get_NOX_only_this();
}

const tdouble RBRV_set_Nataf::get_pdf_x_eval_log()
{
  // transform x2y -> we will need it later
    transform_x2y();
  // calculate the determinant of the Jacobi matrix
    tdouble detJ = (L?(L->det_log()):ZERO);        // log-transform!
    for (tuint i=0;i<Nentries;++i) {
      const tdouble xi = entries[i]->get_value();
      if (std::isinf(xi)) continue;
      detJ += rv_phi_log(entries[i]->transform_x2y(xi)) - entries[i]->calc_pdf_x_log(xi);
    }
  detJ = get_pdf_y_eval_log() - detJ;
  // transform y2x -> so we can check the error
    set_is_valid(false);
    transform_y2x();
  return detJ;
}

const std::string RBRV_set_Nataf::get_rv_name(const tuint index)
{
  if (index>=get_NOX()) throw FlxException_Crude("RBRV_set_base::get_rv_name");
  return entries[index]->name;
}

void RBRV_set_Nataf::find_dependent_sets(std::vector< RBRV_set_base* >& setvec)
{
  // make sure that it is not already in the list
    for (tuint i=0;i<setvec.size();++i) {
      if (setvec[i]==this) return;
    }
  setvec.push_back(this);
}


RBRV_set_noise::RBRV_set_noise(const bool internal, const tuint sRV, const std::string& name, const bool noID, RBRV_entry*const transf, const tuint Nparents, RBRV_set_base**const parents)
 : RBRV_set_parents(internal,sRV,name,Nparents,parents,noID), x_of_set(sRV), transf(transf), is_stdN(dynamic_cast<RBRV_entry_RV_stdN*>(transf))
  #if FLX_DEBUG
    , valid(false)
  #endif
{
  
}

void RBRV_set_noise::transform_y2x()
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_noise::transform_y2x");
    }
  #endif
  transf->eval_para();
  for (tuint i=0;i<sRV;++i) {
    #if FLX_DEBUG
      transf->set_is_valid(false);
    #endif
    transf->transform_y2x(y_of_set.get_tmp_vptr_const()+i);
    x_of_set[i] = transf->get_value();
  }
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_noise::transform_x2y()
{
  for (tuint i=0;i<sRV;++i) {
    y_of_set[i] = transf->transform_x2y(x_of_set[i]);
  }
}

void RBRV_set_noise::transform_y2w(const tdouble*const y_vec, tdouble*const w_vec)
{
  const flxVec yv(y_vec,get_NRV());
  flxVec wv(w_vec,get_NRV());
  wv = yv;
}

const bool RBRV_set_noise::allow_x2y() const
{
  return (Nparents==0);
}

const tdouble RBRV_set_noise::get_pdf_x_eval_log()
{
  if (is_stdN) {
    return -(sRV*log(2*PI)+y_of_set.get_Norm2_NOroot())/2;
  } else {
    tdouble s = ZERO;
    for (tuint i=0;i<sRV;++i) {
      s += transf->calc_pdf_x_log(x_of_set[i]);
    }
    return s;
  }
}

void RBRV_set_noise::set_x(const tdouble*const x_vec)
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_noise::set_x");
    }
  #endif
  transf->eval_para();
  flxVec tv(x_vec,sRV);
  x_of_set = tv;
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_noise::get_x(tdouble*const x_vec)
{
  #if FLX_DEBUG
  if (valid) {
  #endif
  flxVec tv(x_vec,sRV);
  tv = x_of_set;
  #if FLX_DEBUG
  } else {
    std::ostringstream ssV;
    ssV << "The set '" << name << "' does not have valid realizations.";
    throw FlxException("RBRV_set_noise::get_x", ssV.str() );
  }
  #endif
}

const bool RBRV_set_noise::check_xVec(const tdouble* xp)
{
  for (tuint i=0;i<sRV;++i) {
    if ( transf->check_x(xp[i]) == false ) return false;
  }
  return true;
}

void RBRV_set_noise::get_mean(tdouble*const m_vec)
{
  flxVec m(m_vec,sRV);
  m = transf->get_mean_current_config();
}

void RBRV_set_noise::get_sd(tdouble*const s_vec)
{
  flxVec s(s_vec,sRV);
  s = transf->get_sd_current_config();
}

const FlxMtxSparsLTri* RBRV_set_noise::calc_Jinv_1()
{
  return NULL;
}

void RBRV_set_noise::calc_Jinv_2(tdouble* dxdw)
{
  // compute w-vector
    flxVec wv(dxdw,sRV);
    wv = y_of_set;
  // calc dx/dw (diagonal matrix!)
    for (tuint i=0;i<sRV;++i) {
      dxdw[i] = rv_phi(dxdw[i])/(transf->calc_pdf_x(x_of_set[i]));
    }
}

void RBRV_set_noise::add_covMTX(FlxMtxSym& cm)
{
  const tdouble sd2 = pow2(transf->get_sd_current_config());
  for (tuint i=0;i<sRV;++i) {
    cm(i,i) += sd2;
  }
}

void RBRV_set_noise::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << "uncorrelated noise";
  if (printID) {
    sout << "  ( RV-ID: [" << counter << ";" << counter+get_NOX_only_this() << "[ )";
  }
  sout << std::endl;  
  transf->eval_para();
  transf->print(sout,prelim+"  ",counter,false);
  counter += get_NOX_only_this();
}


RBRV_set_proc::RBRV_set_proc(const bool internal, const tuint Nv, const tuint Mv, const std::string& name, const bool noID, RBRV_entry*const transf, FlxFunction* rho, const tdouble dx, const tuint Nparents, RBRV_set_base**const parents, const int ev_solver, const bool only_once, const bool rho_Gauss)
 : RBRV_set_parents(internal,((Mv==0)?Nv:Mv),name,Nparents,parents,noID), x_of_set(Nv), transf(transf), rho(rho), dx(dx), N(Nv), M(Mv), ev_solver(ev_solver), rho_Gauss(rho_Gauss), eole_err(-ONE), Jdet(ZERO), only_once(only_once)
  #if FLX_DEBUG
    , valid(false)
  #endif
  , Eigenvalues(NULL), Lt(NULL),xhelp(NULL)
{
  
}

RBRV_set_proc::~RBRV_set_proc()
{
  delete transf;
  delete rho;
  if (Eigenvalues) delete Eigenvalues;
  if (Lt) delete Lt;
  if (xhelp) delete xhelp;
}

void RBRV_set_proc::assemble_rhoPrime(FlxMtxSym& rhoPrime)
{
  // get global constants
    tdouble &gx = *(FlxDataBase::get_data().ConstantBox.get("gx"));
    tdouble &gx2 = *(FlxDataBase::get_data().ConstantBox.get("gx2"));
    tdouble &deltax = *(FlxDataBase::get_data().ConstantBox.get("deltax"));
    const tdouble old_gx = gx;
    const tdouble old_gx2 = gx2;
    const tdouble old_deltax = deltax;
  RBRV_entry_RV_normal *pn = dynamic_cast<RBRV_entry_RV_normal*>(transf);
  if (pn || rho_Gauss) {        // Gaussian process
    for (tuint i=0;i<N;++i) {
      gx = dx*i;
      for (tuint j=0;j<i;++j) {
        gx2 = dx*j;
        deltax = gx - gx2;
        const tdouble rhot = rho->calc();
        if (rhot>=ONE||rhot<=-ONE) {
          std::ostringstream ssV;
          ssV << "Invalid correlation (" << GlobalVar.Double2String(rhot) << ") at position (" << i << "; " << j << "), "
            << "gx=" << GlobalVar.Double2String(gx) << ", gx2=" << GlobalVar.Double2String(gx2) << ", deltax=" << GlobalVar.Double2String(deltax) << ".";
          throw FlxException("RBRV_set_proc::get_x", ssV.str() );
        }
        rhoPrime(i,j) = rhot;
      }
      rhoPrime(i,i) = ONE;
    }
  } else {
    throw FlxException_NotImplemented("RBRV_set_proc::transform_y2x_3");
  }
  // restore global constants
    gx = old_gx;
    gx2 = old_gx2;
    deltax = old_deltax;
}

void RBRV_set_proc::assemble_system()
{
  if (M==0) {        // make sure that the Cholesky decomposition was performed
    if (Lt==NULL || only_once==false) {
      FlxMtxSym rhoPrime(N);
      assemble_rhoPrime(rhoPrime);
      // perform the Cholesky decomposition
        if (Lt==NULL) Lt = new FlxMtxLTri(N);
        try {
          Lt->CholeskyDec(rhoPrime);
        } catch (FlxException& e) {
          FLXMSG("RBRV_set_proc::assemble_system_1",1);
          delete Lt; Lt = NULL;
          throw;
        }
        Jdet = Lt->det_log();
    }
  } else {        // apply the EOLE method
    if (Eigenvalues==NULL || only_once==false) {
      // check if M is valid
        if (M>N) {
          std::ostringstream ssV;
          ssV << "M=" << M << " must be smaller than N=" << N << ".";
          throw FlxException("RRBRV_set_proc::assemble_system_2", ssV.str() );
        }
      // set up the correlation matrix
        FlxMtxSym rhoPrime(N);
        assemble_rhoPrime(rhoPrime);
      // compute the eigenvalues and eigenvectors
        if (Eigenvalues==NULL) {
          Eigenvalues = new flxVec(M);
          for (tuint i = 0; i < M; i++) {
            flxVec Vh(N);
            Eigenvectors.push_back(Vh);
          }
        }
        MtxEigenValue(rhoPrime,M,*Eigenvalues,Eigenvectors,ev_solver);
      // Normalize the Eigenvectors
        eole_err = ZERO;
        Jdet = ZERO;
        for (tuint i=0;i<M;i++) {
          flxVec& EV = Eigenvectors[i];
          EV/=EV.get_Norm2();
          EV.check_TOL();
          eole_err += (*Eigenvalues)[i];
          Jdet += log((*Eigenvalues)[i]);
        }
        Jdet /= 2;
        eole_err = ONE - eole_err/N;
      if( xhelp==NULL ) xhelp = new flxVec(N);
    }
  }
}

void RBRV_set_proc::transform_y2x()
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_proc::transform_y2x_1");
    }
  #endif
  assemble_system();
  if (M==0) {
    // transform from uncorrelated to correlated standard normal space
      Lt->MultMv(y_of_set,x_of_set);
  } else {
    // transform from uncorrelated to correlated standard normal space
      x_of_set.set_zero();
      for (tuint i=0;i<M;++i) {
        (*xhelp) = Eigenvectors[i];
        (*xhelp) *= (sqrt((*Eigenvalues)[i])*y_of_set[i]);
        x_of_set += *xhelp;
      }
  }
  // transform from correlated standard normal to original space
    transf->eval_para();
    for (tuint i=0;i<N;++i) {
      #if FLX_DEBUG
        transf->set_is_valid(false);
      #endif
      transf->transform_y2x(x_of_set.get_tmp_vptr_const()+i);
      x_of_set[i] = transf->get_value();
    }
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_proc::transform_x2y()
{
  assemble_system();
  if (M==0) {
    // transform from original space to correlated standard normal
      for (tuint i=0;i<N;++i) {
        y_of_set[i] = transf->transform_x2y(x_of_set[i]);
      }
    const tdouble th = y_of_set.get_sum();
    bool has_nan = (th!=th);
    Lt->MultInv(y_of_set,y_of_set);
    if (!has_nan) {
      tdouble *tht = y_of_set.get_tmp_vptr();
      for (tuint i=0;i<y_of_set.get_N();++i) {
        if (tht[i]!=tht[i]) tht[i] = log(ZERO);
      }
    }
  } else {
    for (tuint i=0;i<N;++i) {
      (*xhelp)[i] = transf->transform_x2y(x_of_set[i]);
    }
    for (tuint i=0;i<M;++i) {
      y_of_set[i] = (xhelp->operator*(Eigenvectors[i])) / sqrt((*Eigenvalues)[i]);
    }
  }
}

const bool RBRV_set_proc::allow_x2y() const
{
  return (Nparents==0);
}

void RBRV_set_proc::set_x(const tdouble*const x_vec)
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_proc::set_x");
    }
  #endif
  transf->eval_para();
  flxVec tv(x_vec,N);
  x_of_set = tv;
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_proc::get_x(tdouble*const x_vec)
{
  #if FLX_DEBUG
  if (valid) {
  #endif
  flxVec tv(x_vec,N);
  tv = x_of_set;
  #if FLX_DEBUG
  } else {
    std::ostringstream ssV;
    ssV << "The set '" << name << "' does not have valid realizations.";
    throw FlxException("RBRV_set_proc::get_x", ssV.str() );
  }
  #endif
}

const bool RBRV_set_proc::check_xVec(const tdouble* xp)
{
  for (tuint i=0;i<N;++i) {
    if ( transf->check_x(xp[i])==false ) return false;
  }
  return true;
}

void RBRV_set_proc::get_mean(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("RBRV_set_proc::get_mean");
}

void RBRV_set_proc::get_sd(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("RBRV_set_proc::get_sd");
}

void RBRV_set_proc::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  if (M>0) assemble_system();
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << "correlated variables with rho=" << rho->write() << "; dx=" << GlobalVar.Double2String(dx) << ";";
  if (M>0) {
    sout << " M=" << M << "; EOLE-err=" << GlobalVar.Double2String(eole_err);
  }
  sout << std::endl;
  if (printID) {
    sout << prelim << "  ( RV-ID: [" << counter << ";" << counter+get_NOX_only_this() << "[ )";
  }
  sout << std::endl;  
  transf->eval_para();
  transf->print(sout,prelim+"  ",counter,false);
  counter += get_NOX_only_this();
}

const tdouble RBRV_set_proc::get_pdf_x_eval_log()
{
  // transform x2y -> we will need it later
    transform_x2y();
  // calculate the determinant of the Jacobi matrix
    tdouble detJ = Jdet;        // log-transform!
    if (M>0 && M<N && only_once==false) {
      GlobalVar.alert.alert("RBRV_set_proc::get_pdf_x_eval_log","M<N and only_once=false should be avoided.");
    }
    for (tuint i=0;i<N;++i) {
      if (std::isinf(x_of_set[i])) continue;
      detJ += rv_phi_log(transf->transform_x2y(x_of_set[i])) - transf->calc_pdf_x_log(x_of_set[i]);
    }
  detJ = get_pdf_y_eval_log() - detJ;
  // transform y2x -> so we can check the error
    if (N>0&&M<N) {
      set_is_valid(false);
      transform_y2x();
    }
  return detJ;
}

void RBRV_set_proc::add_covMTX(FlxMtxSym& cm)
{
  // set up the correlation matrix
    FlxMtxSym rhoPrime(N);
    assemble_rhoPrime(rhoPrime);
  // get the standard deviations
    flxVec sdV(N);
    tdouble &gx = *(FlxDataBase::get_data().ConstantBox.get("gx"));
    for (tuint i=0;i<N;++i) {
      gx = dx*i;
      sdV[i] = transf->get_sd_current_config();
    }
  // perform the operation
    for (tuint i=0;i<N;++i) {
      for (tuint j=0;j<=i;++j) {
        rhoPrime(i,j) = rhoPrime(i,j)*sdV[i]*sdV[j];
      }
    }
    cm += rhoPrime;
}


RBRV_set_MVN::RBRV_set_MVN(const bool internal, const tuint Nv, const tuint Mv, const std::string& name, const bool noID, flxVec* mu, FlxMtxSym* CovM, const int ev_solver)
: RBRV_set_parents(internal,((Mv==0)?Nv:Mv),name,0,NULL,noID), x_of_set(Nv), N(Nv), M(Mv), mu(mu), CovM(CovM), ev_solver(ev_solver), eole_err(-ONE), Jdet(ZERO)
  #if FLX_DEBUG
    , valid(false)
  #endif
  , Eigenvalues(NULL), Lt(NULL),xhelp(NULL)
{
  try {
    assemble_system();
  } catch (FlxException& e) {
    FLXMSG("RBRV_set_MVN::RBRV_set_MVN",1);
    deallocate();
    throw;
  }
}

void RBRV_set_MVN::deallocate()
{
  delete mu;
  delete CovM;
  if (Eigenvalues) delete Eigenvalues;
  if (Lt) delete Lt;
  if (xhelp) delete xhelp;
}


void RBRV_set_MVN::update_EVP()
{
  assemble_system();
}

void RBRV_set_MVN::assemble_system()
{
  if (mu->get_N()!=N) throw FlxException_Crude("RBRV_set_MVN::assemble_system_1");
  if (CovM->nrows()!=N) throw FlxException_Crude("RBRV_set_MVN::assemble_system_2");
  if (M==0) {        // make sure that the Cholesky decomposition was performed
    // perform the Cholesky decomposition
      if (Lt==NULL) Lt = new FlxMtxLTri(N);
      try {
        Lt->CholeskyDec(*CovM);
      } catch (FlxException& e) {
        FLXMSG("RBRV_set_MVN::assemble_system_3",1);
        delete Lt; Lt = NULL;
        throw;
      }
      Jdet = Lt->det_log();
  } else {        // apply the EOLE method
    if (Eigenvalues==NULL) {
      // check if M is valid
        if (M!=N) {
          std::ostringstream ssV;
          ssV << "M=" << M << " must be equal to N=" << N << ".";
          throw FlxException("RBRV_set_MVN::assemble_system_5", ssV.str() );
        }
      // compute the eigenvalues and eigenvectors
        if (Eigenvalues==NULL) {
          Eigenvalues = new flxVec(M);
          for (tuint i = 0; i < M; i++) {
            flxVec Vh(N);
            Eigenvectors.push_back(Vh);
          }
        }
        MtxEigenValue(*CovM,M,*Eigenvalues,Eigenvectors,ev_solver);
      // Normalize the Eigenvectors
        eole_err = ZERO;
        Jdet = ZERO;
        for (tuint i=0;i<M;i++) {
          flxVec& EV = Eigenvectors[i];
          EV/=EV.get_Norm2();
          EV.check_TOL();
          eole_err += (*Eigenvalues)[i];
          Jdet += log((*Eigenvalues)[i]);
        }
        Jdet /= 2;
        eole_err = ONE - eole_err/N;
      if( xhelp==NULL ) xhelp = new flxVec(N);
    } else {
      throw FlxException_Crude("RBRV_set_MVN::assemble_system_6");
    }
  }
}

void RBRV_set_MVN::transform_y2x()
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_MVN::transform_y2x_1");
    }
  #endif
  if (M==0) {
    // transform from uncorrelated to correlated standard normal space
      Lt->MultMv(y_of_set,x_of_set);
  } else {
    // transform from uncorrelated to correlated standard normal space
      x_of_set.set_zero();
      for (tuint i=0;i<M;++i) {
        (*xhelp) = Eigenvectors[i];
        (*xhelp) *= (sqrt((*Eigenvalues)[i])*y_of_set[i]);
        x_of_set += *xhelp;
      }
  }
  x_of_set += (*mu);
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_MVN::transform_x2y()
{
  if (M==0) {
    y_of_set = x_of_set;
    y_of_set -= (*mu);
    Lt->MultInv(y_of_set,y_of_set);
  } else {
    (*xhelp) = x_of_set;
    (*xhelp) -= (*mu);
    for (tuint i=0;i<M;++i) {
      y_of_set[i] = (xhelp->operator*(Eigenvectors[i])) / sqrt((*Eigenvalues)[i]);
    }
  }
}

const bool RBRV_set_MVN::allow_x2y() const
{
  return (Nparents==0);
}

void RBRV_set_MVN::set_x(const tdouble*const x_vec)
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_MVN::set_x");
    }
  #endif
  flxVec tv(x_vec,N);
  x_of_set = tv;
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_MVN::get_x(tdouble*const x_vec)
{
  #if FLX_DEBUG
  if (valid) {
  #endif
  flxVec tv(x_vec,N);
  tv = x_of_set;
  #if FLX_DEBUG
  } else {
    std::ostringstream ssV;
    ssV << "The set '" << name << "' does not have valid realizations.";
    throw FlxException("RBRV_set_MVN::get_x", ssV.str() );
  }
  #endif
}

void RBRV_set_MVN::get_mean(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("RBRV_set_MVN::get_mean");
}

void RBRV_set_MVN::get_sd(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("RBRV_set_MVN::get_sd");
}

void RBRV_set_MVN::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << "multivariate normal distribution;";
  if (M>0) {
    sout << " M=" << M << "; EOLE-err=" << GlobalVar.Double2String(eole_err);
  }
  sout << std::endl;
  if (printID) {
    sout << prelim << "  ( RV-ID: [" << counter << ";" << counter+get_NOX_only_this() << "[ )";
  }
  sout << std::endl;  
  counter += get_NOX_only_this();
}

const tdouble RBRV_set_MVN::get_pdf_x_eval_log()
{
  // transform x2y -> we will need it later
    transform_x2y();
  // calculate the determinant of the Jacobi matrix
    tdouble detJ = Jdet;        // log-transform!
  detJ = get_pdf_y_eval_log() - detJ;
  // transform y2x -> so we can check the error
    if (N>0&&M<N) {
      set_is_valid(false);
      transform_y2x();
    }
  return detJ;
}

void RBRV_set_MVN::add_covMTX(FlxMtxSym& cm)
{
  cm += (*CovM);
}


RBRV_set_MVN_cond::RBRV_set_MVN_cond(const bool internal, const tuint NrndV, const tuint NobsvV, const std::string& name, const bool noID, flxVec* mu, FlxMtxSym* CovM, const flxVec& x_obsvV)
: RBRV_set_parents(internal,NrndV,name,0,NULL,noID), x_of_set(NrndV), x_obsv(x_obsvV), y_obsv(NobsvV), 
  Nrnd(NrndV), Nobsv(NobsvV), Ntotal(Nrnd+Nobsv), mu(mu), CovM(CovM), Jdet(ZERO)
  #if FLX_DEBUG
    , valid(false)
  #endif
  , Lt(NULL), xhelp(Ntotal), yhelp(Ntotal)
{
  try {
    assemble_system();
  } catch (FlxException& e) {
    FLXMSG("RBRV_set_MVN_cond::RBRV_set_MVN_cond",1);
    deallocate();
    throw;
  }
}

void RBRV_set_MVN_cond::assemble_system()
{
  if (mu->get_N()!=Ntotal) throw FlxException_Crude("RBRV_set_MVN_cond::assemble_system_1");
  if (CovM->nrows()!=Ntotal) throw FlxException_Crude("RBRV_set_MVN_cond::assemble_system_2");
  // perform the Cholesky decomposition
    if (Lt==NULL) Lt = new FlxMtxLTri(Ntotal);
    try {
      Lt->CholeskyDec(*CovM);
    } catch (FlxException& e) {
      FLXMSG("RBRV_set_MVN_cond::assemble_system_3",1);
      delete Lt; Lt = NULL;
      throw;
    }
    Jdet = Lt->det_log();
  comp_yobsv();
}

void RBRV_set_MVN_cond::deallocate()
{
  delete mu;
  delete CovM;
  if (Lt) delete Lt;
}

void RBRV_set_MVN_cond::transform_y2x()
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_MVN_cond::transform_y2x_1");
    }
  #endif
  // perpare the y-vector
    {
      tdouble* hvt = yhelp.get_tmp_vptr();
      flxVec y1(hvt,Nobsv);
      y1 = y_obsv;
      flxVec y2(hvt+Nobsv,Nrnd);
      y2 = y_of_set;
    }
  // transform from uncorrelated to correlated standard normal space
    Lt->MultMv(yhelp,xhelp);
    xhelp += (*mu);
  // prepare the x-vector
    {
      tdouble* hvt = xhelp.get_tmp_vptr();
      flxVec x2(hvt+Nobsv,Nrnd);
      x_of_set = x2;
    }
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_MVN_cond::transform_x2y()
{
  // perpare the x-vector
    {
      tdouble* hvt = yhelp.get_tmp_vptr();
      flxVec x1(hvt,Nobsv);
      x1 = x_obsv;
      flxVec x2(hvt+Nobsv,Nrnd);
      x2 = x_of_set;
    }
  // perform the transformation
    yhelp -= (*mu);
    Lt->MultInv(yhelp,yhelp);
  // prepare the y-vector
    {
      tdouble* hvt = yhelp.get_tmp_vptr();
      flxVec y2(hvt+Nobsv,Nrnd);
      y_of_set = y2;
    }
}

const bool RBRV_set_MVN_cond::allow_x2y() const
{
  return (Nparents==0);
}

void RBRV_set_MVN_cond::set_x(const tdouble*const x_vec)
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_MVN_cond::set_x");
    }
  #endif
  flxVec tv(x_vec,Nrnd);
  x_of_set = tv;
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_MVN_cond::get_x(tdouble*const x_vec)
{
  #if FLX_DEBUG
  if (valid) {
  #endif
  flxVec tv(x_vec,Nrnd);
  tv = x_of_set;
  #if FLX_DEBUG
  } else {
    std::ostringstream ssV;
    ssV << "The set '" << name << "' does not have valid realizations.";
    throw FlxException("RBRV_set_MVN_cond::get_x", ssV.str() );
  }
  #endif
}

void RBRV_set_MVN_cond::get_mean(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("RBRV_set_MVN_cond::get_mean");
}

void RBRV_set_MVN_cond::get_sd(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("RBRV_set_MVN_cond::get_sd");
}

void RBRV_set_MVN_cond::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << "multivariate normal distribution conditioned on " << Nobsv << " observations;" << std::endl;
  if (printID) {
    sout << prelim << "  ( RV-ID: [" << counter << ";" << counter+get_NOX_only_this() << "[ )";
  }
  sout << std::endl;  
  counter += get_NOX_only_this();
}

const tdouble RBRV_set_MVN_cond::get_pdf_x_eval_log()
{
  // transform x2y -> we will need it later
    transform_x2y();
  // calculate the determinant of the Jacobi matrix
    tdouble detJ = Jdet;        // log-transform!
  detJ = get_pdf_y_eval_log() - detJ;
  return detJ;
}

void RBRV_set_MVN_cond::comp_yobsv()
{
  // perpare the x-vector
    {
      yhelp.set_zero();
      tdouble* hvt = yhelp.get_tmp_vptr();
      flxVec x1(hvt,Nobsv);
      x1 = x_obsv;
      flxVec m1(mu->get_tmp_vptr(),Nobsv);
      x1 -= m1;
    }
  // perform the transformation
    Lt->MultInv(yhelp,yhelp);
  // prepare the y-vector
    {
      tdouble* hvt = yhelp.get_tmp_vptr();
      flxVec y2(hvt,Nobsv);
      y_obsv = y2;
    }
}

void RBRV_set_MVN_cond::set_x_obsv(const flxVec& x_obsvV)
{
  x_obsv = x_obsvV;
  comp_yobsv();
}

void RBRV_set_MVN_cond::update_EVP()
{
  assemble_system();
}


RBRV_set_psd::RBRV_set_psd(const bool internal, const std::string& name, const tuint Nv, FlxFunction* psd_fun, const tdouble lb, const tdouble ub, const tuint Nparents, RBRV_set_base**const parents, tdouble& wp)
: RBRV_set_parents(internal, 3*Nv, name, Nparents, parents, false), N(Nv), psd_fun(psd_fun), lb(lb), ub(ub), wp(wp)
  #if FLX_DEBUG
    , valid(false)
  #endif
{
  if (ub<=lb) {
    std::ostringstream ssV;
    ssV << "The upper bound (" << GlobalVar.Double2String(ub) << ") must be larger than the lower bound (" << GlobalVar.Double2String(lb) << ").";
    throw FlxException("RBRV_set_psd::RBRV_set_psd", ssV.str() );
  }
}

void RBRV_set_psd::get_mean(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("RBRV_set_psd::get_mean");
}

void RBRV_set_psd::get_sd(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("RBRV_set_psd::get_sd");
}

void RBRV_set_psd::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << "random process defined through its power spectral density function: " << psd_fun->write() << ";" << std::endl;
  sout << prelim << "  " << "number of discretization intervals: " << N << std::endl; 
  counter += get_NOX_only_this();
}

const tdouble RBRV_set_psd::eval_realization(const tdouble t)
{
  const tdouble prev = wp;
  pdouble res = ZERO;
  const tdouble* const yp = y_of_set.get_tmp_vptr_const();
  const tdouble dw = (ub-lb)/N;
  const tdouble dwh = dw/2;
  tuint c = 0;
  tdouble hd; 
  pdouble wi;
  for (tuint i=0;i<N;++i) {
    wi = lb; wi += dwh; wi += i*dw;
    wp = wi.cast2double();
    wi += (rv_Phi(yp[c++])*2-ONE)*dwh;
    hd = wi.cast2double()*t;
    hd = yp[c]*cos(hd) + yp[c+1]*sin(hd);
    c += 2;
    hd *= sqrt(psd_fun->cast2positive_or0()*2*dw);
    res += hd;
  }
  #if FLX_DEBUG
    if (c!=sRV) throw FlxException_Crude("RBRV_set_psd::eval_realization");
  #endif
  wp = prev;
  return res.cast2double();
}


RBRV_set_sphere::RBRV_set_sphere(const bool internal, const tuint sRV, const std::string& name, const bool noID, const tuint Nparents, RBRV_set_base**const parents, FlxFunction* r)
 : RBRV_set_parents(internal,sRV,name,Nparents,parents,noID), x_of_set(sRV), r(r)
  #if FLX_DEBUG
    , valid(false)
  #endif
{
  
}

void RBRV_set_sphere::set_x(const tdouble*const x_vec)
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_sphere::set_x");
    }
  #endif
  flxVec tv(x_vec,sRV);
  x_of_set = tv;
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_sphere::get_x(tdouble*const x_vec)
{
  #if FLX_DEBUG
  if (valid) {
  #endif
  flxVec tv(x_vec,sRV);
  tv = x_of_set;
  #if FLX_DEBUG
  } else {
    std::ostringstream ssV;
    ssV << "The set '" << name << "' does not have valid realizations.";
    throw FlxException("RBRV_set_sphere::get_x", ssV.str() );
  }
  #endif
}

const bool RBRV_set_sphere::check_xVec(const tdouble* xp)
{
  const tdouble rV = r->cast2positive();
  return (x_of_set.get_Norm2()<=rV);
}

void RBRV_set_sphere::transform_y2x()
{
  #if FLX_DEBUG
    if (valid) {
      throw FlxException_Crude("RBRV_set_sphere::transform_y2x");
    }
  #endif
  const tdouble s2 = y_of_set.get_Norm2_NOroot();
  const tdouble rV = r->cast2positive();
  const tdouble a = (rV * pow( flxgamma_rl( tdouble(sRV)/2, s2/2 ), 1/tdouble(sRV) ))/sqrt(s2);
  x_of_set = y_of_set;
  x_of_set *= a;  
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_set_sphere::transform_x2y()
{
  const tdouble r2 = x_of_set.get_Norm2_NOroot();
  const tdouble rV = r->cast2positive();
  const tdouble s2 = flxgamma_rl_inv( tdouble(sRV)/2, pow( sqrt(r2)/rV,tdouble(sRV) ) )*2;
  y_of_set = x_of_set;
  y_of_set *= sqrt(s2/r2);
}

void RBRV_set_sphere::get_mean(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("RBRV_set_sphere::get_mean");
}

void RBRV_set_sphere::get_sd(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("RBRV_set_sphere::get_sd");
}

const tdouble RBRV_set_sphere::get_pdf_x_eval_log()
{
  throw FlxException_NotImplemented("RBRV_set_sphere::get_pdf_x_eval_log");
}

void RBRV_set_sphere::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << "random sample distributed uniformly in " << sRV << "-dimensional hyper-sphere" << std::endl;
  counter += get_NOX_only_this();
}



RBRV_set_box::~RBRV_set_box()
{
  const std::size_t N = set_vec.size();
  for (std::size_t i=0;i<N;++i) delete set_vec[i];
}

void RBRV_set_box::register_set(RBRV_set_base* sta)
{
  const std::string& sn = sta->name;
  // check if 'sn' has already been defined
    if (get_set(sn,false)) {
      throw FlxException_Crude("RBRV_set_box::register_set_1");
    }
  // register in set_vec
    const tuint ID = sta->get_ID();
    if (ID!=set_vec.size()) throw FlxException_Crude("RBRV_set_box::register_set_2");
    set_vec.push_back(sta);
  // register in set_box
    std::pair<std::string,RBRV_set_base*>Element(sn,sta);
    if ( ! set_box.insert(Element).second ) {
      throw FlxException_Crude("RBRV_set_box::register_set_3");
    }
  // register all entries
    // has to be done before this function call
}

RBRV_set_base* RBRV_set_box::get_set(const std::string& name, const bool throwErr) const
{
  std::map<std::string,RBRV_set_base*>::const_iterator pos;
  pos = set_box.find(name);
  if ( pos != set_box.end() ) {
    return pos->second;
  } else {
    if (throwErr) {
      std::ostringstream ssV;
      ssV << "The set '" << name << "' does not exist.";
      throw FlxException("RBRV_set_box::get_set", ssV.str() );
    }
    return NULL;
  }
}

void RBRV_set_box::register_entry(RBRV_entry* eta)
{
  const std::string& sn = eta->name;
  // check if 'sn' has already been defined
    if (get_entry(sn,false)) throw FlxException_Crude("RBRV_set_box::register_entry_1");
  // register in entry_box
    std::pair<std::string,RBRV_entry*>Element(sn,eta);
    if ( ! entry_box.insert(Element).second ) {
      throw FlxException_Crude("RBRV_set_box::register_entry_2");
    }
}

RBRV_entry* RBRV_set_box::get_entry(const std::string& name, const bool throwErr) const
{
  std::map<std::string,RBRV_entry*>::const_iterator pos;
  pos = entry_box.find(name);
  if ( pos != entry_box.end() ) {
    return pos->second;
  } else {
    if (throwErr) {
      std::ostringstream ssV;
      ssV << "The entry '" << name << "' does not exist.";
      throw FlxException("RBRV_set_box::get_entry", ssV.str() );
    }
    return NULL;
  }
}

void RBRV_set_box::print_sets(std::ostream& sout, const std::string prelim)
{
  const tuint Nsets = set_vec.size();
  sout << prelim << "Total number of sets: " << Nsets << std::endl;
  for (tuint i=0;i<set_vec.size();++i) {
    sout << prelim << "- " << set_vec[i]->name;
    RBRV_set_parents* ps = dynamic_cast<RBRV_set_parents*>(set_vec[i]);
    if (ps) {
      ps->print_parents(sout);
    }
    sout << std::endl;
  }
}



RBRV_constructor::RBRV_constructor(const std::vector< RBRV_set_base* >& setvec)
: setvec(setvec), NRV(count_NRV(setvec)),NOX(count_NOX(setvec)), Nsets(setvec.size()), allow_x2y(NRV==NOX)
{
  if (allow_x2y) {
    for (tuint i=0;i<Nsets;++i) {
      allow_x2y = allow_x2y && setvec[i]->allow_x2y();
    }
  }
}

RBRV_constructor::RBRV_constructor(const std::vector<std::string>& set_str_vec, RBRV_set_box &rbrv_box)
{
  RBRV_constructor::find_dependent_sets(set_str_vec,setvec,rbrv_box);
  NRV = count_NRV(setvec);
  NOX = count_NOX(setvec);
  Nsets = setvec.size();
  allow_x2y = (NRV==NOX);
}

void RBRV_constructor::transform_y2x()
{
  set_is_valid(false);
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->transform_y2x();
  }
}

void RBRV_constructor::transform_x2y()
{
  if (allow_x2y==false) {
    throw FlxException("RBRV_constructor::transform_x2y","A x2y-transformation is not allowed for this set of random variables.");
  }
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->transform_x2y();
  }
}

void RBRV_constructor::propose_y(flxVec& yV)
{
  #if FLX_DEBUG
    if (yV.get_N()!=NRV) throw FlxException_Crude("RBRV_constructor::propose_y");
  #endif
  set_is_valid(false);
  tuint rvc = 0;
  tdouble* const yp = yV.get_tmp_vptr();
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    const tuint N = cs->get_NRV();
    flxVec yt(yp+rvc,N);
    yt = cs->propose_y();
    rvc += N;
  }
}

void RBRV_constructor::propose_y()
{
  set_is_valid(false);
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->propose_y();
  }
}

void RBRV_constructor::gen_smp()
{
  propose_y();
  transform_y2x();
}

void RBRV_constructor::set_smp(const flxVec& y_new)
{
  set_smp_y(y_new);
  transform_y2x();
}

void RBRV_constructor::set_smp_y(const flxVec& y_new)
{
  set_is_valid(false);
  tuint rvc = 0;
  const tdouble* const yp = y_new.get_tmp_vptr_const();
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->set_y(yp+rvc);
    rvc += cs->get_NRV();
  }
}

void RBRV_constructor::set_smp_x(const flxVec& x_new)
{
  set_is_valid(false);
  tuint rvc = 0;
  const tdouble* const xp = x_new.get_tmp_vptr_const();
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->set_x(xp+rvc);
    rvc += cs->get_NOX();
  }
}

void RBRV_constructor::set_smp_x_transform(const flxVec& x_new)
{
  set_smp_x(x_new);
  transform_x2y();
}

void RBRV_constructor::set_is_valid(const bool is_valid)
{
  #if FLX_DEBUG
    for (std::vector<RBRV_set_base*>::const_iterator it = setvec.begin() ; it != setvec.end(); ++it) {
      (*it)->set_is_valid(is_valid);
    }
  #endif
}

void RBRV_constructor::get_y_Vec(tdouble*const y_vec) const
{
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->get_y(y_vec+rvc);
    rvc += cs->get_NRV();
  }
}

void RBRV_constructor::get_x_Vec(tdouble*const x_vec) const
{
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->get_x(x_vec+rvc);
    rvc += cs->get_NOX();
  }
}

void RBRV_constructor::get_mean_Vec(tdouble*const m_vec) const
{
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->get_mean(m_vec+rvc);
    rvc += cs->get_NOX();
  }
}

void RBRV_constructor::get_sd_Vec(tdouble*const s_vec) const
{
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->get_sd(s_vec+rvc);
    rvc += cs->get_NOX();
  }
}

const bool RBRV_constructor::check_xVec(const flxVec& xV) const
{
  const tdouble* xp = xV.get_tmp_vptr_const();
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    if (cs->check_xVec(xp+rvc)==false) return false;
    rvc += cs->get_NOX();
  }
  return true;
}

void RBRV_constructor::calc_Jinv(FlxMtxLTri& dxdy)
{
  if (get_NRV()!=get_NOX()) throw FlxException_NotImplemented("RBRV_constructor::calc_Jinv_01");
  if (get_NRV()!=dxdy.ncols()) throw FlxException_Crude("RBRV_constructor::calc_Jinv_02");
  dxdy.set_zeroMtx();
  flxVec dxdw(get_NRV());
  tdouble* dxdwp = dxdw.get_tmp_vptr();
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    const tuint csNOX = cs->get_NOX();
    const FlxMtxSparsLTri* L = cs->calc_Jinv_1();
    cs->calc_Jinv_2(dxdwp+rvc);
    if (L) {
      for (tuint k=0;k<csNOX;++k) {
        for (tuint l=0;l<=k;++l) {
          dxdy.add_value(k+rvc,l+rvc,L->operator()(k,l)*dxdwp[rvc+k]);
        }
      }
    } else {
      for (tuint k=0;k<csNOX;++k) {
        dxdy.add_value(k+rvc,k+rvc,dxdwp[rvc+k]);
      }
    }
    rvc += csNOX;
  }
}

void RBRV_constructor::transform_y2w(const tdouble* const y_vec, tdouble* const w_vec)
{
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->transform_y2w(y_vec+rvc,w_vec+rvc);
    rvc += cs->get_NRV();
  }
}


void RBRV_constructor::find_dependent_sets(const std::vector<std::string>& set_str_vec, std::vector< RBRV_set_base* >& setvec, RBRV_set_box& RBRVbox)
{
  if (set_str_vec.empty()) {
    std::ostringstream ssV;
    ssV << "An empty list of sets is not allowed.";
    throw FlxException("RBRV_constructor::find_dependent_sets_1", ssV.str() );
  }
  // collect all relevant sets of random variables
    for (size_t i=0;i<set_str_vec.size();++i) {
      RBRV_set_base* sp = RBRVbox.get_set(set_str_vec[i],true);
      sp->find_dependent_sets(setvec);
    }
  // count the number of random variables
    const tuint NRVt = count_NRV_long(setvec);
    if (NRVt==0) {
      std::ostringstream ssV;
      ssV << "The set does not contain any random variables.";
      throw FlxException("RBRV_constructor::find_dependent_sets_2", ssV.str() );
    }
  tuint i;
  tuint nS = setvec.size();
  for (i=0;i<nS;++i) {
    tuint rm = setvec[nS-i-1]->group_dependent_sets(setvec,nS-i-1);
    nS -= rm;
  }
  #if FLX_DEBUG
    if (NRVt!=count_NRV(setvec)) throw FlxException_Crude("RBRV_constructor::find_dependent_sets");
  #endif
}

const tuint RBRV_constructor::count_NRV(const std::vector< RBRV_set_base* >& setvec)
{
  tuint tN = 0;
  for (tuint i=0;i<setvec.size();++i) {
    RBRV_set_base* srb = setvec[i];
    tN += srb->get_NRV();
  }
  return tN;
}

const tuint RBRV_constructor::count_NRV_long(const std::vector< RBRV_set_base* >& setvec)
{
  tuint tN = 0;
  for (tuint i=0;i<setvec.size();++i) {
    RBRV_set_base* srb = setvec[i];
    tN += srb->get_NRV_only_this();
  }
  return tN;
}

const tuint RBRV_constructor::count_NOX(const std::vector< RBRV_set_base* >& setvec)
{
  tuint tN = 0;
  for (tuint i=0;i<setvec.size();++i) {
    RBRV_set_base* srb = setvec[i];
    tN += srb->get_NOX();
  }
  return tN;
}

const tuint RBRV_constructor::count_NOX_long(const std::vector< RBRV_set_base* >& setvec)
{
  tuint tN = 0;
  for (tuint i=0;i<setvec.size();++i) {
    RBRV_set_base* srb = setvec[i];
    tN += srb->get_NOX_only_this();
  }
  return tN;
}

void RBRV_constructor::print_info(std::ostream& sout, const std::string prelim)
{
  sout << prelim << "Number of random variables in standard normal space: " << get_NRV() << std::endl;
  sout << prelim << "Number of random variables in original space:        " << get_NOX() << std::endl;
  sout << prelim << "Number of sets in the constructor:                   " << Nsets << std::endl;
  sout << prelim << "Sets in the constructor:" << std::endl;
  tuint counter = 0;
  for (tuint i=0;i<Nsets;++i) {
    setvec[i]->print(sout,prelim,counter,true);
  }
  #if FLX_DEBUG
    if (get_NOX()!=counter) throw FlxException_Crude("RBRV_constructor::print_info");
  #endif
}

const std::string RBRV_constructor::get_rv_name(const tuint index)
{
  if (index>=NOX) throw FlxException_Crude("RBRV_constructor::get_rv_name_1");
  tuint c = 0;
  for (tuint j=0;j<Nsets;++j) {
    RBRV_set_base* sj = setvec[j];
    if (c+sj->get_NOX_only_this()>index) {
      return sj->get_rv_name(index-c);
    }
    c+=sj->get_NOX_only_this();
  }
  throw FlxException_Crude("RBRV_constructor::get_rv_name_2");
}



