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

#include "flxBayUp.h"
#include "flxobjrandom.h"
#include "flxfunction_ope_calc.h"
#include "flxobjrbrv.h"
#include "flxfunction_fun_calc.h"

#include <algorithm>
#include <fstream>
#include <limits>

#if FLX_DEBUG
  int FlxBayUpBox::Cinst = 0;
#endif
  

// ====================================================================================================
// 1D - Optimization
// ====================================================================================================

tdouble perfFun1D_flxBayUp_adaptive_ctrl_opti_jump(const tdouble x, void *params) {
  flxBayUp_adaptive_ctrl_opti_jump *p = (flxBayUp_adaptive_ctrl_opti_jump *)params;
  const tdouble sd = rv_Phi(x);
  return -(p->perfFun(sd));
}
  
 
// ====================================================================================================
// GSL - Optimization
// ====================================================================================================

#include <gsl/gsl_multimin.h>
  

double LSF_f (const gsl_vector *v, void *params)
{
  flxBayUp_adaptive_ctrl_dcs *p = (flxBayUp_adaptive_ctrl_dcs *)params;
  
  const tdouble x = rv_Phi(gsl_vector_get(v, 0));
  const tdouble y = rv_Phi(gsl_vector_get(v, 1));

  const tdouble res = -(p->LSF(x,y,0));
  
  return res;
}

double LSF_f2 (const gsl_vector *v, void *params)
{
  flxBayUp_adaptive_ctrl_dcs *p = (flxBayUp_adaptive_ctrl_dcs *)params;
  
  const tdouble x = rv_Phi(gsl_vector_get(v, 0));

  const tdouble res = -(p->LSF(x,ONE,1));
  
  return res;
}

double LSF_f3 (const gsl_vector *v, void *params)
{
  flxBayUp_adaptive_ctrl_dcs *p = (flxBayUp_adaptive_ctrl_dcs *)params;
  
  const tdouble x = rv_Phi(gsl_vector_get(v, 0));

  const tdouble res = -(p->LSF(ONE,x,2));
  
  return res;
}

double LSF_f4 (const gsl_vector *v, void *params)
{
  flxBayUp_adaptive_ctrl_dcs *p = (flxBayUp_adaptive_ctrl_dcs *)params;
  
  const tdouble x = rv_Phi(gsl_vector_get(v, 0));

  const tdouble res = -(p->LSF(ONE,x,3));
  
  return res;
}

double LSF_f5 (const gsl_vector *v, void *params)
{
  flxBayUp_adaptive_ctrl_dcs *p = (flxBayUp_adaptive_ctrl_dcs *)params;
  
  const tdouble x = rv_Phi(gsl_vector_get(v, 0));

  const tdouble res = -(p->LSF(ONE,x,4));
  
  return res;
}

 
// ====================================================================================================

  
FlxBayUpBox::FlxBayUpBox()
{
  #if FLX_DEBUG
    ++Cinst;
    if (Cinst > 1) {
      std::ostringstream ssV;
      ssV << "More than one instance of 'FlxBayUpBox' created ...";
      throw FlxException("FlxBayUpBox::FlxBayUpBox", ssV.str() );
    }
  #endif
}

FlxBayUpBox::~FlxBayUpBox()
{
  for (std::map<std::string, flxBayUp*>::iterator pos = box.begin(); pos != box.end(); ++pos) {
    delete pos->second;
  }
  for (std::map<std::string, flxBayDA*>::iterator pos = box_DA.begin(); pos != box_DA.end(); ++pos) {
    delete pos->second;
  }
}

void FlxBayUpBox::insert(const std::string& nameID, flxBayUp* obj, const bool errSerious)
{
  std::pair<std::string, flxBayUp*> Element(nameID, obj);
  std::map<std::string, flxBayUp*>::iterator pos = box.find(nameID);
  if ( pos != box.end() ) {
    delete pos->second; 
    pos->second = obj;
    return;
  }
  box.insert(Element);
}

flxBayUp& FlxBayUpBox::get(const std::string& nameID)
{
  std::map<std::string, flxBayUp*>::iterator pos = box.find(nameID);
  if ( pos != box.end() ) {
    flxBayUp* res = pos->second;
    if (res==NULL) {
      std::ostringstream ssV;
      ssV << "The BayUp-object '" << nameID << "' does not exist.";
      throw FlxException("FlxBayUpBox::get_1", ssV.str(), "In oder to use a BayUp-object, you have to define it first."); 
    }
    return (*res);
  } else {
    std::ostringstream ssV;
    ssV << "The BayUp-object '" << nameID << "' does not exist.";
    throw FlxException("FlxBayUpBox::get_2", ssV.str(), "In oder to use a BayUp-object, you have to define it first."); 
  }
}

void FlxBayUpBox::insert_DA(const std::string& nameID, flxBayDA* obj, const bool errSerious)
{
  std::pair<std::string, flxBayDA*> Element(nameID, obj);
  std::map<std::string, flxBayDA*>::iterator pos = box_DA.find(nameID);
  if ( pos != box_DA.end() ) {
    std::ostringstream ssV;
    ssV << "The BayDA-object '" << nameID << "' exists already.";
    throw FlxException("FlxBayUpBox::insert_DA", ssV.str());
  }
  box_DA.insert(Element);
}

flxBayDA& FlxBayUpBox::get_DA(const std::string& nameID)
{
  std::map<std::string, flxBayDA*>::iterator pos = box_DA.find(nameID);
  if ( pos != box_DA.end() ) {
    flxBayDA* res = pos->second;
    if (res==NULL) {
      std::ostringstream ssV;
      ssV << "The BayDA-object '" << nameID << "' does not exist.";
      throw FlxException("FlxBayUpBox::get_DA_1", ssV.str(), "In oder to use a BayDA-object, you have to define it first.");
    }
    return (*res);
  } else {
    std::ostringstream ssV;
    ssV << "The BayDA-object '" << nameID << "' does not exist.";
    throw FlxException("FlxBayUpBox::get_DA_2", ssV.str(), "In oder to use a BayDA-object, you have to define it first.");
  }
}


flxBayUp_adaptive_ctrl_base::flxBayUp_adaptive_ctrl_base(FlxFunction* updatesAfterNsamples, const tuint smpl_order)
: updatesAfterNsamples(updatesAfterNsamples), smpl_order(smpl_order)
{
  if (smpl_order>3) {
    std::ostringstream ssV;
    ssV << "ID of adaptive_smpl_order '" << smpl_order << "' is not valid.";
    // free memory
      delete updatesAfterNsamples;
    throw FlxException("flxBayUp_adaptive_ctrl_base::flxBayUp_adaptive_ctrl",ssV.str());
  }
}

flxBayUp_adaptive_ctrl_base::~flxBayUp_adaptive_ctrl_base()
{
  delete updatesAfterNsamples;
}

const tuint flxBayUp_adaptive_ctrl_base::get_updatesAfterNsamples()
{
  return updatesAfterNsamples->cast2tuintW0();
}

const tuint flxBayUp_adaptive_ctrl_base::get_smpl_order() const
{
  if (is_adaptive()==false) return 0;
  return smpl_order;
}

void flxBayUp_adaptive_ctrl_base::print_info(std::ostream& sout) const
{
  sout << "  perform adaptive step after:  " << updatesAfterNsamples->write() << " model calls" << std::endl;
  sout << "  order to process seeds:       ";
  switch (smpl_order) {
    case 0:
      sout << "no reordering";
      break;
    case 1:
      break;
    case 2:
      sout << "random reordering";
    case 3:
      break;
  }
  sout << std::endl;
}

flxBayUp_adaptive_ctrl_bounds::flxBayUp_adaptive_ctrl_bounds(FlxFunction* factor, FlxFunction* lower, FlxFunction* upper, FlxFunction* maxUpdatesPerCStep, const tuint smpl_order)
: flxBayUp_adaptive_ctrl_base(maxUpdatesPerCStep,smpl_order), factor_d(ZERO),lower_d(ZERO),upper_d(ZERO), factor(factor),lower(lower),upper(upper)
{
  
}

flxBayUp_adaptive_ctrl_bounds::~flxBayUp_adaptive_ctrl_bounds()
{
  delete upper;
  delete lower;
  delete factor;
}

flxBayUp_adaptive_ctrl_bounds* flxBayUp_adaptive_ctrl_bounds::copy()
{
  return new flxBayUp_adaptive_ctrl_bounds(new FlxFunction(*factor),new FlxFunction(*lower),new FlxFunction(*upper),new FlxFunction(*updatesAfterNsamples),smpl_order);
}

void flxBayUp_adaptive_ctrl_bounds::eval()
{
  factor_d = factor->cast2positive_or0();
  lower_d = lower->cast2positive_or0();
  upper_d = upper->cast2positive_or0();
  if (upper_d < lower_d) {
    std::ostringstream ssV;
    ssV << "Lower bound '" << lower_d << "' must be smaller than the upper '" << upper_d << "' bound.";
    throw FlxException("SSS_adaptive_ctrl::eval",ssV.str());
  }
}

void flxBayUp_adaptive_ctrl_bounds::requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm)
{  
  if (acr<lower_d) {
    csm.adptv_spread_multiply(ONE-factor_d);
  } else if (acr>upper_d) {
    csm.adptv_spread_multiply(ONE+factor_d);
  }
}

void flxBayUp_adaptive_ctrl_bounds::print_info(std::ostream& sout) const
{
  sout << "  adaptive factor:              " << factor->write() << std::endl;
  sout << "  adaptive acr-bounds:          [" << lower->write() << ";" << upper->write() << "]" << std::endl;
  flxBayUp_adaptive_ctrl_base::print_info(sout);
}

flxBayUp_adaptive_ctrl_log::flxBayUp_adaptive_ctrl_log(FlxFunction* f1, FlxFunction* f2, FlxFunction* ftacr, FlxFunction* maxUpdatesPerCStep, const tuint smpl_order)
: flxBayUp_adaptive_ctrl_base(maxUpdatesPerCStep,smpl_order), f1(f1),f2(f2),ftacr(ftacr)
{

}

flxBayUp_adaptive_ctrl_log::~flxBayUp_adaptive_ctrl_log()
{
  delete f1;
  delete f2;
  delete ftacr;
}

flxBayUp_adaptive_ctrl_log* flxBayUp_adaptive_ctrl_log::copy()
{
  return new flxBayUp_adaptive_ctrl_log(new FlxFunction(*f1),new FlxFunction(*f2),new FlxFunction(*ftacr),new FlxFunction(*updatesAfterNsamples),smpl_order);
}

void flxBayUp_adaptive_ctrl_log::requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm)
{
  const tdouble ac1d = csm.get_ac1d();
  const tdouble d1 = f1->cast2positive_or0();
  const tdouble d2 = f2->cast2positive();
  const tdouble tacr = ftacr->cast2positive();
  const tdouble acr_ = (ac1d<acr && ac1d<tacr)?ac1d:acr;
  const tdouble f = d2*exp(d1*(acr_-tacr));
  csm.adptv_spread_multiply(f);
}

void flxBayUp_adaptive_ctrl_log::print_info(std::ostream& sout) const
{
  sout << "  adaptive factors:             f1=" << f1->write() << "; f2=" << f2->write() << "; target_acr=" << ftacr->write() << std::endl;
  flxBayUp_adaptive_ctrl_base::print_info(sout);
}

flxBayUp_adaptive_ctrl_dcs::flxBayUp_adaptive_ctrl_dcs(FlxFunction* maxUpdatesPerCStep, FlxFunction* factor, FlxFunction* psd_max, const tuint smpl_order)
: flxBayUp_adaptive_ctrl_base(maxUpdatesPerCStep, smpl_order), csm_dcs(NULL), csm_csus(NULL), factor_d(ZERO), factor(factor), psd_max_d(0.9), psd_max(psd_max),
  smpl_i(0), smpl_N(0), smpl_Nmax(0), smpl_list(NULL), smpl_acc(NULL), 
  smpl_weightsP(NULL), smpl_resP(NULL),
  cur_sdR_ut(ZERO), cur_sdW_ut(ZERO), cur_sdWS_ut(ZERO), cur_pSD(ZERO), shift(8),
  omega_sum(ZERO), omega_N(0), sdR_sum(ZERO), sdR_N(0), sdW_sum(ZERO), sdW_N(0), sdWS_sum(ZERO), sdWS_N(0), pSD_sum(ZERO), pSD_N(0)
{
  
}

flxBayUp_adaptive_ctrl_dcs::~flxBayUp_adaptive_ctrl_dcs()
{
  if (factor) delete factor;
  if (psd_max) delete psd_max;
  if (smpl_list) delete [] smpl_list;
  if (smpl_acc) delete [] smpl_acc;
  if (smpl_weightsP) delete [] smpl_weightsP;
  if (smpl_resP) delete [] smpl_resP;
}

flxBayUp_adaptive_ctrl_dcs* flxBayUp_adaptive_ctrl_dcs::copy()
{
  return new flxBayUp_adaptive_ctrl_dcs(new FlxFunction(*updatesAfterNsamples),new FlxFunction(*factor),new FlxFunction(*psd_max),smpl_order);
}

void flxBayUp_adaptive_ctrl_dcs::eval()
{
  smpl_i = 0;
  smpl_N = 0;
  omega_sum = ZERO;
  omega_N = 0;
  // remember optimal spread from last level
    if (sdR_N>0) {
      sdR_sum /= sdR_N;
      sdR_N = 1;
    }
    if (sdW_N>0) {
      sdW_sum /= sdW_N;
      sdW_N = 1;
    }
    if (sdWS_N>0) {
      sdWS_sum /= sdWS_N;
      sdWS_N = 1;
    }
    if (pSD_N>0) {
      pSD_sum /= pSD_N;
      pSD_N = 1;
    }
  factor_d = factor->cast2positive_or0(false);
  psd_max_d = psd_max->cast2positive(false);
}

void flxBayUp_adaptive_ctrl_dcs::plot_shift()
{
  const tuint sold = shift;
  // plot 
    std::ostringstream ssV;
    ssV << "plot_" << callN << ".txt";
    std::ofstream ofs (ssV.str().c_str(), std::ofstream::out);
    for (tdouble y1=-4;y1<=4;y1+=tdouble(8)/100) {
      for (tdouble y2=-6;y2<=6;y2+=tdouble(12)/100) {
        const tdouble x = rv_Phi(y1);
        const tdouble y = rv_Phi(y2);
        shift = 6;
          const tdouble s6 = LSF(x,y,true);
        shift = 7;
          const tdouble s7 = LSF(x,y,true);
        shift = 8;
          const tdouble s8 = LSF(x,y,true);
        shift = 9;
          const tdouble s9 = LSF(x,y,true);
        shift = 10;
          const tdouble s10 = LSF(x,y,true);
        shift = 11;
          const tdouble s11 = LSF(x,y,true);
        shift = 12;
          const tdouble s12 = LSF(x,y,true);
        shift = 13;
          const tdouble s13 = LSF(x,y,true);
        shift = 14;
          const tdouble s14 = LSF(x,y,true);
        shift = 15;
          const tdouble s15 = LSF(x,y,true);
        shift = 16;
          const tdouble s16 = LSF(x,y,true);

        //             1               2            3             4             5            6             7              8               9              10              11           12             13
        ofs << "\t" << y1 << "\t" << y2 << "\t" << s6 << "\t" << s7 << "\t" << s8 << "\t" << s9 << "\t" << s10 << "\t" << s11 << "\t" << s12 << "\t" << s13 << "\t" << s14 << "\t" << s15 << "\t" << s16 << "\t" << std::endl;
      }
      ofs << std::endl;
    }
  
  shift = sold;
  ++callN;
}

void flxBayUp_adaptive_ctrl_dcs::plot_smpls()
{
  const tdouble sdR = rv_Phi(cur_sdR_ut);
  
  std::ostringstream ssV;
  ssV << "smpls_" << callN << ".txt";
  std::ofstream ofs (ssV.str().c_str(), std::ofstream::out);
  tdouble* c = smpl_list;
  bool* acc = smpl_acc;
  for (tuint i=0;i<smpl_N;++i) {
    const tdouble yR = c[0]*c[2]+c[5]*sqrt(ONE-pow2(c[2]));
        const tdouble sratio = c[2]/sdR;
        const tdouble sdRrec = 1/pow2(sdR);
        const tdouble irR = (rv_phi(sratio*c[0]+c[5]*(sqrt(sdRrec-pow2(sratio))-sqrt(sdRrec-ONE)))/sdR)
                            /(rv_phi(c[0])/c[2]);
    //     1              2               3                4              5                6                          7               8
    ofs << yR << "\t" << c[1] << "\t" << c[5] << "\t" << c[6] << "\t" << c[7] << "\t" << c[8] << "\t" << ((*acc)?ONE:ZERO) << "\t" << irR << std::endl;
    c += acdcs_dim;
    ++acc;
  }
}


const tdouble flxBayUp_adaptive_ctrl_dcs::smpl_mean(const tuint shiftV, const bool consider_acc, const bool dois)
{
  const tdouble sdR = rv_Phi(cur_sdR_ut);
  const tdouble sdW = rv_Phi(cur_sdW_ut);
  
  tdouble res0 = ZERO;
  size_t count0 = 0;
  tdouble res = ZERO;
  tdouble irs = ZERO;
  tdouble* c = smpl_list;
  bool* acc = smpl_acc;
  for (tuint i=0;i<smpl_N;++i) {
    // importance sampling ratio
      // sdR
        const tdouble sratio = c[2]/sdR;
        const tdouble sdRrec = 1/pow2(sdR);
        const tdouble irR = (rv_phi(sratio*c[0]+c[5]*(sqrt(sdRrec-pow2(sratio))-sqrt(sdRrec-ONE)))/sdR)
                            /(rv_phi(c[0])/c[2]);
      // sdW
        const tdouble irW = (rv_phi(c[1]/sdW)/sdW)
                            /(rv_phi(c[1]/c[3])/c[3]);
      const tdouble ir = irR*irW;
      
    if ( (consider_acc==false) || *acc) {
      tdouble q = ZERO;
      if (shiftV<acdcs_dim) {
        q = c[shiftV];
      } else {
        switch (shift-acdcs_dim) {
          case 0:
            q += (c[6]*c[8]);
            break;
          default:
            throw FlxException_Crude("flxBayUp_adaptive_ctrl_dcs::LSF");
        }
      }
      res0 += q;
      ++count0;
      res += q*ir;
      irs += ir;
    }
    c += acdcs_dim;
    ++acc;
  }
  res0 /= count0;
  res /= irs;
  if (dois) return res;
  else return res0;
}

void flxBayUp_adaptive_ctrl_dcs::do_gsl_opti(const tuint mode)
{
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  /* Starting point */
    x = gsl_vector_alloc(mode?1:2);
    switch (mode) {
      case 0:
        gsl_vector_set (x, 0, cur_sdR_ut);
        gsl_vector_set (x, 1, cur_sdW_ut);
        break;
      case 1:
        gsl_vector_set (x, 0, cur_sdR_ut);
        break;
      case 2:
        gsl_vector_set (x, 0, cur_sdW_ut);
        break;
      case 3:
        gsl_vector_set (x, 0, cur_sdWS_ut);
        break;
      case 4:
        gsl_vector_set (x, 0, cur_sdR_ut);
        break;
    };
  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc(mode?1:2);
  gsl_vector_set_all(ss, ONE);

  /* Initialize method and iterate */
  minex_func.n = (mode?1:2);
  switch (mode) {
    case 0:
      minex_func.f = LSF_f;
      break;
    case 1:
      minex_func.f = LSF_f2;
      break;
    case 2:
      minex_func.f = LSF_f3;
      break;
    case 3:
      minex_func.f = LSF_f4;
      break;
    case 4:
      minex_func.f = LSF_f5;
      break;
  };
  minex_func.params = this;

  s = gsl_multimin_fminimizer_alloc (T, (mode?1:2));
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status) 
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS)
      {
//         printf ("converged to minimum at\n");
        switch (mode) {
          case 0:
            cur_sdR_ut = gsl_vector_get (s->x, 0);
            cur_sdW_ut = gsl_vector_get (s->x, 1);
            break;
          case 1:
            cur_sdR_ut = gsl_vector_get (s->x, 0);
            break;
          case 2:
            cur_sdW_ut = gsl_vector_get (s->x, 0);
            break;
          case 3:
            cur_sdWS_ut = gsl_vector_get (s->x, 0);
            break;
          case 4:
            cur_sdR_ut = gsl_vector_get (s->x, 0);
            break;
        };
    
        if (cur_sdR_ut>3.) cur_sdR_ut = 3.;
        if (cur_sdW_ut>3.) cur_sdW_ut = 3.;
        if (cur_sdWS_ut>3.) cur_sdWS_ut = 3.;
      }

//       printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
//               iter,
//               gsl_vector_get (s->x, 0), 
//               gsl_vector_get (s->x, 1), 
//               s->fval, size);

    }
  while (status == GSL_CONTINUE && iter < 100);
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

const tdouble flxBayUp_adaptive_ctrl_dcs::adopt_to_acr(const tdouble acr, const tdouble cur_sd_ut, const tdouble prev_sd)
{
  const tdouble d1 = ONE;
  const tdouble d2 = ONE;
  const tdouble tacr = 0.44;
  const tdouble f = d2*exp(d1*(acr-tacr));
  if (f*prev_sd<rv_Phi(cur_sd_ut)) {
    return rv_InvPhi_noAlert(f*prev_sd);
  } else {
    return cur_sd_ut;
  }
}

void flxBayUp_adaptive_ctrl_dcs::requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm)
{
  if (smpl_list==NULL) {
    throw FlxException_Crude("flxBayUp_adaptive_ctrl_dcs_so::requires_adptv_step_02");
  }
  if (smpl_N<=2) return;
  
  if (csm_csus) requires_adptv_step_csus(acr,csm);
  else requires_adptv_step_dcs(acr,csm);
}

void flxBayUp_adaptive_ctrl_dcs::requires_adptv_step_csus(const tdouble acr, FlxBayUP_csm_base& csm)
{
  if (csm_csus==NULL) {
    throw FlxException_Crude("flxBayUp_adaptive_ctrl_dcs::requires_adptv_step_csus");
  }
  
  const tdouble limit_acr = 0.2;        // try 'desperately' to increase the acceptance rate if acr smaller than this value
  
  const tdouble prev_sd = rv_Phi(cur_sdR_ut);
  shift = 8;
  do_gsl_opti(4); 
  // average out optimal spread
    sdR_sum += cur_sdR_ut;
    ++sdR_N;
    cur_sdR_ut = sdR_sum/sdR_N;
  if (acr<limit_acr) {
    cur_sdR_ut = adopt_to_acr(acr,cur_sdR_ut,prev_sd);
  }
  
  csm_csus->set_cur_spread(rv_Phi(cur_sdR_ut));
  
  smpl_i = 0;
  smpl_N = 0;
}

void flxBayUp_adaptive_ctrl_dcs::requires_adptv_step_dcs(const tdouble acr, FlxBayUP_csm_base& csm)
{
  if (csm_dcs==NULL) {
    throw FlxException("flxBayUp_adaptive_ctrl_dcs::requires_adptv_step_dcs","This adaptive strategy must be used in combination with the MCMCM algorithm 'dcs'.");
  }
  
  const tdouble limit_acr = 0.2;        // try 'desperately' to increase the acceptance rate if acr smaller than this value
  
  // optimize spread of radius 
    const tdouble prev_sdR = rv_Phi(cur_sdR_ut);
    shift = 6;
    do_gsl_opti(1); 
    // average out optimal spread
      sdR_sum += cur_sdR_ut;
      ++sdR_N;
      cur_sdR_ut = sdR_sum/sdR_N;
    if (acr<limit_acr) {
      cur_sdR_ut = adopt_to_acr(acr,cur_sdR_ut,prev_sdR);
    }
    
  if (factor_d<GlobalVar.TOL()) {
    // optimize spread of angle (isotropic)
      const tdouble prev_sdW = rv_Phi(cur_sdW_ut);
      shift = 8;
      do_gsl_opti(2); 
      // average out optimal spread
        sdW_sum += cur_sdW_ut;
        ++sdW_N;
        cur_sdW_ut = sdW_sum/sdW_N;
      if (acr<limit_acr) {
        cur_sdW_ut = adopt_to_acr(acr,cur_sdW_ut,prev_sdW);
      }
  } else {
    // acceptance rate of angle-move
      tdouble acr1=ZERO;
      tdouble acr2=ZERO;
      {
        tuint acr1N = 0;
        tuint acr2N = 0;
        tdouble* c = smpl_list;
        for (tuint i=0;i<smpl_N;++i) {
          if (c[9]==ZERO) {
            if (smpl_acc[i]) acr1 += ONE;
            ++acr1N;
          } else {
            if (smpl_acc[i]) acr2 += ONE;
            ++acr2N;
          }
          c += acdcs_dim;
        }
        if (acr1N>2) acr1 /= acr1N;
        else acr1 = acr;
        if (acr2N>2) acr2 /= acr2N;
        else acr2 = acr;
      }
      
    // optimize spread of angle (isotropic)
      const tdouble prev_sdW = rv_Phi(cur_sdW_ut);
      shift = 8;
      do_gsl_opti(2); 
      // average out optimal spread
        sdW_sum += cur_sdW_ut;
        ++sdW_N;
        cur_sdW_ut = sdW_sum/sdW_N;
      if (acr1<limit_acr) {
        cur_sdW_ut = adopt_to_acr(acr1,cur_sdW_ut,prev_sdW);
      }
      
    // optimize spread of angle (seed based)
      const tdouble prev_sdWS = rv_Phi(cur_sdWS_ut);
      shift = 8;
      do_gsl_opti(3); 
      // average out optimal spread
        sdWS_sum += cur_sdWS_ut;
        ++sdWS_N;
        cur_sdWS_ut = sdWS_sum/sdWS_N;
      if (acr2<limit_acr) {
        cur_sdWS_ut = adopt_to_acr(acr2,cur_sdWS_ut,prev_sdWS);
      }
      
    // adapt pDS
      tdouble velo1 = ZERO;
      tdouble velo2 = ZERO;
      tuint velo1N = 0;
      tuint velo2N = 0;
      tdouble* c = smpl_list;
      for (tuint i=0;i<smpl_N;++i) {
        if (c[9]==ZERO) {
          if (smpl_acc[i]) {
            velo1 += c[8];
          }
          ++velo1N;
        } else {
          if (smpl_acc[i]) {
            velo2 += c[8];
          }
          ++velo2N;
        }
        c += acdcs_dim;
      }
      if (velo1N>2 && velo2N>2) {
        const tdouble prev_pSD = cur_pSD;
        velo1 /= velo1N;
        velo2 /= velo2N;
        if (velo1>=velo2) {
          cur_pSD = ZERO;
        } else {
          cur_pSD = ONE-velo1/velo2;
        }
        if (cur_pSD>psd_max_d) {
          cur_pSD = (prev_pSD>psd_max_d)?prev_pSD:psd_max_d;
        }
        if (prev_pSD>ZERO && cur_pSD==ZERO) {
          if (prev_pSD<0.1) cur_pSD = prev_pSD;
          else cur_pSD = 0.1;
        }
        // average out optimal value
          pSD_sum += cur_pSD;
          ++pSD_N;
          cur_pSD = pSD_sum/pSD_N;
          cur_pSD = factor_d*cur_pSD+(ONE-factor_d)*prev_pSD;
      }
  }
    
  const tdouble cur_sdR = rv_Phi(cur_sdR_ut);
  const tdouble cur_sdW = rv_Phi(cur_sdW_ut);
  const tdouble cur_sdWS = rv_Phi(cur_sdWS_ut);
  csm_dcs->set_cur_spread(cur_sdR,cur_sdW,cur_sdWS,cur_pSD);
  
  // remember average omega
    omega_sum += smpl_mean(7,true,false);
    ++omega_N;
  
  smpl_i = 0;
  smpl_N = 0;
}

const tdouble flxBayUp_adaptive_ctrl_dcs::LSF(const tdouble sdR, const tdouble sdW, const tuint mode)
{
  flxVec smpl_weights(smpl_weightsP,smpl_N);
  flxVec smpl_res(smpl_resP,smpl_N);
  tuint irl1N = 0;
  tdouble irl1s = ZERO;                // sum of importance ratios larger than one
  tdouble* c = smpl_list;
  tuint smpl_Nc = 0;
  for (tuint i=0;i<smpl_N;++i) {
    // importance sampling ratio
      // sdR
        const tdouble sratio = c[2]/sdR;
        const tdouble sdRrec = 1/pow2(sdR);
        const tdouble irR = (mode>=2)?(ONE):((rv_phi(sratio*c[0]+c[5]*(sqrt(sdRrec-pow2(sratio))-sqrt(sdRrec-ONE)))/sdR)
                            /(rv_phi(c[0])/c[2]));
      // sdW
        tdouble irW;
        if (mode==4) {
          const tdouble l2n = pow2(c[7]*sqrt(c[8])-c[4]*(ONE-sqrt(ONE-pow2(sdW))))+(ONE-pow2(c[7]))*c[8];
          irW = (rv_pdf_ChiSquare(c[1],l2n/pow2(sdW))/pow2(sdW))/(rv_pdf_ChiSquare(c[1],c[0])/pow2(c[2]));
        } else {
          irW = (mode==1)?(ONE):((rv_phi(c[1]/sdW)/sdW)/(rv_phi(c[1]/c[3])/c[3]));
          if ( (mode==2 && c[9]==ONE) || (mode==3 && c[9]==ZERO) ) {
            smpl_weights[i] = ZERO;
            smpl_res[i] = ZERO;
            c += acdcs_dim;
            continue;
          }
        }
      const tdouble ir = irR*irW;
        smpl_weights[i] = ir;
        if (ir>ONE) {
          irl1s += ir;
          ++irl1N;
        }
    // quantity of interest
      if (smpl_acc[i]) {
        if (shift<acdcs_dim) {
          smpl_res[i] = c[shift];
        } else {
          switch (shift-acdcs_dim) {
            case 0:
            case 1:
              smpl_res[i] = (c[6]*c[8]);
              break;
            default:
              smpl_res[i] = ZERO;
          }
        }
      } else {
        smpl_res[i] = ZERO;
      }
    c += acdcs_dim;
    ++smpl_Nc;
  }
  if (smpl_Nc<=2) return ZERO;
  // obtain C.o.V. of sample weights
    const tdouble ir_mean = smpl_weights.get_sum()/smpl_Nc;
    if (ir_mean<GlobalVar.TOL()) return ZERO;
    const tdouble ir_sd = sqrt((pow2(smpl_weights.get_sd(ir_mean))*smpl_N)/smpl_Nc);
    tdouble ir_CoV = ir_sd/ir_mean;
      if (ir_CoV!=ir_CoV || (ir_sd==ZERO&&ir_mean!=ONE)) ir_CoV = std::numeric_limits<tdouble>::infinity();
    // make sure that C.o.V. is large enough
      const tdouble ir_CoV_low = sqrt(tdouble(smpl_Nc))*fabs(ONE-ir_mean);
    // return some values
      if (ir_CoV<ir_CoV_low) ir_CoV = ir_CoV_low;
      if (shift==acdcs_dim+3) return ir_mean;
      if (shift==acdcs_dim+4) return ir_CoV;
  // evaluate correction factors
    const tdouble ir_CoV_weight = ONE+pow(ir_CoV/5,4);
      if (shift==acdcs_dim+5) return ONE/ir_CoV_weight;
  // compute averaged (and weighted) value
    // try to correct weights
      if (irl1N>0) {
        const tdouble excess = smpl_Nc * (ir_mean-ONE);
        const tdouble target_sum = irl1s-excess;
        const tdouble target_sum1 = target_sum - irl1N;
        const tdouble base_sum1 = irl1s-irl1N;
        const tdouble tf = target_sum1/base_sum1;
        for (tuint i=0;i<smpl_N;++i) {
          if (smpl_weights[i]>ONE) {
            smpl_weights[i] = (smpl_weights[i]-ONE)*tf+ONE;
          }
        }
      } 
      const tdouble ir_mean2 = smpl_weights.get_sum()/smpl_Nc;
    // evaluate acceptance rate
      tdouble acr_weight = ONE;
      {
        tdouble acr = ZERO;
        for (tuint i=0;i<smpl_N;++i) {
          if (smpl_acc[i]) {
            acr += smpl_weights[i];
          }
        }
        acr /= smpl_Nc;
        if (ir_mean2>ONE) acr /= ir_mean2;
        if (shift==acdcs_dim+2) return acr;
        // obtain a weight
          if (acr<0.3) {
            acr_weight = ONE+pow(fabs(log(acr/0.3))/0.4,4);
          }
        if (shift==acdcs_dim+6) return ONE/acr_weight;
        if (shift==acdcs_dim+7) return ONE/(ir_CoV_weight*acr_weight);
      }
    tdouble res = smpl_res.operator*(smpl_weights) / smpl_Nc;
    // decrease estimate if weight mean is larger than ONE
      if (ir_mean2>ONE) {
        res /= ir_mean2;
      }
  // punish estimate
    if (shift==acdcs_dim) {
      res /= ir_CoV_weight*acr_weight;
    }
  return res;
}

void flxBayUp_adaptive_ctrl_dcs::register_csm(FlxBayUP_csm_dcs_MCMC* csm_dcsV)
{
  csm_dcs = csm_dcsV;
  tdouble cur_sdR, cur_sdW, cur_sdWS;
  csm_dcs->get_cur_spread(cur_sdR,cur_sdW,cur_sdWS,cur_pSD);
  cur_sdR_ut = rv_InvPhi_noAlert(cur_sdR);
  cur_sdW_ut = rv_InvPhi_noAlert(cur_sdW);
  cur_sdWS_ut = rv_InvPhi_noAlert(cur_sdWS);
  // make sure values are not too large
    if (cur_sdR_ut>3.) cur_sdR_ut = 3.;
    if (cur_sdW_ut>3.) cur_sdW_ut = 3.;
    if (cur_sdWS_ut>3.) cur_sdWS_ut = 3.;
    cur_sdR = rv_Phi(cur_sdR_ut);
    cur_sdW = rv_Phi(cur_sdW_ut);
    cur_sdWS = rv_Phi(cur_sdWS_ut);
    csm_dcs->set_cur_spread(cur_sdR,cur_sdW,cur_sdWS,cur_pSD);
    
  // allocate memory
    if (smpl_list) throw FlxException_Crude("flxBayUp_adaptive_ctrl_dcs::register_csm_01");
    smpl_Nmax = get_updatesAfterNsamples();
    smpl_list = new tdouble[smpl_Nmax*acdcs_dim];
    smpl_acc = new bool[smpl_Nmax*acdcs_dim];
    smpl_weightsP = new tdouble[smpl_Nmax];
    smpl_resP = new tdouble[smpl_Nmax];
}

void flxBayUp_adaptive_ctrl_dcs::register_csm(FlxBayUP_csm_csus_MCMC* csm_csusV)
{
  csm_csus = csm_csusV;
  tdouble cur_sd;
  csm_csus->get_cur_spread(cur_sd);
  cur_sdR_ut = rv_InvPhi_noAlert(cur_sd);
  // make sure values are not too large
    if (cur_sdR_ut>3.) cur_sdR_ut = 3.;
    cur_sd = rv_Phi(cur_sdR_ut);
    csm_csus->set_cur_spread(cur_sd);
    
  // allocate memory
    if (smpl_list) throw FlxException_Crude("flxBayUp_adaptive_ctrl_dcs::register_csm_02");
    smpl_Nmax = get_updatesAfterNsamples();
    smpl_list = new tdouble[smpl_Nmax*acdcs_dim];
    smpl_acc = new bool[smpl_Nmax*acdcs_dim];
    smpl_weightsP = new tdouble[smpl_Nmax];
    smpl_resP = new tdouble[smpl_Nmax];
}

void flxBayUp_adaptive_ctrl_dcs::append_smpl(const flxVec& lastS, const bool accpetedV)
{
  if (smpl_list==NULL) {
    throw FlxException_Crude("flxBayUp_adaptive_ctrl_dcs::append_smpl_02");
  }
  const size_t pos = acdcs_dim*smpl_i;
  flxVec sample(smpl_list+pos,acdcs_dim);  
  sample = lastS;
  smpl_acc[smpl_i] = accpetedV;
  ++smpl_i;
  if (smpl_N<smpl_i) smpl_N = smpl_i;
  if (smpl_i==smpl_Nmax) smpl_i=0;
}

void flxBayUp_adaptive_ctrl_dcs::write_adaptive_info(std::ostream& sout)
{
  if (csm_dcs==NULL) throw FlxException_Crude("flxBayUp_adaptive_ctrl_dcs::write_adaptive_info");
  sout << std::format("  sdR={:6.2e}  sdW={:6.2e}  ", rv_Phi(cur_sdR_ut), rv_Phi(cur_sdW_ut));
  if (cur_pSD>ZERO) {
    sout << std::format("sdS={:4.2f}  ", rv_Phi(cur_sdWS_ut));
    sout << std::format("pSD={:4.2f}  ", (cur_pSD));
  }
  if (omega_N>0) {
    sout << std::format("cosw={:4.2f}  ", (omega_sum/omega_N));
  }
}

flxBayUp_adaptive_ctrl_velo::flxBayUp_adaptive_ctrl_velo(FlxRndCreator& RndCreator, FlxFunction* vspread, FlxFunction* maxUpdatesPerCStep, const tuint smpl_order)
: flxBayUp_adaptive_ctrl_base(maxUpdatesPerCStep, smpl_order), vspreadh_d(ZERO), RndCreator(RndCreator), vspread(vspread),
  smpl_i(0), smpl_N(0), smpl_Nmax(0), smpl_list(NULL), cur_sd_ut(ONE),
  sd_sum(ZERO), sd_N(0), velo_seq_count(0)
{
  
}

flxBayUp_adaptive_ctrl_velo::~flxBayUp_adaptive_ctrl_velo()
{
  delete vspread;
  if (smpl_list) delete [] smpl_list;
}

flxBayUp_adaptive_ctrl_velo* flxBayUp_adaptive_ctrl_velo::copy()
{
  return new flxBayUp_adaptive_ctrl_velo(RndCreator, new FlxFunction(*vspread), new FlxFunction(*updatesAfterNsamples), smpl_order);
}

void flxBayUp_adaptive_ctrl_velo::eval()
{
  smpl_i = 0;
  smpl_N = 0;
  // remember optimal spread from last level
    if (sd_N>0) {
      sd_sum /= sd_N;
      sd_N = 1;
    }
  vspreadh_d = (vspread->cast2positive_or0(false))/2;
  // allocate memory
    if (smpl_list==NULL) {
      smpl_Nmax = get_updatesAfterNsamples();
      smpl_list = new tdouble[smpl_Nmax*acvelo_dim];
    }
}

void flxBayUp_adaptive_ctrl_velo::requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm)
{
  if (smpl_N<2) return;
  // compute upper and lower velocity
    tdouble velo1 = ZERO;
    tdouble velo2 = ZERO;
    tuint velo1N = 0;
    tuint velo2N = 0;
    const tdouble sd = rv_Phi(cur_sd_ut);
    tdouble* c = smpl_list;
    tdouble acc = ZERO;
    for (tuint i=0;i<smpl_N;++i) {
      if (c[0]<=sd) {
        ++velo1N;
        velo1 += c[1]*c[2];
      } else {
        ++velo2N;
        velo2 += c[1]*c[2];
      }
      acc += c[2];
      c += acvelo_dim;
    }
    acc /= smpl_N;
    if (acc<0.15) {
      const tdouble tacr = 0.44;
      const tdouble f = exp(acr-tacr);
      cur_sd_ut = rv_InvPhi_noAlert(f*rv_Phi(cur_sd_ut));
      return;
    }
    if (velo1N<2||velo2N<2) return;
    velo1 /= velo1N;
    velo2 /= velo2N;
    if ( velo1<GlobalVar.TOL() && velo2<GlobalVar.TOL() ) return;
  // compute spread weight
    const tdouble r = (velo2-velo1)/(velo2+velo1);
    if (r>0) {
      if (velo_seq_count>=0) ++velo_seq_count;
      else velo_seq_count /= 2;
    } else { 
      if (velo_seq_count<=0) --velo_seq_count;
      else velo_seq_count /= 2;
    }
  // adopt spread
    if (fabs(r)<GlobalVar.TOL()) return;
    const tdouble vs = get_dynamic_spread();
    ++sd_N;
    const tdouble f = pow(fabs(r),ONE/((velo_seq_count<0)?tdouble(pow2(velo_seq_count)):ONE));
    if (r>ZERO) {
      cur_sd_ut += f*vs;
    } else {
      cur_sd_ut -= f*vs;
    }
    if (cur_sd_ut>3.) cur_sd_ut = 3.;
    sd_sum += cur_sd_ut;
    cur_sd_ut = sd_sum/sd_N;
}

void flxBayUp_adaptive_ctrl_velo::append_smpl(const flxVec& lastS)
{
  if (smpl_list==NULL) {
    throw FlxException_Crude("flxBayUp_adaptive_ctrl_velo::append_smpl");
  }
  const size_t pos = acvelo_dim*smpl_i;
  flxVec sample(smpl_list+pos,acvelo_dim);  
  sample = lastS;
  ++smpl_i;
  if (smpl_N<smpl_i) smpl_N = smpl_i;
  if (smpl_i==smpl_Nmax) smpl_i=0;
}

const tdouble flxBayUp_adaptive_ctrl_velo::get_dynamic_spread() const
{
  return vspreadh_d*(ONE/pow(sd_N,0.6));
}

const tdouble flxBayUp_adaptive_ctrl_velo::get_working_sd() const
{
  const tdouble vs = get_dynamic_spread();
  const tdouble step = (RndCreator.gen_smp_binary(tdouble(0.5)))?(-vs):vs;
  return rv_Phi(cur_sd_ut+step);
}

void flxBayUp_adaptive_ctrl_velo::write_adaptive_info(std::ostream& sout)
{
  sout << std::format("  h={:4.2f}   ", rv_Phi(cur_sd_ut));
}

void flxBayUp_adaptive_ctrl_velo::print_info(std::ostream& sout) const
{
  sout << "  adaptive factors:             vspread=" << vspread->write() << std::endl;
  flxBayUp_adaptive_ctrl_base::print_info(sout);
}

flxBayUp_adaptive_ctrl_opti_jump::flxBayUp_adaptive_ctrl_opti_jump(FlxRndCreator& RndCreator, FlxFunction* acr_min, FlxFunction* esjd_scale, FlxFunction* pw_p1, FlxFunction* pw_p2, FlxFunction* aeps, FlxFunction* Nmax_fun, FlxFunction* maxUpdatesPerCStep, const tuint smpl_order)
: flxBayUp_adaptive_ctrl_base(maxUpdatesPerCStep, smpl_order), RndCreator(RndCreator), acr_min(acr_min), esjd_scale(esjd_scale), pw_p1(pw_p1), pw_p2(pw_p2), aeps(aeps), Nmax_fun(Nmax_fun),
  M_dim(0), sd_prev_min(-ONE),
  smpl_i(0), smpl_N(0), smpl_Nmax(0), smpl_list(NULL), swl(NULL), twl(NULL),spread_hist_vec(spread_hist_N_max),spread_hist_N(0),
  deactivated(false), limit_acr(ZERO), ipds(100)
{
  
}

flxBayUp_adaptive_ctrl_opti_jump::~flxBayUp_adaptive_ctrl_opti_jump()
{
  delete Nmax_fun;
  delete acr_min;
  delete esjd_scale;
  if (smpl_list) delete [] smpl_list;
  if (swl) delete [] swl;
  if (twl) delete [] twl;
  if (pw_p1) delete pw_p1;
  if (pw_p2) delete pw_p2;
  if (aeps) delete aeps;
}

flxBayUp_adaptive_ctrl_opti_jump* flxBayUp_adaptive_ctrl_opti_jump::copy()
{
  return new flxBayUp_adaptive_ctrl_opti_jump(RndCreator, new FlxFunction(*acr_min), new FlxFunction(*esjd_scale), new FlxFunction(*pw_p1), new FlxFunction(*pw_p2), new FlxFunction(*aeps), new FlxFunction(*Nmax_fun), new FlxFunction(*updatesAfterNsamples), smpl_order);
}

void flxBayUp_adaptive_ctrl_opti_jump::initialize(const tuint M_dim_v, const tuint Nfinal)
{
  M_dim = M_dim_v;
  // allocate memory
    if (smpl_list==NULL) {
      smpl_Nmax = Nmax_fun->cast2tuint();
      if (smpl_Nmax>Nfinal) smpl_Nmax = Nfinal;
      smpl_list = new tdouble[smpl_Nmax*acopti_jump_dim];
      swl = new tdouble[smpl_Nmax];
      twl = new tdouble[smpl_Nmax];
    } else {
      throw FlxException_Crude("flxBayUp_adaptive_ctrl_opti_jump::initialize");
    }
}

void flxBayUp_adaptive_ctrl_opti_jump::eval()
{
  smpl_i = 0;
  smpl_N = 0;
  spread_hist_N = 0;
  deactivated = false;
}

const tdouble flxBayUp_adaptive_ctrl_opti_jump::get_pweight()
{
  // determine step-number
    const tuint istep = spread_hist_N+1;
    const tdouble x = tdouble(smpl_N)/(spread_hist_N+1);
  // determine 1st parameter
    const tdouble b = 500./pow(1.2,2.9*pow(x-9.,ONE/2))+0.7;
  // determine 2nd parameter
    const tdouble a = pow((log(b+0.9)/log(5.)),0.52);
  tdouble res = pow(istep+a,6.)/pow(pow(istep,3)+b,2.)/5;
  res = pow(res,pw_p1->cast2positive_or0())*pw_p2->cast2positive();
  if (res>ONE) res = ONE;
  if (res<ZERO) res = ZERO;
  return res;
}

const tdouble flxBayUp_adaptive_ctrl_opti_jump::compute_overall_acr()
{
  const tdouble* c = smpl_list+4;
  tuint ac = 0;
  for (tuint i=0;i<smpl_N;++i) {
    if (*c>ZERO) ++ac;
    c += acopti_jump_dim;
  }
  const tdouble acr_avg = tdouble(ac)/smpl_N;
  return acr_avg;
}

const bool flxBayUp_adaptive_ctrl_opti_jump::skip_adaptive_step(const tdouble acr)
{
  // make sure that enough data points are available
    if (spread_hist_N<3) return false;
  // evaluate the average acceptance rate
    // perform the adaptive step for sure if average acceptance rate is too small
    if (deactivated) {
      if (acr>limit_acr) return true;
      if (compute_overall_acr()>limit_acr) return true;
      deactivated = false;
      return false;
    }
  // estimate sample statistics
    const tuint N = (spread_hist_N<spread_hist_N_max)?spread_hist_N:spread_hist_N_max; 
    tdouble w = ONE;
    tdouble wms = ZERO;
    pdouble wssp;
    tdouble ws = ZERO;
    const tdouble wf = 0.6;        // factor by which the weight is multiplied in each step
    for (tuint i=0;i<N;++i) {
      wms += w*spread_hist_vec[i];
      wssp += w*pow2(spread_hist_vec[i]);
      ws += w;
      w *= wf;
    }
    wms /= ws;
    wssp -= ws*pow2(wms);
    const tdouble wss = sqrt(wssp.cast2double()/((ws*(N-1))/N));
    const tdouble wCoV = wss/wms;
    if (wCoV<aeps->cast2positive_or0()) {
      deactivated = true;
      return true;
    }
  return false;
}

// const tdouble flxBayUp_adaptive_ctrl_opti_jump::proposal_pdf_ln(const tdouble* c, const tdouble sd, tdouble &ej) const
// {
//   const tdouble ls = c[0];
//   const tdouble a = c[1];
//   const tdouble b2 = c[2];
//   const tdouble sd2 = pow2(sd);
//   const bool acpt = (c[4]>ZERO);
//   const tdouble lrho = sqrt(ONE-sd2)*ls;
//   const tdouble k2 = pow2(ls-lrho);
//   
//   const tdouble d2 = pow2(ls-a);
//   const tdouble esjd = d2+b2;
//   const tdouble g2 = d2*k2/esjd;
//   const tdouble h2 = b2*k2/esjd;
//     
//   tdouble log_q;
//   if (acpt) {                // point is inside failure domain
//     tdouble alpha, beta, m;
//     if (a<=ls) {                // point is 'before' seed
//       alpha = (sqrt(g2)-sqrt(esjd))/sd;
//       beta = sqrt(g2)/sd;
//       m = sqrt(g2);
//     } else {                        // point is 'behind' seed
//       alpha = sqrt(g2)/sd;
//       beta = (sqrt(g2)+sqrt(esjd))/sd;
//       m = -sqrt(g2);
//     }
//     const tdouble q = (alpha<=ZERO)?(rv_Phi(beta) - rv_Phi(alpha)):(rv_Phi(-alpha)-rv_Phi(-beta));
//     log_q = log(q);
//     const tdouble tv = (rv_phi(alpha)-rv_phi(beta))/q;
//     const tdouble mean = tv*sd+m;
//     const tdouble var = sd2*(ONE+(alpha*rv_phi(alpha)-beta*rv_phi(beta))/q-pow2(tv));
//     ej = var+pow2(mean);
//   } else {                // point is outside failure domain
//     ej = ZERO;
//     if (a<=ls) {                // point is 'before' seed
//       const tdouble beta = (sqrt(g2)-sqrt(esjd))/sd;
//       const tdouble q = rv_Phi(beta);
//       log_q = log(q);
//     } else {                        // point is 'behind' seed
//       const tdouble alpha = (sqrt(g2)+sqrt(esjd))/sd;
//       const tdouble q = rv_Phi(-alpha);
//       log_q = log(q);
//     }
//   }
//   #if FLX_DEBUG
//     if (ej<ZERO) throw FlxException_Crude("flxBayUp_adaptive_ctrl_opti_jump::proposal_pdf_ln");
//   #endif
//   return log_q - h2/(sd2*2) - log(sd)*M_dim;
// }
  
const tdouble flxBayUp_adaptive_ctrl_opti_jump::proposal_pdf_ln(const tdouble* c, const tdouble sd) const 
{
  const tdouble sd2 = pow2(sd);
  const tdouble l2 = pow2(c[0]*sqrt(ONE-sd2)-c[1])+c[2];
//   return exp(-l2/(sd2*2))/pow(sd2,tdouble(M_dim)/2);
  return -l2/(sd2*2) - log(sd)*M_dim;
//   return exp(-l2/(sd2*2))/pow_int<tdouble>(sd,M_dim);
}

void flxBayUp_adaptive_ctrl_opti_jump::compute_seed_weights()
{
  // find different rho* employed during the sampling process
    std::vector<tdouble> arv;
    std::vector<tdouble> anv;
    const tdouble* c = smpl_list+3;
    tdouble r_prev = *c;
    tuint n_prev = 1;
    arv.push_back(r_prev);
    for (tuint i=1;i<smpl_N;++i) {
      c += acopti_jump_dim;
      if (*c==r_prev) {
        ++n_prev;
      } else {
        anv.push_back(tdouble(n_prev)/smpl_N);
        n_prev = 1;
        r_prev = *c;
        arv.push_back(r_prev);
      }
    }
    anv.push_back(tdouble(n_prev)/smpl_N);
    const tuint Ndr = arv.size();
    #if FLX_DEBUG
      if (Ndr!=anv.size()) throw FlxException_Crude("flxBayUp_adaptive_ctrl_opti_jump::compute_seed_weights_01");
      tdouble s = ZERO;
      for (tuint i=0;i<Ndr;++i) {
        s += anv[i];
      }
      if (fabs(s-ONE)>1e-8) throw FlxException_Crude("flxBayUp_adaptive_ctrl_opti_jump::compute_seed_weights_02");
    #endif
  // compute the contribution of the seeds to the weights
    c = smpl_list;
    for (tuint i=0;i<smpl_N;++i) {
      tdouble s = ZERO;
      if (Ndr==1) {
        s = proposal_pdf_ln(c,arv[0]);
      } else {
        for (tuint j=0;j<Ndr;++j) {
          s += anv[j]*exp(proposal_pdf_ln(c,arv[j]));
        }
        s = log(s);
      }
      swl[i] = s;
      c += acopti_jump_dim;
    }
}

void flxBayUp_adaptive_ctrl_opti_jump::requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm)
{
  if (smpl_list==NULL || M_dim==0) {
    throw FlxException_Crude("flxBayUp_adaptive_ctrl_opti_jump::requires_adptv_step_01");
  }
  if (smpl_N<=2) return;
  FlxBayUP_csm_csus_MCMC* csm_csus = dynamic_cast<FlxBayUP_csm_csus_MCMC*>(&csm);
  if (csm_csus==NULL) throw FlxException_NotImplemented("flxBayUp_adaptive_ctrl_opti_jump::requires_adptv_step_02");
  limit_acr = acr_min->cast2positive_or0();        // try 'desperately' to increase the acceptance rate if acr smaller than this value
  const tdouble escale = esjd_scale->cast2positive();
  const tdouble pweight_val = get_pweight();
  
  // check if the spread needs to be adapted
    if (skip_adaptive_step(acr)) return;
  
  // initialize before starting the optimization
    compute_seed_weights();
    
  // perform the optimization
    if (sd_prev_min<=ZERO) {
      sd_prev_min = 0.8;
    }
    tdouble x_min = rv_InvPhi(sd_prev_min);
    tdouble x_lb = x_min - ONE/2;
    tdouble x_ub = x_min + ONE/2;
    ipds.reset();
    try {
      // find the spread that optimizes the ESJD
        flx_optim(x_lb,x_ub,x_min,&perfFun1D_flxBayUp_adaptive_ctrl_opti_jump,this,true,true,100,10,1e-2,1e-2,NULL);
    } catch (FlxException_Crude& e) {
      #if FLX_DEBUG
        GlobalVar.slogcout(3) << std::endl << e.what();
      #endif
      throw;
    } catch (FlxException& e) { 
      #if FLX_DEBUG
        GlobalVar.slogcout(3) << std::endl << e.what();
      #endif
      try {
        x_min = rv_InvPhi_noAlert(ipds.find_1st_x_before_xs_smaller_than_f(sd_prev_min,ZERO,false));
      } catch (FlxException &e) {
        #if FLX_DEBUG
          GlobalVar.slogcout(3) << std::endl << e.what();
        #endif
      }
      // in case of an error, just ignore the optimization
    }
    tdouble sd_proposed = rv_Phi(x_min);
    // find the spread that belongs to a reduced ESJD
      const tdouble e_opt = ipds.interpolate(sd_proposed);
      try {
        sd_proposed = ipds.find_1st_x_before_xs_smaller_than_f(sd_proposed,e_opt*escale);
      } catch (FlxException_Crude& e) {
        #if FLX_DEBUG
          GlobalVar.slogcout(3) << std::endl << e.what();
        #endif
        throw;
      } catch (FlxException& e) {
        #if FLX_DEBUG
          GlobalVar.slogcout(3) << std::endl << e.what();
        #endif
        // in case of an error, just ignore it ...
      }
    // make sure the spread is not too close to one (for numerical reasons)
      if (sd_proposed>rv_Phi(3.)) sd_proposed = rv_Phi(3);
    // remember employed spread
      {
        ++spread_hist_N;
        // move existing entries one down
          if (spread_hist_N>1) {
            const tuint N = (spread_hist_N<spread_hist_N_max)?spread_hist_N:spread_hist_N_max; 
            memmove(spread_hist_vec.get_tmp_vptr()+1,spread_hist_vec.get_tmp_vptr(),(N-1)*sizeof(tdouble));
          }
        spread_hist_vec[0] = sd_proposed;
      }
    // decide how to modify the spread based on some heuristic rules
      bool modify_spread = true;
      if ( sd_proposed>=sd_prev_min ) {
        const tdouble acr_oa = compute_overall_acr();
        if (spread_hist_N<=1 || smpl_N<=100 || acr_oa < limit_acr) {
          modify_spread=false; 
          if (acr<limit_acr) {
            sd_prev_min *= 0.75;
          }
        }
      }
      if (modify_spread) {
        sd_prev_min = (ONE-pweight_val)*sd_prev_min+pweight_val*sd_proposed;
      }
    // actually set the new spread
        csm_csus->set_cur_spread(sd_prev_min);
    #if FLX_DEBUG
      if (sd_proposed>rv_Phi(x_min)) throw FlxException_Crude("flxBayUp_adaptive_ctrl_opti_jump::requires_adptv_step_03");
    #endif
}

const tdouble flxBayUp_adaptive_ctrl_opti_jump::perfFun(const tdouble sd)
{
  const tdouble WMEAN_PUNISH = 0.05;        // punish weights if mean is smaller than this value
  const tdouble ACRDD = 0.05;
  // evaluate the weights
    const tdouble* c = smpl_list;
    for (tuint i=0;i<smpl_N;++i) {
      twl[i] = exp(proposal_pdf_ln(c,sd)-swl[i]);
      #if FLX_DEBUG
        if (twl[i]!=twl[i]) {
          throw FlxException_Crude("flxBayUp_adaptive_ctrl_opti_jump::perfFun_01");
        }
      #endif
      c += acopti_jump_dim;
    }
  // regression estimate for the weights [MATH_IS_Hesterberg_1995]
    flxVec wvec(twl,smpl_N);
    tdouble b_reg;
    const tdouble w_mean = wvec.get_Mean();
    tdouble spot = ONE;
      if (w_mean<GlobalVar.TOL()) {
        ipds.append(sd,-1e2);
        return -1e2;
      } else if (fabs(w_mean-ONE)<GlobalVar.TOL()) {
        b_reg = ZERO;
      } else {
        const tdouble w_sd = wvec.get_sd(w_mean);
        if (w_mean<WMEAN_PUNISH) {
          spot = pow2(w_mean/WMEAN_PUNISH);
          b_reg = (ONE-w_mean)/pow(w_sd,spot*2)*spot;
        } else {
          b_reg = (ONE-w_mean)/pow2(w_sd);
        }
      }
    for (tuint i=0;i<smpl_N;++i) {
      twl[i] *= (ONE+b_reg*(twl[i]-w_mean))/smpl_N;
      if (twl[i]<ZERO) twl[i] = ZERO;
    }
    const tdouble weight_sum = wvec.get_sum();
  // Perform the Importance Sampling
    c = smpl_list;
    tdouble esjd = ZERO;
    tdouble aacr = ZERO;
    for (tuint i=0;i<smpl_N;++i) {
      const tdouble d = c[4]*twl[i];
      aacr += d;
      esjd += c[5]*d;
      c += acopti_jump_dim;
    }
    aacr /= weight_sum;
    esjd /= weight_sum;
    #if FLX_DEBUG
      if (aacr>ONE) throw FlxException_Crude("flxBayUp_adaptive_ctrl_opti_jump::perfFun_02");
    #endif
  // check acceptance rate
    if (aacr<limit_acr) {
      const tdouble dd = ACRDD;
      const tdouble d =(limit_acr-aacr)/dd;
      esjd *= (ONE-rv_Phi(sd)*d);
    }
    if (esjd>ZERO) esjd*= spot;        // punish ESJD for weights < WMEAN_PUNISH
  ipds.append(sd,esjd);
  return esjd;
}

void flxBayUp_adaptive_ctrl_opti_jump::append_smpl(const flxVec& lastS)
{
  if (smpl_list==NULL) {
    throw FlxException_Crude("flxBayUp_adaptive_ctrl_opti_jump::append_smpl");
  }
  const size_t pos = acopti_jump_dim*smpl_i;
  flxVec sample(smpl_list+pos,acopti_jump_dim);  
  sample = lastS;
  ++smpl_i;
  if (smpl_N<smpl_i) smpl_N = smpl_i;
  if (smpl_i==smpl_Nmax) smpl_i=0;
}

void flxBayUp_adaptive_ctrl_opti_jump::print_info(std::ostream& sout) const
{
  sout << "  adaptive factors:             acr_min=" << acr_min->write() << ", escl=" << esjd_scale->write() << ", Nmax=" << Nmax_fun->write() << std::endl;
  flxBayUp_adaptive_ctrl_base::print_info(sout);
}

void flxBayUp_adaptive_ctrl_opti_jump::write_adaptive_info(std::ostream& sout)
{
  if (spread_hist_N==0) return;
  sout << "Na=" << spread_hist_N << " ";
}


flxBayUp::flxBayUp(const std::string& nameID, const tdouble& scaleconst, const tdouble& cStartV, const std::string& parentsetstr)
: is_subsetRel(false), scaleconst(log(scaleconst)),cStart(cStartV),cStart_init(cStart),pa_maxL(ONE),N_LklVec(0),global_Likelihood(NULL),global_Likelihood_is_log(false),methCat(UNDEFINED),LklSet(NULL),p_rv(NULL),RndBox(NULL),nameID(nameID),
  updater(data->RndCreator)
{
  // resolve dependencies of the parents
    const std::vector<std::string> set_str_vec = parse_strseq_as_vec(parentsetstr);
    RBRV_constructor::find_dependent_sets(set_str_vec,setvec,data->rbrv_box);
}

flxBayUp::flxBayUp(const std::string& parentsetstr)
: is_subsetRel(true), scaleconst(ZERO),cStart(ZERO),cStart_init(cStart),pa_maxL(ONE),N_LklVec(0),global_Likelihood(NULL),global_Likelihood_is_log(false),methCat(UNDEFINED),LklSet(NULL),p_rv(NULL),RndBox(NULL),nameID("dummy_for_sus"),
  updater(data->RndCreator)
{
  // resolve dependencies of the parents
    const std::vector<std::string> set_str_vec = parse_strseq_as_vec(parentsetstr);
    RBRV_constructor::find_dependent_sets(set_str_vec,setvec,data->rbrv_box);
}

flxBayUp::~flxBayUp()
{
  if (global_Likelihood) delete global_Likelihood;
  if (RndBox) delete RndBox;
  if (LklSet) { 
    if (is_subsetRel) {
      delete LklSet;
    } else {
      //     delete LklSet; (done by RBRV-management)
    }
  } else {
    for (std::size_t i=0;i<LklVec.size();++i) delete LklVec[i];
  }
}

void flxBayUp::set_TMCMC()
{
  if (methCat == TMCMC) return;
  if (LklSet || methCat!=UNDEFINED) {
    throw FlxException_Crude("flxBayUp::set_TMCMC");
  }
  methCat = TMCMC;
}

void flxBayUp::freeze()
{
  if (LklSet) return;        // object was already frozen
  if (methCat==UNDEFINED) methCat=BUS;
  // check whether there is a problem
    if (global_Likelihood==NULL&&N_LklVec==0) {
      std::ostringstream ssV;
      switch (methCat) {
        case (BUS):
          ssV << "At least one (global or local) likelihood function must be defined in '" << nameID << "'.";
          break;
        case (ABC):
          ssV << "A metric is required but not defined in '" << nameID << "'.";
          break;
        case (RA):
          ssV << "A limit-state function is required but not defined in '" << nameID << "'.";
          break;
        default:
          throw FlxException_Crude("flxBayUp::freeze_1");
      };
      throw FlxException("flxBayUp::freeze_2", ssV.str() ); 
    }
  // allocate LklSet
    switch (methCat) {
      case (BUS):
      {
        RBRV_entry** entries = new RBRV_entry*[N_LklVec+1];
        for (tuint i=0;i<N_LklVec;++i) entries[i] = LklVec[i];
        p_rv = new RBRV_entry_RV_stdN("p",0);
        entries[N_LklVec] = p_rv;
        LklSet = new RBRV_set(true,1,nameID + "::local",false,N_LklVec+1,entries,0,NULL,true);
        break;
      }
      case (TMCMC):
      {
        RBRV_entry** entries = new RBRV_entry*[N_LklVec];
        for (tuint i=0;i<N_LklVec;++i) entries[i] = LklVec[i];
        LklSet = new RBRV_set(true,0,nameID + "::local",false,N_LklVec,entries,0,NULL,true);
        break;
      }
      case (ABC):
      case (RA):
        if (global_Likelihood==NULL || N_LklVec>0) {
          throw FlxException_Crude("flxBayUp::freeze_3");
        }
        if (global_Likelihood_is_log && methCat==RA ) throw FlxException_Crude("flxBayUp::freeze_4");
        // define an empty set
          LklSet = new RBRV_set(true,0,nameID + "::local",is_subsetRel,0,new RBRV_entry*[0],0,NULL,true);
        break;
      default:
        throw FlxException_Crude("flxBayUp::freeze_4");
    };
  setvec.push_back(LklSet);
  if (is_subsetRel==false) {
    data->rbrv_box.register_set(LklSet);
  }
  RndBox = new RBRV_constructor(setvec);
}

void flxBayUp::add_localLkl(RBRV_entry* lklEntry)
{
  if (LklSet) {
    std::ostringstream ssV;
    ssV << "It is not possible anymore to add a lokal Likelihood functions to '" << nameID << "'.";
    throw FlxException_NeglectInInteractive("flxBayUp::add_localLkl_1", ssV.str() ); 
  }
  #if FLX_DEBUG
    if (lklEntry==NULL) throw FlxException_Crude("flxBayUp::add_localLkl_2");
    // make sure the name is fine
      std::ostringstream ssV;
      ssV << nameID << "::" << N_LklVec;
      if (lklEntry->name != ssV.str() ) throw FlxException_Crude("flxBayUp::add_localLkl_3");
  #endif
  data->rbrv_box.register_entry(lklEntry);
  LklVec.push_back(lklEntry);
  ++N_LklVec;
}

void flxBayUp::add_localLkl(flxBayUp_uncertobsv_set* ts)
{
  if (LklSet) {
    std::ostringstream ssV;
    ssV << "It is not possible anymore to add a lokal Likelihood functions to '" << nameID << "'.";
    throw FlxException_NeglectInInteractive("flxBayUp::add_localLkl_90", ssV.str() ); 
  }
  // create the name of the entry
    std::ostringstream ssV;
    ssV << nameID << "::" << N_LklVec;
    const std::string nameE = ssV.str();
  RBRV_entry_ref_log* lklEntry = new RBRV_entry_ref_log(nameE,ts->get_logLike());
  try {
    add_localLkl(lklEntry);
  } catch (FlxException &e) {
    FLXMSG("flxBayUp::add_localLkl_95",1);
    delete lklEntry;
    throw;
  }
  setvec.push_back(ts);
}

void flxBayUp::set_globalLkl(FlxFunction& global_LikelihoodV, const bool is_log, const flxBayUp::MethCategory methCatV)
{
  if (global_Likelihood) {
    std::ostringstream ssV;
    ssV << "A 'global likelihood'/'metric'/'limit-state' function has already been defined for '" << nameID << "'.";
    throw FlxException_NeglectInInteractive("flxBayUp::set_globalLkl", ssV.str() ); 
  }
  global_Likelihood = new FlxFunction(global_LikelihoodV);
  global_Likelihood_is_log = is_log;
  if (methCatV!=flxBayUp::UNDEFINED) {
    methCat = methCatV;
    freeze();
  }
}

const tdouble flxBayUp::eval_Likelihood()
{
  #if FLX_DEBUG
    if (LklSet==NULL) throw FlxException_Crude("flxBayUp::eval_Likelihood_1");
  #endif
  tdouble res;
  switch (methCat) {
    case (BUS):
    case (TMCMC):
    {
      if (global_Likelihood) {
        if (global_Likelihood_is_log) res = global_Likelihood->calc();
        else res = log(global_Likelihood->cast2positive_or0());
      } else {
        #if FLX_DEBUG
          if (N_LklVec==0) {
            std::ostringstream ssV;
            ssV << "There is no likelihood function defined.";
            throw FlxException_NeglectInInteractive("flxBayUp::eval_Likelihood_2", ssV.str() ); 
          }
        #endif
        res = scaleconst;
        for (tuint i=0;i<N_LklVec;++i) res += LklVec[i]->get_value_log();
      }
      break;
    }
    case (ABC):
      if (global_Likelihood_is_log) res = global_Likelihood->calc();
      else res = log(global_Likelihood->cast2positive_or0());
      break;
    default:
      throw FlxException_Crude("flxBayUp::eval_Likelihood_3");
  };
  if (res!=res) {        // check for NaN
    throw FlxException("flxBayUp::eval_Likelihood_3", "Likelihood value is 'NaN'." );
  }
  return res;
}

const tdouble flxBayUp::eval_RAlsf()
{
  #if FLX_DEBUG
    if (LklSet==NULL) throw FlxException_Crude("flxBayUp::eval_RAlsf_1");
    if (methCat!=RA) throw FlxException_Crude("flxBayUp::eval_RAlsf_2");
  #endif
  return global_Likelihood->calc();
}

const tdouble flxBayUp::get_p()
{
  #if FLX_DEBUG
    if (methCat!=BUS) throw FlxException_Crude("flxBayUp::get_p");
  #endif
  return p_rv->get_value();
}

RBRV_constructor& flxBayUp::get_RndBox()
{
  if (RndBox==NULL) freeze();
  return *RndBox;
}

FlxBayUp_Update_List::FlxBayUp_Update_List(flxBayUp& parent, const tuint Nc, const tuint Ncl, const tuint Ns_finalV, const tuint Nburn, const FlxBayUp_Update_List::randomizeTec randomizeID, flxBayUp_adaptive_ctrl_base* adpt_ctrl, const tuint max_runs, const bool use_cStart, const FlxBayUp_Update_List::MethType meth_id, const bool log_LSF, const bool find_multiple)
: parent(parent), RndBoxRB(parent.get_RndBox()), N_RV(RndBoxRB.get_NRV()), N_OX(RndBoxRB.get_NOX()), Nc(Nc), Ncl(Ncl), Ns_final((Nc*Ncl>Ns_finalV)?(Nc*Ncl):Ns_finalV), Nburn(Nburn), meth_id(meth_id), c(parent.get_cStart()),
max_L(use_cStart?parent.get_cStart_init():c), s_thr(99.), fullList(false), y_list(new tdouble[N_RV*Ns_final]), yh_vec(N_RV), x_list(new tdouble[N_OX*Ns_final]), xh_vec(N_OX), p_list(NULL), s_list(new tdouble[Ns_final]), sh_list(new tdouble[Ns_final]), L_list(NULL), multpleVec_help(find_multiple?(new bool[Ns_final]):NULL), i_list(new int[Ns_final]), seed_idx(NULL), chain_length(NULL),
max_runs(max_runs),log_LSF(log_LSF),oldC_N(0),oldC_list(NULL),oldS_list(NULL),
adpt_ctrl(adpt_ctrl),randomizeID(randomizeID),finalized(false),curID(0),curSID(0)
{
  for (tuint i=0;i<Ns_final;++i) i_list[i] = 0;
  switch (meth_id) {
    case (FlxBayUp_Update_List::BUS):
    case (FlxBayUp_Update_List::UBUS):
    case (FlxBayUp_Update_List::BUST):
      oldC_list = new tdouble[max_runs];
      oldS_list = new tdouble[max_runs];
      p_list = new tdouble[Ns_final];
      L_list = new tdouble[Ns_final];
    case (FlxBayUp_Update_List::ABCSUBSIM):
    case (FlxBayUp_Update_List::RASUBSIM):
      for (tuint i=Nc*Ncl;i<Ns_final;++i) i_list[i] = -1;
      seed_idx = new tuint[Ns_final];
      for (tuint i=0;i<Ns_final;++i) seed_idx[i] = Ns_final+1;
      chain_length = new tuint[Ns_final];
      for (tuint i=0;i<Ns_final;++i) chain_length[i] = 0;
      break;
    case (FlxBayUp_Update_List::TMCMC):
      L_list = new tdouble[Ns_final];
      seed_idx = new tuint[Ns_final];
      break;
    case (FlxBayUp_Update_List::MHRS):
      seed_idx = new tuint[Ns_final];
    case (FlxBayUp_Update_List::RS):
    case (FlxBayUp_Update_List::LS):
      L_list = new tdouble[Ns_final];
    case (FlxBayUp_Update_List::ABCRS):
    case (FlxBayUp_Update_List::RAMCI):
      break;
  };
}

FlxBayUp_Update_List::~FlxBayUp_Update_List()
{
  if (adpt_ctrl) delete adpt_ctrl;
  if (oldS_list) delete [] oldS_list;
  if (oldC_list) delete [] oldC_list;
  delete [] y_list;
  delete [] x_list;
  if (p_list) delete [] p_list;
  delete [] s_list;
  if (sh_list) delete [] sh_list;
  if (L_list) delete [] L_list;
  delete [] i_list;
  if (seed_idx) delete [] seed_idx;
  if (chain_length) delete [] chain_length;
  if (multpleVec_help) delete [] multpleVec_help;
}

const tdouble FlxBayUp_Update_List::eval_LSF(const tdouble p, const tdouble L) const
{
  switch (meth_id) {
    case (BUS):
    {
      const tdouble r = p-rv_InvPhi_noAlert(exp(L-c));
      if (oldC_N>0) {
        if (r > s_thr) return r;
        for (tuint i=0;i<oldC_N;++i) {
          if ((p-rv_InvPhi_noAlert(exp(L-oldC_list[i]))) > oldS_list[i]) {
            return s_thr+54.4e6;
          }
        }
      }
      return r;
    }
    case (UBUS):
    {
      const tdouble r = log(rv_Phi(p))-L+c;
      #if FLX_DEBUG
        if (oldC_N>0) throw FlxException_Crude("FlxBayUp_Update_List::eval_LSF_1");      
      #endif
      return r;
    }
    case (BUST):
    {
      const tdouble r = p-rv_InvPhi_noAlert(exp(s_thr*(L-c)));
      if (oldC_N>0) {
        if (r > ZERO) return r;
        for (tuint i=0;i<oldC_N;++i) {
          if ((p-rv_InvPhi_noAlert(exp(oldS_list[i]*(L-oldC_list[i])))) > ZERO) {
            return 54.4e6;
          }
        }
      }
      return r;
    }
    case ABCSUBSIM:
    case RASUBSIM:
    case MHRS:
    case RS:
    case ABCRS:
    case RAMCI:
    case LS:
      return p-rv_InvPhi_noAlert(exp(L-c));
    default:
    {
      throw FlxException_Crude("FlxBayUp_Update_List::eval_LSF_2");
    }
  };
//   return p-rv_InvPhi_noAlert(exp(L-c));
//   return p-exp(L-c);
}

const int FlxBayUp_Update_List::get_cur_i_list() const
{
  if (curID>=Ns_final) return -1;
  return i_list[curID];
}

tuint& FlxBayUp_Update_List::get_curSID_ref()
{
  return seed_idx[curSID];
}

tdouble* FlxBayUp_Update_List::get_seed_y_list()
{
  return &(y_list[seed_idx[curSID]*N_RV]);
}

#if FLX_DEBUG

  const tdouble FlxBayUp_Update_List::get_seed_L_list() const
  {
    return L_list[seed_idx[curSID]];
  }

  const tdouble FlxBayUp_Update_List::get_seed_p_list() const
  {
    return p_list[seed_idx[curSID]];
  }

  const tdouble FlxBayUp_Update_List::get_seed_s_list() const
  {
    return s_list[seed_idx[curSID]];
  }

#endif

const tdouble* FlxBayUp_Update_List::get_cur_y_list()
{
  return &(y_list[curID*N_RV]);
}

const tdouble* FlxBayUp_Update_List::get_cur_x_list()
{
  return &(x_list[curID*N_OX]);
}

const tdouble FlxBayUp_Update_List::get_cur_L()
{
  return L_list[curID];
}

void FlxBayUp_Update_List::set_next(const bool init)
{
  if (init) {
    if (get_cur_i_list()==0) return;
  } else {
    i_list[curID] = 1;
  }
  while (true) {
    ++curID;
    if (get_cur_i_list()<=0) return;
  }
}

void FlxBayUp_Update_List::write_smpl(const tuint ID2write, std::ofstream& os_smpl)
{
  os_smpl << ID2write << std::endl;
  // write some other stuff
    switch (meth_id) {
      // ***************************************************************************************************************
      case (FlxBayUp_Update_List::BUS):
      case (FlxBayUp_Update_List::UBUS):      
      case (FlxBayUp_Update_List::BUST):
      {
        os_smpl << "  ";
        os_smpl << GlobalVar.D2S_totalPrec(L_list[ID2write]) << ", ";
        os_smpl << GlobalVar.D2S_totalPrec(p_list[ID2write]) << ", ";
        os_smpl << GlobalVar.D2S_totalPrec(s_list[ID2write]);
        os_smpl << std::endl;
        break;
      }
      // ***************************************************************************************************************
      case (FlxBayUp_Update_List::ABCSUBSIM):
      {
        os_smpl << "  ";
        os_smpl << GlobalVar.D2S_totalPrec(s_list[ID2write]);
        os_smpl << std::endl;
        break;
      }
      // ***************************************************************************************************************
      case (FlxBayUp_Update_List::RASUBSIM):
      {
        os_smpl << "  ";
        os_smpl << GlobalVar.D2S_totalPrec(s_list[ID2write]);
        os_smpl << std::endl;
        break;
      }
      // ***************************************************************************************************************
      case (FlxBayUp_Update_List::TMCMC):
      {
        os_smpl << "  ";
        os_smpl << GlobalVar.D2S_totalPrec(L_list[ID2write]);
        os_smpl << std::endl;
        break;
      }
      // ***************************************************************************************************************
      case (FlxBayUp_Update_List::RS):
      case (FlxBayUp_Update_List::ABCRS):
      case (FlxBayUp_Update_List::RAMCI):
      {
        break;
      }
      default:
        throw FlxException_Crude("FlxBayUp_Update_List::write_smpl_1");
    };
  // write the samples
    const flxVec sy(&(y_list[ID2write*N_RV]),N_RV);        // x-vector to write
      os_smpl << "  ";
      flxVec_totalPrec_plot(os_smpl,sy);
      os_smpl << std::endl;
    const flxVec sx(&(x_list[ID2write*N_OX]),N_OX);        // y-vector to write
      os_smpl << "  ";
      flxVec_totalPrec_plot(os_smpl,sx);
      os_smpl << std::endl;
}

const bool FlxBayUp_Update_List::insert_entry(const bool is_first, bool rejectIt, const bool is_posterior, const bool just_accpt_chk, std::ofstream* os_smpl, const tdouble tLik, tdouble* acr)
{
  #if FLX_DEBUG
    if (curID>=Ns_final) {
      throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_0");
    }
    if (!is_posterior) {
      if (i_list[curID]!=0) throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_1");
    }
    if (is_first && rejectIt) throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_2");
    if (is_first && is_posterior) throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_3");
    if (rejectIt && is_posterior) throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_4");
    if (meth_id!=FlxBayUp_Update_List::TMCMC && is_posterior && is_gt_zero()==false ) throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_5");
    if (just_accpt_chk && is_posterior==false ) throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_6");
  #endif
  flxVec sy(&(y_list[curID*N_RV]),N_RV);        // current x-vector
  flxVec sx(&(x_list[curID*N_OX]),N_OX);        // current y-vector
  switch (meth_id) {
    // ***************************************************************************************************************
    case (FlxBayUp_Update_List::BUS):
    case (FlxBayUp_Update_List::UBUS):
    case (FlxBayUp_Update_List::BUST):
    {
      i_list[curID] = 1;        // mark entry as occupied
      tdouble tL = ZERO;
      if (!rejectIt) {
        try {
          tL = parent.eval_Likelihood();
        } catch (FlxException_NeglectInInteractive& e) {
          FLXMSG("FlxBayUp_Update_List::insert_entry_BUS_1",1);
          rejectIt = true;
          if (is_first) throw;
        }
      }
      const tdouble p_here = parent.get_p();
      const tdouble ts = (is_first||rejectIt)?1.9876:eval_LSF(p_here,tL);
      if (log_LSF && rejectIt==false) {
        GlobalVar.slog(4) << std::endl << GlobalVar.Double2String(ts,false,-1,0) << GlobalVar.Double2String(p_here,false,-1,0) << GlobalVar.Double2String(tL,false,-1,0);
      }
      bool b;
      if (acr && rejectIt==false) {
        if (meth_id==FlxBayUp_Update_List::BUST) {
          if (ts!=54.4e6) {
            tdouble tt1 = exp(s_thr*(tL-c));
            if (tt1>ONE) tt1=ONE;
            *acr += tt1;
          }
        } else {
          if (ts!=s_thr+54.4e6) {
            if (tL-c>=ZERO) {
              *acr += ONE;
            } else {
              *acr += rv_Phi(rv_InvPhi_noAlert(exp(tL-c))+s_thr);
            }
          }
        }
      }
      if ((is_first||((meth_id==BUST&&ts<=ZERO) || ((meth_id!=BUST)&&ts<=s_thr)))&&(rejectIt==false)) {        // accept current proposal
        if (just_accpt_chk) {
          if (tL>c) {        // update max. observed likelihood value
            update_c_posterior(false,tL);
          }
          return true;
        }
        L_list[curID] = tL;
        p_list[curID] = p_here;
        s_list[curID] = ts;
        RndBoxRB.get_y_Vec(sy.get_tmp_vptr());
        RndBoxRB.get_x_Vec(sx.get_tmp_vptr());
        b = true;
      } else {                        // reject it ... take the last one
        if (is_first) throw FlxException_NeglectInInteractive("FlxBayUp_Update_List::insert_entry_BUS_2","Sample had to be rejected.");
        if (just_accpt_chk) return false;
        if (!is_posterior) {
          const tuint seed_pos = seed_idx[curSID];
          L_list[curID] = L_list[seed_pos];
          p_list[curID] = p_list[seed_pos];
          s_list[curID] = s_list[seed_pos];
          sy = flxVec(&(y_list[seed_pos*N_RV]),N_RV);
          sx = flxVec(&(x_list[seed_pos*N_OX]),N_OX);
        }
        b = false;
      }
      if (os_smpl) write_smpl(curID,*os_smpl);
      if (!rejectIt) {
        if (is_posterior) {
          // check if scaling constant changed
          if (tL>c) {        // we have to adopt cStart
            #if FLX_DEBUG
              if (!b) throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_BUS_3");
            #endif
            // split up the current seed (i.e. our sample)
              b = !update_c_posterior();
          }
        } else {
          if (tL>max_L) max_L = tL;
        }
      }
      if (!is_posterior) {
        seed_idx[curSID] = curID;
        set_next();
      }
      return b;
    }
    // ***************************************************************************************************************
    case (FlxBayUp_Update_List::ABCSUBSIM):
    {
      tdouble tM = ONE;
      if (!rejectIt) {
        try {
          tM = parent.eval_Likelihood();
          if (log_LSF) {
            GlobalVar.slog(4) << std::endl << GlobalVar.Double2String(tM,false,-1,0);
          }
        } catch (FlxException_NeglectInInteractive& e) {
          FLXMSG("FlxBayUp_Update_List::insert_entry_ABC_1",1);
          rejectIt = true;
          if (is_first) throw;
        }
      }
      bool b;
      if ( is_first ||( (tM<=s_thr)&&(rejectIt==false) ) ) {        // accept current proposal
        if (just_accpt_chk) return true;
        s_list[curID] = tM;
        RndBoxRB.get_y_Vec(sy.get_tmp_vptr());
        RndBoxRB.get_x_Vec(sx.get_tmp_vptr());
        b = true;
      } else {                        // reject it ... take the last one
        if (is_first) throw FlxException_NeglectInInteractive("FlxBayUp_Update_List::insert_entry_ABC_2","Sample had to be rejected.");
        if (just_accpt_chk) return false;
        if (!is_posterior) {
          const tuint seed_pos = seed_idx[curSID];
          s_list[curID] = s_list[seed_pos];
          sy = flxVec(&(y_list[seed_pos*N_RV]),N_RV);
          sx = flxVec(&(x_list[seed_pos*N_OX]),N_OX);
        }
        b = false;
      }
      if (os_smpl) write_smpl(curID,*os_smpl);
      i_list[curID] = 1;        // mark entry as occupied
      if (!is_posterior) {
        seed_idx[curSID] = curID;
        set_next();
      }
      return b;
    }
    // ***************************************************************************************************************
    case (FlxBayUp_Update_List::RASUBSIM):
    {
      tdouble tLSF = ONE;
      if (!rejectIt) {
        try {
          tLSF = parent.eval_RAlsf();
          if (log_LSF) {
            GlobalVar.slog(4) << std::endl << GlobalVar.Double2String(tLSF,false,-1,0);
          }
        } catch (FlxException_NeglectInInteractive& e) {
          FLXMSG("FlxBayUp_Update_List::insert_entry_RA_1",1);
          rejectIt = true;
          if (is_first) throw;
        }
      }
      bool b;
      if ( is_first ||( (tLSF<=s_thr)&&(rejectIt==false) ) ) {        // accept current proposal
        if (just_accpt_chk) return true;
        s_list[curID] = tLSF;
        RndBoxRB.get_y_Vec(sy.get_tmp_vptr());
        RndBoxRB.get_x_Vec(sx.get_tmp_vptr());
        b = true;
      } else {                        // reject it ... take the last one
        if (is_first) throw FlxException_NeglectInInteractive("FlxBayUp_Update_List::insert_entry_RA_2","Sample had to be rejected.");
        if (just_accpt_chk) return false;
        if (!is_posterior) {
          const tuint seed_pos = seed_idx[curSID];
          s_list[curID] = s_list[seed_pos];
          sy = flxVec(&(y_list[seed_pos*N_RV]),N_RV);
          sx = flxVec(&(x_list[seed_pos*N_OX]),N_OX);
        }
        b = false;
      }
      if (os_smpl) write_smpl(curID,*os_smpl);
      i_list[curID] = 1;        // mark entry as occupied
      if (!is_posterior) {
        seed_idx[curSID] = curID;
        set_next();
      }
      return b;
    }
    // ***************************************************************************************************************
    case (FlxBayUp_Update_List::TMCMC):
    {
      if (!is_posterior || just_accpt_chk) throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_TMCMC_1");
      i_list[curID] = 1;        // mark entry as occupied
      // evaluate likelihood
        tdouble tL = ZERO;
        if (!rejectIt) {
          try {
            tL = parent.eval_Likelihood();
          } catch (FlxException_NeglectInInteractive& e) {
            FLXMSG("FlxBayUp_Update_List::insert_entry_BUS_1",1);
            rejectIt = true;
            if (is_first) throw;
          }
        }
      // evaluate likelihood ratio
        if (!rejectIt) {
          #ifndef FLX_CV_1
            tdouble y_prop[N_RV];
          #else
            tdouble* y_prop = new tdouble[N_RV];
          #endif
          RndBoxRB.get_y_Vec(y_prop);
          double alpha = ONE;
          for (tuint j=0;j<N_RV;++j) {
            alpha *= rv_phi(y_prop[j])/rv_phi(sy[j]);
          }
          alpha *= exp(tL-L_list[curID]);
          rejectIt = (RndCreator->gen_smp_uniform()>alpha);
          #ifdef FLX_CV_1
            delete [] y_prop;
          #endif
        }
      bool b;
      if (!rejectIt) {        // accept current proposal
        L_list[curID] = tL;
        RndBoxRB.get_y_Vec(sy.get_tmp_vptr());
        RndBoxRB.get_x_Vec(sx.get_tmp_vptr());
        b = true;
      } else {                        // reject it
        b = false;
      }
      if (os_smpl) write_smpl(curID,*os_smpl);
      return b;
    }
    // ***************************************************************************************************************
    case (FlxBayUp_Update_List::MHRS):
    case (FlxBayUp_Update_List::RS):
    case (FlxBayUp_Update_List::ABCRS):
    case (FlxBayUp_Update_List::RAMCI):
    case (FlxBayUp_Update_List::LS):
    {
      if (just_accpt_chk) {
        if (meth_id==FlxBayUp_Update_List::RS) {
          tdouble tL = ZERO;
          if (!rejectIt) {
            try {
              tL = parent.eval_Likelihood();
            } catch (FlxException_NeglectInInteractive& e) {
              FLXMSG("FlxBayUp_Update_List::insert_entry_RS_1",1);
              rejectIt = true;
            }
          }
          if (rejectIt) return false;
          const tdouble p_here = parent.get_p();
          const tdouble ts = eval_LSF(p_here,tL);
          return ts<=ZERO;
        } else {
          throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_RS_2");
        }
      }
      if (!rejectIt) {
        if (tLik>c) {
          if (meth_id==FlxBayUp_Update_List::RS) {
            std::ostringstream ssV;
            ssV << "cStart is not large enough: ln(cStart)=" << GlobalVar.Double2String(c) << "; "
                << "observed log-likelihood: " << GlobalVar.Double2String(tLik);
            throw FlxException("FlxBayUp_Update_List::insert_entry_RS_3", ssV.str() ); 
          } else if (meth_id==FlxBayUp_Update_List::MHRS) {
            if (tLik>max_L) max_L = tLik;
          }
        }
        RndBoxRB.get_y_Vec(sy.get_tmp_vptr());
        RndBoxRB.get_x_Vec(sx.get_tmp_vptr());
        if (meth_id==FlxBayUp_Update_List::MHRS||meth_id==FlxBayUp_Update_List::RS||meth_id==FlxBayUp_Update_List::LS) {
          L_list[curID] = tLik;
        }
      }
      if (os_smpl) write_smpl(curID,*os_smpl);
      if (!is_posterior) set_next();
      return !rejectIt;
    }
  };
  throw FlxException_Crude("FlxBayUp_Update_List::insert_entry_3");
}

const tdouble FlxBayUp_Update_List::expectation_likelihood() const
{
  if (finalized) {
    pdouble E_L;
    for (tuint i=0;i<Ns_final;++i) {
      E_L += L_list[i];
    }
    E_L /= Ns_final;
    return E_L.cast2double();
  } else {
    pdouble E_L;
    tuint i;
    tuint Ntot = 0;
    for (i=0;i<Ns_final;++i) {
      if (i_list[i]<=0) {
        if (i_list[i]==-1) break;
        else continue;
      }
      E_L += L_list[i];
      ++Ntot;
    }
    E_L /= Ntot;
    return E_L.cast2double();
  }
}

const tdouble FlxBayUp_Update_List::propose_qlnL(const tdouble pa_maxL)
{
  if (pa_maxL>=ONE) return max_L;
  // define some variables
    const tuint Ns = Nc*Ncl;
    const tuint N_in_list = (fullList?Ns_final:Ns);
  // sort Likelihood vector (temporary)
    {
      flxVec Lh_vec(sh_list,N_in_list);
      flxVec L_vec(L_list,N_in_list);
      Lh_vec = L_vec;
    }
    std::sort(sh_list,sh_list+N_in_list);
  // retrieve the percentile
    return flx_percentile(sh_list,N_in_list,pa_maxL);
}

void FlxBayUp_Update_List::swap_smpls(const tuint id1, const tuint id2)
{
  if (id1==id2) return;
  // y-list
  {
    flxVec y1(y_list+N_RV*id1,N_RV);
    flxVec y2(y_list+N_RV*id2,N_RV);
    yh_vec = y1;
    y1 = y2;
    y2 = yh_vec;
  }
  // x-list
  {
    flxVec x1(x_list+N_OX*id1,N_OX);
    flxVec x2(x_list+N_OX*id2,N_OX);
    xh_vec = x1;
    x1 = x2;
    x2 = xh_vec;
  }
  // p-list
    tdouble pt = p_list[id1];
    p_list[id1] = p_list[id2];
    p_list[id2] = pt;
  // s-list
    pt = s_list[id1];
    s_list[id1] = s_list[id2];
    s_list[id2] = pt;
  // L-list
    pt = L_list[id1];
    L_list[id1] = L_list[id2];
    L_list[id2] = pt;
}

void FlxBayUp_Update_List::update_LSF_vals(const tuint N_in_list, const bool do_check)
{
  for (tuint i=0;i<N_in_list;++i) {
    if (i_list[i]<=0) break;
    #if FLX_DEBUG
      const tdouble s_prev = s_list[i];
    #endif
    s_list[i] = eval_LSF(p_list[i],L_list[i]);
    #if FLX_DEBUG
      if (do_check && ((s_prev-s_thr)*(s_list[i]-s_thr)<ZERO)) {
        throw FlxException_Crude("FlxBayUp_Update_List::update_LSF_vals_1");
      }
      if ( (meth_id==UBUS || meth_id==BUS) && s_thr<=ZERO && s_list[i]>ZERO ) throw FlxException_Crude("FlxBayUp_Update_List::update_LSF_vals_2");
    #endif
  }
}

void FlxBayUp_Update_List::assign_new_p_vals()
{
  // NOTE: only for BUS (NOT for BUST)
  const tuint Ns = Nc*Ncl;
  const tuint N_in_list = (fullList?Ns_final:Ns);
  for (tuint i=0;i<N_in_list;++i) {
    const tdouble ht = exp(eval_L_thresh(L_list[i],c,s_thr)-c);
    if (ht>ZERO) {
      p_list[i] = rv_InvPhi_noAlert(parent.get_data().RndCreator.gen_smp_uniform()*ht);
      y_list[N_RV*i+(N_RV-1)] = p_list[i];
      x_list[N_OX*i+(N_OX-1)] = p_list[i];
    }
  }
}

const tdouble FlxBayUp_Update_List::MHRS_uBUS_help(const tdouble maxLc, const flxVec& ut_vec, flxVec& Lt_vec)
{
  const tuint N_in_list = (fullList?Ns_final:Nc*Ncl);
  // compute the weight vector
    flxVec w_vec(sh_list,N_in_list);
    tdouble w_sum = ZERO;
    for (tuint i=0;i<N_in_list;++i) {
      sh_list[i] = exp(eval_L_thresh(L_list[i],maxLc,s_thr)-eval_L_thresh(L_list[i],max_L,s_thr));
      w_sum += sh_list[i];
    }
    if (isinf(w_sum)) return ZERO;
    w_vec /= w_sum;        // normalize w_vec
    w_sum /= N_in_list;
  // assign temporary vectors
    for (tuint i=0;i<ut_vec.get_N();++i) {
      tuint sid = parent.get_data().RndCreator.gen_smp_index2_help(ut_vec[i],w_vec);
      Lt_vec[i] = L_list[sid];
    }
  // estimate acceptance rate
    tuint count_cand_less = 0;
    tdouble pr_fins_less = ZERO;
    tdouble e_3 = ZERO;
    tdouble e_4 = ZERO;
    tuint loop_count = 0;
    for (tuint i=0;i<N_in_list;++i) {
      if (L_list[i]+s_thr<=max_L) {
        pr_fins_less += sh_list[i];
        ++count_cand_less;
      } else {
        const tdouble lt = (L_list[i]+s_thr>maxLc)?maxLc:(L_list[i]+s_thr);
        e_3 += exp(max_L-lt);
        tdouble e_4h = ZERO;
        // approximate e_4h
          const tuint e_4haN = 3;
          for (tuint i=0;i<e_4haN;++i) {
            const tdouble lt2 = (Lt_vec[loop_count]+s_thr>maxLc)?maxLc:(Lt_vec[loop_count]+s_thr);
            e_4h += ((lt>lt2)?ONE:exp(lt-lt2));
            loop_count = (loop_count+1>=ut_vec.get_N())?0:(loop_count+1);
          }
          e_4 += e_4h/e_4haN;
        // do a full analysis to get e_4h
          // for (tuint j=0;j<N_in_list;++j) {
          //   if (L_list[j]+s_thr>max_L) {
          //     const tdouble lt2 = (L_list[j]+s_thr>maxLc)?maxLc:(L_list[j]+s_thr);
          //     e_4h += ((lt>lt2)?ONE:exp(lt-lt2))*sh_list[j];
          //   }
          // }
          // e_4 += e_4h;
      }
    }
    const tdouble pr_cand_less = tdouble(count_cand_less)/tdouble(N_in_list);
    const tdouble pr_12 = pr_fins_less;
    const tdouble pr_3 = pr_cand_less * (ONE-pr_12);
    const tdouble pr_4 = (ONE-pr_cand_less)*(ONE-pr_12);
    // e_4 /= (ONE-pr_12);   // only if full analysis is used to get e_4h!!!
    e_4 /= tdouble(N_in_list-count_cand_less);
    e_3 /= tdouble(N_in_list-count_cand_less);
    const tdouble acr = pr_12+pr_3*e_3+pr_4*e_4;
  return acr;
}

void FlxBayUp_Update_List::MHRS_uBUS(tdouble& p_mod)
{
  const tdouble pa_maxL = parent.get_pa_maxL();
  if (pa_maxL >= ONE) return;
  if (meth_id!=UBUS) {
    throw FlxException_NotImplemented("FlxBayUp_Update_List::MHRS_uBUS_1");
  }
  if (pa_maxL<=ZERO) throw FlxException_Crude("FlxBayUp_Update_List::MHRS_uBUS_2");
  if (s_thr<=ZERO) return;
  #if FLX_DEBUG
    if (c!=max_L) throw FlxException_Crude("FlxBayUp_Update_List::MHRS_uBUS_3");
  #endif
  const tuint Ns = Nc*Ncl;
  const tuint N_in_list = (fullList?Ns_final:Ns);
  // compute the weights (w.r.t. maxL)
    const tdouble maxLc = max_L + s_thr;        // in log-transform
    flxVec w_vec(sh_list,N_in_list);
    // prepare 'random' vector
      flxVec ut_vec(100);
      for (tuint i=0;i<ut_vec.get_N();++i) {
        ut_vec[i] = parent.get_data().RndCreator.gen_smp_uniform();
      }
      flxVec Lt_vec(ut_vec.get_N());
    // what is max. acr?
      const tdouble acr_max = MHRS_uBUS_help(maxLc,ut_vec,Lt_vec);
    // select the threshold value for c
      tdouble c_next;
      if (acr_max > pa_maxL) {
        c_next = maxLc;
      } else {
        // root search on interval [max_L,maxLc]
          tdouble c_upper = maxLc;
          tdouble f_upper = acr_max - pa_maxL;
          tdouble c_lower = max_L;
          tdouble f_lower = ONE - pa_maxL;
          tdouble f_mid;
          tuint loop_count = 0;
          const tuint iter_max = 20;
          const tdouble tol = 1e-4;
          while (loop_count<iter_max) {
            #if FLX_DEBUG
              if (f_upper*f_lower>ZERO) throw FlxException_Crude("FlxBayUp_Update_List::MHRS_uBUS_4");
            #endif
            c_next = (c_upper*f_lower-c_lower*f_upper)/(f_lower-f_upper);
            f_mid = MHRS_uBUS_help(c_next,ut_vec,Lt_vec) - pa_maxL;
            if (f_lower*f_mid<ZERO) {
              c_upper = c_lower;
              f_upper = f_lower;
              c_lower = c_next;
              f_lower = f_mid;
            } else {
              const tdouble m = f_lower/(f_lower+f_mid);
              f_upper = m*f_upper;
              c_lower = c_next;
              f_lower = f_mid;
            }
            // check for convergence
              if (fabs(f_mid)<=tol) {
                break;
              }
              if (fabs(c_upper-c_lower)<=tol) {
                c_next = (c_upper+c_lower)/2;
                f_mid = MHRS_uBUS_help(c_next,ut_vec,Lt_vec) - pa_maxL;
                break;
              }
            ++loop_count;
          }
      }
    // evaluate the weight vector
      pdouble pw_i;
      for (tuint i=0;i<N_in_list;++i) {
        pw_i += exp(eval_L_thresh(L_list[i],c_next,s_thr)-eval_L_thresh(L_list[i],max_L,s_thr));
      }
      tdouble w_i = pw_i.cast2double();
      w_i /= N_in_list;
      // evaluate C.o.V. of weights
        tdouble w_cov = ZERO;
        for (tuint i=0;i<N_in_list;++i) {
          w_cov += pow2(exp(eval_L_thresh(L_list[i],c_next,s_thr)-eval_L_thresh(L_list[i],max_L,s_thr))-w_i);
        }
        w_cov /= N_in_list;
        w_cov = sqrt(w_cov)/w_i;
      w_i = log(w_i);
      w_i += (max_L - c_next);
      p_mod += w_i;
      ext_out << std::format("qw={:4.2f} ", exp(w_i));
  // perform a re-sampling to alter the sampling distribution
    // pick the seed randomly
      tuint sid = parent.get_data().RndCreator.gen_smp_index2(w_vec);
      swap_smpls(0,sid);
    // perform the resampling
      tdouble acr = ZERO;
      for (tuint i=1;i<N_in_list;++i) {
        // make sure the samples are not ordered
          sid = i + parent.get_data().RndCreator.gen_smp_index(N_in_list-i);
          swap_smpls(i,sid);
        // evaluate the acceptance probability
          const tdouble l_prev = eval_L_thresh(L_list[i-1],c_next,s_thr);
          const tdouble l_cand = eval_L_thresh(L_list[i],c_next,s_thr);
          tdouble ut = exp(((l_cand<max_L)?max_L:l_cand)-l_prev);
          ut = (ut>ONE)?ONE:ut;
          acr += ut;
        // reject/accept the sample
          const tdouble u = parent.get_data().RndCreator.gen_smp_uniform();
          if (u>ut) {        // reject it ...
            // y-list
              flxVec y1(y_list+N_RV*(i-1),N_RV);
              flxVec y2(y_list+N_RV*i,N_RV);
              y2 = y1;
            // x-list
              flxVec x1(x_list+N_OX*(i-1),N_OX);
              flxVec x2(x_list+N_OX*i,N_OX);
              x2 = x1;
            // s-list
              s_list[i] = s_list[i-1];
            // L-list
              L_list[i] = L_list[i-1];
          }
      }
      acr /= (N_in_list-1);
      ext_out << std::format("qacr={:4.2f} qCoV={:4.2f} ", acr, w_cov);
  // correct s_thr
    s_thr += max_L - c_next;
    ext_out << std::format("gt2={:9.2e} ", s_thr);
}

const bool FlxBayUp_Update_List::update_c(tdouble& p_mod, const bool is_first)
{
  bool zeroLSF = is_gt_zero();
  #if FLX_DEBUG
    if (!(meth_id==BUS||meth_id==UBUS||meth_id==BUST)) throw FlxException_Crude("FlxBayUp_Update_List::update_c_0");
    if (is_first&&zeroLSF) throw FlxException_Crude("FlxBayUp_Update_List::update_c_1");
  #endif
  if (max_L==c && !is_first) {
    #if FLX_DEBUG
      for (tuint i=0;i<Ns_final;++i) {
        if (i_list[i]<=0) break;
        if ( fabs(s_list[i] - eval_LSF (p_list[i],L_list[i]))>GlobalVar.TOL() ) {
          throw FlxException_Crude("FlxBayUp_Update_List::update_c_4");
        }
      }
    #endif
    if (meth_id==UBUS) {
      MHRS_uBUS(p_mod);
      assign_new_p_vals();
      update_LSF_vals(Ns_final,false);
    }
    return true;
  }
  // remember this c
    if (!is_first) {
      if (meth_id==UBUS) {
        const tdouble delta = max_L - c;
        s_thr += delta;
        ext_out << std::format("gt1={:9.2e} ", s_thr);
      } else {
        if (oldC_N>=max_runs) throw FlxException_Crude("FlxBayUp_Update_List::update_c_5");
        oldC_list[oldC_N] = c;                // this is the max. Likelihood of the previous step
        oldS_list[oldC_N] = s_thr;        // this is the threshold of the previous step
        ++oldC_N;
      }
    }
  c = max_L;
  if (meth_id==BUST) return false;
  if (meth_id==UBUS && !is_first) {
    MHRS_uBUS(p_mod);
    assign_new_p_vals();
    zeroLSF = is_gt_zero();
  }
  if (zeroLSF && fullList) {        // we have to kill some samples
    bool fixed = false;
    for (tuint i=0;i<Ns_final;++i) {
      if (i_list[i]==-1) break;        // end of list
      s_list[i] = eval_LSF(p_list[i],L_list[i]);
      if (s_list[i]>ZERO) {
        #if FLX_DEBUG
          if (i_list[i]<0) throw FlxException_Crude("FlxBayUp_Update_List::update_c_7");
        #endif
        fixed = true;
      }
    }
    if (!fixed && (meth_id==UBUS) ) {
      s_thr = ZERO;
    }
    return !fixed;
  } else {                // no samples are killed !!! -> LSF is just updated
    update_LSF_vals(Ns_final,!(is_first||meth_id==BUS||meth_id==UBUS));
    return false;
  }
}

void FlxBayUp_Update_List::nested_c()
{
  #if FLX_DEBUG
    if (!(meth_id==BUS||meth_id==BUST)) throw FlxException_Crude("FlxBayUp_Update_List::nested_c_0");
  #endif
  for (tuint i=oldC_N;i>0;--i) {
    tuint j = i-1;
    // check if entry is obsolete
      bool erase = false;
      if (meth_id==BUS) {
        erase = (oldS_list[j]>=s_thr);
      } else {
        erase = (oldS_list[j]<=s_thr);
      }
    // erase the entry
      if (erase) {
        for (tuint k=j+1;k<oldC_N;++k) {
          oldS_list[k-1] = oldS_list[k];
          oldC_list[k-1] = oldC_list[k];
        }
        --oldC_N;
      }
  }
}

const bool FlxBayUp_Update_List::update_c_posterior(const bool smpl_part_of_posterior, const tdouble lMax_)
{
  #if FLX_DEBUG
    if (meth_id!=BUS && meth_id!=UBUS && meth_id!=BUST) {
      throw FlxException_Crude("FlxBayUp_Update_List::update_c_posterior_1");
    }
  #endif
  const tdouble p_true = (smpl_part_of_posterior?(p_list[curID]):ZERO);        // stdN-transform
  // rescale the probabilities in the p_list (-> scale p-values to lower region)
    max_L = (smpl_part_of_posterior?(L_list[curID]):lMax_);
    #if FLX_DEBUG
      if (max_L<=c) throw FlxException_Crude("FlxBayUp_Update_List::update_c_posterior_2");
    #endif
    const tdouble c_old = c;        // log-transform
    c = max_L;                        // log-transform
    const tdouble scale = exp(c_old-c);        // < 1.0
    for (tuint i=0;i<Ns_final;++i) {
      const tdouble pt = rv_InvPhi(rv_Phi(p_list[i])*scale);                // scale the p down !!!
      p_list[i] = pt;
      y_list[i*N_RV+N_RV-1] = pt;
      s_list[i] = eval_LSF(pt,L_list[i]);
      #if FLX_DEBUG
        if (s_list[i]>ZERO) throw FlxException_Crude("FlxBayUp_Update_List::update_c_posterior_3");
      #endif
    }
    if (!smpl_part_of_posterior) return false;
  // rescale the new probability (the p-value of the sample in the upper region)
    const tdouble p_new = rv_InvPhi(rv_Phi(p_true)*(ONE-scale)+scale);
  // how many samples to draw from upper region
    // -> lets gather all samples in the lower region in the set bulk_1
    // -> lets gather all samples in the upper region in the set bulk_2
    const tdouble pr_bulk_2 = (ONE-scale)/((ONE-scale)+(Ns_final*scale));  // probability that sample is drawn from bulk_2
    RBRV_entry_RV_Binomial bino_t("dummy",0,new FlxFunction(new FunNumber(pr_bulk_2)),new FlxFunction(new FunNumber(tdouble(Ns_final))),true);
    // number of samples to draw from upper region
      const tuint k = round_flx(bino_t.transform_y2x(RndCreator->gen_smp()));
  // which indices to replace?
    #ifndef FLX_CV_1
      tuint k_id[k];
    #else
      tuint* k_id = new tuint[k];
    #endif
    for (tuint i=0;i<k;++i) {
      k_id[i] = RndCreator->gen_smp_index(Ns_final-i);
      bool b;
      do {
        b=false;
        for (tuint j=0;j<i;++j) {
          if (k_id[j]==k_id[i]) {
            b = true;
            ++(k_id[i]);
            break;
          }
        }
      } while (b);
    }
  // loop over all entries to replace
    bool b = false;        // true: the current sample was replaced
    for (tuint i=0;i<k;++i) {
      if (k_id[i]==curID) {
        #if FLX_DEBUG
          if (b) throw FlxException_Crude("FlxBayUp_Update_List::update_c_posterior_4");;
        #endif
        b = true;
      } else {
        L_list[k_id[i]] = L_list[curID];
        flxVec sy_curID(&(y_list[curID*N_RV]),N_RV);
        flxVec sx_curID(&(x_list[curID*N_OX]),N_OX);
        flxVec sy(&(y_list[k_id[i]*N_RV]),N_RV);
        flxVec sx(&(x_list[k_id[i]*N_OX]),N_OX);
        sy = sy_curID;
        sx = sx_curID;
        i_list[k_id[i]] = 0;
      }
      p_list[k_id[i]] = p_new;
      y_list[k_id[i]*N_RV+N_RV-1] = p_new;
      s_list[k_id[i]] = eval_LSF(p_new,L_list[k_id[i]]);
      #if FLX_DEBUG
        if (s_list[k_id[i]]>ZERO) throw FlxException_Crude("FlxBayUp_Update_List::update_c_posterior_5");
      #endif
    }
    #ifdef FLX_CV_1
      delete [] k_id;
    #endif
  return b;
}

void FlxBayUp_Update_List::compute_gamma(Flx_SuS_CLevelStat& curLevelStat, const Flx_SuS_Control& sctrl) const
{
  const tuint Nc_now = curLevelStat.Nchains;
  if (Nc_now==0) return;
  const tuint max_nCl = curLevelStat.g_Ncl_max;
  if (max_nCl<=1) return;
  const bool consider_seed_corr = ((curLevelStat.get_level()>=2)?sctrl.consider_seed_corr:false);
  // some initializations
    tuint N_00 = 0;
    tuint N_11 = 0;
    // allocate local vector of lsf-relizations
      #ifndef FLX_CV_1
        tdouble sl[max_nCl];
        tuint sl_h[max_nCl*2];
      #else
        tdouble* sl = new tdouble[max_nCl];
        tuint* sl_h = new tuint[2*max_nCl];
      #endif
      memset(sl_h,0, sizeof(tuint)*(2*max_nCl));
    // allocate level-statistics vector
      curLevelStat.g_N_entries = max_nCl-1;
      curLevelStat.g_N = new tuint[curLevelStat.g_N_entries*3];
        for (tuint i=0;i<curLevelStat.g_N_entries*3;++i) {
          curLevelStat.g_N[i] = 0;
        }
    // allocate memory for seed history
      if (curLevelStat.get_level()>0) {        // level=0: perfect Monte Carlo samples (uncorrelated)
        curLevelStat.seed_chainID = new tuint[curLevelStat.Nfailures];
        curLevelStat.seed_chainPos = new tuint[curLevelStat.Nfailures];
      }
      tuint N_seeds_found_so_far = 0;        // number of seeds found so far
    // variables required for the Gelman-estimate
      tdouble sum_sj2 = ZERO;
      tdouble sum_r = ZERO;
    const tdouble st = ((meth_id==BUST)?ZERO:s_thr);        // the current threshold
    const tuint Ns = curLevelStat.Nsamples;
    const tdouble p_i = tdouble(curLevelStat.Nfailures)/tdouble(Ns);
    const tdouble R0 = p_i*(ONE-p_i);
    tuint i = 0;                // counts the currend index in the sample list
    tuint i_sp = 0;        // index of the next seed index that nees to be skipped
  // find the initial i
    while (true) {        // set entry to first non-seed
      if (i!=curLevelStat.g_seed_ID_original[i_sp]) break;
      ++i;
      ++i_sp;
      if (i>=Ns) break;
      if (i_sp>=Nc_now) break;
    }
  for (tuint j=0;j<Nc_now;++j) {        // loop over all seeds
    const tuint nCl = curLevelStat.g_chain_length[j];
    tuint nh_count = 0;        // count the number of hits in the current chain
    // gather all lsf-values in one vector
      sl[0] = s_list[curLevelStat.g_seed_ID[j]];
      // statistics of the ith element in the chain
        if (sl[0]<=st) {
          ++(sl_h[0]);
          ++nh_count;
        }
        ++(sl_h[1]);
      for (tuint k=1;k<nCl;++k) {
        #if FLX_DEBUG
          if (i>=Ns) {
            throw FlxException_Crude("FlxBayUp_Update_List::compute_Rk_c1");
          }
        #endif
        sl[k] = s_list[i];
        // staistics of the ith element in the chain
          if (sl[k]<=st) {
            ++(sl_h[2*k]);
            ++nh_count;
          }
          ++(sl_h[2*k+1]);
        // find the next entry that is not a seed
          ++i;
          while (true) {
            if (i>=Ns) break;
            if (i_sp>=Nc_now) break;
            if (i!=curLevelStat.g_seed_ID_original[i_sp]) break;
            ++i;
            ++i_sp;
          }
      }                // end: loop over all seeds
    // store seed-history
      if (curLevelStat.get_level()>0) {
        for (tuint k=0;k<nCl;++k) {
          if (sl[k]<=st) {
            curLevelStat.seed_chainID[N_seeds_found_so_far] = j;
            curLevelStat.seed_chainPos[N_seeds_found_so_far] = k;
            ++N_seeds_found_so_far;
          }
        }
      }
    // store values in slc
    if (consider_seed_corr) {
      // find seed it in sorted list (original order of seeds in previous list)
        tuint sid = curLevelStat.g_seed_ID[j];
        {
          bool found = false;
          for (tuint k=0;k<Nc_now;++k) {
            if (curLevelStat.g_seed_ID_original[k] == sid) {
              sid = k;
              found = true;
              break;
            }
          }
          if (found==false) throw FlxException_Crude("FlxBayUp_Update_List::compute_Rk_c2");
        }
      // change meaning of g_seed_ID
        curLevelStat.g_seed_ID[j] = sid;
    }
    // compute Rk's
      for (tuint k=1;k<nCl;++k) {
        tuint N_hits1 = 0;        // number of successes
        tuint N_hits2 = 0;        // number of successes
        const tuint N_t = nCl-k;
        tuint *g_N_p = &(curLevelStat.g_N[3*(k-1)]);
        g_N_p[0] += N_t;
        for (tuint l=0;l<N_t;++l) {
          if ( sl[l]<=st ) ++N_hits1;
          if ( sl[l]<=st && sl[l+k]<=st ) ++N_hits2;
        }
        g_N_p[1] += N_hits1;
        g_N_p[2] += N_hits2;
      }
    // compute p_00 and p_11
      for (tuint k=1;k<nCl;++k) {
        if (sl[k-1]>st&&sl[k]>st) ++N_00;
        if (sl[k-1]<=st&&sl[k]<=st) ++N_11;
      }
    // compute values for Gelman-estimator
      // for sj2 we divide by nCl instead of (nCl-1) ... this is accounted for when the final estimate is evaluated
      tdouble sj2 = tdouble(nh_count)/nCl;        // sample mean of chain
      sum_r += pow2(sj2-p_i);                        // ... for the variance of the sample mean
      sj2 = (sj2)*(ONE-sj2);                        // biased estimate of sample variance
      sum_sj2 += sj2;
  }
  // evaluate the Gelman-estimate
    sum_sj2 /= Nc_now;
    sum_r /= (Nc_now-1);
    curLevelStat.eff_Gelman = (sum_sj2 + sum_r)/((curLevelStat.Nsamples/Nc_now)*sum_r);        // Gelman-estimate
  // compute statistics of the ith element in the chain
    curLevelStat.pi_chain = new tdouble [max_nCl];
    for (tuint k=0;k<max_nCl;++k) {
      if (sl_h[2*k+1]>0) {
        curLevelStat.pi_chain[k] = tdouble(sl_h[2*k])/tdouble(sl_h[2*k+1]);
      } else {
        curLevelStat.pi_chain[k] = tdouble(curLevelStat.Nfailures)/tdouble(curLevelStat.Nsamples);
      }
    }
  #ifdef FLX_CV_1
    delete [] sl;
    delete [] sl_h;
  #endif
  #if FLX_DEBUG
    if (N_seeds_found_so_far!=curLevelStat.Nfailures) throw FlxException_Crude("FlxBayUp_Update_List::compute_Rk_c3");
    if (i_sp != Nc_now) throw FlxException_Crude("FlxBayUp_Update_List::compute_Rk_2");
    if (i!= Ns) throw FlxException_Crude("FlxBayUp_Update_List::compute_Rk_3");
    if (curLevelStat.g_N_entries==0) throw FlxException_Crude("FlxBayUp_Update_List::compute_Rk_4");
  #endif
  // finally, compute gamma (& thus, the efficiency)
    // evaluate gamma
      // contribution chains
        #if FLX_DEBUG
          if (curLevelStat.g_Ncl_max-1!=curLevelStat.g_N_entries) throw FlxException_Crude("FlxBayUp_Update_List::compute_Rk_5");
        #endif
        const tdouble p_i_2 = pow2(p_i);
        pdouble gam_res;
        // compute the lag-1 autocorrelation
          curLevelStat.lag1_corr = tdouble(curLevelStat.g_N[2])/tdouble(curLevelStat.g_N[1]);
          curLevelStat.lag1_corr = (curLevelStat.lag1_corr-curLevelStat.pi)/(ONE-curLevelStat.pi);
        // compute p_00 and p_11
          if (curLevelStat.g_N[2]!=N_11) throw FlxException_Crude("FlxBayUp_Update_List::compute_Rk_6");
          curLevelStat.p_00 = tdouble(N_00)/tdouble(curLevelStat.g_N[0]);
          curLevelStat.p_11 = tdouble(N_11)/tdouble(curLevelStat.g_N[0]);
        for (tuint k=0;k<curLevelStat.g_N_entries;++k) {
          tuint *g_N_p = &(curLevelStat.g_N[3*k]);
          // determin rho
            pdouble rho = tdouble(g_N_p[2]);
            rho -= g_N_p[0]*p_i_2;
            rho /= R0;
            // we do not allow for negative entries
              if (rho.cast2double()<ZERO) rho = ZERO;
          rho *= 2;
          rho /= tdouble(curLevelStat.Nsamples);
          gam_res += rho;
        }
        curLevelStat.gamma_chain = gam_res.cast2double();
      // evluate seed correlation
        if (consider_seed_corr) {
          if (sctrl.empirical_corr) {
            curLevelStat.empirical_Corr();
          } else {
            curLevelStat.add2seedCorr();
          }
        }
      // contribution seed correlation
        curLevelStat.gamma = curLevelStat.gamma_chain + curLevelStat.gamma_from_seed;
        #if FLX_DEBUG
          if (!consider_seed_corr && curLevelStat.gamma_from_seed!=ZERO) throw FlxException_Crude("FlxBayUp_Update_List::compute_Rk_7");
        #endif
    // calculate efficiency
      curLevelStat.eff = ONE/(ONE+curLevelStat.gamma);
  // correlation of the pi
    if ( consider_seed_corr && sctrl.consider_pi_corr && sctrl.empirical_corr==false ) {
      curLevelStat.add2piCorr();
    }
}

void FlxBayUp_Update_List::prepare_for_gamma_comp_1(Flx_SuS_CLevelStat& curLevelStat_next) const
{
  const tuint Nc_now = curLevelStat_next.Nchains;
  if (Nc_now==0) return;
  curLevelStat_next.g_seed_ID_original = new tuint[Nc_now];
  for (tuint i=0;i<Nc_now;++i) {
    curLevelStat_next.g_seed_ID_original[i] = seed_idx[i];
  }
  curLevelStat_next.g_chain_length = new tuint[Nc_now];
  for (tuint i=0;i<Nc_now;++i) {
    curLevelStat_next.g_chain_length[i] = chain_length[i];
  }
  // determin the maximum chain length
    tuint max_nCl = 0;
    for (tuint i=0;i<Nc_now;++i) {
      if (max_nCl<chain_length[i]) max_nCl = chain_length[i];
    }
    curLevelStat_next.g_Ncl_max = max_nCl;
}

void FlxBayUp_Update_List::prepare_for_gamma_comp_2(Flx_SuS_CLevelStat& curLevelStat_next) const
{
  const tuint Nc_now = curLevelStat_next.Nchains;
  if (Nc_now==0) return;
  curLevelStat_next.g_seed_ID = new tuint[Nc_now];
  for (tuint i=0;i<Nc_now;++i) {
    curLevelStat_next.g_seed_ID[i] = seed_idx[i];
  }
}

const tdouble FlxBayUp_Update_List::get_perc_BUST()
{
  if (s_thr<=ZERO) return ONE;
  const tuint Ns = fullList?Ns_final:(Nc*Ncl);
  // update LSF values
    for (tuint i=0;i<Ns;++i) {
      #if FLX_DEBUG
        if (i_list[i]<=0) throw FlxException_Crude("FlxBayUp_Update_List::get_perc_BUST_1");
      #endif
      s_list[i] = eval_LSF(p_list[i],L_list[i]);
    }
  // count number of negative entries
    tuint Nneg = 0;
    tdouble mneg = ONE;
    tdouble mpos = -ONE;
    for (tuint i=0;i<Ns;++i) {
      if (s_list[i]<=ZERO) {
        ++Nneg;
        if (mneg>ZERO || mneg<s_list[i]) {
          mneg = s_list[i];
        }
      } else {
        if (mpos<ZERO || mpos>s_list[i]) {
          mpos = s_list[i];
        }
      }
    }
  if (Nneg==0) return ZERO;
  if (Nneg>=Ns) return ONE;
  const tdouble dp = ONE/Ns;
  const tdouble pa = (ZERO-mneg)/(mpos-mneg);
  return (Nneg-ONE/2+pa)*dp;
}

const tuint FlxBayUp_Update_List::update_thr_BUST(tdouble& pnow, Flx_SuS_CLevelStat& curLevelStat_next, Flx_SuS_CLevelStat& curLevelStat)
{
  #if FLX_DEBUG
    if (meth_id!=BUST) throw FlxException_Crude("FlxBayUp_Update_List::update_thr_BUST");
  #endif
  const tdouble p0 = ONE/Ncl;
  // count number of entries different from -inf
    const tuint Ns = fullList?Ns_final:(Nc*Ncl);
    tuint NnotInf = 0;
    for (tuint i=0;i<Ns;++i) {
      if (L_list[i]!=log(ZERO)) ++NnotInf;
    }
  if ( tdouble(NnotInf)/tdouble(Ns) < p0 ) {
    s_thr = GlobalVar.TOL();
    pnow = get_perc_BUST();
  } else {
    // determin pnow
      tdouble perc_prev = ONE;
      tdouble perc_next;
      tuint i; 
      const tuint Nin = 5;
      const int Npot = 8;
      for (i=1;i<=Nin;++i) {
        s_thr =  ((i==Nin)?ONE:pow(i/(ONE*Nin),Npot));
        perc_next =get_perc_BUST();
        if (perc_next<=p0) break;
        else perc_prev = perc_next;
      }
      if (i<=Nin) {                // get the root
        tdouble lb = pow((i-1)/(ONE*Nin),Npot);
        tdouble ub = s_thr;
        const tuint NMAX = tuint(1e2);
        const tdouble dx = 1e-5;
        const tdouble dy = 1e-4;
        tuint c = 0;
        tdouble f1 = perc_prev-p0;
        tdouble f2 = perc_next-p0;
        bool found = false;
        if (fabs(f1)<=GlobalVar.TOL()) {
          s_thr = lb;
          get_perc_BUST();
          found = true;
        }
        if (!found && fabs(f2)<=GlobalVar.TOL()) {
          found = true;
        }
        if (found==false) {
          while (fabs(ub-lb)*2/fabs(ub+lb)>dx && c<NMAX) {
            ++c;
            const tdouble res = (lb*f2-ub*f1)/(f2-f1);
            s_thr = res;
            perc_next = get_perc_BUST();
            const tdouble f3 = perc_next-p0;
            if (fabs(f3/p0)<=dy) {
              found = true;
              break;
            }
            if ( f2*f3<ZERO ) {
              lb = ub;
              f1 = f2;
              ub = res;
              f2 = f3;
            } else {
              const tdouble m = f2/(f2+f3);
              f1 = m*f1;
              ub = res;
              f2 = f3;
            }
          }
        }
        if (found==false) {
          if (ub>lb) {
            s_thr = lb;
          } else {
            s_thr = ub;
          }
          perc_next = get_perc_BUST();
        }
      } else {
        s_thr = ONE;        // just to make sure
      }
    pnow = perc_next;
  }
  if (log_LSF) {
    GlobalVar.slog(4) << std::endl << "BUST: " << GlobalVar.Double2String(pnow,false,-1,0) << GlobalVar.Double2String(s_thr,false,-1,0);
  }
  // call update_thr
    tuint tmpri = update_thr(curLevelStat_next, curLevelStat);
    if (log_LSF) {
      GlobalVar.slog(4) << " " << tmpri;
      tuint cti = 0;
      for (tuint i=0;i<Ns_final;++i) {
        if (i_list[i]<0) break;
        else if (i_list[i]==2) {
          GlobalVar.slog(4) << std::endl << GlobalVar.Double2String(cti,false,-1,0) << GlobalVar.Double2String(s_list[i],false,-1,0) << GlobalVar.Double2String(p_list[i],false,-1,0)  << GlobalVar.Double2String(L_list[i],false,-1,0);
          ++cti;
        }
      }
    }
    return tmpri;
}

const bool FlxBayUp_Update_List::is_gt_zero() const
{
  if (parent.get_methCat()==flxBayUp::ABC) { 
    return (s_thr==parent.get_cStart());
  }
  if (meth_id==BUST || meth_id==TMCMC) {
    return (fabs(s_thr-ONE)<=GlobalVar.TOL());
  }
  if (meth_id==RS) {
    return true;
  }
  return (fabs(s_thr)<=GlobalVar.TOL());
}

const tuint FlxBayUp_Update_List::update_thr(Flx_SuS_CLevelStat& curLevelStat_next, Flx_SuS_CLevelStat& curLevelStat)
{
  curID = 0;        // reset counter
  curSID = 0;        // reset seed-position
  // ====================================================
  // Obtain the requested percentile
  //   Nse: number of seeds
  //   s_thr: threshold value
  //   (sets also fullList ...)
  // ====================================================
  const tuint Ns = Nc*Ncl;
  const tuint N_in_list = (fullList?Ns_final:Ns);
  // count number of samples that already follow posterior
    const tdouble z_zero = ((parent.get_methCat()==flxBayUp::ABC)?(parent.get_cStart()):ZERO);
    tuint Nse = 0;        // Number of seeds
    for (tuint i=0;i<N_in_list;++i) if (s_list[i]<=z_zero) ++Nse;
  if (meth_id==BUST) {                // this is a special case -> seeds are already known
    if (fabs(s_thr-ONE)<=GlobalVar.TOL()) {        // we reached the last level with BUST
      fullList = true;
      if (Nse<Nc) throw FlxException_Crude("FlxBayUp_Update_List::update_thr_a1");
    } 
  }
  {
    // sort the list
      // copy to temporary vector
      {
        flxVec sh_vec(sh_list,N_in_list);
        flxVec s_vec(s_list,N_in_list);
        sh_vec = s_vec;
      }
      std::sort(sh_list,sh_list+N_in_list);
      const tdouble mins = sh_list[Nc-1];
      tuint cP = Nc;        // number of seeds in list
      for (tuint i=Nc;i<N_in_list;++i) {
        if (sh_list[i]==mins) {                // MULTIPLE!
          ++cP;
        } else {
          break;
        }
      }
      // cP is the number of seeds!!!
        if (cP>=Ns) {        // only seeds in the list!!!
          std::ostringstream ssV;
          ssV << "There is a problem with the correlation of the chain: in the list are too many entries that are the same.";
          throw FlxException("FlxBayUp_Update_List::update_thr_a2",ssV.str());
        }
    // adapt threshold for BUS
      if (meth_id!=BUST) {        // Method is NOT BUST
        if (Nse<Nc) {                // we have to sort
          s_thr = (mins+sh_list[cP])/2;
          #if FLX_DEBUG
            Nse = 0;
            for (tuint i=0;i<N_in_list;++i) {
              if (s_list[i]<=s_thr) {
                ++Nse;
              }
            }
            if (Nse!=cP) throw FlxException_Crude("FlxBayUp_Update_List::update_thr_a4");
          #endif
          Nse = cP;
        } else {
          s_thr = z_zero;
          fullList = true;
        }
      }
    // find multiple samples
      if (multpleVec_help) {
        // empty vector
          curLevelStat.find_multiples = new tuint[N_in_list];
          tuint* multpleVec = curLevelStat.find_multiples;
          for (tuint i=0;i<N_in_list;++i) {
            multpleVec[i] = 0;
          }
          for (tuint i=0;i<N_in_list;++i) {
            multpleVec_help[i] = false;
          }
        for (tuint i=0;i<N_in_list;++i) {
          if (multpleVec_help[i]) continue;
          tuint tmp_n_rep = 0;
          flxVec vti(y_list+i*N_RV,N_RV);
          for (tuint j=i+1;j<N_in_list;++j) {
            if (s_list[i]==s_list[j]) {              
              flxVec vtj(y_list+j*N_RV,N_RV);
              if (vti==vtj) {
                ++tmp_n_rep;
                multpleVec_help[j] = true;
              }
            }
          }
          ++(multpleVec[tmp_n_rep]);
        }
      }
    // generate extented output
      {
        std::ostringstream ext_out2;
        tdouble s_prev = sh_list[0];
        tuint mult_seq_count = 0;
        tuint mult_seq_this = 0;
        tuint longest_rep = 0;
        tuint single_count = 0;
        for (tuint i=1;i<Nse;++i) {
          if (sh_list[i]==s_prev) {
            if (mult_seq_this==0) ++mult_seq_count;
            ++mult_seq_this;
          } else {
            s_prev = sh_list[i];
            if (mult_seq_this>0) {
              if (mult_seq_this>longest_rep) longest_rep = mult_seq_this;
              mult_seq_this = 0;
            } else {
              ++single_count;
            }
          }
        }
        if (mult_seq_this==0) ++single_count;;
        if (mult_seq_this>longest_rep) longest_rep = mult_seq_this;
        if (mult_seq_count>0) {
          ext_out2 << "multiple:";
          ext_out2 << round_flx(tdouble(mult_seq_count)/tdouble(single_count+mult_seq_count),2) 
                  << '/' << round_flx(tdouble(longest_rep+1)/tdouble(single_count+mult_seq_count),2)
                  << '/' << round_flx(tdouble(Nse-(single_count+mult_seq_count))/tdouble(Nse),2);
          if (mult_seq_this>0 ) {
            ext_out2 << '/' << mult_seq_this+1;
          }
        }
        if (Nse<Nc) {
          if (ext_out2.str().empty()==false) ext_out2 << '/';
          ext_out2 << "pt=" << round_flx(tdouble(Nse)/tdouble(N_in_list)*100) << '%';
        }
        if (ext_out2.str().empty()==false && ext_out2.str().size()<27) {
          ext_out2 << std::string(27-ext_out2.str().size(),' ');
        }
        // seed statistics
          const tdouble s_thr_tmp = (meth_id==BUST?ZERO:s_thr);
          tuint cP2 = 0; tuint cP3 = 0;
          for (tuint i=0;i<N_in_list;++i) {
            if (i_list[i] != 2) continue;
            ++cP2;
            if (s_list[i]<=s_thr_tmp) {        // it is a seed value
              ++cP3;
            }
          }
          if (cP2>0) {
            const tdouble pts = (tdouble(cP3)/tdouble(cP2))/(tdouble(Nse)/tdouble(N_in_list));
            ext_out2 << std::format("pts={:4.2f}  ", pts);
          }
        // append output to ext_out2
          if (ext_out.str().empty()==false && ext_out2.str().empty()==false) {
            ext_out << std::endl << "            ";
          }
          ext_out << ext_out2.str();
      }   // end of extented output
  }
  // fill up the list
    if (fullList&&(N_in_list!=Ns_final)) {
      for (tuint i=Ns;i<Ns_final;++i) s_list[i] = 12345.6;
      for (tuint i=Ns;i<Ns_final;++i) i_list[i] = 0;
    }
  const tuint Ns2 = (fullList?Ns_final:Ns);
  // ====================================================
  // Find the seeds
  //   i_list: (0:no seed; 1:seed; -1:end)
  //   seed_idx: indices of seed-positions (ends with Ns_final+1)
  // ====================================================
  {
    tuint cP = 0;        // current seed
    // for BUST s_thr has a different meaning -> account for it:
      const tdouble s_thr_tmp = (meth_id==BUST?ZERO:s_thr);
    for (tuint i=0;i<N_in_list;++i) {
      if (s_list[i]<=s_thr_tmp) {        // it is a seed value
        seed_idx[cP] = i;
        ++cP;
        i_list[i] = 2;
      } else {                        // it is not a seed value
        i_list[i] = 0;
      }
    }
    if (Nse!=cP) throw FlxException_Crude("FlxBayUp_Update_List::update_thr_b1");
    if (cP<Ns_final) seed_idx[cP] = Ns_final+1;
    for (tuint i=N_in_list;i<Ns2;++i) {
      i_list[i] = 0;
    }
  }
  // ====================================================
  // compute length of individual chains
  // ====================================================
  {
    tuint Ngiven = Ns2/Nse;
    #if FLX_DEBUG
      if (Ngiven==0) throw FlxException_Crude("FlxBayUp_Update_List::update_thr_c5");
    #endif
    if (adpt_ctrl->get_smpl_order()==2) {
      const tuint Ndist = Ns2 - (Ngiven*Nse);
      for (tuint i=0;i<Ndist;++i) {
        chain_length[i] = Ngiven + 1;
      }
      for (tuint i=Ndist;i<Nse;++i) {
        chain_length[i] = Ngiven;
      }
    } else {
      for (tuint i=0;i<Nse;++i) {
        chain_length[i] = Ngiven;
      }
      Ngiven *= Nse;
      Ngiven = Ns2 - Ngiven;
      for (tuint i=0;i<Ngiven;++i) {
        ++chain_length[RndCreator->gen_smp_index(Nse)];
      }
    }
    #if FLX_DEBUG
      tuint t_c = 0;
      for (tuint i=0;i<Nse;++i) {
        t_c += chain_length[i];
      }
      if (t_c!=Ns2) throw FlxException_Crude("FlxBayUp_Update_List::update_thr_c6");
    #endif
  }
  // ====================================================
  // Store results for statistical purposes
  // ====================================================
  curLevelStat_next.Nchains = Nse;
  curLevelStat_next.Nsamples = Ns2;
  prepare_for_gamma_comp_1(curLevelStat_next);
  // ====================================================
  // Shuffle seeds?
  // ====================================================
  {
    const tuint so = adpt_ctrl->get_smpl_order();
    switch (so) {
      case 0:
        break;
      case 1:        // large to small g
        throw FlxException_NotImplemented("FlxBayUp_Update_List::update_thr_c1");
        break;
      case 2:        // random order (random shuffling)
        std::shuffle(seed_idx,seed_idx+Nse,get_rng());
        break;
      case 3:        // length
        throw FlxException_NotImplemented("FlxBayUp_Update_List::update_thr_c3");
        break;
      default:
        throw FlxException_Crude("FlxBayUp_Update_List::update_thr_c4");
    }
  }
  prepare_for_gamma_comp_2(curLevelStat_next);
  // ====================================================
  // for BUS & BUST -> can we forget the history???
  // ====================================================
  if (meth_id==BUS||meth_id==BUST) {
    nested_c();
    if ( 
      ((meth_id==BUS)&&(s_thr==ZERO))
      ||((meth_id==BUST)&&(fabs(s_thr-ONE)<=GlobalVar.TOL()))
    ) {
      oldC_N = 0;
    }
  }
  // ====================================================
  // return number of seeds
  // ====================================================
  set_next(true);
  return Nse;
}

const tuint FlxBayUp_Update_List::finalize()
{
  if (finalized) throw FlxException_Crude("FlxBayUp_Update_List::finalize_1");
  switch (meth_id) {
    case (FlxBayUp_Update_List::BUS):
    case (FlxBayUp_Update_List::UBUS):
    case (FlxBayUp_Update_List::BUST):
    case (FlxBayUp_Update_List::ABCSUBSIM):
    case (FlxBayUp_Update_List::RASUBSIM):
    case (FlxBayUp_Update_List::TMCMC):
    case (FlxBayUp_Update_List::MHRS):
    {
      #if FLX_DEBUG
      tuint Nneg = 0;
      const tdouble z_zero = ((parent.get_methCat()==flxBayUp::ABC)?(parent.get_cStart()):ZERO);
      #endif
      for (tuint i=0;i<Ns_final;++i) {
        #if FLX_DEBUG
        if (meth_id==FlxBayUp_Update_List::MHRS || s_list[i]<=z_zero) {
        #endif
          i_list[i] = 0;
        #if FLX_DEBUG
          ++Nneg;
        } else {
          i_list[i] = -2;
          throw FlxException_Crude("FlxBayUp_Update_List::finalize_2a");
        }
        #endif
      }
      #if FLX_DEBUG
      if (Nneg!=Ns_final && meth_id!=FlxBayUp_Update_List::MHRS) throw FlxException_Crude("FlxBayUp_Update_List::finalize_2");
      #endif
      // randomize list (relevant for randomizeID==INIT)
        for (tuint i=0;i<Ns_final;++i) {
          seed_idx[i] = i;
        }
        std::shuffle(seed_idx,seed_idx+Ns_final,get_rng());
      break;
    }
    case (FlxBayUp_Update_List::RS):
    case (FlxBayUp_Update_List::ABCRS):
    case (FlxBayUp_Update_List::RAMCI):
    case (FlxBayUp_Update_List::LS):
    {
      for (tuint i=0;i<Ns_final;++i) {
        i_list[i] = 0;
      }
      break;
    }
    default:
      throw FlxException_Crude("FlxBayUp_Update_List::finalize_3");
  };
  curSID = 0;
  if (randomizeID==INIT) {
    curID = seed_idx[curSID];
  } else {
    curID = 0;
  }
  finalized = true;
  return Ns_final;
}

void FlxBayUp_Update_List::reset_finalized()
{
  for (tuint i=0;i<Ns_final;++i) {
    if (i_list[i]==1) {
      i_list[i] = 0;
    }
  }
  curSID = 0;
  if (randomizeID==INIT) {
    curID = seed_idx[curSID];
  } else {
    curID = 0;
  }
}

void FlxBayUp_Update_List::set_next_draw()
{
  #if FLX_DEBUG
    if (i_list[curID]<0) throw FlxException_Crude("FlxBayUp_Update_List::set_next_draw_1");
  #endif
  i_list[curID] = 1;
  switch (randomizeID) {
    case (NONE):
      ++curID;
      if (curID>=Ns_final) curID = 0;
      break;
    case (INIT):
      ++curSID;
      if (curSID>=Ns_final) curSID = 0;
      curID = seed_idx[curSID];
      break;
    case (RNDPICK):
      curID = RndCreator->gen_smp_index(Ns_final);
      break;
    default:
      throw FlxException_Crude("FlxBayUp_Update_List::set_next_draw_2");
  };
  #if FLX_DEBUG
    if (i_list[curID]<0) throw FlxException_Crude("FlxBayUp_Update_List::set_next_draw_3");
  #endif
}

void FlxBayUp_Update_List::fill_sbox(FlxStatBox& sbox)
{
  sbox.clear();
  if (finalized) {
    for (tuint i=0;i<Ns_final;++i) {
      if (i_list[i]>=0) {
        flxVec vt(y_list+i*N_RV,N_RV);
        sbox.add(vt);
      }
    }
  } else {
    for (tuint i=0;i<Ns_final;++i) {
      if (i_list[i]==-1) break;
      if (i_list[i]==2) {
        flxVec vt(y_list+i*N_RV,N_RV);
        sbox.add(vt);
      }
    }
  }
}

void FlxBayUp_Update_List::fill_slist(std::vector< tdouble* >& slist)
{
  if (slist.empty()) {
    slist.reserve((Nc*11)/10);
  }
  slist.clear();
  for (tuint i=0;i<Ns_final;++i) {
    if (i_list[i]==-1) break;
    if (i_list[i]==2) {
      slist.push_back(y_list+i*N_RV);
    }
  }
}

void FlxBayUp_Update_List::print_ext_out(std::ostream& os)
{
  if (ext_out.str().empty()) return;
  os << std::endl << "            " << ext_out.str();
  ext_out.clear();
  ext_out.str("");
}

const tdouble FlxBayUp_Update_List::get_velo(const tuint NcNow) const
{
  const tuint Ns = (is_fullList()?Ns_final:Ncl*Nc);
  pdouble v;
  tuint n = 0;
  if (NcNow==0) {
    for (tuint i=1;i<Ns;++i) {
      v += calc_distance(&(y_list[(i-1)*N_RV]),&(y_list[i*N_RV]),N_RV);
      ++n;
      #if FLX_DEBUG
        if (i_list[i]!=1) throw FlxException_Crude("FlxBayUp_Update_List::get_velo_1");
      #endif
    }
  } else {
    tuint i = 0;
    while (true) {        // set entry to first non-seed
      if (i_list[i]==1) break;
      ++i;
      if (i>=Ns) break;
    }
    for (tuint j=0;j<NcNow;++j) {
      const tuint nCl = chain_length[j];
      if (nCl<=1) continue;
      tuint i_prev = seed_idx[j];
      for (tuint k=1;k<nCl;++k) {
        #if FLX_DEBUG
          if (i>=Ns) throw FlxException_Crude("FlxBayUp_Update_List::get_velo_2");
          if (i_list[i]!=1) throw FlxException_Crude("FlxBayUp_Update_List::get_velo_3");
        #endif
        v += calc_distance(&(y_list[i_prev*N_RV]),&(y_list[i*N_RV]),N_RV);
        ++n;
        // find the next entry that is not a seed
        i_prev = i;
        while (true) {
          ++i;
          if (i>=Ns) break;
          if (i_list[i]==1) break;
        }
      }
    }
  }
  #if FLX_DEBUG
    if (NcNow>0 && n!=(Ns-NcNow)) {
      throw FlxException_Crude("FlxBayUp_Update_List::get_velo_4");
    }
  #endif
  return v.cast2double()/n/sqrt(tdouble(N_RV));
}

FlxBayUP_csm_base::FlxBayUP_csm_base(FlxRndCreator& RndCreator, FlxFunction* h_fun)
 : RndCreator(RndCreator), h_fun(h_fun), adpt_velo(NULL), lastVeloInfo(acvelo_dim)
{
  
}

void FlxBayUP_csm_base::register_adpt_ctrl(flxBayUp_adaptive_ctrl_base* adpt_ctrl)
{
  adpt_velo = dynamic_cast<flxBayUp_adaptive_ctrl_velo*>(adpt_ctrl);
}

void FlxBayUP_csm_base::acceptance_feedback(const bool was_accepted)
{
  if (adpt_velo) {
    lastVeloInfo[2] = was_accepted?ONE:ZERO;
    adpt_velo->append_smpl(lastVeloInfo);
  }
}

FlxBayUP_csm_cwmh_MCMC::FlxBayUP_csm_cwmh_MCMC(FlxRndCreator& RndCreator, const std::string& kernelName, const tdouble h_value, FlxFunction* h_fun)
: FlxBayUP_csm_base(RndCreator,h_fun), kernel(FlxRndKernel_base::read(kernelName,false)), N1Dacc(0), N1Dtotal(0)
{
  kernel->set_h(h_value);
}

void FlxBayUP_csm_cwmh_MCMC::prepare()
{
  N1Dacc = 0;
  N1Dtotal = 0;
  if (h_fun) {
    const tdouble hval = h_fun->cast2positive();
    (*data->ConstantBox.get("sus_kernel_h",true)) = hval;
    kernel->set_h(hval);
  }
}

const bool FlxBayUP_csm_cwmh_MCMC::propose(flxVec& y_prop, const flxVec& y_prev)
{
  if (adpt_velo) {
    const tdouble cur_sd = adpt_velo->get_working_sd();
    kernel->set_h(cur_sd);
    lastVeloInfo[0] = cur_sd;
  }
  const tdouble* yp_prev_p = y_prev.get_tmp_vptr_const();
  tdouble* yp_prop_p = y_prop.get_tmp_vptr();
  const tuint N_RV = y_prev.get_N();
  RndCreator.gen_smp(y_prop);
  bool diff = false;        // true if sample is different from seed
  tdouble velo2 = ZERO;
  // for each coordinate ...
    for (tuint i=0;i<N_RV;++i) {
      const tdouble xi = kernel->transform_y2x(yp_prop_p[i])+yp_prev_p[i];
      const tdouble r = rv_phi(xi)/rv_phi(yp_prev_p[i]);
      if (r>=ONE) {
        yp_prop_p[i] = xi;
        diff = true;
        velo2 += pow2(yp_prop_p[i]-yp_prev_p[i]);
        ++N1Dacc;
      } else {
        const tdouble pr = RndCreator.gen_smp_uniform();
        if (pr<=r) {
          yp_prop_p[i] = xi;
          diff = true;
          velo2 += pow2(yp_prop_p[i]-yp_prev_p[i]);
          ++N1Dacc;
        } else {
          yp_prop_p[i] = yp_prev_p[i];
        }
      }
    }
  lastVeloInfo[1] = velo2;
  N1Dtotal += N_RV;
  return diff;
}

void FlxBayUP_csm_cwmh_MCMC::adptv_spread_multiply(const tdouble f)
{
  const tdouble hval = kernel->get_h()*f;
  (*data->ConstantBox.get("sus_kernel_h",true)) = hval;
  kernel->set_h(hval);
}

const std::string FlxBayUP_csm_cwmh_MCMC::print_info()
{
  std::ostringstream ssV;
  ssV << "component-wise Metropolis-Hastings" << std::endl;
  ssV << "  Type of the sampling kernel:  "; 
  kernel->print_info(ssV);
  if (h_fun) {
    if (h_fun->dependOn_Const(data->ConstantBox.get("sus_iter"))) {
      ssV << std::endl << "                     kernel_h:  " << h_fun->write();
    }
  }
  return ssV.str();
}

void FlxBayUP_csm_cwmh_MCMC::write_adaptive_info(std::ostream& sout, const bool is_adaptive)
{
  if (!is_adaptive && h_fun==NULL) return;
  if (adpt_velo) {
    adpt_velo->write_adaptive_info(sout);
    return;
  }
  sout << std::format("  h={:4.2f}   ", kernel->get_h());
}

FlxBayUP_csm_cov_MCMC::FlxBayUP_csm_cov_MCMC(FlxRndCreator& RndCreator, const tuint M, const std::string& kernelName, const tdouble h_value, FlxFunction* h_fun, const tdouble p, const tuint Nmax, const tdouble p_single, const tuint Nmax_single, FlxBayUp_Update_List& listV)
: FlxBayUP_csm_base(RndCreator,h_fun), M(M), h(h_value), p(p), Nmax(Nmax), p_single(p_single), Nmax_single(Nmax_single),
  sde(M), hvN1(listV.get_Ns_final()), hvN2(listV.get_Ns_final()), hvM(M), ivN(listV.get_Ns_final()),
  covmtx(M), Tinv(M,M), sbox(listV.get_Ns_final(),M), list(listV), kernel(FlxRndKernel_base::read(kernelName,false)),
  N1Dacc(0), N1Dtotal(0)
{
  for (tuint i = 0; i < M; i++) {
    flxVec Vh(M);
    Eigenvectors.push_back(Vh);
  }
}

void FlxBayUP_csm_cov_MCMC::prepare()
{
  N1Dacc = 0;
  N1Dtotal = 0;
  if (h_fun) {
    h = h_fun->cast2positive();
    (*data->ConstantBox.get("sus_kernel_h",true)) = h;
  }
  list.fill_sbox(sbox);  
  if (p<=ONE) return;
  // solve the cov-problem only once (at each conditioning step) for all samples
    sbox.get_smpl_eigen(covmtx,sde,Eigenvectors,Tinv);
}

void FlxBayUP_csm_cov_MCMC::prepareS(const flxVec& y)
{
  if (p<=ONE) {
    sde = y;
    sbox.get_smpl_eigen(p,Nmax,sde,covmtx,hvN1,hvN2,hvM,ivN,Eigenvectors,Tinv,p_single,Nmax_single);
  }
}

const bool FlxBayUP_csm_cov_MCMC::propose(flxVec& y_prop, const flxVec& y_prev)
{
  if (adpt_velo) throw FlxException_NotImplemented("FlxBayUP_csm_cov_MCMC::propose");
  tdouble* yp_prop_p = y_prop.get_tmp_vptr();
  tdouble* hvMp = hvM.get_tmp_vptr();
  const tdouble* sdep = sde.get_tmp_vptr_const();
  // rotate y_prev to eigen-space
    for (tuint i=0;i<M;++i) {
      hvMp[i] = Eigenvectors[i].operator*(y_prev);
    }
  RndCreator.gen_smp(y_prop);
  bool diff = false;        // true if sample is different from seed
  // for each coordinate ... (we work in eigen-space)
    for (tuint i=0;i<M;++i) {
      kernel->set_h_fast(sdep[i]*h);
      const tdouble xi = kernel->transform_y2x(yp_prop_p[i])+hvMp[i];
      const tdouble r = rv_phi(xi)/rv_phi(hvMp[i]);
      if (r>=ONE) {
        yp_prop_p[i] = xi;
        diff = true;
        ++N1Dacc;
      } else {
        const tdouble pr = RndCreator.gen_smp_uniform();
        if (pr<=r) {
          yp_prop_p[i] = xi;
          diff = true;
          ++N1Dacc;
        } else {
          yp_prop_p[i] = hvMp[i];
        }
      }
    }
    N1Dtotal += M;
  // rotate y_prop back to standard normal space (from eigen-space)
    hvM = y_prop;
    Tinv.MultMv(hvM,y_prop);
  return diff;
}

void FlxBayUP_csm_cov_MCMC::adptv_spread_multiply(const tdouble f)
{
  h *= f;
  (*data->ConstantBox.get("sus_kernel_h",true)) = h;
}

const std::string FlxBayUP_csm_cov_MCMC::print_info()
{
  std::ostringstream ssV;
  ssV << "sample-covariance Metropolis-Hastings" << std::endl;
  ssV << "  Type of the sampling kernel:  "; 
  kernel->set_h(h);
  kernel->print_info(ssV);
  if (h_fun) {
    if (h_fun->dependOn_Const(data->ConstantBox.get("sus_iter"))) {
      ssV << std::endl << "                     kernel_h:  " << h_fun->write();
    }
  }
  ssV << std::endl;
  ssV << "  fraction of samples for COV:  ";
  tuint Nd = 0;
  if (p<=ONE) {
    ssV << "p=" << GlobalVar.Double2String(p) << "; Npmax=" << GlobalVar.Double2String(Nmax)
        << "; p_single=" << GlobalVar.Double2String(p_single) << "; Nmax_single=" << GlobalVar.Double2String(Nmax_single);
      Nd = list.get_Nc()*p;
      if (Nd>Nmax) Nd = Nmax;
  } else {
    ssV << "  all samples are considered.";
    Nd = list.get_Nc();
  }
  // ratio given dimension
    const tdouble fdim = tdouble(Nd)/list.get_Nrv();
      ssV << " (fdim=" << GlobalVar.Double2String(fdim) << ")";
    if (fdim<=ONE) {
      GlobalVar.alert.alert("FlxBayUP_csm_cov_MCMC::print_info","The number of cov-samples should be chosen such that 'fdim' becomes larger than 1.0");
    }
  return ssV.str();
}

void FlxBayUP_csm_cov_MCMC::write_adaptive_info(std::ostream& sout, const bool is_adaptive)
{
  if (!is_adaptive && h_fun==NULL) return;
  sout << std::format("  h=%4.2f   ", h);
}

FlxBayUP_csm_csus_MCMC::FlxBayUP_csm_csus_MCMC(FlxRndCreator& RndCreator, const tdouble sD, FlxFunction* h_fun)
: FlxBayUP_csm_base(RndCreator,h_fun), rho(ZERO), sD(sD), adpt_ctrl(NULL), lastS(acopti_jump_dim)
{
  if (sD > ONE) {
    std::ostringstream ssV;
    ssV << "'kernel_h' must be within the interval ]0;1]; and not '" << GlobalVar.Double2String(sD) << "'.";
    throw FlxException_NeglectInInteractive("FlxBayUP_csm_csus_MCMC::FlxBayUP_csm_csus_MCMC", ssV.str() ); 
  }
  rho = sqrt(ONE-pow2(sD));
}

void FlxBayUP_csm_csus_MCMC::register_adpt_ctrl(flxBayUp_adaptive_ctrl_base* adpt_ctrlV)
{
  FlxBayUP_csm_base::register_adpt_ctrl(adpt_ctrlV);
  adpt_ctrl = dynamic_cast<flxBayUp_adaptive_ctrl_opti_jump*>(adpt_ctrlV);
}

void FlxBayUP_csm_csus_MCMC::prepare()
{
  if (h_fun) {
    sD = h_fun->cast2positive();
    (*data->ConstantBox.get("sus_kernel_h",true)) = sD;
    if (sD > ONE) {
      std::ostringstream ssV;
      ssV << "'kernel_h' must be within the interval ]0;1]; and not '" << GlobalVar.Double2String(sD) << "'.";
      throw FlxException_NeglectInInteractive("FlxBayUP_csm_csus_MCMC::prepare", ssV.str() ); 
    }
    rho = sqrt(ONE-pow2(sD));
  }
}

const bool FlxBayUP_csm_csus_MCMC::propose(flxVec& y_prop, const flxVec& y_prev)
{
  const tdouble cur_sd = (adpt_velo)?(adpt_velo->get_working_sd()):sD;
  const tdouble cur_rho = (adpt_velo)?(sqrt(ONE-pow2(cur_sd))):rho;
  const tdouble* yp_prev_p = y_prev.get_tmp_vptr_const();
  tdouble* yp_prop_p = y_prop.get_tmp_vptr();
  const tuint N_RV = y_prev.get_N();
  RndCreator.gen_smp(y_prop);
  tdouble ss = ZERO;
  tdouble vecdot = ZERO;
  tdouble vecl2 = ZERO;
  for (tuint i=0;i<N_RV;++i) {
    const tdouble hv = yp_prop_p[i]*cur_sd;
    vecl2 += pow2(hv);
    vecdot += hv*yp_prev_p[i];
    yp_prop_p[i] = hv+yp_prev_p[i]*cur_rho;
    ss += pow2(yp_prop_p[i]-yp_prev_p[i]);
  }
  if (adpt_ctrl) {
    const tdouble ls2 = y_prev.get_Norm2_NOroot();
    lastS[0] = sqrt(ls2);                                        // length of seed
    lastS[1] = vecdot/lastS[0] + cur_rho*lastS[0];                // a
    lastS[2] = vecl2-pow2(vecdot)/ls2;        // b^2
    lastS[3] = cur_sd;                                                // spread*
    lastS[4] = ZERO;
    lastS[5] = y_prev.comp_dist_NOroot(y_prop);                        // ESJD
    #if FLX_DEBUG
      const tdouble l2 = pow2(lastS[0]*sqrt(ONE-pow2(cur_sd))-lastS[1])+lastS[2];
      if (fabs(y_prop.get_Norm2_NOroot()-(pow2(lastS[1])+lastS[2]))>1e-6) throw FlxException_Crude("FlxBayUP_csm_csus_MCMC::propose_1");
      if (fabs(l2-vecl2)>1e-6) throw FlxException_Crude("FlxBayUP_csm_csus_MCMC::propose_2");
    #endif
  }
  lastVeloInfo[0] = cur_sd;
  lastVeloInfo[1] = ss;
  return true;
}

void FlxBayUP_csm_csus_MCMC::acceptance_feedback(const bool was_accepted)
{
  if (adpt_ctrl) {
    if (was_accepted) lastS[4] = ONE;
    adpt_ctrl->append_smpl(lastS);
  }
  else FlxBayUP_csm_base::acceptance_feedback(was_accepted);
}

void FlxBayUP_csm_csus_MCMC::adptv_spread_multiply(const tdouble f)
{
  sD *= f;
  if (sD>ONE) sD = ONE;
  (*data->ConstantBox.get("sus_kernel_h",true)) = sD;
  rho = sqrt(ONE-pow2(sD));
}

const std::string FlxBayUP_csm_csus_MCMC::print_info()
{
  std::ostringstream ssV;
  ssV << "simulation of conditional samples in U-space; sd=";
  if (h_fun) {
    if (h_fun->dependOn_Const(data->ConstantBox.get("sus_iter"))) {
      ssV << h_fun->write();
    } else {
      ssV << GlobalVar.Double2String(sD);
    }
  } else {
    ssV << GlobalVar.Double2String(sD);
  }
  return ssV.str();
}

void FlxBayUP_csm_csus_MCMC::write_adaptive_info(std::ostream& sout, const bool is_adaptive)
{
  if (!is_adaptive && h_fun==NULL) return;
  sout << std::format("  h={:4.2f}   ", sD);
  if (adpt_ctrl) adpt_ctrl->write_adaptive_info(sout);
}

void FlxBayUP_csm_csus_MCMC::get_cur_spread(tdouble& sdV) const
{
  sdV = sD;
}

void FlxBayUP_csm_csus_MCMC::set_cur_spread(const tdouble& sdV)
{
  sD = sdV;
    rho = sqrt(ONE-pow2(sD));
}


FlxBayUP_csm_dcs_MCMC::FlxBayUP_csm_dcs_MCMC(FlxRndCreator& RndCreator, const tdouble sdV, const tdouble pSD, FlxFunction* h_fun, FlxBayUp_Update_List& list)
: FlxBayUP_csm_base(RndCreator,h_fun), rhoR(ZERO), sdR(sdV), sdW(sdR), sdWS(sdR), pSD(pSD), adpt_ctrl(NULL), lastS(acdcs_dim), list(list)
{
  if (sdR > ONE) {
    std::ostringstream ssV;
    ssV << "'kernel_h' must be within the interval ]0;1]; and not '" << GlobalVar.Double2String(sdR) << "'.";
    throw FlxException_NeglectInInteractive("FlxBayUP_csm_dcs_MCMC::FlxBayUP_csm_dcs_MCMC", ssV.str() ); 
  }
  rhoR = sqrt(ONE-pow2(sdR));
}

void FlxBayUP_csm_dcs_MCMC::register_adpt_ctrl(flxBayUp_adaptive_ctrl_base* adpt_ctrlV)
{
  FlxBayUP_csm_base::register_adpt_ctrl(adpt_ctrlV);
  adpt_ctrl = dynamic_cast<flxBayUp_adaptive_ctrl_dcs*>(adpt_ctrlV);
}

void FlxBayUP_csm_dcs_MCMC::prepare()
{
  list.fill_slist(seedVec);
  if (h_fun) {
    sdR = h_fun->cast2positive();
    (*data->ConstantBox.get("sus_kernel_h",true)) = sdR;
    if (sdR > ONE) {
      std::ostringstream ssV;
      ssV << "'kernel_h' must be within the interval ]0;1]; and not '" << GlobalVar.Double2String(sdR) << "'.";
      throw FlxException_NeglectInInteractive("FlxBayUP_csm_dcs_MCMC::prepare", ssV.str() ); 
    }
    rhoR = sqrt(ONE-pow2(sdR));
    sdW = sdR;
  }
}

const bool FlxBayUP_csm_dcs_MCMC::propose(flxVec& y_prop, const flxVec& y_prev)
{
  if (adpt_velo) throw FlxException_NotImplemented("FlxBayUP_csm_dcs_MCMC::propose");
  const tdouble* yp_prev_p = y_prev.get_tmp_vptr_const();
  tdouble* yp_prop_p = y_prop.get_tmp_vptr();
  const tuint N_RV = y_prev.get_N();
  const tdouble ls2 = y_prev.get_Norm2_NOroot();        // length^2 of seed-vector
  // generate randomly oriented hyperplane
    tdouble cosw = ZERO;
    tdouble l2 = ZERO;                // length^2 of generated sample
    tdouble sdW_help = sdW;
    if (RndCreator.gen_smp_uniform()<pSD) {
      do {
        const tuint si = RndCreator.gen_smp_index(seedVec.size());
        y_prop = flxVec(seedVec[si],N_RV);
        l2 = y_prop.get_Norm2_NOroot();
        cosw = y_prop.operator*(y_prev)/sqrt(l2*ls2);        // compute angle between y_prop and y_prev
      } while (fabs(ONE-fabs(cosw))<=GlobalVar.TOL() || l2<=GlobalVar.TOL() );
      lastS[9] = ONE;
      sdW_help = sdWS;
    } else {
      do {
        RndCreator.gen_smp(y_prop);                        // generate N_RV independent standard Normal realizations
        l2 = y_prop.get_Norm2_NOroot();
        cosw = y_prop.operator*(y_prev)/sqrt(l2*ls2);        // compute angle between y_prop and y_prev
      } while (fabs(ONE-fabs(cosw))<=GlobalVar.TOL() || l2<=GlobalVar.TOL() );
      lastS[9] = ZERO;
    }
  // find vector in hyperplane that is perpendicular to y_prev
    {
      const tdouble bs = y_prop.operator*(y_prev)/ls2;
      for (tuint i=0;i<N_RV;++i) {
        yp_prop_p[i] -= bs*yp_prev_p[i];
      }
      l2 = y_prop.get_Norm2();
      y_prop /= l2;
    }
  // sample vector in plane that has conditional angle 
    {
      const tdouble uW = RndCreator.gen_smp();
      const tdouble yW = sdW_help*uW;
      lastS[1] = yW;
      lastS[2] = sdR;
      lastS[3] = sdW_help;
      const tdouble omega = flxerf(yW/sqrt(2*ONE))*PI;                // target angle
      const tdouble cosw = cos(omega);
      lastS[7] = cosw;
      const tdouble sinw = sin(omega);
      y_prop *= sinw;
      const tdouble bs = cosw/sqrt(ls2);
      for (tuint i=0;i<N_RV;++i) {
        yp_prop_p[i] += bs*yp_prev_p[i];
      }
    }
  // scale vector to target length
    {
      const tdouble N_RVh = tdouble(N_RV)/2;
      const tdouble uR = RndCreator.gen_smp();
      const tdouble uRs = (ls2<=N_RV)?(rv_InvPhi(flxgamma_rl(N_RVh,ls2/2))):(-rv_InvPhi(flxgamma_ru(N_RVh,ls2/2)));
      const tdouble yR = rhoR*uRs + sdR*uR;
      const tdouble lt = sqrt((yR<=ZERO)?(flxgamma_rl_inv(N_RVh,rv_Phi(yR))*2):(flxgamma_ru_inv(N_RVh,rv_Phi(-yR))*2));
      lastS[0] = uR;                // u-transform of yR
      lastS[4] = sqrt(ls2);
      lastS[5] = uRs;
      lastS[6] = pow2(lt-lastS[4]);
      y_prop *= lt;
    }
  // obtain step-size
    tdouble ss = ZERO;
    for (tuint i=0;i<N_RV;++i) {
      ss += pow2(yp_prop_p[i]-yp_prev_p[i]);
    }
    lastS[8] = ss;
  return true;
}

void FlxBayUP_csm_dcs_MCMC::acceptance_feedback(const bool was_accepted)
{
  if (adpt_ctrl==NULL) return;
  adpt_ctrl->append_smpl(lastS,was_accepted);
}

void FlxBayUP_csm_dcs_MCMC::adptv_spread_multiply(const tdouble f)
{
  if (adpt_ctrl==NULL) {
    sdR *= f;
    if (sdR>ONE) sdR = ONE;
    (*data->ConstantBox.get("sus_kernel_h",true)) = sdR;
    rhoR = sqrt(ONE-pow2(sdR));
    sdW *= f;
    if (sdW>ONE) sdW = ONE;
  }
}

const std::string FlxBayUP_csm_dcs_MCMC::print_info()
{
  if (adpt_ctrl) {
    adpt_ctrl->register_csm(this);
  }
  std::ostringstream ssV;
  ssV << "directional conditional sampling; sdR=";
  if (h_fun) {
    if (h_fun->dependOn_Const(data->ConstantBox.get("sus_iter"))) {
      ssV << h_fun->write();
    } else {
      ssV << GlobalVar.Double2String(sdR);
    }
  } else {
    ssV << GlobalVar.Double2String(sdR);
  }
  ssV << "; sdW=";
  if (h_fun) {
    ssV << "sdR";
  } else {
    ssV << GlobalVar.Double2String(sdW);
  }
  ssV << "; pSD=" << GlobalVar.Double2String(pSD);
  return ssV.str();
}

void FlxBayUP_csm_dcs_MCMC::write_adaptive_info(std::ostream& sout, const bool is_adaptive)
{
  if (adpt_ctrl==NULL) throw FlxException_Crude("FlxBayUP_csm_dcs_MCMC::write_adaptive_info");
  if (!is_adaptive && h_fun==NULL) return;
  adpt_ctrl->write_adaptive_info(sout);
}

void FlxBayUP_csm_dcs_MCMC::get_cur_spread(tdouble& sdRV, tdouble& sdWV, tdouble& sdWSV, tdouble& pSDV) const
{
  sdRV = sdR;
  sdWV = sdW;
  sdWSV = sdWS;
  pSDV = pSD;
}

void FlxBayUP_csm_dcs_MCMC::set_cur_spread(const tdouble& sdRV, const tdouble& sdWV, const tdouble& sdWSV, const tdouble& pSDV)
{
  sdR = sdRV;
    rhoR = sqrt(ONE-pow2(sdR));
  sdW = sdWV;
  sdWS = sdWSV;
  pSD = pSDV;
}

FlxBayUP_csm_TMCMC::FlxBayUP_csm_TMCMC(FlxRndCreator& RndCreator, const tuint M, const tdouble beta, FlxFunction* h_fun)
: FlxBayUP_csm_base(RndCreator, h_fun), beta(beta), Acov(M), is_set(false)
{

}

void FlxBayUP_csm_TMCMC::prepare(const bool is_posterior)
{
  if (h_fun) {
    beta = h_fun->cast2positive();
  }
}

void FlxBayUP_csm_TMCMC::prepare(const flxVec& covMtx)
{
  prepare();
  Acov.CholeskyDec(covMtx);
  is_set = true;
}

const bool FlxBayUP_csm_TMCMC::propose(flxVec& y_prop, const flxVec& y_prev)
{
  if (adpt_velo) throw FlxException_NotImplemented("FlxBayUP_csm_TMCMC::propose");
  RndCreator.gen_smp(y_prop);
  Acov.MultMv(y_prop,y_prop);
  y_prop *= beta;
  y_prop += y_prev;
  return true;
}

void FlxBayUP_csm_TMCMC::adptv_spread_multiply(const tdouble f)
{
  beta *= f;
}

const std::string FlxBayUP_csm_TMCMC::print_info()
{
  std::ostringstream ssV;
  ssV << "sample-covariance proposal for TMCMC (beta=" << (h_fun?(h_fun->write()):(GlobalVar.Double2String(beta))) << ")" << std::endl;
  return ssV.str();
}

void FlxBayUP_csm_TMCMC::write_adaptive_info(std::ostream& sout, const bool is_adaptive)
{
  if (!is_adaptive && h_fun==NULL) return;
  sout << std::format("  h={:4.2f}   ", beta);
}

Flx_SuS_CLevelStat::Flx_SuS_CLevelStat(const tuint level, Flx_SuS_CLevelStat* prev)
 : level(level),maxLevelDepth(get_MaxLevelDepth()), prev(prev),
    pi(ZERO), Nsamples(0), Nchains(0), Nfailures(0), g_t(ZERO), 
    gamma(ZERO), gamma_chain(ZERO), gamma_from_seed(ZERO), g_Ncl_max(0), g_seed_ID_original(NULL), g_seed_ID(NULL), g_chain_length(NULL), g_N(NULL), g_N_entries(0), pi_chain(NULL), p_00(ZERO), p_11(ZERO), lag1_corr(ZERO), g_seed_corrE(NULL), eff(ONE), eff_Gelman(ONE),
    corr_pi_prev(ZERO), corr_pi_prev_negFrac(ZERO), frac_rel_chains(ZERO), seed_chainID(NULL), seed_chainPos(NULL), find_multiples(NULL)
{

}

Flx_SuS_CLevelStat::~Flx_SuS_CLevelStat()
{
  if (g_N) delete [] g_N;
  if (g_seed_ID_original) delete [] g_seed_ID_original;
  if (g_seed_ID) delete [] g_seed_ID;
  if (pi_chain) delete [] pi_chain;
  if (g_chain_length) delete [] g_chain_length;
  if (seed_chainID) delete [] seed_chainID;
  if (seed_chainPos) delete [] seed_chainPos;
  if (find_multiples) delete [] find_multiples;
}

Flx_SuS_CLevelStat::Flx_SuS_CLevelStat(const Flx_SuS_CLevelStat& rhs)
: level(rhs.level),maxLevelDepth(get_MaxLevelDepth()), prev(rhs.prev),
  pi(rhs.pi), Nsamples(rhs.Nsamples), Nchains(rhs.Nchains), Nfailures(rhs.Nfailures),g_t(rhs.g_t), gamma(rhs.gamma),gamma_chain(rhs.gamma_chain), gamma_from_seed(rhs.gamma_from_seed), g_Ncl_max(rhs.g_Ncl_max),
  g_seed_ID_original(NULL), g_seed_ID(NULL), g_chain_length(NULL), g_N(NULL), g_N_entries(0), pi_chain(NULL),p_00(rhs.p_00),p_11(rhs.p_11),lag1_corr(rhs.lag1_corr), g_seed_corrE(NULL),
  eff(rhs.eff), eff_Gelman(rhs.eff_Gelman),
  corr_pi_prev(rhs.corr_pi_prev), corr_pi_prev_negFrac(rhs.corr_pi_prev_negFrac), frac_rel_chains(rhs.frac_rel_chains), seed_chainID(NULL), seed_chainPos(NULL), find_multiples(NULL)
{
  if (rhs.g_N) throw FlxException_Crude("Flx_SuS_CLevelStat::Flx_SuS_CLevelStat");
}

Flx_SuS_CLevelStat& Flx_SuS_CLevelStat::operator=(const Flx_SuS_CLevelStat& rhs)
{
  if (rhs.g_N) throw FlxException_Crude("Flx_SuS_CLevelStat::operator=_1");
  level = rhs.level;
  maxLevelDepth = rhs.maxLevelDepth;
  #if FLX_DEBUG
    if (maxLevelDepth!=get_MaxLevelDepth()) throw FlxException_Crude("Flx_SuS_CLevelStat::operator=_2");
  #endif
  prev = rhs.prev;
  pi = rhs.pi;
  Nsamples = rhs.Nsamples;
  Nchains = rhs.Nchains;
  Nfailures = rhs.Nfailures;
  g_t = rhs.g_t;
  gamma = rhs.gamma;
  gamma_chain = rhs.gamma_chain;
  gamma_from_seed = rhs.gamma_from_seed;
  g_Ncl_max = rhs.g_Ncl_max;
  if (g_seed_ID_original) delete [] g_seed_ID_original;
    g_seed_ID_original = NULL;
  if (g_seed_ID) delete [] g_seed_ID;
    g_seed_ID = NULL;
  if (g_chain_length) delete [] g_chain_length;
    g_chain_length = NULL;
  if (g_N) delete [] g_N;
    g_N = NULL;
  g_N_entries = 0;
  if (pi_chain) delete [] pi_chain;
    pi_chain = NULL;
  p_00 = rhs.p_00;
  p_11 = rhs.p_11;
  lag1_corr = rhs.lag1_corr;
  #if FLX_DEBUG
    if (g_seed_corrE) throw FlxException_Crude("Flx_SuS_CLevelStat::operator=_3");
  #endif
  eff = rhs.eff;
  eff_Gelman = rhs.eff_Gelman;
  corr_pi_prev = rhs.corr_pi_prev;
  corr_pi_prev_negFrac = rhs.corr_pi_prev_negFrac;
  frac_rel_chains = rhs.frac_rel_chains;
  if (seed_chainID) delete [] seed_chainID;
    seed_chainID = NULL;
  if (seed_chainPos) delete [] seed_chainPos;
    seed_chainPos = NULL;
  return *this;
}

const tuint Flx_SuS_CLevelStat::get_MaxLevelDepth() const
{
  const tuint mld = 10;                // the maximum depth
  if (level<=1) return 0;
  const tuint ld = level-1;
  if (ld<=mld) return ld;
  return mld;
}

const tuint Flx_SuS_CLevelStat::get_seed_N_groups() const
{
  const tuint mld = get_MaxLevelDepth();
  if (mld==0) return 0;
  return mld+2;
}

const tuint Flx_SuS_CLevelStat::get_seed_group(const tuint delta_level, const tuint delta_pos) const
{
  #if FLX_DEBUG
    if (delta_level<1) throw FlxException_Crude("Flx_SuS_CLevelStat::get_seed_group");
  #endif
  if (delta_level==1) {
    if (delta_pos<=1) return 0;
    else if (delta_pos<=3) return 1;
    else return 2;
  } else {
    return delta_level+1;
  }
}

const tuint Flx_SuS_CLevelStat::get_seed_group_size(tuint depth) const
{
  return depth*depth;
}

void Flx_SuS_CLevelStat::allocate_g_seed_corrE(const bool isSeedCorr)
{
  const tuint gS = isSeedCorr?get_seed_N_groups():get_pic_N_groups();
  g_seed_corrE = new tuint*[gS];
  for (tuint i=0;i<gS;++i) {
    const tuint mf = (!isSeedCorr&&i==0)?3:2;
    const tuint d = isSeedCorr?get_seed_group_depth(i):get_pic_group_depth(i);
    const tuint S = get_seed_group_size(d)*mf + d*4;
    g_seed_corrE[i] = new tuint[S];
    memset(g_seed_corrE[i],0, sizeof(tuint)*S);
  }
}

void Flx_SuS_CLevelStat::deallocate_g_seed_corrE(const bool isSeedCorr)
{
  if (g_seed_corrE) {
    const tuint Ni = isSeedCorr?get_seed_N_groups():get_pic_N_groups();
    for (tuint i=0;i<Ni;++i) {
      delete [] g_seed_corrE[i];
    }
    delete [] g_seed_corrE;
    g_seed_corrE = NULL;
  }
}

const tuint Flx_SuS_CLevelStat::get_pic_N_groups() const
{
  const tuint mld = get_MaxLevelDepth();
//   if (mld>=5) return 5;
  return mld;
}

const tuint Flx_SuS_CLevelStat::get_seed_group_depth(const tuint seed_group) const
{
  const tuint res = (g_Ncl_max>prev->g_Ncl_max)?(prev->g_Ncl_max):g_Ncl_max;
  if (res>8) return 8;
  return res;
}

const tuint Flx_SuS_CLevelStat::get_pic_group_depth(const tuint pic_group) const
{
  const tuint res = (g_Ncl_max>prev->g_Ncl_max)?(prev->g_Ncl_max):g_Ncl_max;
  if (res>6) return 6;
  return res;
}


const bool Flx_SuS_CLevelStat::find_common_seed(const tuint chainID_this, const tuint posInChain_this, const tuint level_rhs, const tuint chainID_rhs, const tuint posInChain_rhs, tuint& delta_level, tuint& delta_pos, tuint mLdp)
{
  if (level_rhs<1) return false;        // perfect MCS!!!
  if (level>level_rhs) {
    #if FLX_DEBUG
      if (prev==NULL) throw FlxException_Crude("Flx_SuS_CLevelStat::find_common_seed_1");
      if (level!=prev->level+1) throw FlxException_Crude("Flx_SuS_CLevelStat::find_common_seed_2");
      if (prev->Nfailures!=Nchains) throw FlxException_Crude("Flx_SuS_CLevelStat::find_common_seed_3");
      if (chainID_this>=Nchains) throw FlxException_Crude("Flx_SuS_CLevelStat::find_common_seed_4");
      if (delta_level!=0) throw FlxException_Crude("Flx_SuS_CLevelStat::find_common_seed_5");
      if (mLdp!=0) throw FlxException_Crude("Flx_SuS_CLevelStat::find_common_seed_6");
    #endif
    const tuint seedID_orig = g_seed_ID[chainID_this];                // the chain stems from the 'i'th seed in previous list
    const tuint chainID_prev = prev->seed_chainID[seedID_orig];
    const tuint posInChain_prev = prev->seed_chainPos[seedID_orig];
    return prev->find_common_seed(chainID_prev,posInChain_prev,level_rhs,chainID_rhs,posInChain_rhs,delta_level,delta_pos);
  } else if (level==level_rhs) {        // level==level_rhs
    if (chainID_this==chainID_rhs) {
      if (posInChain_this>=posInChain_rhs) {
        delta_pos = posInChain_this - posInChain_rhs;
      } else {
        delta_pos = posInChain_rhs - posInChain_this;
      }
      return true;
    } else {
      if (level<=1) return false;
      if (mLdp==0) mLdp = maxLevelDepth;
      if (delta_level>=mLdp) return false;
      ++delta_level;
      // this
        const tuint seedID_orig = g_seed_ID[chainID_this];                // the chain stems from the 'i'th seed in previous list
        const tuint chainID_prev = prev->seed_chainID[seedID_orig];
        const tuint posInChain_prev = prev->seed_chainPos[seedID_orig];
      // rhs
        const tuint rhs_seedID_orig = g_seed_ID[chainID_rhs];                // the chain stems from the 'i'th seed in previous list
        const tuint rhs_chainID_prev = prev->seed_chainID[rhs_seedID_orig];
        const tuint rhs_posInChain_prev = prev->seed_chainPos[rhs_seedID_orig];
      return prev->find_common_seed(chainID_prev,posInChain_prev,level-1,rhs_chainID_prev,rhs_posInChain_prev,delta_level,delta_pos,mLdp);
    }
  } else {
    throw FlxException_Crude("Flx_SuS_CLevelStat::find_common_seed_99");
  }
}

const tuint Flx_SuS_CLevelStat::find_start_in_seed_chainID(const tuint cID) const
{
  tuint ip = 0;
  tuint l = Nfailures;
  while (true) {
    if (l<=1) return ip;
    const tuint mid = ip + l/2;
    const tuint mid_f = seed_chainID[mid];
    if (mid_f==cID) {
      ip = mid;
      while (ip>0) {
        if (cID!=seed_chainID[ip-1]) {
          return ip;
        }
        --ip;
      }
      return ip;
    } else if (mid_f>cID) {
      l = mid-ip;
    } else {
      const tuint end = ip + l;
      ip = mid + 1;
      if (ip>=end) return 0;
      l = end-ip;
    }    
  }
}

void Flx_SuS_CLevelStat::add2seedCorr_2group(const tuint cID1, const tuint cID2, const tuint group, const bool isSeedCorr)
{
  const tuint d = isSeedCorr?(get_seed_group_depth(group)):(get_pic_group_depth(group));
  Flx_SuS_CLevelStat* refp = isSeedCorr?this:prev;
  const tuint cl_cID1 = g_chain_length[cID1];
  const tuint cl_cID2 = refp->g_chain_length[cID2];
  // increase the total sample size
    tuint c = 0;
    tuint* gcE = g_seed_corrE[group];
    for (tuint i=0;i<d;++i) {
      for (tuint j=0;j<d;++j) {
        if (i<cl_cID1&&j<cl_cID2) {
          ++( gcE[c] );
        }
        c += 2;
      }
    }
  const tuint beg1 = find_start_in_seed_chainID(cID1);
  const tuint beg2 = refp->find_start_in_seed_chainID(cID2);
  const tuint l2Nf = refp->Nfailures;
  // add 'samples' (1 hit)
    {
      // count total samples
        tuint a_pos = d*d*2;
        const tuint l1 = (d>cl_cID1)?cl_cID1:d;
        for (tuint i=0;i<l1;++i) {
          ++( gcE[a_pos] );
          a_pos += 2;
        }
      // count hits
        a_pos = d*d*2 + 1;
        tuint i1 = beg1;
        while (cID1==seed_chainID[i1]) {
          const tuint l1p = seed_chainPos[i1];
          if (l1p>=d) break;
          ++( gcE[ a_pos + l1p*2 ] );
          ++i1;
          if (i1>=Nfailures) break;
        }
    }
    {
      // count total samples
        tuint a_pos = d*d*2+d*2;
        const tuint l2 = (d>cl_cID2)?cl_cID2:d;
        for (tuint i=0;i<l2;++i) {
          ++( gcE[a_pos] );
          a_pos += 2;
        }
      // count hits
        a_pos = d*d*2+2*d + 1;
        tuint i2 = beg2;
        while (cID2==refp->seed_chainID[i2]) {
          const tuint l2p = refp->seed_chainPos[i2];
          if (l2p>=d) break;
          ++( gcE[ a_pos + l2p*2 ] );
          ++i2;
          if (i2>=l2Nf) break;
        }
    }
  // add 'samples' (2 hits)
    tuint i1 = beg1;
    while (cID1==seed_chainID[i1]) {
      const tuint l1p = seed_chainPos[i1];
      if (l1p>=d) break;
      const tuint l2 = (d>cl_cID2)?cl_cID2:d;
      tuint a_pos = l1p*d*2;
      tuint i2 = beg2;
      while (cID2==refp->seed_chainID[i2]) {
        const tuint l2p = refp->seed_chainPos[i2];
        if (l2p>=l2) break;
        ++( gcE[ a_pos + l2p*2+1 ] );
        ++i2;
        if (i2>=l2Nf) break;
      }
      ++i1;
      if (i1>=Nfailures) break;
    }
}

void Flx_SuS_CLevelStat::add2piCorr_2group0(const tuint cID1, const tuint cID2, const tuint delta_pos)
{
  const tuint d = get_pic_group_depth(0);
  tuint* gcE = g_seed_corrE[0];
  tuint gcE_size = (d*d+d)/2;
  const tuint cl_cID1 = g_chain_length[cID1];
  const tuint begin_1 = find_start_in_seed_chainID(cID1);
  tuint beg2 = prev->find_start_in_seed_chainID(cID2);
  const tuint p2_before = (delta_pos>=d)?(delta_pos+1-d):0;
  const tuint p2_after = (prev->g_chain_length[cID2]-delta_pos>d)?(delta_pos+d-1):(prev->g_chain_length[cID2]-1);
  for (tuint i2=p2_before;i2<=p2_after;++i2) {
    const tuint l1p = (i2<delta_pos)?(delta_pos-i2):(i2-delta_pos);
    #if FLX_DEBUG
      if (l1p>=d) throw FlxException_Crude("Flx_SuS_CLevelStat::add2piCorr_2group0_1");
    #endif
    const tuint l2 = d-l1p;
    // adjust position of beg2
      bool beg_valid = beg2<prev->Nfailures && cID2==prev->seed_chainID[beg2];
      while (beg_valid && prev->seed_chainPos[beg2]<i2) {
        ++beg2;
        beg_valid = beg2<prev->Nfailures && cID2==prev->seed_chainID[beg2];
      }
    tuint beg1 = begin_1;
    if (beg_valid) {
      beg_valid = (prev->seed_chainPos[beg2]==i2);
    }
    const bool beg_valid_l1 = beg_valid;        // the first indicator is true!!!
    if (beg_valid) {
      beg_valid = beg1<Nfailures && cID1==seed_chainID[beg1];
    }
    tuint* gcE_p = gcE + (gcE_size - (l2*l2+l2)/2)*3;
    for (tuint i1=0;i1<l2;++i1) {
      if (i1>=cl_cID1) break;
      const tuint l2p = i1;
      ++(*gcE_p);
      if (beg_valid) {
        // adjust position of beg1
          while (beg_valid && seed_chainPos[beg1]<l2p) {
            ++beg1;
            beg_valid = beg1<Nfailures && cID1==seed_chainID[beg1];
          }
          #if FLX_DEBUG
            if ( gcE_p==gcE && prev->seed_chainPos[beg2]!=i2 ) throw FlxException_Crude("Flx_SuS_CLevelStat::add2piCorr_2group0_2");
          #endif
        if (beg_valid && seed_chainPos[beg1]==l2p) {
          ++(gcE_p[1]);
        }
      }
      if (beg_valid_l1) {
        ++(gcE_p[2]);
      }
      gcE_p += 3;
    }
  }
}

void Flx_SuS_CLevelStat::add2seedCorr()
{
  const bool localPi = true;        // use global cond. prob. or cond prob. that depends on chain position
  const tdouble NEG_WEIGHT = 0.9;
  if (level<=1) return;
  // allocate memory
    allocate_g_seed_corrE(true);
  // fill memory
    for (tuint k=0;k<Nchains;++k) {
      for (tuint l=k+1;l<Nchains;++l) {
        tuint delta_level = 0;
        tuint delta_pos = 0;
        if ( find_common_seed(k,0,level,l,0,delta_level,delta_pos) ) {
          const tuint s_group = get_seed_group(delta_level,delta_pos);
          add2seedCorr_2group(k,l,s_group,true);
        }
      }
    }
  // evaluate memory
    const tdouble p_i = tdouble(Nfailures)/tdouble(Nsamples);
    const tdouble R0 = p_i*(ONE-p_i);
    const tdouble p_i_2 = pow2(p_i);
    pdouble rs = ZERO;
    const tuint gS = get_seed_N_groups();
    for (tuint k=0;k<gS;++k) {
      const tuint d = get_seed_group_depth(k);
      tuint* gscE = g_seed_corrE[k];
      tuint* gscE_c1 = gscE + d*d*2;
      tuint* gscE_c2 = gscE_c1 + 2*d;
      const tuint gscN = get_seed_group_size(d);
      tuint c = 0;
      for (tuint l=0;l<d;++l) {
        for (tuint m=0;m<d;++m) {
          if (gscE[0]>0) {
            tdouble rh = tdouble(gscE[1])/tdouble(gscE[0]);
            const tdouble pi_l_1 = (gscE_c1[l*2+0]+gscE_c2[l*2+0]>0)?(tdouble(gscE_c1[l*2+1]+gscE_c2[l*2+1])/tdouble(gscE_c1[l*2+0]+gscE_c2[l*2+0])):p_i;
            const tdouble pi_l_2 = (gscE_c1[m*2+0]+gscE_c2[m*2+0]>0)?(tdouble(gscE_c1[m*2+1]+gscE_c2[m*2+1])/tdouble(gscE_c1[m*2+0]+gscE_c2[m*2+0])):p_i;
            const tdouble pi_tmp = localPi?( pi_l_1 * pi_l_2 ):p_i_2;
            rh -= pi_tmp;
            rh /= R0;
            rh *= 2*(tdouble(gscE[0])/tdouble(Nsamples));
            if (rh<ZERO) {
              rh *= NEG_WEIGHT;
            }
            rs += rh;
          }
          gscE += 2;
          ++c;
        }
      }
      if (c!=gscN) throw FlxException_Crude("Flx_SuS_CLevelStat::add2seedCorr_1");
    }
    if (rs.cast2double()>ZERO) gamma_from_seed = rs.cast2double();
  // deallocate memory
    deallocate_g_seed_corrE(true);
}

void Flx_SuS_CLevelStat::add2piCorr()
{
  const bool localPi = true;        // use global cond. prob. or cond prob. that depends on chain position
  const tdouble NEG_WEIGHT = 0.9;
  if (level<=1) return;                // levels are not correlated
  // the weighting factor (needed later)
    const tdouble pi_1 = tdouble(Nfailures)/tdouble(Nsamples);
    const tdouble pi_2 = tdouble(prev->Nfailures)/tdouble(prev->Nsamples);
    const tdouble pi_12 = pi_1 * pi_2;
    const tdouble var_1 = (pi_1*(ONE-pi_1))/(Nsamples*eff);
    const tdouble var_2 = (pi_2*(ONE-pi_2))/(prev->Nsamples*prev->eff);
    const tdouble COV_weight = tdouble(Nsamples)*tdouble(prev->Nsamples) * sqrt( var_1 * var_2 );
  // allocate memory
    allocate_g_seed_corrE(false);
  // fill memory
    const tuint gS = get_pic_N_groups();
    tulong oc = 0;
    #if FLX_DEBUG
      tuint oc2 = 0;
    #endif
    const tuint Nc_prev = prev->Nchains;
    for (tuint i=0;i<Nchains;++i) {
      for (tuint j=0;j<Nc_prev;++j) {
        tuint delta_level = 0;
        tuint delta_pos = 0;
        if ( find_common_seed(i,0,level-1,j,0,delta_level,delta_pos) ) {
          if (delta_level==0) {
            add2piCorr_2group0(i,j,delta_pos);
            #if FLX_DEBUG
              ++oc2;
            #endif
          } else {
            if ( delta_level>=gS ) continue;
            add2seedCorr_2group(i,j,delta_level,false);
          }
          ++oc;
        }
      }
    }
    frac_rel_chains = (tdouble(oc)/Nc_prev)/Nchains;
    #if FLX_DEBUG
      if (oc2!=Nchains) throw FlxException_Crude("Flx_SuS_CLevelStat::add2piCorr_1");
    #endif
  // evaluate memory
    // first group is special
      pdouble rs_0 = ZERO;
      pdouble rs_pos, rs_neg;
      {
        pdouble rs_0_neg = ZERO;
        tuint* gscE = g_seed_corrE[0];
        const tuint g0_depth = get_pic_group_depth(0);
        const tuint g0_size = (g0_depth*g0_depth+g0_depth)/2;
        tuint c = 0;        // counter
        for (tuint k=0;k<g0_depth;++k) {
          for (tuint l=0;l<g0_depth-k;++l) {
            if (gscE[0]>0) {
              pdouble rh = tdouble(gscE[1]);
              const tdouble pi_tmp = localPi?pi_chain[l]:pi_1;
              rh -= pi_tmp*gscE[2];
              if (rh.cast2double() >= ZERO) {
                rs_0 += rh;
              } else {
                rs_0_neg += rh;
              }
            }
            gscE += 3;
            ++c;
          }
        }
        if (c!=g0_size) throw FlxException_Crude("Flx_SuS_CLevelStat::add2piCorr_2");
        // deal with negative contributions
          rs_pos += rs_0;
          rs_neg += rs_0_neg;
          rs_0_neg *= NEG_WEIGHT;
          rs_0 += rs_0_neg;
          if (rs_0.cast2double()<ZERO) {
            rs_0 = ZERO;
          }
      }
    pdouble rs = ZERO;
    for (tuint m=1;m<gS;++m) {
      pdouble rs_n0;
      pdouble rs_n0_neg;
      const tuint d = get_pic_group_depth(m);
      tuint* gscE = g_seed_corrE[m];
      tuint* gscE_c1 = gscE + d*d*2;
      tuint* gscE_c2 = gscE_c1 + 2*d;
      const tuint gscN = get_seed_group_size(d);
      tuint c=0;
      for (tuint k=0;k<d;++k) {
        for (tuint l=0;l<d;++l) {
          if (gscE[0]>0) {
            pdouble rh = tdouble(gscE[1]);
            const tdouble pi_l_1 = (gscE_c1[l*2+0]>0)?(tdouble(gscE_c1[l*2+1])/tdouble(gscE_c1[l*2+0])):pi_1;
            const tdouble pi_l_2 = (gscE_c2[k*2+0]>0)?(tdouble(gscE_c2[k*2+1])/tdouble(gscE_c2[k*2+0])):pi_2;
            const tdouble pi_tmp = localPi?( pi_l_1 * pi_l_2 ):pi_12;
            rh -= pi_tmp*gscE[0];
            if (rh.cast2double() >= ZERO) {
              rs_n0 += rh;
            } else {
              rs_n0_neg += rh;
            }
          }
          gscE += 2;
          ++c;
        }
      }
      if (c!=gscN) throw FlxException_Crude("Flx_SuS_CLevelStat::add2piCorr_3");
      // deal with negative contributions
        rs_pos += rs_n0;
        rs_neg += rs_n0_neg;
        rs_n0_neg *= NEG_WEIGHT;
        rs_n0 += rs_n0_neg;
        if (rs_n0.cast2double()<ZERO) {
          rs_n0 = ZERO;
        }
        rs += rs_n0;
    }
    rs += rs_0;
    // compute negative fraction
      corr_pi_prev_negFrac = fabs(rs_neg.cast2double()/rs_pos.cast2double());
    corr_pi_prev = rs.cast2double() / COV_weight;
  // deallocate memory
    deallocate_g_seed_corrE(false);
}

void Flx_SuS_CLevelStat::empirical_seedGamma(const tuint lid, tdouble* e_l, const tuint mld, const tdouble NsmpX, const tdouble* para)
{
  if (lid>mld) return;
  empirical_seedGamma(lid+1,e_l,mld,NsmpX,para);
  // a measure of not-efficiency
    // 0: if gamma_chain=ZERO  -> uncorrelated
    // 1: if gamma_chain=infinite        -> totally correlated
    const tdouble p_3 = para[2*3+0]+para[2*3+1]*NsmpX+para[2*3+2]*pow2(NsmpX);
    const tdouble p_4 = para[3*3+0]+para[3*3+1]*NsmpX+para[3*3+2]*pow2(NsmpX);
    const tdouble p_5 = para[4*3+0]+para[4*3+1]*NsmpX+para[4*3+2]*pow2(NsmpX);
    const tdouble p_6 = para[5*3+0]+para[5*3+1]*NsmpX+para[5*3+2]*pow2(NsmpX);
    const tdouble not_eff = (ONE-ONE/(ONE+p_3*(gamma_chain+p_4*pow(gamma_chain,p_5))));
  for (tuint k=lid-1;k<mld;++k) {
    e_l[k] *= not_eff;
  }
  for (tuint k=lid;k<mld;++k) {
    e_l[k] = pow(e_l[k],p_6);
  }
}

void Flx_SuS_CLevelStat::empirical_corrLevel(const tuint lid, tdouble* e_l, const tuint mld, const tdouble NsmpX, const tdouble* para)
{
  if (lid>=mld) return;
  empirical_corrLevel(lid+1,e_l,mld,NsmpX,para);
    const tdouble p_3 = para[8*3+0]+para[8*3+1]*NsmpX+para[8*3+2]*pow2(NsmpX);
    const tdouble p_4 = para[9*3+0]+para[9*3+1]*NsmpX+para[9*3+2]*pow2(NsmpX);
    const tdouble p_5 = para[10*3+0]+para[10*3+1]*NsmpX+para[10*3+2]*pow2(NsmpX);
    const tdouble p_6 = para[11*3+0]+para[11*3+1]*NsmpX+para[11*3+2]*pow2(NsmpX);
  const tdouble not_eff = (ONE-ONE/(ONE+p_3*(gamma_chain+p_4*pow(gamma_chain,p_5))));
  for (tuint k=lid-1;k<mld;++k) {
    e_l[k] *= not_eff;
  }
  for (tuint k=lid;k<mld;++k) {
    e_l[k] = pow(e_l[k],p_6);
  }
}

void Flx_SuS_CLevelStat::empirical_Corr()
{
  if (level<=1) return;
  // **************************************************************
  // Retrieve parameters of empirical model
  // **************************************************************
    const tuint paraSize = 36;
    tdouble para[paraSize];
    {
      tuint i=0;
      // seed correlation
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p1_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p1_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p1_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p2_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p2_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p2_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p3_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p3_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p3_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p4_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p4_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p4_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p5_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p5_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p5_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p6_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p6_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_sc_p6_2");
      // level correlation
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p1_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p1_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p1_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p2_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p2_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p2_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p3_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p3_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p3_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p4_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p4_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p4_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p5_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p5_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p5_2");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p6_0");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p6_1");
        para [i++] = data->ConstantBox.getRef("sus_empi_lc_p6_2");
      if (i!=paraSize) throw FlxException_Crude("Flx_SuS_CLevelStat::empirical_Corr_1");
    }
  // **************************************************************
  // Seed Correlation
  // **************************************************************
  const tdouble NsmpX = log10(tdouble(Nsamples))-3*ONE;
  const tuint mld = get_MaxLevelDepth();
  #ifndef FLX_CV_1
    tuint N_l[mld];
  #else
    tuint* N_l = new tuint[mld];
  #endif
    memset(N_l,0, sizeof(tuint)*mld);
  // evaluate seed history
    for (tuint k=0;k<Nchains;++k) {
      for (tuint l=k+1;l<Nchains;++l) {
        tuint delta_level = 0;
        tuint delta_pos = 0;
        if ( find_common_seed(k,0,level,l,0,delta_level,delta_pos) ) {
          ++(N_l[delta_level-1]);
        }
      }
    }
  // compute the fractions
    #ifndef FLX_CV_1
      tdouble f_l[mld];
    #else
      tdouble* f_l = new tdouble[mld];
    #endif
    for (tuint k=0;k<mld;++k) {
      f_l[k] = tdouble(N_l[k])/tdouble(Nchains);
    }
  // compute the contribution of seed correlation to gamma
    #ifndef FLX_CV_1
      tdouble e_l[mld];
    #else
      tdouble* e_l = new tdouble[mld];
    #endif
    for (tuint k=0;k<mld;++k) e_l[k]=ONE;
    prev->empirical_seedGamma(1,e_l,mld,NsmpX,para);
  // compute the seed-gamma
    const tdouble af = para[0*3+0]+para[0*3+1]*NsmpX+para[0*3+2]*pow2(NsmpX);        // reduction factor
    const tdouble pf = para[1*3+0]+para[1*3+1]*NsmpX+para[1*3+2]*pow2(NsmpX);
    tdouble gsc = ZERO;
    for (tuint k=0;k<mld;++k) {
      gsc += af * gamma_chain*pow(e_l[k],pf)*f_l[k];
    }
  gamma_from_seed = gsc;
  // **************************************************************
  // Seed Correlation
  // **************************************************************
    memset(N_l,0, sizeof(tuint)*mld);
  // evaluate seed history
    const tuint Nc_prev = prev->Nchains;
    for (tuint i=0;i<Nchains;++i) {
      for (tuint j=0;j<Nc_prev;++j) {
        tuint delta_level = 0;
        tuint delta_pos = 0;
        if ( find_common_seed(i,0,level-1,j,0,delta_level,delta_pos) ) {
          ++(N_l[delta_level]);
        }
      }
    }
  // compute the fractions
    for (tuint k=0;k<mld;++k) {
      f_l[k] = tdouble(N_l[k])/tdouble(Nchains);
    }
  // compute the contribution of seed correlation to gamma
    for (tuint k=0;k<mld;++k) e_l[k]=ONE;
    prev->empirical_corrLevel(0,e_l,mld,NsmpX,para);
  // compute the seed-gamma
    const tdouble af_l = para[6*3+0]+para[6*3+1]*NsmpX+para[6*3+2]*pow2(NsmpX);        // reduction factor
    const tdouble pf_l = para[7*3+0]+para[7*3+1]*NsmpX+para[7*3+2]*pow2(NsmpX);
    gsc = ZERO;
    for (tuint k=0;k<mld;++k) {
      gsc += af_l * pow(e_l[k],pf_l)*f_l[k];
    }
  if (gsc>0.99) gsc = 0.99;
  corr_pi_prev = gsc;
  #ifdef FLX_CV_1
    delete [] N_l;
    delete [] e_l;
    delete [] f_l;
  #endif
}

Flx_SuS_Control::Flx_SuS_Control(const Flx_SuS_Control& rhs)
: prt_alert(rhs.prt_alert), verbose(rhs.verbose), credEst(rhs.credEst), pc(NULL), N_cred_smpl(rhs.N_cred_smpl), 
  comp_gamma(rhs.comp_gamma), consider_seed_corr(rhs.consider_seed_corr), consider_pi_corr(rhs.consider_pi_corr),
  empirical_corr(rhs.empirical_corr), find_multiples(rhs.find_multiples), os_samples(NULL), 
  TMCMC_target_COV(NULL), TMCMC_update_weights(rhs.TMCMC_update_weights),
  LS_SPNT(NULL), LS_tol(NULL), LS_max_iter(NULL), pa_maxL(NULL)
{
  if (rhs.pc) pc = new FlxMtxConstFun(*(rhs.pc));
  if (rhs.os_samples) os_samples = new FlxString(*(rhs.os_samples));
  if (rhs.TMCMC_target_COV) TMCMC_target_COV = new FlxFunction(*(rhs.TMCMC_target_COV));
  if (rhs.TMCMC_alpha) TMCMC_alpha = new FlxFunction(*(rhs.TMCMC_alpha));
  if (rhs.LS_SPNT) LS_SPNT = new FlxMtxConstFun(*(rhs.LS_SPNT));
  if (rhs.LS_tol) LS_tol = new FlxFunction(*(rhs.LS_tol));
  if (rhs.LS_max_iter) LS_max_iter = new FlxFunction(*(rhs.LS_max_iter));
  if (rhs.pa_maxL) pa_maxL = new FlxFunction(*(rhs.pa_maxL));
}

Flx_SuS_Control::~Flx_SuS_Control()
{
  if (pc) delete pc;
  if (os_samples) delete os_samples;
  if (TMCMC_target_COV) delete TMCMC_target_COV;
  if (TMCMC_alpha) delete TMCMC_alpha;
  if (LS_SPNT) delete LS_SPNT;
  if (LS_tol) delete LS_tol;
  if (LS_max_iter) delete LS_max_iter;
  if (pa_maxL) delete pa_maxL;
}

Flx_SuS_Control& Flx_SuS_Control::operator=(const Flx_SuS_Control& rhs)
{
  if ( this != &rhs ) {
    prt_alert = rhs.prt_alert;
    verbose = rhs.verbose;
    credEst = rhs.credEst;
    if (pc) delete pc;
      if (rhs.pc) {
        pc = new FlxMtxConstFun(*(rhs.pc));
      } else {
        pc = NULL;
      }
    N_cred_smpl = rhs.N_cred_smpl;
    comp_gamma = rhs.comp_gamma;
    consider_seed_corr = rhs.consider_seed_corr;
    consider_pi_corr = rhs.consider_pi_corr;
    empirical_corr = rhs.empirical_corr;
    find_multiples = rhs.find_multiples;
    if (os_samples) delete os_samples;
      if (rhs.os_samples) {
        os_samples = new FlxString(*(rhs.os_samples));
      } else {
        os_samples = NULL;
      }
    if (TMCMC_target_COV) delete TMCMC_target_COV;
      if (rhs.TMCMC_target_COV) {
        TMCMC_target_COV = new FlxFunction(*(rhs.TMCMC_target_COV));
      } else {
        TMCMC_target_COV = NULL;
      }
    TMCMC_update_weights = rhs.TMCMC_update_weights;
    if (TMCMC_alpha) delete TMCMC_alpha;
      if (rhs.TMCMC_alpha) {
        TMCMC_alpha = new FlxFunction(*(rhs.TMCMC_alpha));
      } else {
        TMCMC_alpha = NULL;
      }
    if (LS_SPNT) delete LS_SPNT;
      if (rhs.LS_SPNT) {
        LS_SPNT = new FlxMtxConstFun(*(rhs.LS_SPNT));
      } else {
        LS_SPNT = NULL;
      }
    if (LS_tol) delete LS_tol;
      if (rhs.LS_tol) {
        LS_tol = new FlxFunction(*(rhs.LS_tol));
      } else {
        LS_tol = NULL;
      }
    if (LS_max_iter) delete LS_max_iter;
      if (rhs.LS_max_iter) {
        LS_max_iter = new FlxFunction(*(rhs.LS_max_iter));
      } else {
        LS_max_iter = NULL;
      }
    if (pa_maxL) delete pa_maxL;
      if (rhs.pa_maxL) {
        pa_maxL = new FlxFunction(*(rhs.pa_maxL));
      } else {
        pa_maxL = NULL;
      }
  }
  return *this;
}

Flx_SuS_Control::credibleEstim Flx_SuS_Control::parse_credibleEstim(const std::string& strCredEst)
{
  if (strCredEst=="none") {
    return none;
  }
  if (strCredEst=="simple") {
    return simpleBayes;
  }
  if (strCredEst=="ccorr") {
    return ccorrBayes;
  }
  if (strCredEst=="fcorr") {
    return fcorrBayes;
  }
  if (strCredEst=="icorr") {
    return icorrBayes;
  }
  std::ostringstream ssV;
  ssV << "Unknown identifier (" << strCredEst << ").";
  throw FlxException_NeglectInInteractive("Flx_SuS_Control::parse_credibleEstim", ssV.str() );
}

const std::string Flx_SuS_Control::get_credibleStr(const Flx_SuS_Control::credibleEstim cID)
{
  switch (cID) {
    case none:
      return "none";
    case simpleBayes:
      return "simple";
    case ccorrBayes:
      return "ccorr";
    case fcorrBayes:
      return "fcorr";
    case icorrBayes:
      return "icorr";
  }
  throw FlxException_Crude("Flx_SuS_Control::get_credibleStr");
}


FlxBayUp_Update::FlxBayUp_Update(FlxRndCreator& RndCreator)
: list(NULL), csm(NULL), RndCreator(RndCreator),burbrvs(NULL),iadpt(*(data->ConstantBox.get("sus_iadpt",true))), 
  PrMod(ZERO), DataFit(ZERO), relEntr(ZERO), N_tot_LSFc_sim(0), 
  post_adpt_calls(0), post_adptcount_N(0), post_adptcount(0), post_adpt_accsmpl(0)
{
}

FlxBayUp_Update::~FlxBayUp_Update()
{
  if (list) delete list;
  if (csm) delete csm;
  for (tuint i=0;i<CLevelStat.size();++i) {
    delete CLevelStat[i];
  }
}

void FlxBayUp_Update::define_constants()
{
  // define constants
    data->ConstantBox.declareC("sus_pr");                // probability of failure / evidence of the model
    data->ConstantBox.declareC("bayup_mlnl");
    data->ConstantBox.declareC("sus_lsf_calls");
    data->ConstantBox.declareC("sus_iadpt");                // number of the adaptive step @ current cond. level (starts with 1)
    data->ConstantBox.declareC("sus_iter");                // number of the conditioning level (MC-step is 0)
    data->ConstantBox.declareC("sus_fwd_coeffofvar");
    data->ConstantBox.declareC("sus_fwd_coeffofvar_fc");
    data->ConstantBox.declareC("sus_pf_post_mean");
    data->ConstantBox.declareC("sus_pf_post_sd");
    data->ConstantBox.declareC("sus_pf_post_coeffofvar");
    data->ConstantBox.declareC("sus_kernel_h");
    data->ConstantBox.declareC("bayup_mhrs_acr");
  // define matrices
    data->ConstMtxBox.declareC("sus_gt_vec");
    data->ConstMtxBox.declareC("sus_eff_vec"); 
  // define constants for empirical estimation of seed- and level-correlation
    // seed correlation
      data->ConstantBox.declareC("sus_empi_sc_p1_0",ONE);
      data->ConstantBox.declareC("sus_empi_sc_p1_1",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p1_2",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p2_0",ONE);
      data->ConstantBox.declareC("sus_empi_sc_p2_1",-0.1);
      data->ConstantBox.declareC("sus_empi_sc_p2_2",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p3_0",0.1);
      data->ConstantBox.declareC("sus_empi_sc_p3_1",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p3_2",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p4_0",0.6);
      data->ConstantBox.declareC("sus_empi_sc_p4_1",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p4_2",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p5_0",2.5);
      data->ConstantBox.declareC("sus_empi_sc_p5_1",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p5_2",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p6_0",1.2);
      data->ConstantBox.declareC("sus_empi_sc_p6_1",ZERO);
      data->ConstantBox.declareC("sus_empi_sc_p6_2",ZERO);
    // level correlation
      data->ConstantBox.declareC("sus_empi_lc_p1_0",0.06);
      data->ConstantBox.declareC("sus_empi_lc_p1_1",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p1_2",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p2_0",ONE);
      data->ConstantBox.declareC("sus_empi_lc_p2_1",-0.1);
      data->ConstantBox.declareC("sus_empi_lc_p2_2",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p3_0",0.3);
      data->ConstantBox.declareC("sus_empi_lc_p3_1",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p3_2",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p4_0",0.8);
      data->ConstantBox.declareC("sus_empi_lc_p4_1",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p4_2",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p5_0",2.0);
      data->ConstantBox.declareC("sus_empi_lc_p5_1",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p5_2",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p6_0",1.2);
      data->ConstantBox.declareC("sus_empi_lc_p6_1",ZERO);
      data->ConstantBox.declareC("sus_empi_lc_p6_2",ZERO);
}

/**
* @brief arithmetic class: evaluates a function -> returns the result if it is larger than zero
*/
class FunEvalHelper_bayup_1 : public FunBase {
  protected:
    FunBase* fun;
  public:
    FunEvalHelper_bayup_1 (FunBase* fun) : fun(fun) {};
    virtual ~FunEvalHelper_bayup_1() { delete fun; }
    virtual const tdouble calc() { tdouble val = fun->calc(); return ((val>ZERO)?val:ZERO); }
    const bool search_circref(FlxFunction* fcr) { return fun->search_circref(fcr); }
    const std::string write() { return "(...)"; }
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return fun->dependOn_Const(thenumber); }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return fun->optimize(optf,foi); }
    virtual const bool evalw() { return fun->evalw(); }
};
 
/**
* @brief arithmetic class: evaluates a function -> stores the result to dpl_ptr
*/
class FunEvalHelper_bayup_2 : public FunEvalHelper_bayup_1 {
  private:
    tdouble* dpl_ptr;
  public:
    FunEvalHelper_bayup_2 (FunBase* fun, tdouble* dpl_ptr) : FunEvalHelper_bayup_1(fun), dpl_ptr(dpl_ptr) {};
    virtual ~FunEvalHelper_bayup_2() { }
    virtual const tdouble calc() { *dpl_ptr = fun->calc(); return *dpl_ptr; }
};
 

void FlxBayUp_Update::output_forwardEstimators()
{
  const tuint N = CLevelStat.size();
  // export the g_t-vector
    tuint N_hlp = N;
    {
      tdouble* gt_vec = data->ConstMtxBox.get_Vec("sus_gt_vec",N_hlp);
      for (tuint i=0;i<N;++i) {
        const Flx_SuS_CLevelStat& cls = *(CLevelStat[i]);
        gt_vec[i] = cls.g_t;
      }
    }
  // export the efficiency-vector
    {
      tdouble* eff_vec = data->ConstMtxBox.get_Vec("sus_eff_vec",N_hlp);
      for (tuint i=0;i<N;++i) {
        const Flx_SuS_CLevelStat& cls = *(CLevelStat[i]);
        eff_vec[i] = cls.eff;
      }
    }
  if (susControl.comp_gamma==false) return;
  // compute delta factors
    #ifndef FLX_CV_1
      tdouble delta[N];
    #else
      tdouble* delta = new tdouble[N];
    #endif
    for (tuint i=0;i<N;++i) {
      const Flx_SuS_CLevelStat& cls = *(CLevelStat[i]);
      delta[i] = sqrt((ONE-cls.pi)/(cls.pi*cls.Nsamples)*(ONE+cls.gamma));
//       delta[i] = sqrt((tdouble(cls.Nsamples-cls.Nfailures+inv_eff)*inv_eff)/((cls.Nsamples+3*inv_eff)*(cls.Nfailures+inv_eff)));
    }
  // the output stream
    std::ostream& scout3 = GlobalVar.slogcout(3);
  scout3 << "  'Forward' estimators:" << std::endl;
  // some computations
    pdouble p1;
    for (tuint i=0;i<N;++i) {
      p1 += pow2(delta[i]);
    }
    pdouble p2;
    for (tuint i=0;i<N;++i) {
      for (tuint j=i+1;j<N;++j) {
        p2 += delta[i]*delta[j];
      }
    }
    p2 *= 2;
    #ifdef FLX_CV_1
      delete [] delta;
    #endif
  // C.o.V.
    scout3 << "                      C.o.V.: " << std::format("{:9.2e}", sqrt(p1.cast2double()))
            << "          (assuming that the p_i are uncorrelated)" << std::endl;
    data->ConstantBox.insert("sus_fwd_coeffofvar",sqrt(p1.cast2double()));
    p1 += p2;
    scout3 << "                              " << std::format("{:9.2e}", sqrt(p1.cast2double()))
            << " + o(1/N) (assuming that the p_i are fully correlated)" << std::endl;
    data->ConstantBox.insert("sus_fwd_coeffofvar_fc",sqrt(p1.cast2double()));
  // bias
    scout3 << "        upper bound for bias: " << std::format("{:9.2e}", (p2.cast2double()/2))
            << " + o(1/N) (equality for full correlation of the p_i)" << std::endl;
    
}

void FlxBayUp_Update::output_credibleIntervals()
{
  if (susControl.credEst==Flx_SuS_Control::none) return;
  // the output stream
    std::ostream& scout3 = GlobalVar.slogcout(3);
  FlxSMtx* pcM = data->ConstMtxBox.get(susControl.pc->eval(),true);        // matrix that stores the credible limits to evaluate
  RBRV_set* post_PfSim = NULL;                        // set needed to draw from in order to get pr_ref
  const tdouble* pr_ref = NULL;                // reference to realization of failure probability
  const tuint N = CLevelStat.size();
  // variables for case Flx_SuS_Control::fcorrBayes
    #ifndef FLX_CV_1
      tdouble rev_pi_vec_orig[N];
      tdouble rev_pi_vec[N];
      const tdouble* rev_eff_vec[N];
    #else
      tdouble* rev_pi_vec_orig = new tdouble[N];
      tdouble* rev_pi_vec = new tdouble[N];
      const tdouble** rev_eff_vec = new const tdouble*[N];
    #endif
  try {
    if (pcM->get_nrows()*pcM->get_ncols()==0) {
      #ifdef FLX_CV_1
        delete [] rev_pi_vec_orig;
        delete [] rev_pi_vec;
        delete [] rev_eff_vec;
      #endif
      return;
    }
    // create the rbrv-set used for simulation
    switch (susControl.credEst) {
      case Flx_SuS_Control::none:
        throw FlxException_Crude("FlxBayUp_Update::output_credibleIntervals_a1");
      case Flx_SuS_Control::simpleBayes:
      case Flx_SuS_Control::ccorrBayes:
      {
        RBRV_entry** rev = new RBRV_entry*[N+1];
        FunBase* f_pr = NULL;
        for (tuint i=0;i<N;++i) {
          const Flx_SuS_CLevelStat& cls = *(CLevelStat[i]);
          std::ostringstream enam;
          enam << "dummy_" << (i+1);
          // get the parameters of the prior
            const tdouble prior_alpha = ONE;
            const tdouble prior_beta = ONE;
          const tdouble eff = ((susControl.credEst==Flx_SuS_Control::simpleBayes)?ONE:cls.eff);
          const tdouble post_alpha = prior_alpha+tdouble(cls.Nfailures)*eff;
          const tdouble post_beta  = prior_beta+tdouble(cls.Nsamples-cls.Nfailures)*eff;
          FlxFunction* alpha = new FlxFunction(new FunNumber(post_alpha));
          FlxFunction* beta = new FlxFunction(new FunNumber(post_beta));
          RBRV_entry_RV_base* rev_i = new RBRV_entry_RV_beta(enam.str(),i,false,alpha,beta,NULL,NULL,true);
          rev[i] = rev_i;
          // account for level correlation
            if (cls.corr_pi_prev > 0.01 && i>0) {
              tdouble cpi_ = cls.corr_pi_prev;
              if (cpi_>0.99) cpi_ = 0.99;
              rev_i->set_corr(dynamic_cast<RBRV_entry_RV_base*>(rev[i-1]),new FlxFunction(new FunNumber(cpi_)),true,false);
            }
          if (i==0) {
            f_pr = new FunConst(rev[i]->get_value_addr());
          } else {
            f_pr = new FunMult(f_pr,new FunConst(rev[i]->get_value_addr()));
          }
        }
        #if FLX_DEBUG
          if (f_pr==NULL) throw FlxException_Crude("FlxBayUp_Update::output_credibleIntervals_b1");
        #endif
        rev[N] = new RBRV_entry_fun("dummy_pr",new FlxFunction(f_pr));
        pr_ref = rev[N]->get_value_addr();
        post_PfSim = new RBRV_set(true,N,"dummy_for_sus",true,N+1,rev,0,NULL,false);
        break;
      }
      case Flx_SuS_Control::fcorrBayes:
      case Flx_SuS_Control::icorrBayes:
      {
        // get the parameters of the prior
          const tdouble prior_alpha = ONE;
          const tdouble prior_beta = ONE;
        std::vector<RBRV_entry*> rev_vec;
        FunBase* f_pr = NULL;
        tuint sRV = 0;                // number of true random variables allocated
        for (tuint i=0;i<N;++i) {
          const Flx_SuS_CLevelStat& cls = *(CLevelStat[i]);
          rev_pi_vec_orig[i] = tdouble(cls.Nfailures)/tdouble(cls.Nsamples);
          rev_pi_vec[i] = rev_pi_vec_orig[i];
          // compute the efficiency of the chain
            FunBase* f_eff = new FunNumber(ONE);
            for (tuint j=0;j<cls.g_N_entries;++j) {
              const tuint* gp = &(cls.g_N[j*3]);
              const tdouble post_alpha = prior_alpha+tdouble(gp[2]);
              const tdouble post_beta  = prior_beta+tdouble(gp[1]-gp[2]);
              FlxFunction* alpha = new FlxFunction(new FunNumber(post_alpha));
              FlxFunction* beta = new FlxFunction(new FunNumber(post_beta));
              // define a beta-random variable
                std::ostringstream enam;
                enam << "dummy_" << (i+1) << "_" << (j+1);
                rev_vec.push_back( new RBRV_entry_RV_beta(enam.str(),sRV,false,alpha,beta,NULL,NULL,true) );
                ++sRV;
              FunBase* f_h1 = new FunConst( rev_vec[rev_vec.size()-1]->get_value_addr() );
              f_h1 = new FunMult(f_h1,new FunConst(rev_pi_vec+i));
              f_h1 = new FunSub(f_h1,new FunPower(new FunConst(rev_pi_vec+i),new FunNumber(tdouble(2.))));                // = R_i
              FunBase* f_R0 = new FunMult(new FunConst(rev_pi_vec+i), new FunSub(new FunNumber(ONE),new FunConst(rev_pi_vec+i)) );
              f_h1 = new FunMult(f_h1,new FunMult_Div(new FunNumber(2*(tdouble(gp[0])/tdouble(cls.Nsamples))),f_R0));
              f_h1 = new FunEvalHelper_bayup_1(f_h1);
              f_eff = new FunAdd(f_eff,f_h1);
            }
            f_eff = new FunMult_Div(new FunNumber(ONE),f_eff);
          // add efficiency to set
            {
              std::ostringstream enam;
              enam << "dummy_" << (i+1) << "_eff";
              rev_vec.push_back( new RBRV_entry_fun(enam.str(),new FlxFunction(f_eff)) );
              rev_eff_vec[i] = rev_vec[rev_vec.size()-1]->get_value_addr();
            }
          // assemble alpha
            FunBase* f_alpha = new FunNumber(tdouble(cls.Nfailures));
            f_alpha = new FunMult(f_alpha,new FunConst(rev_eff_vec[i]));
            f_alpha = new FunAdd(new FunNumber(prior_alpha),f_alpha);
            FlxFunction* alpha = new FlxFunction(f_alpha);
          // assemble beta
            FunBase* f_beta = new FunNumber(tdouble(cls.Nsamples-cls.Nfailures));
            f_beta = new FunMult(f_beta,new FunConst(rev_eff_vec[i]));
            f_beta = new FunAdd(new FunNumber(prior_beta),f_beta);
            FlxFunction* beta = new FlxFunction(f_beta);
          // define the model for p_i
          {
            std::ostringstream enam;
            enam << "dummy_" << (i+1);
            rev_vec.push_back( new RBRV_entry_RV_beta(enam.str(),sRV,false,alpha,beta,NULL,NULL,false) );
            ++sRV;
            FunBase* f_const = new FunConst(rev_vec[rev_vec.size()-1]->get_value_addr());
            f_const = new FunEvalHelper_bayup_2(f_const,rev_pi_vec+i);
            if (f_pr) {
              f_pr = new FunMult(f_pr,f_const);
            } else {
              f_pr = f_const;
            }
          }
        }
        #if FLX_DEBUG
          if (f_pr==NULL) throw FlxException_Crude("FlxBayUp_Update::output_credibleIntervals_b1");
        #endif
        rev_vec.push_back( new RBRV_entry_fun("dummy_pr",new FlxFunction(f_pr)) );
        pr_ref = rev_vec[rev_vec.size()-1]->get_value_addr();
        RBRV_entry** rev = new RBRV_entry*[rev_vec.size()];
        for (tuint i=0;i<rev_vec.size();++i) {
          rev[i] = rev_vec[i];
        }
        post_PfSim = new RBRV_set(true,sRV,"dummy_for_sus",true,rev_vec.size(),rev,0,NULL,false);
        break;
      }
      default:
        throw FlxException_Crude("FlxBayUp_Update::output_credibleIntervals_h1");
    };
    // allocate isteram
      const tuint Ncs = susControl.N_cred_smpl;
      FlxIstream_vector* istrm = data->IstreamBox.get_isVector("sus_pr_bayesian");
      if (istrm==NULL) {        // ... vector-input-stream with that name does not exist
        try {
          istrm = new FlxIstream_vector("sus_pr_bayesian",NULL,false,Ncs);
          data->IstreamBox.insert("sus_pr_bayesian", istrm);
        } catch (FlxException &e) {
          FLXMSG("FlxObjInputVectorStream::task",1);
          if (istrm) delete istrm;
          throw;
        }
      } else {
        istrm->clear();
      }
    scout3 << "  Bayesian post-processor: ";
    // activate the progress bar
      std::ostream &op = *GlobalVar.get_cout();
      FlxProgress prg(op,true);
      prg.start(Ncs);
    // =============================================================================================
    // generate samples
    // =============================================================================================
      switch (susControl.credEst) {
        case Flx_SuS_Control::none:
          throw FlxException_Crude("FlxBayUp_Update::output_credibleIntervals_m1");
        case Flx_SuS_Control::simpleBayes:
        case Flx_SuS_Control::ccorrBayes:
        {
          for (tuint i=0;i<Ncs;++i) {
            prg.tick(i);
            post_PfSim->set_is_valid(false);
            post_PfSim->propose_y();
            post_PfSim->transform_y2x();
            istrm->appendNumber(*pr_ref);
          }
          break;
        }
        case Flx_SuS_Control::fcorrBayes:
        case Flx_SuS_Control::icorrBayes:
        {
          const tuint Niter = (susControl.credEst==Flx_SuS_Control::fcorrBayes)?1:3;
          flxVec v_orig(rev_pi_vec_orig,N);
          flxVec v_iter(rev_pi_vec,N);
          for (tuint i=0;i<Ncs;++i) {
            prg.tick(i);
            v_iter = v_orig;
            post_PfSim->set_is_valid(false);
            post_PfSim->propose_y();
            for (tuint j=0;j<Niter;++j) {
              post_PfSim->transform_y2x();
            }
            istrm->appendNumber(*pr_ref);
          }
          break;
        }
        default:
          throw FlxException_Crude("FlxBayUp_Update::output_credibleIntervals_p1");
      };
      prg.stop();
    // output some statistics
      scout3 << std::endl;
      // run 'FlxObjStatSmp'
      {
        FlxObjStatSmp sts(false,"cout",istrm,"sus_pf_post",Ncs,1,true);
        sts.exec();
      }
      scout3 << "                 estim. mean:  " << GlobalVar.Double2String(*data->ConstantBox.get("sus_pf_post_mean"),false,2)  << std::endl;
      scout3 << "            estim. std. dev.:  " << GlobalVar.Double2String(*data->ConstantBox.get("sus_pf_post_sd"),false,2)  << std::endl;
      scout3 << "               estim. C.o.V.:  " << GlobalVar.Double2String(*data->ConstantBox.get("sus_pf_post_coeffofvar"),false,2)  << std::endl;
      scout3 << "      number of samples used:  " << Ncs << std::endl;
//     scout3 << "        upper bound for bias: " << 
      istrm->sortStream();
      const tdouble* const tp = istrm->get_tmpPtr();
    // where to store the credible intervals (two-sided) to?:
      tuint Nentries = pcM->get_nrows()*pcM->get_ncols();
      tuint Ncols = 3;
      tdouble *crp = data->ConstMtxBox.get_Mtx("sus_credible_2sided",Nentries,Ncols);
      tuint crp_id = 0;
      scout3 << "    Credible intervals (two-sided) of estimate:" << std::endl;
      scout3 << "                          Bayesian" << std::endl;
      for ( tuint i=0; i<pcM->get_nrows(); ++i ) {
        for ( tuint j=0; j<pcM->get_ncols(); ++j ) {
          const tdouble pcT = pcM->operator()(i,j);
          // Calculate the credible interval
            scout3 << "       " << GlobalVar.Double2String(pcT,false,4,6) << ":    ";
            crp[crp_id++] = pcT;
            const tdouble bay_low = FunSmpCDF::inv_cdf((ONE-pcT)/2,tp,Ncs);
            const tdouble bay_up = FunSmpCDF::inv_cdf(ONE-(ONE-pcT)/2,tp,Ncs);
            crp[crp_id++] = bay_low;
            crp[crp_id++] = bay_up;
            scout3 << "[" << std::format("{:10.3e}", bay_low) << ";" << std::format("{:10.3e}", bay_up) << "]" << std::endl;
        }
      }
      // where to store the credible intervals (upper bound) to?:
        Ncols = 2;
        crp = data->ConstMtxBox.get_Mtx("sus_credible_ubound",Nentries,Ncols);
        crp_id = 0;
      scout3 << "    Credible intervals (upper bound) of estimate:" << std::endl;
      scout3 << "                          Bayesian" << std::endl;
      for ( tuint i=0; i<pcM->get_nrows(); ++i ) {
        for ( tuint j=0; j<pcM->get_ncols(); ++j ) {
          const tdouble pcT = pcM->operator()(i,j);
          // Calculate the cerdible interval
            scout3 << "       " << GlobalVar.Double2String(pcT,false,4,6) << ":    ";
            crp[crp_id++] = pcT;
            const tdouble bay_ub = FunSmpCDF::inv_cdf(pcT,tp,Ncs);
            crp[crp_id++] = bay_ub;
            scout3 << "       " << std::format("{:10.3e}", bay_ub) << std::endl;;
        }
      }
  } catch (FlxException& e) {
    if (post_PfSim) delete post_PfSim;
    #ifdef FLX_CV_1
      delete [] rev_pi_vec_orig;
      delete [] rev_pi_vec;
      delete [] rev_eff_vec;
    #endif
    throw;
  }
  #ifdef FLX_CV_1
    delete [] rev_pi_vec_orig;
    delete [] rev_pi_vec;
    delete [] rev_eff_vec;
  #endif
  if (post_PfSim) delete post_PfSim;
}

void FlxBayUp_Update::draw_realization(flxVec& smpl)
{
  const tuint id = list->get_cur_i_list();
  #if FLX_DEBUG
    if (id < 0) throw FlxException_Crude("FlxBayUp_Update::draw_realization_1");
  #endif
  const tuint N_RV = smpl.get_N();
  #if FLX_DEBUG
    if (list->get_Nrv() != N_RV) throw FlxException_Crude("FlxBayUp_Update::draw_realization_2");
  #endif
  RBRV_constructor& RndBox(list->get_RndBox());
  bool proposal_accepted = false;        // must be set to true if the proposed sample is accepted
  // id==0 -> current sample was not used before
  if (id>0) {        // propose a new realization
    switch (list->meth_id) {
      // ***************************************************************************************************************
      case (FlxBayUp_Update_List::BUS):
      case (FlxBayUp_Update_List::UBUS):
      case (FlxBayUp_Update_List::BUST):
      case (FlxBayUp_Update_List::ABCSUBSIM):
      case (FlxBayUp_Update_List::RASUBSIM):
      case (FlxBayUp_Update_List::TMCMC):
      {
        flxVec y_prev(list->get_cur_y_list(),N_RV);
        csm->prepareS(y_prev);
        bool rejectIt = false;
        try {
          rejectIt = !(csm->propose(smpl,y_prev));
          if (rejectIt==false) {
            RndBox.set_smp(smpl);
          }
        } catch (FlxException &e) {
          FLXMSG("FlxBayUp_Update::draw_realization_BUS_1",1);
          rejectIt = true;
        }
        if (!rejectIt) {
          proposal_accepted = list->insert_entry(false,rejectIt,true,false,NULL);        // true: proposal was accepted && max. L did not change!
        }
        // adative feature 
          ++post_adptcount;
          if (proposal_accepted) {
            ++post_adpt_accsmpl;
          }
          if (post_adptcount_N>0 && post_adptcount>=post_adptcount_N) {        // perform the adaptive step
            ++post_adpt_calls;
            const tdouble iadpt_prev = iadpt;
            iadpt = tdouble(post_adpt_calls);
            list->get_adpt_ctrl().requires_adptv_step(tdouble(post_adpt_accsmpl)/tdouble(post_adptcount),*csm);
            iadpt = iadpt_prev;
            post_adpt_accsmpl = 0;
            post_adptcount = 0;
            csm->prepare();        // just to make sure that this is called from time to time ...
          }
        break;
      }
      // ***************************************************************************************************************
      case (FlxBayUp_Update_List::MHRS):
      case (FlxBayUp_Update_List::RS):
      {        
        tdouble ts=ONE, tL=ZERO;
        bool err;
        do {        // loop until we can accept smpl
          try {
            // draw a realization from the prior
              RndBox.gen_smp(); 
            // check if proposed sample is in failure domain
              tL = list->parent.eval_Likelihood();
              RndBox.get_y_Vec(smpl.get_tmp_vptr());
              ts = list->eval_LSF(smpl[N_RV-1],tL);
              err = false;
          } catch (FlxException &e) {
            FLXMSG("FlxBayUp_Update::draw_realization_RS_1",1);
            err = true;
            continue;
          }
        } while (ts>ZERO || err);
        bool acpt = true;
        if (list->meth_id==FlxBayUp_Update_List::MHRS) {        // Metropolis-Hastings step
          const tdouble c = list->parent.get_cStart();
          const tdouble ur = ((tL>c)?tL:c)/list->get_cur_L();
          acpt = (RndCreator.gen_smp_uniform()<=ur);
        }
        list->insert_entry(false,!acpt,true,false,NULL,tL);   
        proposal_accepted = acpt;
        break;
      }
      // ***************************************************************************************************************
      case (FlxBayUp_Update_List::ABCRS):
      {
        // register the progress-bar
        const tdouble thr_m = list->parent.get_cStart();
        tdouble metric=ONE;
        bool err;
        do {        // loop until we can accept smpl
          try {
            // draw a realization from the prior
              RndBox.gen_smp(); 
            // check if proposed sample is in failure domain
              metric = list->parent.eval_Likelihood();        // evaluates the ABC-metric (returned as log-transform)
              err = false;
          } catch (FlxException &e) {
            FLXMSG("FlxBayUp_Update::draw_realization_ABCRS_1",1);
            err = true;
            continue;
          }
        } while (metric>thr_m || err);
        RndBox.get_y_Vec(smpl.get_tmp_vptr());
        list->insert_entry(false,false,true,false,NULL); 
        proposal_accepted = true;
        break;
      }
      // ***************************************************************************************************************
      case (FlxBayUp_Update_List::RAMCI):
      {
        tdouble lsf=ONE;
        bool err;
        do {        // loop until we can accept smpl
          try {
            // draw a realization from the prior
              RndBox.gen_smp(); 
            // check if proposed sample is in failure domain
              lsf = list->parent.eval_RAlsf();        // evaluates the limit-state function
              err = false;
          } catch (FlxException &e) {
            FLXMSG("FlxBayUp_Update::draw_realization_RAMCI_1",1);
            err = true;
            continue;
          }
        } while (lsf>ZERO || err);
        RndBox.get_y_Vec(smpl.get_tmp_vptr());
        list->insert_entry(false,false,true,false,NULL); 
        proposal_accepted = true;
        break;
      }
      default:
        throw FlxException_Crude("FlxBayUp_Update::draw_realization_x3"); 
    }
  }
  if (!proposal_accepted) {        // forget the proposed coordinates!
    flxVec y(list->get_cur_y_list(),list->get_Nrv());
    list->get_RndBox().set_smp_y(y);
    flxVec x(list->get_cur_x_list(),list->get_NOX());
    list->parent.get_RndBox().set_smp_x(x);
    smpl = y;
  }
  // move curID to the next sample
    list->set_next_draw();
}

const bool FlxBayUp_Update::chk_accept_cur_smpl()
{
  return list->insert_entry(false,false,true,true,NULL);
}

const FlxBayUp_Update_List& FlxBayUp_Update::get_list() const
{
  if (list==NULL) throw FlxException_Crude("FlxBayUp_Update::get_list");
  return *list;
}

void FlxBayUp_Update::reset_finalized_smpls()
{
  // check if the updating problem was solved.
    if (list==NULL || list->is_finalized()==false) {
      std::ostringstream ssV;
      ssV << "Please perform a Bayesian updating of the set first.";
      throw FlxException_NeglectInInteractive("FlxObjBayUp_Reset_Smpls::task_1", ssV.str() );
    }
  // reset the status of the posterior samples
    list->reset_finalized();
}

void FlxBayUp_Update::get_sus_level_info(const std::string vecs, const tuint pid, const tuint pid2)
{
  switch (pid) {
    case 1:
    {
      tdouble* gv = data->ConstMtxBox.get_Vec(CLevelStat.size(),vecs);
      for (tuint i=0;i<CLevelStat.size();++i) {
        Flx_SuS_CLevelStat& csi = *(CLevelStat[i]);
        gv[i] = csi.g_t;
      }
      break;
    }
    case 2:
    {
      tdouble* gv = data->ConstMtxBox.get_Vec(CLevelStat.size(),vecs);
      for (tuint i=0;i<CLevelStat.size();++i) {
        Flx_SuS_CLevelStat& csi = *(CLevelStat[i]);
        gv[i] = csi.pi;
      }
      break;
    }
    case 3:
    {
      tdouble* gv = data->ConstMtxBox.get_Vec(CLevelStat.size(),vecs);
      for (tuint i=0;i<CLevelStat.size();++i) {
        Flx_SuS_CLevelStat& csi = *(CLevelStat[i]);
        gv[i] = csi.eff;
      }
      break;
    }
    case 4:
    {
      tdouble* gv = data->ConstMtxBox.get_Vec(CLevelStat.size(),vecs);
      for (tuint i=0;i<CLevelStat.size();++i) {
        Flx_SuS_CLevelStat& csi = *(CLevelStat[i]);
        gv[i] = csi.eff_Gelman;
      }
      break;
    }
    case 5:
    {
      if (pid2>=CLevelStat.size()) {
        std::ostringstream ssV;
        ssV << "Index '" << pid2 << "' must be smaller than " << CLevelStat.size() << ".";
        throw FlxException("FlxBayUp_Update::get_sus_level_info_01",ssV.str());
      }
      Flx_SuS_CLevelStat& csi = *(CLevelStat[pid2]);
      if (csi.find_multiples==NULL) {
        throw FlxException("FlxBayUp_Update::get_sus_level_info_02","Search for multiple files has not been activated.");
      }
      tdouble* gv = data->ConstMtxBox.get_Vec(csi.Nsamples,vecs);
      for (tuint i=0;i<csi.Nsamples;++i) {
        gv[i] = tdouble(csi.find_multiples[i]);
      }
      break;
    }
    default:
      {
        std::ostringstream ssV;
        ssV << "ID '" << pid << "' not defined.";
        throw FlxException("FlxBayUp_Update::get_sus_level_info_03",ssV.str());
      }
  }
}

std::ofstream* FlxBayUp_Update::open_smpl_file4write()
{
  if (susControl.os_samples) {
    const std::string fn = susControl.os_samples->eval(false);
    if (fn.empty()) return NULL;
    std::ofstream* os_smpl = new std::ofstream(fn.c_str());
    if ( os_smpl == NULL || !os_smpl->is_open() ) {
      throw FlxException("FlxBayUp_Update::update_b0", "File '" + fn + "' could not be opened." );
    }
    return os_smpl;
  }
  return NULL;
}

void FlxBayUp_Update::TMCMC_weight_vec(const tdouble q_prev, const tdouble q_now, const flxVec& Lvec, flxVec& Wvec)
{
  Wvec = Lvec;
  tdouble* wp = Wvec.get_tmp_vptr();
  size_t N = Wvec.get_N();
  const tdouble dq = q_now-q_prev;
  if (q_now!=q_now) throw FlxException_Crude("FlxBayUp_Update::TMCMC_weight_vec");
  for (size_t i=0;i<N;++i) {
    wp[i] = exp(wp[i]*dq);
  }
}

const tdouble FlxBayUp_Update::TMCMC_COV(const tdouble q_prev, const tdouble q_now, const flxVec& Lvec, flxVec& Wvec)
{
  TMCMC_weight_vec(q_prev,q_now,Lvec,Wvec);
  const tdouble mu = Wvec.get_Mean();
  const tdouble sd = Wvec.get_sd(mu);
  tdouble res = sd/mu;
  if (res!=res || sd==ZERO) res = std::numeric_limits<tdouble>::infinity();
  return res;
}

const tdouble FlxBayUp_Update::TMCMC_new_q(const tdouble q_prev, const tdouble target_cov, const flxVec& Lvec, flxVec& Wvec)
{
  if (q_prev>=ONE) throw FlxException_Crude("FlxBayUp_Update::TMCMC_new_q_1");
  // initialize variables
    tdouble q_low = q_prev;
    tdouble q_up = ONE;
    tdouble cov_low = ZERO;
    tdouble cov_up = TMCMC_COV(q_prev,q_up,Lvec,Wvec);
  // check if we can already stop
    if (cov_up<=target_cov) {
      return ONE;
    }
  // search a proper starting interval
    tuint i;
    const tuint Nin = 5;
    const tuint Npot = 8;
    for (i=1;i<=Nin;++i) {
      q_up = ((i==Nin)?ONE:pow(i/(ONE*Nin),Npot));
      if (q_up<q_prev) continue;
      cov_up = TMCMC_COV(q_prev,q_up,Lvec,Wvec);
      if (cov_up>=target_cov) break;
      else {
        q_low = q_up;
        cov_low = cov_up;
      }
    }
    if (i>Nin) throw FlxException_Crude("FlxBayUp_Update::TMCMC_new_q_2");
    // to avoid cov_up=inf !!!
      while (cov_up>tdouble(100.)) {
        q_up = (q_up+q_low)/2;
        cov_up = TMCMC_COV(q_prev,q_up,Lvec,Wvec);
      }
      if (cov_up<=cov_low) throw FlxException_Crude("FlxBayUp_Update::TMCMC_new_q_3");
  // start the actual iteration
    const tuint NMAX = tuint(1e2);
    const tdouble dx = 1e-5;
    const tdouble dy = 1e-4;
    tuint count = 0;
    tdouble f1 = cov_low-target_cov;
    tdouble f2 = cov_up-target_cov;
    if (fabs(f1)<=GlobalVar.TOL()) {
      throw FlxException_Crude("FlxBayUp_Update::TMCMC_new_q_4");
    }
    if (fabs(f2)<=GlobalVar.TOL()) {
      return q_up;
    }
    while (fabs(q_up-q_low)*2/fabs(q_up+q_low)>dx && count<NMAX) {
      ++count;
      tdouble res = (q_low*f2-q_up*f1)/(f2-f1);
      tdouble f3 = TMCMC_COV(q_prev,res,Lvec,Wvec)-target_cov;
      tuint c = 0;
      bool b = false;
      while (f3==std::numeric_limits<tdouble>::infinity() || b) {
        if (b) {
          const tdouble help3 = res + (res-((q_up>q_low)?q_up:q_low))/2;
          if (q_up<q_low) {
            std::swap(q_low,q_up);
            std::swap(f1,f2);
          }
          q_low = q_up;
          f1 = f2;
          q_up = res;
          f2 = f3;
          res = help3;
        } else {
          res = (res + ((q_up>q_low)?q_up:q_low))/2;
        }
        f3 = TMCMC_COV(q_prev,res,Lvec,Wvec)-target_cov;
        if (fabs(f3/target_cov)<=dy) {
          return res;
        }
        b = (f3*f1>ZERO && f3*f2>ZERO);
        if (++c>100) break;
      }
      if (fabs(f3/target_cov)<=dy) {
        return res;
      }
      if ( f2*f3<ZERO ) {
        q_low = q_up;
        f1 = f2;
        q_up = res;
        f2 = f3;
      } else {
        const tdouble m = f2/(f2+f3);
        f1 = m*f1;
        q_up = res;
        f2 = f3;
      }
    }
  // return
    const tdouble res = (q_up>q_low)?q_low:q_up;
    TMCMC_COV(q_prev,res,Lvec,Wvec);
    if (res!=res) throw FlxException_Crude("FlxBayUp_Update::TMCMC_new_q_7");
    return res;
}

void FlxBayUp_Update::TMCMC_assemble_smplCOV(flxVec& y_covV, const flxVec& W_vec, const tdouble* y_list_prev, const tuint Nc, const tuint NRV, flxpVec& y_meanP, flxVec& y_mean, flxVec& y_tmp )
{
  const tdouble W_mean = W_vec.get_Mean();
  y_meanP.set_zero();
  for (tuint i=0;i<Nc;++i) {
    const flxVec sy(&(y_list_prev[i*NRV]),NRV);
    y_meanP.add(sy,W_vec[i]);
  }
  y_mean = y_meanP;
  y_mean /= (W_mean*Nc);
  y_covV.set_zero();
  for (tuint i=0;i<Nc;++i) {
    const flxVec sy(&(y_list_prev[i*NRV]),NRV);
    y_tmp = sy;
    y_tmp -= y_mean;
    size_t c1 = 0;
    for (tuint j=0;j<NRV;++j) {
      for (tuint col=0;col<=j;++col) {
        y_covV[c1++] += y_tmp[j]*y_tmp[col]*W_vec[i];
      }
    }
  }
  y_covV /= (W_mean*Nc);
}

void FlxBayUp_Update::update(FlxBayUp_Update_List* listV, FlxBayUP_csm_base* csmV, const bool use_cStart, const Flx_SuS_Control& susControlV)
{
  // --------------------------------------------------------------------------------------------------------------------------
  // reset the algorithm (and initialize)
  // --------------------------------------------------------------------------------------------------------------------------
  const bool first_run = (list==NULL);
  if (list) {
    if (
      use_cStart 
      || list->meth_id==FlxBayUp_Update_List::ABCSUBSIM
      || list->meth_id==FlxBayUp_Update_List::ABCRS
    ) {
      list->parent.get_cStart() = list->parent.get_cStart_init();
    }
    delete list;
    list = NULL;
    if (csm) { 
      delete csm;
      csm = NULL;
    } else {
      if (csmV) throw FlxException_Crude("FlxBayUp_Update::update_a1"); 
    }
  }
  list = listV;
  csm = csmV;
  // clear CLevelStat
    for (tuint i=0;i<CLevelStat.size();++i) {
      delete CLevelStat[i];
    }
    CLevelStat.clear();
  susControl = susControlV;
    list->parent.get_pa_maxL() = susControl.pa_maxL->cast2positive_or0();
  std::ofstream *os_smpl = NULL;
  try {
  // --------------------------------------------------------------------------------------------------------------------------
  // perform some initial checks
  // --------------------------------------------------------------------------------------------------------------------------
  if (list->parent.get_methCat()==flxBayUp::UNDEFINED) throw FlxException_Crude("FlxBayUp_Update::update_a2");
  switch (list->meth_id) {
    case FlxBayUp_Update_List::BUS:
    case FlxBayUp_Update_List::UBUS:
    case FlxBayUp_Update_List::BUST:
    case FlxBayUp_Update_List::MHRS:
    case FlxBayUp_Update_List::RS:
    case FlxBayUp_Update_List::LS:
      if (list->parent.get_methCat()!=flxBayUp::BUS) {
        std::ostringstream ssV;
        ssV << "The likelihood is not defined in a correct way.";
        throw FlxException("FlxBayUp_Update::update_a3", ssV.str() ); 
      }
      break;
    case FlxBayUp_Update_List::TMCMC:
      if (list->parent.get_methCat()!=flxBayUp::TMCMC) {
        std::ostringstream ssV;
        ssV << "The likelihood is not defined in a correct way.";
        throw FlxException("FlxBayUp_Update::update_a4", ssV.str() ); 
      }
      break;
    case FlxBayUp_Update_List::ABCSUBSIM:
    case FlxBayUp_Update_List::ABCRS:
      if (list->parent.get_methCat()!=flxBayUp::ABC) {
        std::ostringstream ssV;
        ssV << "The metric is not defined in a correct way - use object 'bayup_abcmetric'.";
        throw FlxException("FlxBayUp_Update::update_a5", ssV.str() ); 
      }
      break;
    case FlxBayUp_Update_List::RASUBSIM:
    case FlxBayUp_Update_List::RAMCI:
      if (list->parent.get_methCat()!=flxBayUp::RA) {
        std::ostringstream ssV;
        ssV << "The limit-state function is not defined in a correct way - use object 'bayup_ralsf'.";
        throw FlxException("FlxBayUp_Update::update_a6", ssV.str() ); 
      }
      break;
    default:
      throw FlxException_Crude("FlxBayUp_Update::update_a7"); 
  };
  // --------------------------------------------------------------------------------------------------------------------------
  // output properties of the algorithm
  // --------------------------------------------------------------------------------------------------------------------------
  std::ostream& scout = GlobalVar.slogcout(4);
  if (list->parent.is_subsetRel) {
    scout << "Reliability analysis: Subset simulation" << std::endl;
  } else {
    scout << "BayUp: performs the Bayesian updating numerically" << std::endl;
    scout << "  Updating method:              ";
    switch (list->meth_id) {
      case FlxBayUp_Update_List::BUS:
        scout << "BUS  ";
        break;
      case FlxBayUp_Update_List::UBUS:
        scout << "uBUS  ";
        break;
      case FlxBayUp_Update_List::BUST:
        scout << "BUST ";
        break;
      case FlxBayUp_Update_List::ABCSUBSIM:
        scout << "ABCsubSim";
        break;
      case FlxBayUp_Update_List::RASUBSIM:
        scout << "RAsubSim";
        break;
      case FlxBayUp_Update_List::TMCMC:
        scout << "TMCMC ";
        switch (susControl.TMCMC_update_weights) {
          case 0:
            scout << "(original variant of TMCMC)";
            break;
          case 1:
            scout << "(weights are updated after each sample - variant 1)";
            break;
          case 2:
            scout << "(weights are updated after each sample - variant 2)";
            break;
          default:
          {
            std::ostringstream ssV;
            ssV << "The optional parameter 'tmcmc::tmcmc_update_weights' can only take value 0, 1, 2; and not '" << susControl.TMCMC_update_weights << "'.";
            throw FlxException("FlxBayUp_Update::update_t45",ssV.str());
          }
        }
        break;
      case FlxBayUp_Update_List::MHRS:
        scout << "MHRS";
        break;
      case FlxBayUp_Update_List::RS:
        scout << "RS";
        break;
      case FlxBayUp_Update_List::ABCRS:
        scout << "ABCrs";
        break;
      case FlxBayUp_Update_List::RAMCI:
        scout << "RAmci";
        break;
      case FlxBayUp_Update_List::LS:
        scout << "LS";
        break;
      default:
      {
        std::ostringstream ssV;
        ssV << "Unknown method ID '" << list->meth_id << "'.";
        throw FlxException("FlxBayUp_Update::update_a8", ssV.str() ); 
      }
    };
    if ( list->meth_id==FlxBayUp_Update_List::UBUS ) {
      if (list->parent.get_pa_maxL()<ONE) {
        scout << "      (with pa: " << GlobalVar.Double2String(list->parent.get_pa_maxL()) << ")";
      }
    }
    scout << std::endl;
  }
  scout << "  Number of random variables:   " << std::format("{:<10d}", list->get_Nrv());
    if (list->get_Nrv()!=list->get_NOX()) {
      scout << " (standard normal space)"<< std::endl;
      scout << "                                " << std::format("{:<10d}", list->get_NOX()) << " (original space)";
    }
    scout << std::endl;
  if (list->parent.is_subsetRel==false) {
  scout << "  Number of posterior samples:  " << list->get_Ns_final() << std::endl;
  }
  if (csm) {
    scout << "  Number of samples per step:   " << list->get_Nc()*list->get_Ncl() << std::endl;
    if (list->meth_id==FlxBayUp_Update_List::TMCMC) {
    scout << "  Burn-in phase before level:   " << list->get_Nburn() << std::endl;
    } else {
    scout << "  Threshold level:              " << std::format("{:4.2f}", (ONE/tdouble(list->get_Ncl()))) << std::endl;
    }
    scout << "  Type of conditional sampling: " << csm->print_info() << std::endl;
    list->get_adpt_ctrl().eval();
    if (list->get_adpt_ctrl().is_adaptive()) {
      list->get_adpt_ctrl().print_info(scout);
      flxBayUp_adaptive_ctrl_opti_jump* ap = dynamic_cast<flxBayUp_adaptive_ctrl_opti_jump*>(&(list->get_adpt_ctrl()));
      if (ap) {
        ap->initialize(list->get_Nrv(),list->get_Ns_final());
      }
    }
    if (list->meth_id!=FlxBayUp_Update_List::TMCMC) {
      if (susControl.comp_gamma) {
        scout << "  Estimate efficiency of chain: " << "yes" << std::endl;
        scout << "    approx. correl. of seeds:   " << (susControl.consider_seed_corr?"yes":"no") << std::endl;
        if (susControl.consider_pi_corr) {
        scout << "    approx. correl. of levels:  " << (susControl.consider_pi_corr?"yes":"no") << std::endl;
        }
      }
      if (susControl.credEst!=Flx_SuS_Control::none) {
        scout << "  Method to est. credible int.: " << Flx_SuS_Control::get_credibleStr(susControl.credEst) << std::endl;
      }
    }
  }
  if (list->meth_id!=FlxBayUp_Update_List::LS) {
    scout << " --------------------------------------------------------------------" << std::endl;
  }
  // register progress-bar
    std::ostream &op = *GlobalVar.get_cout();
    FlxProgress prg(op,true);
  // some definitions
    tuint N_runs = 0;                // number of conditioning levels
    tulong Nlsf_calls = 0;        // number of successful LSF-calls
    tulong Nlsf_errors = 0;        // number of failed LSF-calls
    tdouble* iterN = flxBayUp::get_data().ConstantBox.get("sus_iter",true);        // user-accessible variant of N_runs
      *iterN = 0;
    RBRV_constructor& RndBox = list->get_RndBox();
    flxVec y_prop(list->get_Nrv());        // proposal
    tdouble p_mod = ZERO;                // probability of the model (log-transform)
  switch (list->meth_id) {
    case (FlxBayUp_Update_List::BUS):
    case (FlxBayUp_Update_List::UBUS):
    case (FlxBayUp_Update_List::BUST):
    case (FlxBayUp_Update_List::ABCSUBSIM):
    case (FlxBayUp_Update_List::RASUBSIM):
    {
    // --------------------------------------------------------------------------------------------------------------------------
    // BUS, BUST, ABC-SubSim, RASUBSIM
    // --------------------------------------------------------------------------------------------------------------------------
      if ( list->meth_id!=FlxBayUp_Update_List::UBUS && list->parent.get_pa_maxL()<ONE ) {
        throw FlxException_NotImplemented("FlxBayUp_Update::update_b0");
      }
      std::ostringstream scout_ext; bool scout_ext_nempty = false;
      bool err_first = susControl.prt_alert;
      // do MCI in a first step
        scout << std::format("  iter {:2d}:", N_runs) << "  ";
        os_smpl = open_smpl_file4write();
        // register the progress-bar
          op.flush();
          prg.start(list->get_Nc());
        {
          tuint i=0; 
          const tuint NsM = list->get_Ncl();
          while (list->get_cur_i_list()==0) {
            try {
              RndBox.gen_smp(); 
              list->insert_entry(true,false,false,false,os_smpl); 
              ++Nlsf_calls;
            } catch (FlxException_NeglectInInteractive &e) {
              FLXMSG("FlxBayUp_Update::update_b1",1);
              if (err_first) {
                GlobalVar.alert.alert("FlxBayUp_Update::update_b2",e.what());
                err_first = false;
              }
              ++Nlsf_errors;
              continue;
            }
            ++i;
            if (i%NsM==0) prg.tick(i/NsM);
          }
        }
        if (os_smpl) {
          delete os_smpl;
          os_smpl = NULL;
        }
        // deactivate the progress-bar
          prg.stop();
          ++N_runs;
        switch (list->meth_id) {
          case (FlxBayUp_Update_List::BUS):
          case (FlxBayUp_Update_List::UBUS):
          case (FlxBayUp_Update_List::BUST):
          {
            list->update_c(p_mod,true);
            const tdouble mhml = list->get_maxL();
            if (mhml>=log(GlobalVar.TOL()) && mhml<=100.) {
              scout << std::format("maxL={:9.2e}  ", exp(mhml)) << "           ";
            } else {
              scout << std::format("mlnL={:9.2e}  ", mhml) << "           ";
            }
            break;
          }
          case (FlxBayUp_Update_List::ABCSUBSIM):
          case (FlxBayUp_Update_List::RASUBSIM):
          {
            scout << "          ";
            break;
          }
          default:
            throw FlxException_Crude("FlxBayUp_Update::update_b3");
        };
        scout << std::format("velo={:4.2f}  ", list->get_velo());
          const tdouble p_0 = ONE/tdouble(list->get_Ncl());
      // stores statistics of the previous level
          Flx_SuS_CLevelStat* curLevelStat = new Flx_SuS_CLevelStat(0,NULL);
          curLevelStat->Nsamples = list->get_Nc()*list->get_Ncl();
          Flx_SuS_CLevelStat* curLevelStat_next = NULL;
          scout_ext << std::endl << "            ";
      // perform a conditioned sampling
        try {
        while (true) {
          // prepare the seeds
            const bool fullList_prev = list->is_fullList();
            tdouble pnow; // needed just for BUST
            curLevelStat_next = new Flx_SuS_CLevelStat(CLevelStat.size()+1,curLevelStat);
            const tuint Nc_now = ((list->meth_id==FlxBayUp_Update_List::BUST)?(list->update_thr_BUST(pnow,*curLevelStat_next,*curLevelStat)):(list->update_thr(*curLevelStat_next,*curLevelStat)));
            curLevelStat->Nfailures = Nc_now;
            curLevelStat->g_t = list->get_thr();
            bool gt_is_zero = list->is_gt_zero();
            switch (list->meth_id) {
              case (FlxBayUp_Update_List::BUS):
              case (FlxBayUp_Update_List::UBUS):
              case (FlxBayUp_Update_List::ABCSUBSIM):
              case (FlxBayUp_Update_List::RASUBSIM):
              {
                if (gt_is_zero) {
                  if (Nc_now<list->get_Ns_final()) {        // only if we are not already done
                    const tdouble rt = tdouble(Nc_now) / tdouble(curLevelStat->Nsamples);
                    curLevelStat->pi = rt;
                    scout << std::format("pt={:9.2e}",rt);
                  }
                } else {
                  if (fullList_prev) {
                    curLevelStat->pi = tdouble(list->get_Nc())/tdouble(list->get_Ns_final());
                  } else {
                    curLevelStat->pi = p_0;
                  }
                  scout << std::format("gt={:9.2e}",list->get_thr());
                }
                break;
              }
              case (FlxBayUp_Update_List::BUST):
              {
                if (Nc_now<list->get_Ns_final()) {        // only if we are not already done
                  curLevelStat->pi = pnow;
                  if (gt_is_zero) {
                    scout << std::format("pt={:9.2e}",pnow);
                  } else {
                    scout << std::format("gt={:9.2e}",list->get_thr());
                  }
                }
                break;
              }
              default:
                throw FlxException_Crude("FlxBayUp_Update::update_b4");
            };
            // compute the effective number of samples
              if (susControl.comp_gamma && Nc_now!=curLevelStat->Nsamples) {        // compute gamma
                list->compute_gamma(*curLevelStat,susControl);
                if (susControl.verbose && N_runs>1) {
                  scout_ext << std::format("Ge={:2.0f}", (curLevelStat->eff_Gelman*100)) << "% ";
                  scout_ext << std::format("eff={:2.0f}", (curLevelStat->eff*100)) << "% ";
                  scout_ext << std::format("lag1c={:2.0f}", (curLevelStat->lag1_corr*100)) << "% ";
                  if (N_runs>2 && susControl.consider_seed_corr) {
                    scout_ext << std::format("sc2eff={:2.0f}", (curLevelStat->gamma_from_seed/curLevelStat->gamma*100)) << "% ";
                    if (curLevelStat->corr_pi_prev>0.01) {
                      scout_ext << std::format("piCorr={:2.0f}", (curLevelStat->corr_pi_prev*100)) << "% ";
                      #if FLX_DEBUG
                        scout_ext << std::format("piCnF={:1.1f}", (curLevelStat->corr_pi_prev_negFrac)) << " ";
                      #endif
                    }
                  }
                  scout_ext_nempty = true;
                }
              }
            if ( Nc_now!=list->get_Ns_final() ) {
              p_mod += log(curLevelStat->pi);
              CLevelStat.push_back(curLevelStat);
              curLevelStat = curLevelStat_next;
              curLevelStat_next = NULL;
            }
          // check if we can stop
            if (Nc_now==list->get_Ns_final() || (list->parent.is_subsetRel && gt_is_zero) ) {
              if (susControl.verbose) {
                list->print_ext_out(scout);
                if (scout_ext_nempty) scout << scout_ext.str();
                scout_ext.str("");
                scout_ext << std::endl;
                scout_ext_nempty = false;
              }
              break;
            }
          #if FLX_DEBUG
            if (Nc_now>list->get_Ns_final()) throw FlxException_Crude("FlxBayUp_Update::update_b5");
          #endif
          // output info about the state of the MCMC method
            *iterN = N_runs;
            if (list->get_adpt_ctrl().is_adaptive() || N_runs>1) {
              csm->write_adaptive_info(scout,list->get_adpt_ctrl().is_adaptive());
            }
            // some preliminarities of the adaptive scheme
              list->get_adpt_ctrl().eval();
          // output/initialize verbose output
            if (susControl.verbose) {
              list->print_ext_out(scout);
              if (scout_ext_nempty) scout << scout_ext.str();
              scout_ext.str("");
              scout_ext << std::endl << "            ";
              scout_ext_nempty = false;
            }
          scout << std::endl;
          // generate a new set of samples
            if (gt_is_zero) {
              scout << "        *:  ";
            } else {
              scout << "  " << std::format("iter {:2d}:  ", N_runs);
            }
            // visibility of progress-bar
              op.flush();
          // ... the actual generation of new samples
            // counts the total number of samples that were not accepted
              tuint nac = 0;
            // counts the rejected samples in the first stage
              tuint nac1 = 0;
            // number of samples in this level
              tuint tc2 = 0;
            prg.start(Nc_now);        // activate progress-bar
            // gives the number of the calls of the adaptive step (exports it for the user to work with)
              tuint adupd_ncalls = 1;
            // counters for the adaptive process
              const tuint adupd_counter_N = list->get_adpt_ctrl().get_updatesAfterNsamples();
              // counts the number of seeds since the last adaptive-step (to decide whether and adaptive step is needed)
                tuint adupd_counter = 0;
              // counts the total number of LSF-calls since the last adaptive step
                tuint adupd_total = 0;
              // number of accepted samples since the last adaptive-step
                tdouble adupd_acr = ZERO;
//                 tdouble* adupt_acrP = (list->meth_id==FlxBayUp_Update_List::BUS||list->meth_id==FlxBayUp_Update_List::BUST)?(&adupd_acr):NULL;
                tdouble* adupt_acrP = NULL;
                tuint adupd_acrN = 0;
            csm->prepare();
            os_smpl = open_smpl_file4write();
            for (tuint i=0;i<Nc_now;++i) {
              tuint& curSID = list->get_curSID_ref();
              if (os_smpl) list->write_smpl(curSID,*os_smpl);
              const tuint curSID_prev = curSID;
              const flxVec y_chain_seed(list->get_seed_y_list(),list->get_Nrv());
              csm->prepareS(y_chain_seed);
              #if FLX_DEBUG
                try {
                  RndBox.set_smp(y_chain_seed);
                    // NOTE: the previous call might result in an error for nested updating, because the observed maxL of the prior (the input) might have changed
                  if (list->meth_id!=FlxBayUp_Update_List::RASUBSIM) {
                    const tdouble L = list->parent.eval_Likelihood();
                    if (fabs((L-list->get_seed_L_list())/L)>GlobalVar.TOL()*100) {
                      GlobalVar.slog(4) << std::endl << L << "\t" << list->get_seed_L_list() << "\t" << GlobalVar.Double2String(L-list->get_seed_L_list())<< std::endl;
                      throw FlxException_Crude("FlxBayUp_Update::update_b6");
                    }
                  } else {
                    if (fabs(list->parent.eval_RAlsf()-list->get_seed_s_list())>GlobalVar.TOL()*100) {
                      throw FlxException_Crude("FlxBayUp_Update::update_b7");
                    }
                  }
                  flxVec cy(list->get_Nrv());
                  list->get_RndBox().get_y_Vec(cy.get_tmp_vptr());
                  flxVec c1(list->get_seed_y_list(),list->get_Nrv());
                  if (c1!=cy) {
                    GlobalVar.slog(4) << std::endl << c1 << std::endl << cy << std::endl;
                    throw FlxException_Crude("FlxBayUp_Update::update_b8");
                  }
                  if (list->parent.get_methCat()==flxBayUp::BUS) {
                    if (cy[list->get_Nrv()-1]!=list->get_seed_p_list()) {
                      throw FlxException_Crude("FlxBayUp_Update::update_b9");
                    }
                    if ( fabs(list->eval_LSF(list->get_seed_p_list(),list->get_seed_L_list())-list->get_seed_s_list())>GlobalVar.TOL() ) {
                      throw FlxException_Crude("FlxBayUp_Update::update_b10");
                    }
                  }
                } catch (FlxException_NeglectInInteractive& e) {
                  FLXMSG("FlxBayUp_Update::update_b11",1);
                  GlobalVar.alert.alert("FlxBayUp_Update::update_b12",e.what());
                }
              #endif
              const tuint chainL = list->get_cur_chain_length();
              #if FLX_DEBUG
                if (chainL==0) throw FlxException_Crude("FlxBayUp_Update::update_b13");
              #endif
              for (tuint j=1;j<chainL;++j) {        // seed is included in list -> i=1 !!!
                // propose a new realization
                  const flxVec y_prev(list->get_seed_y_list(),list->get_Nrv());
                  bool rejectIt = !(csm->propose(y_prop,y_prev));
                  if (rejectIt) ++nac1;
                // test the new realization
                  if (rejectIt==false) {
                    try {
                      RndBox.set_smp(y_prop);
                      ++Nlsf_calls; 
                    } catch (FlxException_NeglectInInteractive &e) {
                      if (err_first) {
                        FLXMSG("FlxBayUp_Update::update_b14",1);
                        GlobalVar.alert.alert("FlxBayUp_Update::update_b15",e.what());
                        err_first = false;
                      }
                      ++Nlsf_errors;
                      rejectIt = true;
                    }
                  }
                const bool bt = list->insert_entry(false,rejectIt,false,false,os_smpl,ZERO,&adupd_acr);
                if (bt) {
                  ++adupd_acrN;
                } else {
                  ++nac;
                }
                csm->acceptance_feedback(bt);
              }  // end for over chain
              // some things that need to be done after the loop closed
                tc2 += chainL-1;
                adupd_total += chainL-1;
                curSID = curSID_prev;
              // adaptive control
                adupd_counter += chainL-1;
                if (adupd_counter_N>0 && adupd_counter>=adupd_counter_N) {
                  const tdouble iadpt_prev = iadpt;
                  iadpt = tdouble(adupd_ncalls);
                  list->get_adpt_ctrl().requires_adptv_step((adupt_acrP?adupd_acr:adupd_acrN)/tdouble(adupd_total),*csm);
                  iadpt = iadpt_prev;
                  adupd_counter = 0;
                  adupd_acrN = 0;
                  adupd_acr = ZERO;
                  adupd_total = 0;
                  ++adupd_ncalls;
                }
              prg.tick(i+1);
              list->set_next_seed();
            }        // end loop over chains
            if (os_smpl) {
              delete os_smpl;
              os_smpl = NULL;
            }
            prg.stop();        // deactivate progress-bar
            #if FLX_DEBUG
              if (tc2+Nc_now!=list->get_Ns_final() && tc2+Nc_now!=list->get_Nc()*list->get_Ncl()) throw FlxException_Crude("FlxBayUp_Update::update_b16");
              if (list->get_cur_i_list()>=0) throw FlxException_Crude("FlxBayUp_Update::update_b17");
            #endif
          const bool uc_res = (list->meth_id==FlxBayUp_Update_List::BUS || list->meth_id==FlxBayUp_Update_List::UBUS || list->meth_id==FlxBayUp_Update_List::BUST )?(list->update_c(p_mod,false)):false;
            if (list->meth_id==FlxBayUp_Update_List::UBUS) gt_is_zero = list->is_gt_zero();
          // output
            switch (list->meth_id) {
              case (FlxBayUp_Update_List::BUS):
              case (FlxBayUp_Update_List::UBUS):
              case (FlxBayUp_Update_List::BUST):
              {
                const tdouble mhml = list->get_maxL();
                if (mhml>=log(GlobalVar.TOL()) && mhml<=100.) {
                  scout << std::format("maxL={:9.2e}  ", exp(mhml)) << " ";
                } else {
                  scout << std::format("mlnL={:9.2e}  ", mhml) << " ";
                }
                break;
              }
              case (FlxBayUp_Update_List::ABCSUBSIM):
              case (FlxBayUp_Update_List::RASUBSIM):
              {
                // TODO
                break;
              }
              default:
                throw FlxException_Crude("FlxBayUp_Update::update_b18");
            };
          // output acceptance rate
            {
              const tdouble acr_val = ONE-tdouble(nac)/tdouble(tc2);
              scout << std::format("acr={:4.2f}  ", acr_val);
              if (susControl.verbose && acr_val>GlobalVar.TOL()) {
                if (list->get_Nrv()>1 && csm->get_ac1d()<0.995) {
                  scout_ext << std::format("ac1d={:4.2f} ", csm->get_ac1d());
                  scout_ext_nempty = true;
                }
                const tdouble acr1 = ONE-tdouble(nac1)/tdouble(tc2);
                if (acr1<0.995) {
                  scout_ext << std::format("acr1={:4.2f} ", (acr1));
                  scout_ext << std::format("acr2={:4.2f} ", (ONE-tdouble((nac-nac1))/tdouble(tc2-nac1)));
                  scout_ext_nempty = true;
                }
              }
            }
          // output velocity
            scout << std::format("velo={:4.2f}  ", list->get_velo(Nc_now));
          ++N_runs;
          {        // check whether we are done
            bool b_done = false;        
            if (gt_is_zero && list->is_fullList()) {
              if (list->meth_id==FlxBayUp_Update_List::BUS || list->meth_id==FlxBayUp_Update_List::UBUS || list->meth_id==FlxBayUp_Update_List::BUST ) {
                if (uc_res) b_done = true; // ... exit only if maxL did not change
              } else {        // otherwise ... exit the loop
                b_done = true;
              }
            } else {
              if (N_runs>=list->get_max_runs()) {
                b_done = true;
              }
            }
            if (b_done) {
              if (susControl.verbose) {
                list->print_ext_out(scout);
                if (scout_ext_nempty) scout << scout_ext.str();
                scout_ext.str("");
                scout_ext << std::endl << "            ";
              }
              break;
            }
          }
        }
        } catch (FlxException &e) {
          delete curLevelStat;
          if (curLevelStat_next) delete curLevelStat_next;
          throw;
        }
        scout << std::endl;
      delete curLevelStat;
      if (curLevelStat_next) delete curLevelStat_next;
      if (list->is_gt_zero()==false) {
        std::ostringstream ssV;
        ssV << "Simulation did not converge within " << N_runs << " runs.";
        throw FlxException_NeglectInInteractive("FlxBayUp_Update::update_b19",ssV.str());
      }
      post_adptcount_N = list->get_adpt_ctrl().get_updatesAfterNsamples();
      break;
    }
    case (FlxBayUp_Update_List::TMCMC):
    {
    // --------------------------------------------------------------------------------------------------------------------------
    // TMCMC
    // --------------------------------------------------------------------------------------------------------------------------
      const tdouble target_cov = susControl.TMCMC_target_COV->cast2positive(false);
      // paramters that belong to the improved version
        const tuint update_weights = susControl.TMCMC_update_weights;
        const tdouble tmcmc_alpha = susControl.TMCMC_alpha->cast2positive_or0(false);
      tdouble& l_max = list->parent.get_cStart();
      FlxBayUP_csm_TMCMC* csm_TMCMC = dynamic_cast<FlxBayUP_csm_TMCMC*>(csm);
        if (csm_TMCMC==NULL) {
          std::ostringstream ssV;
          ssV << "Only MCMC-method 'tmcmc' allowed in combination with the TMCMC method.";
          throw FlxException("FlxBayUp_Update::update_t2",ssV.str());
        }
      bool err_first = susControl.prt_alert;
      const tuint Nc = list->get_Nc();
      const tuint Nb = list->get_Nburn();
      const tuint NOX = list->get_NOX();
      const tuint NRV = list->get_Nrv();
      // allocate temporary memory
        tdouble* L_list1 = new tdouble[Nc];
        tdouble* L_list2 = new tdouble[Nc];
        tdouble* y_list1 = new tdouble[Nc*NRV];
        tdouble* y_list2 = new tdouble[Nc*NRV];
        tdouble* x_list1 = new tdouble[Nc*NOX];
        tdouble* x_list2 = new tdouble[Nc*NOX];
        tuint* n_list = NULL;                // how many MCMC transition has the current seed/sample already untertaken?
        tdouble* tmcmc_fy_list = NULL;
        tdouble *tmcmc_aq_list = NULL;
        tdouble *tmcmc_aw_list = NULL;
        flxVec W_vec(Nc);
        flxpVec y_meanP(NRV);
        flxVec y_mean(NRV);
        flxVec y_covV((NRV*NRV+NRV)/2);
      // iTMCMC: perform a prior iteration
        tuint tmcmc_iN = (update_weights>0&&tmcmc_alpha>GlobalVar.TOL())?(Nc+Nb):0;
        if (tmcmc_iN>0) {
          // determin number of required transitions for convergence to target distribution
            if (ONE/tmcmc_alpha<tmcmc_iN) tmcmc_iN = ONE/tmcmc_alpha;
            if (tmcmc_iN<1) tmcmc_iN = 1;
          // allocate memory
            n_list = new tuint[Nc];
            tmcmc_aq_list = new tdouble[tmcmc_iN+1];
            tmcmc_aw_list = new tdouble[tmcmc_iN+1];
            if (update_weights<2) {
              tmcmc_fy_list = new tdouble[Nc];
            }
        }
      try {
      // do MCI in a first step
        scout << std::format("  iter {:2d}:", N_runs) << "  ";
        // register the progress-bar
          op.flush();
          prg.start(list->get_Nc());
        {
          tuint i=0; 
          const tuint NsM = list->get_Ncl();
          while (i<Nc) {
            try {
              RndBox.gen_smp(); 
              L_list1[i] = list->parent.eval_Likelihood();
              RndBox.get_y_Vec(y_list1+NRV*i);
              RndBox.get_x_Vec(x_list1+NOX*i);
              ++Nlsf_calls;
            } catch (FlxException_NeglectInInteractive &e) {
              FLXMSG("FlxBayUp_Update::update_t1",1);
              if (err_first) {
                GlobalVar.alert.alert("FlxBayUp_Update::update_t3",e.what());
                err_first = false;
              }
              ++Nlsf_errors;
              continue;
            }
            ++i;
            if (i%NsM==0) prg.tick(i/NsM);
          }
        }
        // deactivate the progress-bar
          prg.stop();
        tdouble* L_list_prev = L_list1;
        tdouble* L_list_now = L_list2;
        tdouble* y_list_prev = y_list1; 
        tdouble* y_list_now = y_list2; 
        tdouble* x_list_prev = x_list1; 
        tdouble* x_list_now = x_list2;
        tdouble* final_L_list = list->TMCMC_get_L_list();
        tdouble* final_y_list = list->TMCMC_get_y_list();
        tdouble* final_x_list = list->TMCMC_get_x_list();
        tdouble q_prev = ZERO;
        tdouble q_now = ZERO;                // potence of the likelihood function
      // perform a conditioned sampling
        while (true) {
          flxVec L_vec_prev(L_list_prev,Nc);
          // output samples of the last step
            os_smpl = open_smpl_file4write();
            if (os_smpl) {
              const tuint Nct = (fabs(q_now-ONE)<=GlobalVar.TOL())?(list->get_Ns_final()):Nc;
              for (tuint i=0;i<Nct;++i) {
                *os_smpl << GlobalVar.D2S_totalPrec(tdouble(i)) << GlobalVar.D2S_totalPrec(tdouble(i)) <<  std::endl;
                const flxVec sy(&(y_list_prev[i*NRV]),NRV);        // y-vector to write
                  *os_smpl << "  ";
                  flxVec_totalPrec_plot(*os_smpl,sy);
                  *os_smpl << std::endl;
                const flxVec sx(&(x_list_prev[i*NOX]),NOX);        // x-vector to write
                  *os_smpl << "  ";
                  flxVec_totalPrec_plot(*os_smpl,sx);
                  *os_smpl << std::endl;
              }
              delete os_smpl;
              os_smpl = NULL;
            }
            {
              const tdouble l_max_level = L_vec_prev.get_max();
              if (l_max_level>l_max) l_max = l_max_level;
              if (exp(l_max_level)>=GlobalVar.TOL() && l_max_level<=100.) {
                scout << std::format("maxL={:9.2e}  ", exp(l_max_level));
              } else {
                scout << std::format("mlnL={:9.2e}  ", l_max_level) << " ";
              }
            }
          // check stopping condition
            if (fabs(q_now-ONE)<=GlobalVar.TOL() || N_runs>=list->get_max_runs()) {
              scout << std::endl;
              if (N_runs==0) {        // prior=posterior
                std::ostringstream ssV;
                ssV << "Posterior equal the prior. This case is not handled in TMCMC!";
                throw FlxException("FlxBayUp_Update::update_t10",ssV.str());
              }
              W_vec = ONE;
              TMCMC_assemble_smplCOV(y_covV,W_vec,y_list_prev,list->get_Ns_final(),NRV,y_meanP,y_mean,y_prop);
              csm_TMCMC->prepare(y_covV);
              break;
            }
          // determin the current Likelihood 'scaling' (its potence)
            q_now = TMCMC_new_q(q_prev,target_cov,L_vec_prev,W_vec);
            if (q_now!=q_now) throw FlxException_Crude("FlxBayUp_Update::update_t10a");
            scout << std::format("q={:9.2e}  ", q_now);
            const tuint Nct = (fabs(q_now-ONE)<=GlobalVar.TOL())?(list->get_Ns_final()):Nc;
          // check if this is the last level
            if (fabs(q_now-ONE)<=GlobalVar.TOL()) {
              L_list_now = final_L_list;
              y_list_now = final_y_list;
              x_list_now = final_x_list;
            }
          // compute mean of Likelihoods
            const tdouble W_mean = W_vec.get_Mean();
            p_mod += log(W_mean);
            scout << std::format("Sk={:9.2e}  ", W_mean);
          // Estimate Covariance matrix
            TMCMC_assemble_smplCOV(y_covV,W_vec,y_list_prev,Nc,NRV,y_meanP,y_mean,y_prop);
          // prepare the intermediate weight estimation
            if (tmcmc_iN>0) {
              // tmcmc_aq_list
                const tdouble* Lvpp = L_vec_prev.get_tmp_vptr_const();
                tmcmc_aq_list[0] = q_prev;
                const tdouble dq2 = (q_now - q_prev)/tmcmc_iN;
                for (tuint i=1;i<tmcmc_iN;++i) {
                  tmcmc_aq_list[i] = tmcmc_aq_list[i-1] + dq2;
                }
                tmcmc_aq_list[tmcmc_iN] = q_now;
              if (update_weights<2 || tmcmc_iN==0) {
                // tmcmc_aw_list
                  tmcmc_aw_list[0] = ONE;
                  for (tuint i=1;i<tmcmc_iN;++i) {
                    tdouble w_tmp = ZERO;
                    size_t N_tmp = W_vec.get_N();
                    const tdouble dq = tmcmc_aq_list[i]-q_prev;
                    for (size_t j=0;j<N_tmp;++j) {
                      w_tmp += exp(Lvpp[j]*dq);
                    }
                    tmcmc_aw_list[i] = w_tmp/N_tmp;
                  }
                  tmcmc_aw_list[tmcmc_iN] = W_mean;
                // tmcmc_fy_list
                  for (tuint i=0;i<Nc;++i) {
                    tmcmc_fy_list[i] = exp(Lvpp[i]*q_prev);
                  }
              } else {
                for (tuint i=0;i<=tmcmc_iN;++i) {
                  tmcmc_aw_list[i] = ONE;
                }
              }
              // n_list
                for (tuint i=0;i<Nc;++i) {
                  n_list[i] = 0;
                }
            }
          // some preliminarities 
            list->get_adpt_ctrl().eval();
            ++N_runs;
            *iterN = N_runs;
            if (list->get_adpt_ctrl().is_adaptive() || N_runs>1) {
              csm_TMCMC->write_adaptive_info(scout,list->get_adpt_ctrl().is_adaptive());
            }
            scout << std::endl;
            scout << "  " << std::format("iter {:2d}:  ", N_runs);
            // visibility of progress-bar
              op.flush();
          // ... the actual generation of new samples
            // counts the total number of samples that were not accepted
              tuint nac = 0;
            prg.start(Nb+Nct);        // activate progress-bar
            // gives the number of the calls of the adaptive step (exports it for the user to work with)
              tuint adupd_ncalls = 1;
            // counters for the adaptive process
              const tuint adupd_counter_N = list->get_adpt_ctrl().get_updatesAfterNsamples();
              // counts the number of samples since the last adaptive-step (to decide whether and adaptive step is needed)
                tuint adupd_counter = 0;
              // number of accepted samples since the last adaptive-step
                tdouble acr_sum = ZERO;
            csm_TMCMC->prepare(y_covV);
            for (tuint i=0;i<Nb+Nct;++i) {
              #if FLX_DEBUG
                if (W_vec.get_min()<ZERO) throw FlxException_Crude("FlxBayUp_Update::update_t11a");
                if (W_vec.get_sum()<=ZERO) throw FlxException_Crude("FlxBayUp_Update::update_t11b");
              #endif
              // select a sample at random
                tuint indx;
                tdouble indx_pr = ZERO;
                if ( update_weights<2 || tmcmc_iN==0) {
                  indx = data->RndCreator.gen_smp_index2(W_vec);
                } else {
                  // select the transition group
                    for (indx=0;indx<tmcmc_iN;++indx) {
                      if (tmcmc_aw_list[indx]<=GlobalVar.TOL()) continue;
                      if (data->RndCreator.gen_smp_binary(tmcmc_aw_list[indx])) break;
                    }
                    #if FLX_DEBUG
                      if (indx==tmcmc_iN && tmcmc_aw_list[indx]!=ONE) throw FlxException_Crude("FlxBayUp_Update::update_t12");
                    #endif
                  // select the index (according to its weight) within the selected transition group
                    // compute the normalizing constant
                      tdouble wnt = ZERO;
                      for (tuint j=0;j<Nc;++j) {
                        if (n_list[j]==indx) {
                          wnt += W_vec[j];
                        }
                      }
                    const tdouble wntp = data->RndCreator.gen_smp_uniform()*wnt;
                    tdouble wnt2 = ZERO;
                    for (tuint j=0;j<Nc;++j) {
                      if (n_list[j]==indx) {
                        wnt2 += W_vec[j];
                        if (wntp<=wnt2) {
                          indx = j;
                          indx_pr = W_vec[j]/wnt;
                          break;
                        }
                      }
                    }
                    #if FLX_DEBUG
                      if (wntp>wnt) throw FlxException_Crude("FlxBayUp_Update::update_t13a");
                      if (indx_pr==ZERO) throw FlxException_Crude("FlxBayUp_Update::update_t13b");
                      if (indx>=Nc) throw FlxException_Crude("FlxBayUp_Update::update_t13c");
                    #endif
                }
              // retrieve the selected sample
                flxVec sy0(&(y_list_prev[indx*NRV]),NRV);        // y-seed
                flxVec sx0(&(x_list_prev[indx*NOX]),NOX);        // x-seed
                if (update_weights==2 && i>=Nb) {        // accept this sample as sample from the target distribution
                  flxVec y_now(&(y_list_now[(i-Nb)*NRV]),NRV);
                  flxVec x_now(&(x_list_now[(i-Nb)*NOX]),NOX);
                  L_list_now[i-Nb] = L_list_prev[indx];
                  y_now = sy0;
                  x_now = sx0;
                }
              // propose a new sample
                csm_TMCMC->propose(y_prop,sy0);
                bool rejectIt = false;
                tdouble l=ZERO;
                try {
                  RndBox.set_smp(y_prop);
                  l = list->parent.eval_Likelihood();
                  ++Nlsf_calls; 
                } catch (FlxException_NeglectInInteractive &e) {
                  if (err_first) {
                    FLXMSG("FlxBayUp_Update::update_t14",1);
                    GlobalVar.alert.alert("FlxBayUp_Update::update_t15",e.what());
                    err_first = false;
                  }
                  ++Nlsf_errors;
                  rejectIt = true;
                }
              // compute acceptance ratio
                if (!rejectIt) {
                  double alpha = ONE;
                  for (tuint j=0;j<NRV;++j) {
                    alpha *= rv_phi(y_prop[j])/rv_phi(sy0[j]);
                  }
                  alpha *= exp(q_now*(l-L_list_prev[indx]));
                  acr_sum += (alpha>ONE)?ONE:alpha;
                  rejectIt = (RndCreator.gen_smp_uniform()>alpha);
                }
                ++adupd_counter;
              // accept or reject the proposed sample
                const tdouble l_prev = L_list_prev[indx];
                if (!rejectIt) {
                  L_list_prev[indx] = l;
                  sy0 = y_prop;
                  RndBox.get_x_Vec(sx0.get_tmp_vptr());
                } else {
                  ++nac;
                }
              // update the weights
                if (update_weights>0) {
                  l = L_list_prev[indx];
                  if (tmcmc_iN>0) {
                    if (update_weights<2) {
                      tuint k_prev = (n_list[indx]>=tmcmc_iN)?tmcmc_iN:(n_list[indx]);
                      ++(n_list[indx]);
                      tuint k_now = (n_list[indx]>=tmcmc_iN)?tmcmc_iN:(n_list[indx]);
                      if (k_prev!=k_now) {
                        for (tuint j=0;j<Nc;++j) {
                          if (L_list_prev[j]<ZERO && std::isinf(L_list_prev[j])) {
                            W_vec[j] = ZERO;
                          } else {
                            tmcmc_fy_list[j] += (
                                exp(L_list_prev[j]*tmcmc_aq_list[k_now])/tmcmc_aw_list[k_now]
                                -exp(L_list_prev[j]*tmcmc_aq_list[k_prev])/tmcmc_aw_list[k_prev]
                              )/Nc;
                            W_vec[j] = exp(L_list_prev[j]*q_now)/tmcmc_fy_list[j];
                          }
                        }
                      }
                      if (l_prev!=l) {
                        tmcmc_fy_list[indx] = ZERO;
                        for (tuint j=0;j<Nc;++j) {
                          tuint k = (n_list[j]>=tmcmc_iN)?tmcmc_iN:(n_list[j]);
                          tmcmc_fy_list[indx] += exp(l*tmcmc_aq_list[k])/tmcmc_aw_list[k];
                        }
                        tmcmc_fy_list[indx] /= Nc;
                        W_vec[indx] = exp(l*q_now)/tmcmc_fy_list[indx];
                      }
                    } else {
                      if (n_list[indx]<tmcmc_iN) ++(n_list[indx]);
                      // correct the transition weights
                        const tdouble nxtfull = ONE-tmcmc_aw_list[n_list[indx]-1];
                        const tdouble raus = tmcmc_aw_list[n_list[indx]-1] * indx_pr;
                        const tdouble nxtraus = nxtfull * (ONE-tmcmc_aw_list[n_list[indx]]);
                        // 'correct' the previous bin
                          tmcmc_aw_list[n_list[indx]-1] *= (ONE-indx_pr);
                        // 'correct' the current bin
                          tmcmc_aw_list[n_list[indx]] = ONE-nxtraus/(nxtfull+raus);
                      W_vec[indx]  = exp(l*(q_now-tmcmc_aq_list[n_list[indx]]));
                    }
                  } else {
                    W_vec[indx]  = exp(l*(q_now-q_prev));
                  }
                }
              // in case burn-in is over
                if (update_weights<2 && i>=Nb) {
                  flxVec y_now(&(y_list_now[(i-Nb)*NRV]),NRV);
                  flxVec x_now(&(x_list_now[(i-Nb)*NOX]),NOX);
                  L_list_now[i-Nb] = L_list_prev[indx];
                  y_now = sy0;
                  x_now = sx0;
                }
              // adaptive control
                if (adupd_counter_N>0 && adupd_counter>=adupd_counter_N) {
                  const tdouble iadpt_prev = iadpt;
                  iadpt = tdouble(adupd_ncalls);
                  acr_sum /= adupd_counter;
                  list->get_adpt_ctrl().requires_adptv_step(acr_sum,*csm);
                  iadpt = iadpt_prev;
                  adupd_counter = 0;
                  acr_sum = ZERO;
                  ++adupd_ncalls;
                }
              prg.tick(i+1);
            }        // end loop over chains
            prg.stop();        // deactivate progress-bar
          // output acceptance rate
            {
              const tdouble acr_val = ONE-tdouble(nac)/tdouble(Nct+Nb);
              scout << std::format("acr={:4.2f}  ", acr_val);
            }
          // swap variables
            tdouble* ptrswap;
            ptrswap = L_list_prev;
              L_list_prev = L_list_now;
              L_list_now = ptrswap;
            ptrswap = y_list_prev;
              y_list_prev = y_list_now;
              y_list_now = ptrswap;
            ptrswap = x_list_prev;
              x_list_prev = x_list_now;
              x_list_now = ptrswap;
          q_prev = q_now;
        }
        if (fabs(q_now-ONE)>GlobalVar.TOL()) {
          std::ostringstream ssV;
          ssV << "Simulation did not converge within " << N_runs << " runs.";
          throw FlxException_NeglectInInteractive("FlxBayUp_Update::update_b19",ssV.str());
        }
      } catch (FlxException &e) {
        delete [] L_list1;
        delete [] L_list2;
        delete [] y_list1; 
        delete [] y_list2; 
        delete [] x_list1; 
        delete [] x_list2; 
        if (n_list) delete [] n_list; 
        if (tmcmc_fy_list) delete [] tmcmc_fy_list;
        if (tmcmc_aq_list) delete [] tmcmc_aq_list;
        if (tmcmc_aw_list) delete [] tmcmc_aw_list;
        throw;
      }
      delete [] L_list1;
      delete [] L_list2;
      delete [] y_list1; 
      delete [] y_list2; 
      delete [] x_list1; 
      delete [] x_list2; 
      if (n_list) delete [] n_list;
      if (tmcmc_fy_list) delete [] tmcmc_fy_list;
      if (tmcmc_aq_list) delete [] tmcmc_aq_list;
      if (tmcmc_aw_list) delete [] tmcmc_aw_list;
      post_adptcount_N = list->get_adpt_ctrl().get_updatesAfterNsamples();
      break;
    }
    case FlxBayUp_Update_List::MHRS:
    case (FlxBayUp_Update_List::RS):
    {
    // --------------------------------------------------------------------------------------------------------------------------
    // Rejectional sampling (of BUS problem)
    // --------------------------------------------------------------------------------------------------------------------------
      // register the progress-bar
      op.flush();
      const tuint Nf = list->get_Ns_final();
      prg.start(Nf);
      const tuint N_RV = y_prop.get_N();
      const tuint N_OX = list->get_NOX();
      os_smpl = open_smpl_file4write();
      std::ofstream* os_smpl1 = (list->meth_id==FlxBayUp_Update_List::RS)?os_smpl:NULL;
      const tdouble c = list->parent.get_cStart();
      for (tuint i=0;i<Nf;++i) {
        tdouble ts=ONE, tL=ZERO;
        bool err; bool err_first=susControl.prt_alert;
        do {        // loop until we can accept smpl
          try {
            // draw a realization from the prior
              RndBox.gen_smp(); 
            // check if proposed sample is in failure domain
              tL = list->parent.eval_Likelihood();
              RndBox.get_y_Vec(y_prop.get_tmp_vptr());
              ts = list->eval_LSF(y_prop[N_RV-1],tL);
              ++Nlsf_calls;
              err = false;
          } catch (FlxException &e) {
            FLXMSG("FlxBayUp_Update::update_c1",1);
            if (err_first) {
              GlobalVar.alert.alert("FlxBayUp_Update::update_c2",e.what());
              err_first = false;
            }
            ++Nlsf_errors;
            err = true;
            continue;
          }
        } while (ts>ZERO || err);
        list->insert_entry(false,false,false,false,os_smpl1,tL);
        prg.tick(i+1);
      }
      // deactivate the progress-bar
        prg.stop();
      // --------------------------------------------------------------
      // ... perform the Metropolis-Hastings step ... IF (meth=MHRS)
      // --------------------------------------------------------------
      p_mod = log(tdouble(Nf)/tdouble(Nlsf_calls));
      if (list->meth_id==FlxBayUp_Update_List::MHRS) {
        p_mod += c;
      }
      if (list->meth_id==FlxBayUp_Update_List::MHRS && list->get_maxL()>c) {
        // compute the likelihood ratios
          flxVec wvec(Nf);
          tdouble* L_list = list->TMCMC_get_L_list();
          for (tuint i=0;i<Nf;++i) {
            wvec[i] = (L_list[i]>c)?exp(L_list[i]-c):ONE;
          }
        // correct the estimate of the evidence
          const tdouble wmean = wvec.get_sum()/Nf;
          p_mod += log(wmean);
        // pick the initial seed of the Metropolis-Hastings chain
          tdouble* y_list = list->TMCMC_get_y_list();
          tdouble* x_list = list->TMCMC_get_x_list();
          {
            tuint i = RndCreator.gen_smp_index2(wvec);
            // y-vector
              {
                flxVec yi(y_list+i*N_RV,N_RV);
                y_prop = yi;
                flxVec y0(y_list,N_RV);
                yi = y0;
                y0 = y_prop;
              }
            // x-vector
              {
                flxVec x_prop(N_OX);
                flxVec xi(x_list+i*N_OX,N_OX);
                x_prop = xi;
                flxVec x0(x_list,N_OX);
                xi = x0;
                x0 = x_prop;
              }
            // likelihood
              tdouble Lt = L_list[i];
              Lt = L_list[i];
              L_list[i] = L_list[0];
              L_list[0] = Lt;
              y_list[N_RV-1] = rv_InvPhi_noAlert(RndCreator.gen_smp_uniform()*exp(Lt-c));
              x_list[N_OX-1] = y_list[N_RV-1];
          }
        // perform the Metropolis-Hastings step
          tdouble tc = ZERO;
          for (tuint i=1;i<Nf;++i) {
            flxVec yi(y_list+i*N_RV,N_RV);
            flxVec xi(x_list+i*N_OX,N_OX);
            // accept/reject 'proposed' sample
              const tdouble ur = exp(((L_list[i]>c)?(L_list[i]):c)-L_list[i-1]);
              if (ur<ONE) {
                tc += ur;
                if (RndCreator.gen_smp_uniform()>ur) {        // DO NOT accept the current sample
                  flxVec yp(y_list+(i-1)*N_RV,N_RV);
                  flxVec xp(x_list+(i-1)*N_OX,N_OX);
                  L_list[i] = L_list[i-1];
                  yi = yp;
                  xi = xp;
                }
              } else {
                tc += ONE;
              }
            // resample the auxiliary random variable
              yi[N_RV-1] = rv_InvPhi_noAlert(RndCreator.gen_smp_uniform()*exp(L_list[i]-c));
              xi[N_OX-1] = yi[N_RV-1];
            if (os_smpl) list->write_smpl(i,*os_smpl);
          }
          // output HM-acceptance rate
            tc /= Nf;
            scout << "  E[scaled likelihood]:       " << GlobalVar.Double2String(wmean,false,2) << std::endl;
            scout << "  MH-acceptance rate:         " << GlobalVar.Double2String(tc,false,2) << std::endl;
              (*data->ConstantBox.get("bayup_mhrs_acr",true)) = tc;
      }
      if (os_smpl) {
        delete os_smpl;
        os_smpl = NULL;
      }
      break;
    }
    case (FlxBayUp_Update_List::ABCRS):
    {
    // --------------------------------------------------------------------------------------------------------------------------
    // Rejectional sampling (of ABC problem)
    // --------------------------------------------------------------------------------------------------------------------------
      // register the progress-bar
      op.flush();
      prg.start(list->get_Ns_final());
      const tdouble thr_m = list->parent.get_cStart();
      os_smpl = open_smpl_file4write();
      for (tuint i=0;i<list->get_Ns_final();++i) {
        tdouble metric=ONE;
        bool err; bool err_first=susControl.prt_alert;
        do {        // loop until we can accept smpl
          try {
            // draw a realization from the prior
              RndBox.gen_smp(); 
            // check if proposed sample is in failure domain
              metric = list->parent.eval_Likelihood();        // evaluates the ABC-metric (returned as log-transform)
              ++Nlsf_calls;
              err = false;
          } catch (FlxException &e) {
            FLXMSG("FlxBayUp_Update::update_d1",1);
            if (err_first) {
              GlobalVar.alert.alert("FlxBayUp_Update::update_d2",e.what());
              err_first = false;
            }
            ++Nlsf_errors;
            err = true;
            continue;
          }
        } while (metric>thr_m || err);
        list->insert_entry(false,false,false,false,os_smpl); 
        prg.tick(i+1);
      }
      if (os_smpl) {
        delete os_smpl;
        os_smpl = NULL;
      }
      // deactivate the progress-bar
        prg.stop();
      p_mod = log(tdouble(list->get_Ns_final())/tdouble(Nlsf_calls));      
      break;
    }
    case (FlxBayUp_Update_List::RAMCI):
    {
    // --------------------------------------------------------------------------------------------------------------------------
    // Rejectional sampling (of RA problem)
    // --------------------------------------------------------------------------------------------------------------------------
      // register the progress-bar
      op.flush();
      prg.start(list->get_Ns_final());
      os_smpl = open_smpl_file4write();
      for (tuint i=0;i<list->get_Ns_final();++i) {
        tdouble lsf=ONE;
        bool err; bool err_first=susControl.prt_alert;
        do {        // loop until we can accept smpl
          try {
            // draw a realization from the prior
              RndBox.gen_smp(); 
            // check if proposed sample is in failure domain
              lsf = list->parent.eval_RAlsf();        // evaluates the limit-state function
              ++Nlsf_calls;
              err = false;
          } catch (FlxException &e) {
            FLXMSG("FlxBayUp_Update::update_e1",1);
            if (err_first) {
              GlobalVar.alert.alert("FlxBayUp_Update::update_e2",e.what());
              err_first = false;
            }
            ++Nlsf_errors;
            err = true;
            continue;
          }
        } while (lsf>ZERO || err);
        list->insert_entry(false,false,false,false,os_smpl); 
        prg.tick(i+1);
      }
      if (os_smpl) {
        delete os_smpl;
        os_smpl = NULL;
      }
      // deactivate the progress-bar
        prg.stop();
      p_mod = log(tdouble(list->get_Ns_final())/tdouble(Nlsf_calls));      
      break;
    }
    case (FlxBayUp_Update_List::LS):
    {
    // --------------------------------------------------------------------------------------------------------------------------
    // Line sampling
    // --------------------------------------------------------------------------------------------------------------------------
      const tuint NLS = list->get_Nburn();
      const tuint NOX = list->get_NOX();
      const tuint NRV = list->get_Nrv();
      const tdouble LS_tol = susControl.LS_tol->cast2positive(true);
      const tuint LS_max_iter = susControl.LS_max_iter->cast2tuint(true);
      scout << "  number of line-searches:      " << NLS << std::endl;
      scout << "  tolerance of line-search:     " << LS_tol << std::endl;
      scout << "  max. line-search iterations:  " << LS_max_iter << std::endl;
      scout << " --------------------------------------------------------------------" << std::endl;
      flxVec rvy(NRV);
      flxVec rvy_dummy(NRV);
      flxVec rvx(NOX);
      // optain the beta-vector
        flxVec betaVec(data->ConstMtxBox.get_Vec(NRV,susControl.LS_SPNT->eval(),true),NRV,true);
      // perform an initial line search (to fix the actual tolerance parameter)
        {
          flx_LS_line_prop line_0 = perform_line_search(rvy,rvy_dummy,betaVec,LS_tol,LS_max_iter,Nlsf_calls);
          if (line_0.is_topo_set()) {
            const tdouble lscl = line_0.get_lower_scale();
            if (lscl>ONE) betaVec *= lscl;
          }
        }
        const tdouble betaNorm2 = betaVec.get_Norm2_NOroot();
        const tdouble betaNorm = sqrt(betaNorm2);
        std::valarray<flx_LS_line_prop> lines(NLS);
        flxVec lines_Pr(NLS);
      // start with the iteration (line searches)
        flxVec ycol(NRV*NLS);
        tdouble* ycolp = ycol.get_tmp_vptr();
        // register the progress-bar
          scout << "  preform line searches ...      ";
          prg.start(NLS);
        const tuint N_LSF_calls_1 = Nlsf_calls;
        for (tuint i=0;i<NLS;++i) {
          // propose a random vector from the prior (standard Normal)
            RndCreator.gen_smp(rvy);
          // create in-plane vector
            const tdouble c = -(rvy.operator*(betaVec)/betaNorm2);
            rvy.add(betaVec,c);
            flxVec yt(ycolp+i*NRV,NRV);
            yt = rvy;
          // perform the line search
            flx_LS_line_prop &line_i = lines[i];
            line_i = perform_line_search(rvy,rvy_dummy,betaVec,LS_tol,LS_max_iter,Nlsf_calls);
            lines_Pr[i] = line_i.get_upper_Pr(betaNorm);
            //line_i.print_topo(scout);
          prg.tick(i+1);
        }
        // deactivate the progress-bar
          prg.stop();
          scout << "[" << GlobalVar.Double2String(tdouble(Nlsf_calls-N_LSF_calls_1)/tdouble(NLS),false,2) << " calls per line]" << std::endl;
      // check if any line-searches where successful ..
          {
            const tdouble Not0sum = lines_Pr.get_sum();
            if (Not0sum==ZERO) {
              std::ostringstream ssV;
              ssV << "None of the line-searches returned a probability larger than zero.";
              throw FlxException("FlxBayUp_Update::update_f0", ssV.str() );
            }
            tuint Nnot0 = 0;
            for (tuint i=0;i<NLS;++i) {
              if (lines_Pr[i]/Not0sum>GlobalVar.TOL()) ++Nnot0;
            }
            scout << "                                 [" << Nnot0 << " lines with Pr>0]" << std::endl;
          }
      // generate the posterior samples
        // register the progress-bar
          scout << "  generate posterior samples ... ";
          prg.start(list->get_Ns_final());
        const tuint N_LSF_calls_2 = Nlsf_calls;
        for (tuint i=0;i<list->get_Ns_final();++i) {
          tdouble Lval=ZERO;
          tdouble lsf=ONE;
          bool err; bool err_first=susControl.prt_alert;
          do {        // loop until we can accept smpl
            // select index from list (at random)
              const tuint j = RndCreator.gen_smp_index2(lines_Pr);
              flx_LS_line_prop &line_j = lines[j];
            // propose a point on line
              const tdouble c = line_j.propose(RndCreator.gen_smp_uniform(),betaNorm);
            try {
              // evaluate its lsf
                flxVec yt(ycolp+j*NRV,NRV);
                lsf = line_search_LSF_call(c,yt,rvy_dummy,betaVec,Nlsf_calls,line_j,Lval);
                err = false;
            } catch (FlxException &e) {
              FLXMSG("FlxBayUp_Update::update_f1",1);
              if (err_first) {
                GlobalVar.alert.alert("FlxBayUp_Update::update_f2",e.what());
                err_first = false;
              }
              ++Nlsf_errors;
              err = true;
              continue;
            }
            // re-evaluate weight of line
              lines_Pr[j] = line_j.get_upper_Pr(betaNorm);
          } while (lsf>ZERO || err);
          list->insert_entry(false,false,false,false,os_smpl,Lval); 
          prg.tick(i+1);
        }
        // deactivate the progress-bar
          prg.stop();
          scout << "[" << GlobalVar.Double2String(tdouble(Nlsf_calls-N_LSF_calls_2)/tdouble(list->get_Ns_final()),false,2) << " calls per line]" << std::endl;
        // determin probability
          p_mod = ZERO;
          for (tuint i=0;i<NLS;++i) {
            flx_LS_line_prop &line_i = lines[i];
            p_mod += line_i.get_mean_Pr(betaNorm);
          }
          p_mod = log(p_mod/NLS);  
        break;
    }
    default:
    {
      throw FlxException_NotImplemented("FlxBayUp_Update::update_g1");
    }
  };

  // --------------------------------------------------------------------------------------------------------------------------
  // General post-processing for all methods
  // --------------------------------------------------------------------------------------------------------------------------
  // finalize the list
    const tuint Ns_final_Nneg = (list->parent.is_subsetRel?0:list->finalize());
  // finish with some output
    scout << " --------------------------------------------------------------------" << std::endl;
    std::ostream& scout3 = GlobalVar.slogcout(3);
    N_tot_LSFc_sim = Nlsf_calls+Nlsf_errors;
    scout3 << " Total number of LSF-calls:   " << N_tot_LSFc_sim << std::endl;
    (*data->ConstantBox.get("sus_lsf_calls",true)) = N_tot_LSFc_sim;
    if (Nlsf_errors>0) {
    scout3 << " Unsuccessful LSF-calls:      " << Nlsf_errors << std::endl;
    }
    if (list->parent.is_subsetRel) {
      PrMod = p_mod;
      scout3 << " Probability of failure:      " << GlobalVar.Double2String(exp(PrMod)) << " (" << GlobalVar.Double2String(PrMod) << ")" << std::endl;
    } else {
      scout3 << " Posterior seed values:       " << Ns_final_Nneg << std::endl;
      if (
        list->meth_id!=FlxBayUp_Update_List::ABCSUBSIM 
        && list->meth_id!=FlxBayUp_Update_List::RASUBSIM 
        && list->meth_id!=FlxBayUp_Update_List::ABCRS
        && list->meth_id!=FlxBayUp_Update_List::RAMCI
      ) {
        const tdouble lnc = (list->meth_id==FlxBayUp_Update_List::BUS||list->meth_id==FlxBayUp_Update_List::UBUS||list->meth_id==FlxBayUp_Update_List::BUST||list->meth_id==FlxBayUp_Update_List::MHRS)?(list->get_maxL()):(list->parent.get_cStart());
        if (list->parent.cStart_has_changed()) {
        scout3 << " max. observed Likel. val.:   ";
        } else {
        scout3 << " upper bound for Likel. val.: ";
        }
        scout3 << GlobalVar.Double2String(exp(lnc)) << " (" << GlobalVar.Double2String(lnc) << ")" << std::endl;
        PrMod = p_mod;
        if (list->meth_id!=FlxBayUp_Update_List::TMCMC && list->meth_id!=FlxBayUp_Update_List::MHRS) {
          PrMod += list->parent.get_cStart();
        }
        (*data->ConstantBox.get("bayup_mlnl",true)) = lnc;
      } else {
        PrMod = p_mod;
      }
      scout3 << " f(Observations|Model):       " << GlobalVar.Double2String(exp(PrMod)) << " (" << GlobalVar.Double2String(PrMod) << ")" << std::endl;
      if (
        list->meth_id!=FlxBayUp_Update_List::ABCSUBSIM 
        && list->meth_id!=FlxBayUp_Update_List::RASUBSIM 
        && list->meth_id!=FlxBayUp_Update_List::ABCRS
        && list->meth_id!=FlxBayUp_Update_List::RAMCI
      ) {
        DataFit = list->expectation_likelihood();
        relEntr = DataFit-PrMod;
        scout3 << "   data fit:                  " << GlobalVar.Double2String(DataFit) << std::endl;
        scout3 << "   relative entropy:          " << GlobalVar.Double2String(relEntr) << std::endl;
      }
    }
    (*data->ConstantBox.get("sus_pr",true)) = PrMod;
  // output credible intervals
    if (csm && list->meth_id!=FlxBayUp_Update_List::TMCMC) {
      output_forwardEstimators();
      output_credibleIntervals();
    }
  // create a new RBRV
    if (list->parent.is_subsetRel==false) {
      if (first_run) {
        if (data->rbrv_box.get_set(list->parent.get_name(),false)==NULL) {
          burbrvs = new flxBayUp_RBRV_set(true,list->parent);
          data->rbrv_box.register_set(burbrvs);
          scout << "rbrv_set: created new set '" << burbrvs->name << "' (" << burbrvs->get_NRV() << "/" << burbrvs->get_NOX() << ")." << std::endl;
        }
      }
    // some settings before we start with posterior sampling
      if (csm) csm->prepare();        // for posterior sampling
    }
  } catch (...) {
    FLXMSG("FlxBayUp_Update::update_g1",1);
    delete list; list = NULL;
    if (csm) {
      delete csm; 
      csm = NULL;
    }
    if (os_smpl) delete os_smpl;
    throw;
  }
}

const tdouble FlxBayUp_Update::line_search_LSF_call(const tdouble c, const flxVec& rv_base, flxVec& rv_prop, const flxVec& betaVec, tulong& N_LSF_calls, flx_LS_line_prop& lsp, tdouble& L) const
{
  RBRV_constructor& RndBox = list->get_RndBox();
  const tuint NRV = rv_base.get_N();
  rv_prop = rv_base;
  rv_prop.add(betaVec,c);
  RndBox.set_smp(rv_prop);
  ++N_LSF_calls;
  L = list->parent.eval_Likelihood();
  const tdouble g = list->eval_LSF(rv_prop[NRV-1], L);
  lsp.register_c(c,g);
  return g;
}

const flx_LS_line_prop FlxBayUp_Update::perform_line_search(const flxVec& rv_base, flxVec& rv_prop, const flxVec& betaVec, const tdouble tol, const tuint iter_max, tulong& N_LSF_calls) const
{
  flx_LS_line_prop res;
  const tdouble norm_beta = betaVec.get_Norm2();
  const tdouble tol1 = tol*(ONE+rv_base.get_Norm2())/(ONE+norm_beta);
  tdouble Ldummy=ZERO;
  // evaluate the initial value
    tdouble c_end = tdouble(1.2)*ONE;
    tdouble g_2 = line_search_LSF_call(c_end,rv_base,rv_prop,betaVec,N_LSF_calls,res,Ldummy);
  // evaluate the 2nd starting value
    tdouble c_start = tdouble(0.8)*ONE;
    tdouble g_1 = line_search_LSF_call(c_start,rv_base,rv_prop,betaVec,N_LSF_calls,res,Ldummy);
  // start the iteration
    tuint c = 0;
    while (c<iter_max) {
      // compute current tolerance
        tdouble dx;
        if (res.is_topo_set()) {
          dx = -rv_InvPhi_noAlert(res.get_upper_Pr(norm_beta));
          if (dx<ZERO) dx=ZERO;
          dx += ONE;
          if (dx>100*ONE) dx = 100*ONE;
        } else {
          dx = ONE;
        }
        dx *= tol1;
        const tdouble& dy = dx;
      tdouble g_3;
      if (g_1*g_2<=ZERO) {        // Pegasus-algorithm
        // evaluate a new point
          const tdouble c_res = (c_start*g_2-c_end*g_1)/(g_2-g_1);
          g_3 = line_search_LSF_call(c_res,rv_base,rv_prop,betaVec,N_LSF_calls,res,Ldummy);
        if ( g_2*g_3<ZERO ) {
          c_start = c_end;
          g_1 = g_2;
          c_end = c_res;
          g_2 = g_3;
        } else {
          const tdouble m = g_2/(g_2+g_3);
          g_1 = m*g_1;
          c_end = c_res;
          g_2 = g_3;
        }
      } else {                        // Sekantenverfahren
        tdouble c_res = c_end - ((c_end-c_start)/(g_2-g_1))*g_2;
        if (fabs((c_res-(c_end+c_start)/2)/(c_end-c_start)) > tdouble(10.)) {
          const tdouble cmtmp = (c_end+c_start)/2;
          const tdouble cdtmp = fabs(c_end-c_start)*10;
          c_res = cmtmp + ((c_res<cmtmp)?(-cdtmp):(cdtmp));
        }
        g_3 = line_search_LSF_call(c_res,rv_base,rv_prop,betaVec,N_LSF_calls,res,Ldummy);
        g_1 = g_2;
        c_start = c_end;
        g_2 = g_3;
        c_end = c_res;
      }
      // check for convergence
        if (fabs(c_end-c_start)<=dx) break;
        if (fabs(g_3)<=dy) break;
      ++c;
    }
    if (c>=iter_max) {
      return res;
    }
    // if only failures were observed up to now
      if (res.is_topo_set()==false || res.is_only_in()) {
        // propose a point that is very likely outside the failure region
          const tdouble t = -((g_2>g_1)?g_2:g_1);
          const tdouble m = (g_2-g_1)/(c_end-c_start);
          const tdouble c_res = c_start + (t-g_1)/m;
          if (res.is_topo_set()==false && c_res>tdouble(20.)) {
            res.force_topo(c_start,g_1,c_end,g_2);
            return res;
          }
        // evaluate the LSF for this point
          line_search_LSF_call(c_res,rv_base,rv_prop,betaVec,N_LSF_calls,res,Ldummy);
        // final check
//           if (res.is_only_in()) {
//             throw FlxException_NeglectInInteractive("FlxBayUp_Update::perform_line_search_1");
//           }
      }
    return res; 
}

flx_LS_line_prop::flx_LS_line_prop(const flx_LS_line_prop& rhs)
: topo_set(rhs.topo_set), l_out(rhs.l_out), l_in(rhs.l_in), u_out(rhs.u_out), u_in(rhs.u_in), ostack(NULL)
{
  if (rhs.ostack) {
    ostack = new std::stack<tdouble>(*(rhs.ostack));
  }
}

flx_LS_line_prop& flx_LS_line_prop::operator=(const flx_LS_line_prop& rhs)
{
  if (this!=&rhs) {
    topo_set = rhs.topo_set;
    l_out = rhs.l_out;
    l_in = rhs.l_in;
    u_out = rhs.u_out;
    u_in = rhs.u_in;
    if (ostack) {
      delete ostack;
      ostack = NULL;
    }
    if (rhs.ostack) {
      ostack = new std::stack<tdouble>(*(rhs.ostack));
    }
  }
  return *this;
}

flx_LS_line_prop::~flx_LS_line_prop()
{
  if (ostack) delete ostack;
}

void flx_LS_line_prop::set_topo()
{
  if (!topo_set) {
    topo_set = true;
    if (ostack) {
      while (!ostack->empty()) {
        register_out(ostack->top());
        ostack->pop();
      }
      delete ostack;
      ostack = NULL;
    }
  }
}

void flx_LS_line_prop::register_in(const tdouble c)
{
  if (l_in>c) l_in = c;
  if (u_in<c) u_in = c;
  set_topo();
}

void flx_LS_line_prop::register_out(const tdouble c)
{
  only_in = false;
  if (!topo_set) {
    if (ostack==NULL) {
      ostack = new std::stack<tdouble>;
    }
    ostack->push(c);
    return;
  }
  if (c<l_in) {
    if (l_out<c) l_out = c;
    #if FLX_DEBUG
      if (l_in<l_out) throw FlxException_Crude("flx_LS_line_prop::register_out_1");
    #endif
  } else if (c>u_in) {
    if (u_out>c) u_out = c;
    #if FLX_DEBUG
      if (u_in>u_out) throw FlxException_Crude("flx_LS_line_prop::register_out_2");
    #endif
  } else {
    throw FlxException_Crude("flx_LS_line_prop::register_out_3");
  }
}

void flx_LS_line_prop::register_c(const tdouble c, const tdouble g)
{
  if (g<=ZERO) register_in(c);
  else register_out(c);
}

const tdouble flx_LS_line_prop::get_upper_Pr(const tdouble betaNorm) const
{
  if (topo_set) {
    return rv_Phi(-l_out*betaNorm)-rv_Phi(-u_out*betaNorm);
  } else {
    return ZERO;
  }
}

const tdouble flx_LS_line_prop::get_mean_Pr(const tdouble betaNorm) const
{
  if (topo_set) {
    tdouble res = rv_Phi(-((l_out>-1e4)?(l_out+l_in)/2:l_out)*betaNorm);
    res -= rv_Phi(-((u_out<1e4)?(u_out+u_in)/2:u_out)*betaNorm);
    return res;
  } else {
    return ZERO;
  }
}

const tdouble flx_LS_line_prop::propose(const tdouble pr, const tdouble betaNorm) const
{
  if (!topo_set) {
    throw FlxException_Crude("flx_LS_line_prop::propose");
  } 
  const tdouble p1 = rv_Phi(-u_out*betaNorm);
  const tdouble p2 = rv_Phi(-l_out*betaNorm);
  const tdouble dp = pr*(p2-p1);
  const tdouble p3 = p1+dp;
  return -rv_InvPhi(p3)/betaNorm;
}

void flx_LS_line_prop::force_topo(const tdouble c1, const tdouble g1, const tdouble c2, const tdouble g2)
{
  if (topo_set) return;
  // find upper and lower values
    tdouble c_l, c_u, g_l, g_u;
    if (c1<c2) {
      c_l = c1;
      c_u = c2;
      g_l = g1;
      g_u = g2;
    } else if (c1>c2) {
      c_l = c2;
      c_u = c1;
      g_l = g2;
      g_u = g1;
    } else {
      throw FlxException_Crude("flx_LS_line_prop::force_topo_1");
    }
  if (g_l>g_u) {
    l_out = c_u;
  } else if (g_u>g_l) {
    u_out = c_l;
  } else {
    throw FlxException_Crude("flx_LS_line_prop::force_topo_1");
  }
  set_topo();
}

const tdouble flx_LS_line_prop::get_lower_scale() const
{
  if (topo_set==false) throw FlxException_Crude("flx_LS_line_prop::get_lower_scale_1");
  bool bo = (l_out>-1e3);
  bool bi = (l_in<1e3);
  if (bo && bi) return (l_out+l_in)/2;
  if (bo) return l_out;
  if (bi) return l_in;
  throw FlxException_Crude("flx_LS_line_prop::get_lower_scale_2");
}

void flx_LS_line_prop::print_topo(std::ostream& os) const
{
  if ( is_topo_set() ) {
    os << "    lower: ]" << l_out << "," << l_in << "]\t[" << u_in << "," << u_out << "["  << std::endl; 
  } else  {
    os << "    (topo not set)" << std::endl; 
  }
}



flxBayUp_RBRV_set::flxBayUp_RBRV_set(const bool internal, flxBayUp& parentV)
: RBRV_set_base(internal,RBRV_constructor::count_NRV(parentV.get_setvec()),parentV.get_name(),false),parent(parentV),setvec(parent.get_setvec()),
  NRV(sRV), NOX(RBRV_constructor::count_NOX(setvec)),Nsets(setvec.size()),proposed(false)
{

}

void flxBayUp_RBRV_set::set_is_valid(const bool is_valid)
{
  #if FLX_DEBUG
  for (std::vector<RBRV_set_base*>::const_iterator it = setvec.begin() ; it != setvec.end(); ++it) {
    (*it)->set_is_valid(is_valid);
  }
  #endif
}

void flxBayUp_RBRV_set::set_y(const tdouble*const y_vec)
{
  proposed = false;
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->set_y(y_vec+rvc);
    rvc += cs->get_NRV();
  }
  #if FLX_DEBUG
    if (rvc!=get_NRV()) throw FlxException_Crude("flxBayUp_RBRV_set::set_y");
  #endif
}

void flxBayUp_RBRV_set::set_y_only_this(const tdouble*const y_vec)
{
  proposed = false;
}

void flxBayUp_RBRV_set::get_y(tdouble*const y_vec)
{
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->get_y(y_vec+rvc);
    rvc += cs->get_NRV();
  }
}

void flxBayUp_RBRV_set::get_y_only_this(tdouble*const y_vec)
{

}

flxVec& flxBayUp_RBRV_set::get_y()
{
  get_y(y_of_set.get_tmp_vptr());
  return y_of_set;
}

void flxBayUp_RBRV_set::get_x(tdouble*const x_vec)
{
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->get_x(x_vec+rvc);
    rvc += cs->get_NOX();
  }
}

void flxBayUp_RBRV_set::get_x_only_this(tdouble*const x_vec)
{

}

void flxBayUp_RBRV_set::set_x(const tdouble*const x_vec)
{
  proposed = false;
  tuint rvc = 0;
  for (tuint i=0;i<Nsets;++i) {
    RBRV_set_base* cs = setvec[i];
    cs->set_x(x_vec+rvc);
    rvc += cs->get_NOX();
  }
}

void flxBayUp_RBRV_set::set_x_only_this(const tdouble*const x_vec)
{
  proposed = false;
}

const flxVec& flxBayUp_RBRV_set::propose_y()
{
  parent.updater.draw_realization(y_of_set);
  proposed = true;
  return y_of_set;
}

void flxBayUp_RBRV_set::transform_y2x()
{
  if (proposed) {
    #if FLX_DEBUG
      set_is_valid(true);
    #endif
    return;
  }
  // do the transformation
    for (tuint i=0;i<Nsets;++i) {
      RBRV_set_base* cs = setvec[i];
      cs->transform_y2x();
    }
  // check whether sample is in failure domain
    if (parent.updater.chk_accept_cur_smpl()==false) {
      std::ostringstream ssV;
      ssV << "The current sample must be rejected.";
      throw FlxException_NeglectInInteractive("flxBayUp_RBRV_set::transform_y2x", ssV.str() );
    }
}

const bool flxBayUp_RBRV_set::check_xVec(const tdouble* xp)
{
  throw FlxException_NotImplemented("flxBayUp_RBRV_set::check_xVec");
}

void flxBayUp_RBRV_set::get_mean(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("flxBayUp_RBRV_set::get_mean");
}

void flxBayUp_RBRV_set::get_mean_only_this(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("flxBayUp_RBRV_set::get_mean_only_this");
}

void flxBayUp_RBRV_set::get_sd(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("flxBayUp_RBRV_set::get_sd");
}

void flxBayUp_RBRV_set::get_sd_only_this(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("flxBayUp_RBRV_set::get_sd_only_this");
}

void flxBayUp_RBRV_set::find_dependent_sets(std::vector< RBRV_set_base* >& setvecV)
{
  // make sure that it is not already in the list
    for (tuint i=0;i<setvecV.size();++i) {
      if (setvecV[i]==this) return;
    }
  // add the parents of this set first
    const tuint n1 = setvec.size();
    for (tuint i=0;i<n1;++i) {
      setvec[i]->find_dependent_sets(setvecV);
    }
  setvecV.push_back(this);
}

const tuint flxBayUp_RBRV_set::group_dependent_sets(std::vector< RBRV_set_base* >& setvecV, const tuint pos_this)
{
  // find position in list
    tuint pos=pos_this;
    tuint removed = 0;
  // remove selected entries
    const tuint n1 = setvec.size();
    for (tuint i=0;i<n1;++i) {
      tuint j;
      for (j=0;j<pos;++j) {
        if (setvecV[j]==setvec[i]) break;
      }
      if (j>=pos) {
        std::ostringstream ssV;
        ssV << "There is a conflict with the set '" << setvec[i]->name << "' (current set: '" << name << "').";
        throw FlxException("flxBayUp_RBRV_set::group_dependent_sets_1", ssV.str() );
      } else {
        tuint ir = setvecV[j]->group_dependent_sets(setvecV,j);
        j -= ir;
        pos -= ir;
        removed += ir;
        setvecV.erase(setvecV.begin()+j);
        --pos;
        ++removed;
      }
    }
  return removed;
}

void flxBayUp_RBRV_set::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  for (tuint i=0;i<Nsets;++i) {
    setvec[i]->print(sout,prelim+"  ",counter,printID);
  }
}


RBRV_entry_fun_data::RBRV_entry_fun_data(const std::string& name, FlxFunction* fun, const tuint paraN, FlxIstream_vector* isv, const tdouble is_log)
: RBRV_entry_fun(name,fun), paraN(paraN), dataN(isv->get_total_size()/paraN), dataPtr(NULL), is_log(is_log)
{
  try {
    // check the input stream
      isv->reset_stream();
      const tulong N = isv->get_total_size();
      if (N==0) throw FlxException_Crude("RBRV_entry_fun_data::RBRV_entry_fun_data_1");
      if (N%paraN!=0) throw FlxException_Crude("RBRV_entry_fun_data::RBRV_entry_fun_data_2");
    // allocate and assign dataPtr
      dataPtr = new tdouble[N];
      const tdouble* cp = isv->get_tmpPtr();
      memcpy(dataPtr,cp,N*sizeof(tdouble));
  } catch (FlxException &e) {
    FLXMSG("RBRV_entry_fun_data::RBRV_entry_fun_data_3",1);
    delete  fun;
    throw;
  }
}

RBRV_entry_fun_data::~RBRV_entry_fun_data()
{
  if (dataPtr) delete [] dataPtr;
}

void RBRV_entry_fun_data::transform_y2x(const tdouble*const y_vec)
{
  #if FLX_DEBUG
    if (valid) throw FlxException_Crude("RBRV_entry_fun_data::transform_y2x");
  #endif
  const tdouble* const tT = FunPara::ParaList;
  const tuint tTS = FunPara::ParaListSize;
  FunPara::ParaListSize = paraN;
  tdouble res = 0.;
  if (is_log) {
    for (tuint i=0;i<dataN;++i) {
      FunPara::ParaList = dataPtr+i*paraN;
      res += fun->calc();
    }
  } else {
    for (tuint i=0;i<dataN;++i) {
      FunPara::ParaList = dataPtr+i*paraN;
      res += log(fun->cast2positive_or0());
    }
  }
  value = res;
  fun_val = res;
  FunPara::ParaList = tT;
  FunPara::ParaListSize = tTS;
  #if FLX_DEBUG
    valid = true;
  #endif
}

void RBRV_entry_ref_log::transform_y2x(const tdouble*const y_vec)
{
  #if FLX_DEBUG
    if (valid) throw FlxException_Crude("RBRV_entry_ref_log::transform_y2x");
    valid = true;
  #endif
  value = exp(reflog);
}

const tdouble RBRV_entry_ref_log::get_mean_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_ref_log::get_mean_current_config");
}

const tdouble RBRV_entry_ref_log::get_sd_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_ref_log::get_sd_current_config");
}

const bool RBRV_entry_ref_log::check_x(const tdouble xV)
{
  throw FlxException_NotImplemented("RBRV_entry_ref_log::check_x");
}


flxBayUp_uncertobsv_set::flxBayUp_uncertobsv_set(const std::string& name, RBRV_set* singleObsvV, FlxFunction* flxLike, const tuint NobservV, const tuint paraN, FlxIstream_vector* isv, const bool is_log)
: RBRV_set_base(true,singleObsvV->get_NRV()*NobservV,name,false),
  singleObsv(singleObsvV), flxLike(flxLike), Nobserv(NobservV), paraN(paraN), dataPtr(NULL), Likelihood(ZERO), is_log(is_log)
{
  // allocate and assign dataPtr
    const tulong N = Nobserv*paraN;
    #if FLX_DEBUG
      if (N!=isv->get_total_size()) throw FlxException_Crude("flxBayUp_uncertobsv_set::flxBayUp_uncertobsv_set");
    #endif
    dataPtr = new tdouble[N];
    const tdouble* cp = isv->get_tmpPtr();
    memcpy(dataPtr,cp,N*sizeof(tdouble));
}

flxBayUp_uncertobsv_set::~flxBayUp_uncertobsv_set()
{
  delete flxLike;
  if (dataPtr) delete [] dataPtr;
}

void flxBayUp_uncertobsv_set::set_is_valid(const bool is_valid)
{
  #if FLX_DEBUG
  singleObsv->set_is_valid(is_valid);
  #endif
}

void flxBayUp_uncertobsv_set::transform_y2x()
{
  // prepare the parameter list
    const tdouble* const tT = FunPara::ParaList;
    const tuint tTS = FunPara::ParaListSize;
    FunPara::ParaListSize = paraN;
    const tdouble* yt = y_of_set.get_tmp_vptr_const();
  Likelihood = 0.;
  tuint idx_count = 0;
  tuint nrv_count = 0;
  for (tuint i=0;i<Nobserv;++i) {
    FunPara::ParaList = dataPtr+idx_count;
    singleObsv->set_y(yt+nrv_count);
    singleObsv->transform_y2x();
    if (is_log) {
      Likelihood += flxLike->calc();
    } else {
      Likelihood += log(flxLike->cast2positive_or0());
    }
    idx_count += paraN;
    nrv_count += singleObsv->get_NRV();
    #if FLX_DEBUG
      singleObsv->set_is_valid(false);
    #endif
  }
  // reset the parameter list
    FunPara::ParaList = tT;
    FunPara::ParaListSize = tTS;
}

const bool flxBayUp_uncertobsv_set::check_xVec(const tdouble* xp)
{
  throw FlxException_NotImplemented("flxBayUp_uncertobsv_set::check_xVec");
}

void flxBayUp_uncertobsv_set::get_mean(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("flxBayUp_uncertobsv_set::get_mean");
}

void flxBayUp_uncertobsv_set::get_mean_only_this(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("flxBayUp_uncertobsv_set::get_mean_only_this");
}

void flxBayUp_uncertobsv_set::get_sd(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("flxBayUp_uncertobsv_set::get_sd");
}

void flxBayUp_uncertobsv_set::get_sd_only_this(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("flxBayUp_uncertobsv_set::get_sd_only_this");
}

void flxBayUp_uncertobsv_set::find_dependent_sets(std::vector< RBRV_set_base* >& setvec)
{
  // make sure that it is not already in the list
    for (tuint i=0;i<setvec.size();++i) {
      if (setvec[i]==this) return;
    }
  setvec.push_back(this);
}

void flxBayUp_uncertobsv_set::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << " Number of observations:               " << Nobserv << std::endl;
  sout << prelim << "  " << " Number of parameters per observation: " << paraN << std::endl;
  sout << prelim << "  " << " RBRV-set of an observation: " << std::endl;
  tuint c2 = 0;
  singleObsv->print(sout,prelim+"    ",c2,false);
  counter += get_NOX_only_this();
}


void RBRV_entry_value::transform_y2x(const tdouble*const y_vec)
{
  throw FlxException_NotImplemented("RBRV_entry_value::transform_y2x");
}

const tdouble RBRV_entry_value::get_mean_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_value::get_mean_current_config");
}

const tdouble RBRV_entry_value::get_sd_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_value::get_sd_current_config");
}

const bool RBRV_entry_value::check_x(const tdouble xV)
{
  throw FlxException_NotImplemented("RBRV_entry_value::check_x");
}

const bool RBRV_entry_value::search_circref(FlxFunction* fcr)
{
  throw FlxException_NotImplemented("RBRV_entry_value::search_circref");
}



flxBayUp_mProb_set::flxBayUp_mProb_set(const bool internal, const std::string& name, const tuint Nmodels, flxBayUp_RBRV_set** modelVec, flxVec priorPrVec, const tuint N_model_res, std::vector< std::string >& model_res_list_Str, FlxFunction** model_res_map)
: RBRV_set_base(internal, 0, name,false), Nmodels(Nmodels), modelVec(modelVec), postPrVec(priorPrVec),
  sumPrVec(ZERO), p(RBRV_entry_RV_uniform("p",0,new FlxFunction(new FunNumber(ZERO)),new FlxFunction(new FunNumber(ONE)),true)),
  y_total(NULL),
  N_model_res(N_model_res), model_res_list(new RBRV_entry_value*[N_model_res]), model_res_map(model_res_map)
{
  for (tuint i=0;i<N_model_res;++i) model_res_list[i] = NULL;
  try {
    if (Nmodels<=1) {
      throw FlxException("flxBayUp_mProb_set::flxBayUp_mProb_set_1","You have to take at least 2 models into account.");
    }
    for (tuint i=0;i<Nmodels;++i) {
      postPrVec[i] *= modelVec[i]->get_parent().updater.get_PrMod();
      sumPrVec += postPrVec[i];
    }
    for (tuint i=0;i<Nmodels;++i) {
      modelVec[i]->find_dependent_sets(setvec);
    }
    NRV_total = RBRV_constructor::count_NRV_long(setvec);
    NOX_total = RBRV_constructor::count_NOX_long(setvec);
    y_total = new flxVec(get_NRV());
    for (tuint i=0;i<N_model_res;++i) {
      const std::string s = name + "::" + model_res_list_Str[i];
      model_res_list[i] = new RBRV_entry_value(s);
      data->rbrv_box.register_entry(model_res_list[i]);
    }
  } catch (FlxException &e) {
    FLXMSG("flxBayUp_mProb_set::flxBayUp_mProb_set_2",1);
    free_mem();
    throw;
  }
}

flxBayUp_mProb_set::~flxBayUp_mProb_set()
{
  free_mem();
}


void flxBayUp_mProb_set::free_mem()
{
  delete [] modelVec;
  if (y_total) delete y_total;
  for (tuint i=0;i<N_model_res;++i) {
    if (model_res_list[i]) delete model_res_list[i];
  }
  delete [] model_res_list;
  for (tuint i=0;i<N_model_res*Nmodels;++i) {
    delete model_res_map[i];
  }
  delete [] model_res_map;
}

const tuint flxBayUp_mProb_set::get_model_ID() const
{
  const tdouble t = p.get_value()*sumPrVec;
  tdouble s = 0.;
  for (tuint i=0;i<Nmodels;++i) {                // the efficiency of this can be improved!
    s += postPrVec[i];
    if (t<=s) return i;
  }
  return Nmodels-1;
}

void flxBayUp_mProb_set::update_model_res(const tuint id)
{
  tuint c = id;
  for (tuint i=0;i<N_model_res;++i) {
    model_res_list[i]->eval_para();
    model_res_list[i]->set_x( model_res_map[c]->calc() );
    c += Nmodels;
  }
}

void flxBayUp_mProb_set::set_is_valid(const bool is_valid)
{
  #if FLX_DEBUG
  for (tuint i=0;i<Nmodels;++i) {
    modelVec[i]->set_is_valid(is_valid);
  }
  p.set_is_valid(is_valid);
  for (tuint i=0;i<N_model_res;++i) {
    model_res_list[i]->set_is_valid(is_valid);
  }
  #endif
}

const flxVec& flxBayUp_mProb_set::propose_y()
{
  y_total->operator[](NRV_total) = RndCreator->gen_smp();
  p.RBRV_entry_RV_base::transform_y2x(y_total->get_tmp_vptr()+NRV_total); 
  const tuint id = get_model_ID();
  modelVec[id]->propose_y();
  return get_y();
}

void flxBayUp_mProb_set::transform_y2x()
{
  p.eval_para();
  p.RBRV_entry_RV_base::transform_y2x(y_total->get_tmp_vptr()+NRV_total); 
  const tuint id = get_model_ID();
  modelVec[id]->transform_y2x();
  update_model_res(id);
  #if FLX_DEBUG
    set_is_valid(true);
  #endif
}

void flxBayUp_mProb_set::set_y(const tdouble*const y_vec)
{
  tuint rvc = 0;
  for (tuint i=0;i<setvec.size();++i) {
    RBRV_set_base* cs = setvec[i];
    cs->set_y_only_this(y_vec+rvc);
    rvc += cs->get_NRV_only_this();
  }
  #if FLX_DEBUG
    if (rvc!=NRV_total) throw FlxException_Crude("flxBayUp_mProb_set::set_y");
  #endif
  y_total->operator[](NRV_total) = y_vec[NRV_total];
}

void flxBayUp_mProb_set::get_y(tdouble*const y_vec)
{
  tuint rvc = 0;
  for (tuint i=0;i<setvec.size();++i) {
    RBRV_set_base* cs = setvec[i];
    cs->get_y_only_this(y_vec+rvc);
    rvc += cs->get_NRV_only_this();
  }
  #if FLX_DEBUG
    if (rvc!=NRV_total) throw FlxException_Crude("flxBayUp_mProb_set::get_y");
  #endif
  y_vec[NRV_total] = y_total->operator[](NRV_total);
}

flxVec& flxBayUp_mProb_set::get_y()
{
  get_y(y_total->get_tmp_vptr());
  return *y_total;
}

void flxBayUp_mProb_set::set_x(const tdouble*const x_vec)
{
  tuint rvc = 0;
  for (tuint i=0;i<setvec.size();++i) {
    RBRV_set_base* cs = setvec[i];
    cs->set_x_only_this(x_vec+rvc);
    rvc += cs->get_NOX_only_this();
  }
  #if FLX_DEBUG
    if (rvc!=NOX_total) throw FlxException_Crude("flxBayUp_mProb_set::set_x");
  #endif
  p.eval_para();
  p.set_x(x_vec[rvc]);
  ++rvc;
  for (tuint i=0;i<N_model_res;++i) {
    model_res_list[i]->eval_para();
    model_res_list[i]->set_x(x_vec[rvc]);
    ++rvc;
  }
}

void flxBayUp_mProb_set::get_x(tdouble*const x_vec)
{
  tuint rvc = 0;
  for (tuint i=0;i<setvec.size();++i) {
    RBRV_set_base* cs = setvec[i];
    cs->get_x(x_vec+rvc);
    rvc += cs->get_NOX();
  }
  #if FLX_DEBUG
    if (rvc!=NOX_total) throw FlxException_Crude("flxBayUp_mProb_set::get_x");
  #endif
  x_vec[rvc] = p.get_value();
  ++rvc;
  for (tuint i=0;i<N_model_res;++i) {
    x_vec[rvc] = model_res_list[i]->get_value();
    ++rvc;
  }  
}

const bool flxBayUp_mProb_set::check_xVec(const tdouble* xp)
{
  throw FlxException_NotImplemented("flxBayUp_mProb_set::check_xVec");
}

void flxBayUp_mProb_set::get_x_only_this(tdouble*const x_vec)
{
  throw FlxException_NotImplemented("flxBayUp_mProb_set::get_x_only_this");
}

void flxBayUp_mProb_set::get_y_only_this(tdouble*const y_vec)
{
  throw FlxException_NotImplemented("flxBayUp_mProb_set::get_y_only_this");
}

void flxBayUp_mProb_set::set_x_only_this(const tdouble*const x_vec)
{
  throw FlxException_NotImplemented("flxBayUp_mProb_set::set_x_only_this");
}

void flxBayUp_mProb_set::set_y_only_this(const tdouble*const y_vec)
{
  throw FlxException_NotImplemented("flxBayUp_mProb_set::set_y_only_this");
}

void flxBayUp_mProb_set::get_mean(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("flxBayUp_mProb_set::get_mean");
}

void flxBayUp_mProb_set::get_mean_only_this(tdouble*const m_vec)
{
  throw FlxException_NotImplemented("flxBayUp_mProb_set::get_mean_only_this");
}

void flxBayUp_mProb_set::get_sd(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("flxBayUp_mProb_set::get_sd");
}

void flxBayUp_mProb_set::get_sd_only_this(tdouble*const s_vec)
{
  throw FlxException_NotImplemented("flxBayUp_mProb_set::get_sd_only_this");
}

void flxBayUp_mProb_set::print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID)
{
  sout << prelim << "- " << name << " (" << get_NRV_only_this() << "/" <<  get_NOX_only_this() << ")" << std::endl;
  sout << prelim << "  " << " Number of models: " << Nmodels << std::endl;
  for (tuint i=0;i<Nmodels;++i) {
    sout << prelim << "  " << " - " << modelVec[i]->name << "  (" << GlobalVar.Double2String((postPrVec[i]/sumPrVec)*100) << "%)" << std::endl;
  }
  sout << prelim << "  " << " All sets: " << std::endl;
  tuint c2 = 0;
  for (tuint i=0;i<setvec.size();++i) {
    setvec[i]->print(sout,prelim+"    ",c2,false);
  }
  counter += get_NOX();
}

void flxBayUp_mProb_set::find_dependent_sets(std::vector< RBRV_set_base* >& setvecV)
{
  // make sure that it is not already in the list
    for (tuint i=0;i<setvecV.size();++i) {
      if (setvecV[i]==this) return;
    }
  // add the parents of this set first
    const tuint n1 = setvec.size();
    for (tuint i=0;i<n1;++i) {
      setvec[i]->find_dependent_sets(setvecV);                // this could be done faster ... 
    }
  setvecV.push_back(this);
}

const tuint flxBayUp_mProb_set::group_dependent_sets(std::vector< RBRV_set_base* >& setvecV, const tuint pos_this)
{
  // find position in list
    tuint pos=pos_this;
    tuint removed = 0;
  // remove selected entries
    const tuint n1 = setvec.size();
    for (tuint i=0;i<n1;++i) {
      tuint j;
      for (j=0;j<pos;++j) {
        if (setvecV[j]==setvec[i]) break;
      }
      if (j>=pos) {
        std::ostringstream ssV;
        ssV << "There is a conflict with the set '" << setvec[i]->name << "'.";
        throw FlxException("flxBayUp_mProb_set::group_dependent_sets_1", ssV.str() );
      } else {
        setvecV.erase(setvecV.begin()+j);
        --pos;
        ++removed;
      }
    }
  return removed;
}












