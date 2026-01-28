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

#include "flxrbrv_rvs.h"
#include "flxobjrandom.h"


const tdouble RBRV_entry_RV_stdN::transform_y2x(const tdouble y_val)
{
  return y_val;
}

const tdouble RBRV_entry_RV_stdN::transform_x2y(const tdouble& x_val)
{
  return x_val;
}

const tdouble RBRV_entry_RV_stdN::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  return rv_phi(x_val);
}

const tdouble RBRV_entry_RV_stdN::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  return rv_phi_log(x_val);
}

const tdouble RBRV_entry_RV_stdN::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return rv_Phi(x_val);
}

const tdouble RBRV_entry_RV_stdN::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  return rv_Phi(-x_val);
}

const tdouble RBRV_entry_RV_stdN::calc_entropy()
{
  return log(2*PI*exp(ONE))/2;
}

py::dict RBRV_entry_RV_stdN::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["mean"] = get_mean_current_config();
  res["sd"] = get_sd_current_config();
  res["entropy"] = calc_entropy();
  return res;
}

const tdouble RBRV_entry_RV_stdN::get_HPD(const tdouble p)
{
  return (ONE-p)/2;
}


RBRV_entry_RV_normal::RBRV_entry_RV_normal(const std::string& name, const tuint iID, const int pid, FlxFunction* p1v, FlxFunction* p2v, FlxFunction* p3v, FlxFunction* p4v, const bool eval_once)
: RBRV_entry_RV_base(name,iID), pid(pid), p1(p1v), p2(p2v), p3(p3v), p4(p4v), eval_once(eval_once), mu(ZERO), sigma(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_normal::RBRV_entry_RV_normal_100",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (p3) delete p3;
    if (p4) delete p4;
    throw;
  }
}

RBRV_entry_RV_normal::RBRV_entry_RV_normal(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), pid(0), p1(nullptr), p2(nullptr), p3(nullptr), p4(nullptr), eval_once(false), mu(ZERO), sigma(ZERO)
{
  try {
    if (config.contains("mu")) {          // mean, standard deviation
      pid = 0;
      p1 = parse_py_para("mu", config);
      p2 = parse_py_para("sd", config);
    } else if (config.contains("cov")) {   // C.o.V. and quantile value
      pid = 2;
      p1 = parse_py_para("cov", config);
      p2 = parse_py_para("val_1", config);
      p3 = parse_py_para("pr_1", config);
    } else if (config.contains("sd")) {   // std.dev. and quantile value
      pid = 3;
      p1 = parse_py_para("sd", config);
      p2 = parse_py_para("val_1", config);
      p3 = parse_py_para("pr_1", config);
    } else if (config.contains("pr_1")) {   // quantile values
      pid = 1;
      p2 = parse_py_para("pr_1", config);
      p1 = parse_py_para("val_1", config);
      p3 = parse_py_para("val_2", config);
      p4 = parse_py_para("pr_2", config);
    } else {
      throw FlxException_NeglectInInteractive("RBRV_entry_RV_normal::RBRV_entry_RV_normal_70", "Required parameters to define distribution not found in Python <dict>.");
    }

    eval_once = parse_py_para_as_bool("eval_once", config, false, false);

    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_normal::RBRV_entry_RV_normal_99",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (p3) delete p3;
    if (p4) delete p4;
    throw;
  }
}

RBRV_entry_RV_normal::~RBRV_entry_RV_normal()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (p3) delete p3;
  if (p4) delete p4;
}

void RBRV_entry_RV_normal::eval_para()
{
  if (!eval_once || (eval_once&&p1) ) {
    switch (pid) {
      case 0:                        // defined in terms of mean and standard deviation
        mu = p1->calc();
        sigma = p2->cast2positive();
        break;
      case 1:                        // defined in terms of quantiles
        {
          const tdouble x1 = p1->calc();
          const tdouble pr1 = p2->cast2positive();
          const tdouble x2 = p3->calc();
          const tdouble pr2 = p4->cast2positive();
          get_para_from_quantile(mu,sigma,x1,pr1,x2,pr2);
        }
        break;
      case 2:                        // defined in terms of C.o.V and quantile
        {
          const tdouble delta = p1->cast2positive();
          const tdouble x1 = p2->calc();
          const tdouble pr1 = p3->cast2positive();
          get_para_from_quantile2(mu,sigma,x1,pr1,delta);
        }
        break;
      case 3:                        // defined in terms of C.o.V and quantile
        {
          sigma = p1->cast2positive();
          const tdouble x1 = p2->calc();
          const tdouble pr1 = p3->cast2positive();
          get_para_from_quantile3(mu,x1,pr1,sigma);
        }
        break;
      default:
        throw FlxException_Crude("RBRV_entry_RV_normal::get_paras_1");
    }
    if (eval_once) {
      delete p1; p1=NULL;
      delete p2; p2=NULL;
      if (p3) { delete p3; p3=NULL; }
      if (p4) { delete p4; p4=NULL; }
    }
  }
}

void RBRV_entry_RV_normal::get_para_from_quantile(tdouble& meanV, tdouble& sdV, const tdouble x1, const tdouble p1, const tdouble x2, const tdouble p2)
{
  // check validity of parameter values
    if (p1>=ONE) {
      std::ostringstream ssV;
      ssV << "Expected a probability, which must not be larger than one (" << GlobalVar.Double2String(p1) << ").";
      throw FlxException("RBRV_entry_RV_normal::get_para_from_quantile_01", ssV.str() );
    }
    if (p2>=ONE) {
      std::ostringstream ssV;
      ssV << "Expected a probability, which must not be larger than one (" << GlobalVar.Double2String(p2) << ").";
      throw FlxException("RBRV_entry_RV_normal::get_para_from_quantile_02", ssV.str() );
    }
  // get standardized values
    const tdouble z1 = rv_InvPhi(p1);
    const tdouble z2 = rv_InvPhi(p2);
  // solve linear system of equations
    meanV = (z2*x1-z1*x2)/(z2-z1);    
    sdV = (x2-x1)/(z2-z1);
  // check for validity of standard deviation
    if (sdV <= ZERO) {
      std::ostringstream ssV;
      ssV << "Standard deviation must not become negative or zero (" << GlobalVar.Double2String(sdV) << ").";
      throw FlxException("RBRV_entry_RV_normal::get_para_from_quantile_03", ssV.str() );            
    }
}

void RBRV_entry_RV_normal::get_para_from_quantile2(tdouble& meanV, tdouble& sdV, const tdouble x1, const tdouble p1, const tdouble delta)
{
  // check validity of parameter values
    if (p1>=ONE) {
      std::ostringstream ssV;
      ssV << "Expected a probability, which must not be larger than one (" << GlobalVar.Double2String(p1) << ").";
      throw FlxException("RBRV_entry_RV_normal::get_para_from_quantile2_01", ssV.str() );
    }
    if (delta<=ZERO) {
      std::ostringstream ssV;
      ssV << "Expected a coefficient of variation, which must not be smaller than zero (" << GlobalVar.Double2String(delta) << ").";
      throw FlxException("RBRV_entry_RV_normal::get_para_from_quantile2_02", ssV.str() );
    }
  // get standardized values
    const tdouble z1 = rv_InvPhi(p1);
  // compute standard deviation   
    sdV = x1/(ONE/delta+z1);
    meanV = sdV/delta;
  // check for validity of standard deviation
    if (sdV <= ZERO) {
      std::ostringstream ssV;
      ssV << "Standard deviation must not become negative or zero (" << GlobalVar.Double2String(sdV) << ").";
      throw FlxException("RBRV_entry_RV_normal::get_para_from_quantile2_03", ssV.str() );            
    }
}

void RBRV_entry_RV_normal::get_para_from_quantile3(tdouble& meanV, const tdouble x1, const tdouble p1, const tdouble sdV)
{
  // check validity of parameter values
    if (p1>=ONE) {
      std::ostringstream ssV;
      ssV << "Expected a probability, which must not be larger than one (" << GlobalVar.Double2String(p1) << ").";
      throw FlxException("RBRV_entry_RV_normal::get_para_from_quantile3_01", ssV.str() );
    }
    if (sdV<=ZERO) {
      std::ostringstream ssV;
      ssV << "Expected a standard deviation, which must not be smaller than zero (" << GlobalVar.Double2String(sdV) << ").";
      throw FlxException("RBRV_entry_RV_normal::get_para_from_quantile3_02", ssV.str() );
    }
  // get standardized values
    const tdouble z1 = rv_InvPhi(p1);
  // compute standard deviation  
    meanV = x1 - sdV*z1;
}

const tdouble RBRV_entry_RV_normal::transform_y2x(const tdouble y_val)
{
  return mu+sigma*y_val;
}

const tdouble RBRV_entry_RV_normal::transform_x2y(const tdouble& x_val)
{
  return (x_val-mu)/sigma;
}

const tdouble RBRV_entry_RV_normal::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  return rv_phi((x_val-mu)/sigma)/sigma;
}

const tdouble RBRV_entry_RV_normal::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  return rv_phi_log((x_val-mu)/sigma)-log(sigma);
}

const tdouble RBRV_entry_RV_normal::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return rv_Phi((x_val-mu)/sigma);
}

const tdouble RBRV_entry_RV_normal::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  return rv_Phi((mu-x_val)/sigma);
}

const tdouble RBRV_entry_RV_normal::calc_entropy()
{
  const tdouble s = sigma;
  return log(2*PI*exp(ONE)*pow2(s))/2;
}

const bool RBRV_entry_RV_normal::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || p1->search_circref(fcr) || p2->search_circref(fcr) || (p3?(p3->search_circref(fcr)):false) || (p4?(p4->search_circref(fcr)):false);
}

py::dict RBRV_entry_RV_normal::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["mean"] = get_mean_current_config();
  res["sd"] = get_sd_current_config();
  res["entropy"] = calc_entropy();
  return res;
}

const tdouble RBRV_entry_RV_normal::get_HPD(const tdouble p)
{
  return (ONE-p)/2;
}

RBRV_entry_RV_lognormal::RBRV_entry_RV_lognormal(const std::string& name, const tuint iID, const int pid, FlxFunction* p1, FlxFunction* p2, FlxFunction* p3, FlxFunction* p4, FlxFunction* epsilon, const bool eval_once)
: RBRV_entry_RV_base(name,iID), pid(pid), p1(p1), p2(p2), p3(p3), p4(p4), epsilon(epsilon), eval_once(eval_once), lambda(ZERO), zeta(ZERO), eps(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_lognormal::RBRV_entry_RV_lognormal_100",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (p3) delete p3;
    if (p4) delete p4;
    if (epsilon) delete epsilon;
    throw;
  }
}

RBRV_entry_RV_lognormal::RBRV_entry_RV_lognormal(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), pid(0), p1(nullptr), p2(nullptr), p3(nullptr), p4(nullptr), epsilon(nullptr), eval_once(false), lambda(ZERO), zeta(ZERO), eps(ZERO)
{
  try {
    if (config.contains("lambda")) {          // 0:lambda,zeta;
      pid = 0;
      p1 = parse_py_para("lambda", config);
      p2 = parse_py_para("zeta", config);
    } else if (config.contains("mu")) {   // 1:mean,sd;
      pid = 1;
      p1 = parse_py_para("mu", config);
      p2 = parse_py_para("sd", config);
    } else if (config.contains("mode") && config.contains("sd")) {   // 2:mode,sd;
      pid = 2;
      p1 = parse_py_para("mode", config);
      p2 = parse_py_para("sd", config);
    } else if (config.contains("median") && config.contains("sd")) {   // 3:median,sd;
      pid = 3;
      p1 = parse_py_para("median", config);
      p2 = parse_py_para("sd", config);
    } else if (config.contains("pr_2")) {   // 4: quantile values;
      pid = 4;
      p1 = parse_py_para("val_1", config);
      p2 = parse_py_para("pr_1", config);
      p3 = parse_py_para("val_2", config);
      p4 = parse_py_para("pr_2", config);
    } else if (config.contains("median") && config.contains("cov")) {   // 5: median,C.o.V.;
      pid = 5;
      p1 = parse_py_para("median", config);
      p2 = parse_py_para("cov", config);
    } else if (config.contains("cov") && config.contains("pr_1")) {   // 6: C.o.V,quantile
      pid = 6;
      p1 = parse_py_para("cov", config);
      p2 = parse_py_para("val_1", config);
      p3 = parse_py_para("pr_1", config);
    } else if (config.contains("mode") && config.contains("cov")) {   // 6: mode, C.o.V.
      pid = 7;
      p1 = parse_py_para("mode", config);
      p2 = parse_py_para("cov", config);
    } else {
      throw FlxException_NeglectInInteractive("RBRV_entry_RV_lognormal::RBRV_entry_RV_lognormal_70", "Required parameters to define distribution not found in Python <dict>.");
    }

    epsilon = parse_py_para("epsilon", config, false);

    eval_once = parse_py_para_as_bool("eval_once", config, false, false);

    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_lognormal::RBRV_entry_RV_lognormal_99",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (p3) delete p3;
    if (p4) delete p4;
    if (epsilon) delete epsilon;
    throw;
  }
}

RBRV_entry_RV_lognormal::~RBRV_entry_RV_lognormal()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (p3) delete p3;
  if (p4) delete p4;
  if (epsilon) delete epsilon;
}

void RBRV_entry_RV_lognormal::eval_para()
{
  if (!eval_once || (eval_once&&p1) ) {
    const tdouble p1d = ((pid==6)?p1->cast2positive():p1->calc());
    const tdouble p2d = ((pid==1 || pid==4)?p2->cast2positive():p2->calc());
    eps = epsilon?(epsilon->calc()):ZERO;
    switch (pid) {
      case 0:        // lambda, zeta
        lambda = p1d;
        zeta = p2d;
        break;
      case 1:        // mean, sd
        {
          const tdouble& mu = p1d;
          const tdouble& sigma = p2d;
          if (mu<=eps) {
            std::ostringstream ssV;
            ssV << "The mean (" << GlobalVar.Double2String(mu) << ") must not be smaller than epsilon (" << GlobalVar.Double2String(eps) << ").";
            throw FlxException("RBRV_entry_RV_lognormal::get_paras_0", ssV.str() );
          }
          const tdouble t = std::log(pow2(sigma)/pow2(mu-eps)+ONE);
          lambda = std::log(mu-eps) - t/2;
          zeta = std::sqrt(t);
        }
        break;
      case 2:        // mode ,sd
        {
          const tdouble& mode = p1d;
          const tdouble& sigma = p2d;
          if (mode<=eps) {
            std::ostringstream ssV;
            ssV << "The mode (" << GlobalVar.Double2String(mode) << ") must not be smaller than epsilon (" << GlobalVar.Double2String(eps) << ").";
            throw FlxException("RBRV_entry_RV_lognormal::get_paras_1", ssV.str() );
          }
          // get zeta
            const tdouble c = pow2(sigma/(mode-eps));
            tdouble y = ONE;
            tdouble fy;
            tuint N = 0;
            do {
              const tdouble y2 = pow2(y);
              fy = (pow2(y2)-y*y2-c);
              y -= fy/(4*y*y2-3*y2);
              ++N;
            } while (fabs(fy)/c>=1e-8 && N<100 );
            const tdouble zeta2 = log(y);
            zeta = sqrt(zeta2);
          // get lambda
            lambda = log(mode-eps) + zeta2;
        }
        break;
      case 3:        // median, sd
        {
          const tdouble& median = p1d;
          if (median<=eps) {
            std::ostringstream ssV;
            ssV << "The median (" << GlobalVar.Double2String(median) << ") must not be smaller than epsilon (" << GlobalVar.Double2String(eps) << ").";
            throw FlxException("RBRV_entry_RV_lognormal::get_paras_2", ssV.str() );
          }
          lambda = std::log(median-eps);
          // get zeta
            const tdouble& sigma = p2d;
            const tdouble c = pow2(sigma/(median-eps));
            const tdouble y1 = (ONE+sqrt(ONE+4*c))/2;
            zeta = std::sqrt(std::log(y1));
        }
        break;
      case 4:        // fraktiles
        {
          const tdouble x1 = log(p1d-eps);
          const tdouble p1 = p2d;
          const tdouble x2 = log(p3->calc()-eps);
          const tdouble p2 = p4->cast2positive();
          RBRV_entry_RV_normal::get_para_from_quantile(lambda,zeta,x1,p1,x2,p2);
          break;
        }
      case 5:        // median, C.o.V.
        {
          const tdouble& median = p1d;
          if (median<=eps) {
            std::ostringstream ssV;
            ssV << "The median (" << GlobalVar.Double2String(median) << ") must not be smaller than epsilon (" << GlobalVar.Double2String(eps) << ").";
            throw FlxException("RBRV_entry_RV_lognormal::get_paras_3", ssV.str() );
          }
          lambda = std::log(median-eps);
          // get zeta
            const tdouble& CoV = p2d;
            zeta = std::sqrt(log(pow2(CoV)+ONE));
        }
        break;
      case 6:        // C.o.V., quantile
        {
          const tdouble CoV = p1d;
          const tdouble x1 = log(p2d-eps);
          const tdouble p1 = p3->cast2positive();
          zeta = std::sqrt(log(pow2(CoV)+ONE));
          RBRV_entry_RV_normal::get_para_from_quantile3(lambda,x1,p1,zeta);
        }
        break;
      case 7:        // mode, C.o.V.
        {
          const tdouble& mode = p1d;
          if (mode<=eps) {
            std::ostringstream ssV;
            ssV << "The mode (" << GlobalVar.Double2String(mode) << ") must not be smaller than epsilon (" << GlobalVar.Double2String(eps) << ").";
            throw FlxException("RBRV_entry_RV_lognormal::get_paras_4", ssV.str() );
          }
          // get zeta
            const tdouble& CoV = p2d;
            const tdouble t = log(pow2(CoV)+ONE);
            zeta = std::sqrt(t);
          lambda = std::log(mode-eps)+t;
        }
        break;
      default:
        throw FlxException_Crude("RBRV_entry_RV_lognormal::get_paras_5");
    }
    if (eval_once) {
      delete p1; p1=NULL;
      delete p2; p2=NULL;
      if (p3) { delete p3; p3=NULL; }
      if (p4) { delete p4; p4=NULL; }
      delete epsilon; epsilon = NULL;
    }
  }
}

const tdouble RBRV_entry_RV_lognormal::transform_y2x(const tdouble y_val)
{
  return exp(y_val*zeta+lambda)+eps;
}

const tdouble RBRV_entry_RV_lognormal::transform_x2y(const tdouble& x_val)
{
  if (x_val<=eps) {
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller or equal than '" << GlobalVar.Double2String(eps) << "' is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_lognormal::transform_x2y", ssV.str() );
  }
  return (log(x_val-eps)-lambda)/zeta;
}

const tdouble RBRV_entry_RV_lognormal::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller or equal than '" << GlobalVar.Double2String(eps) << "' is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_lognormal::calc_pdf_x", ssV.str() );
  }
  return exp(-(ONE/2)*pow2((log(x_val-eps)-lambda)/zeta))/((x_val-eps)*zeta*sqrt(2*PI));
}

const tdouble RBRV_entry_RV_lognormal::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return log(ZERO);
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller or equal than '" << GlobalVar.Double2String(eps) << "' is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_lognormal::calc_pdf_x_log", ssV.str() );
  }
  return (-(ONE/2)*pow2((log(x_val-eps)-lambda)/zeta))-(log(x_val-eps)+log(zeta)+log(sqrt(2*PI)));
}

const tdouble RBRV_entry_RV_lognormal::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller or equal than '" << GlobalVar.Double2String(eps) << "' is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_lognormal::calc_cdf_x", ssV.str() );
  }
  return rv_Phi((log(x_val-eps)-lambda)/zeta);
}

const tdouble RBRV_entry_RV_lognormal::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return ONE;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller or equal than '" << GlobalVar.Double2String(eps) << "' is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_lognormal::calc_sf_x", ssV.str() );
  }
  return rv_Phi((lambda-log(x_val-eps))/zeta);
}

const tdouble RBRV_entry_RV_lognormal::calc_entropy()
{
  return (ONE+log(2*PI*pow2(zeta)))/2+lambda;
}

const tdouble RBRV_entry_RV_lognormal::get_mean_current_config()
{
  return exp(lambda+pow2(zeta)/(ONE*2))+eps;
}

const tdouble RBRV_entry_RV_lognormal::get_sd_current_config()
{
  return sqrt(exp(pow2(zeta))-ONE)*exp(lambda+pow2(zeta)/(ONE*2));
}

const tdouble RBRV_entry_RV_lognormal::get_median_current_config()
{
  return exp(lambda)+eps;
}

const tdouble RBRV_entry_RV_lognormal::get_mode_current_config()
{
  return exp(lambda-pow2(zeta))+eps;
}

const bool RBRV_entry_RV_lognormal::check_x(const tdouble xV)
{
  return (xV>eps);
}

const bool RBRV_entry_RV_lognormal::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (p3?(p3->search_circref(fcr)):false) || (p4?(p4->search_circref(fcr)):false) || (epsilon?(epsilon->search_circref(fcr)):false);
}

py::dict RBRV_entry_RV_lognormal::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["lambda"] = lambda;
  res["zeta"] = zeta;
  res["epsilon"] = eps;
  res["mean"] = get_mean_current_config();
  res["sd"] = get_sd_current_config();
  res["entropy"] = calc_entropy();
  return res;
}

const tdouble RBRV_entry_RV_lognormal::get_CoeffOfVar_withoutEpsilon()
{
  return get_sd_current_config()/(get_mean_current_config()-eps);
}


RBRV_entry_RV_uniform::RBRV_entry_RV_uniform(const std::string& name, const tuint iID, FlxFunction* a, FlxFunction* b, const bool eval_once)
: RBRV_entry_RV_base(name,iID), a(a), b(b), eval_once(eval_once), av(ZERO), bv(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_uniform::RBRV_entry_RV_uniform",1);
    if (a) delete a;
    if (b) delete b;
    throw;
  }
}

RBRV_entry_RV_uniform::RBRV_entry_RV_uniform(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), a(nullptr), b(nullptr), eval_once(false), av(ZERO), bv(ZERO)
{

  try {
    a = parse_py_para("a", config);
    b = parse_py_para("b", config);
    eval_once = parse_py_para_as_bool("eval_once", config, false, false);

    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_uniform::RBRV_entry_RV_uniform_99",1);
    if (a) delete a;
    if (b) delete b;
    throw;
  }
}

RBRV_entry_RV_uniform::~RBRV_entry_RV_uniform()
{
  if (a) delete a;
  if (b) delete b;
}

void RBRV_entry_RV_uniform::eval_para()
{
  if (!eval_once || (eval_once && a) ) {
    av = a->calc();
    bv = b->calc();
    if (bv<=av) {
      std::ostringstream ssV;
      ssV << "Upper bound of uniform distribution (" << GlobalVar.Double2String(bv) << ") must not be smaller than lower bound (" << GlobalVar.Double2String(av) << ").";
      throw FlxException("RBRV_entry_RV_uniform::transform_y2x_2", ssV.str() );
    }
    if (eval_once) {
      delete a; a = NULL;
      delete b; b = NULL;
    }
  }
}

const tdouble RBRV_entry_RV_uniform::transform_y2x(const tdouble y_val)
{
  if (y_val<=-Y_INFTY_2) {
    return av;
  }
  if (y_val>=Y_INFTY_2) {
    return bv;
  }
  return (y_val>ZERO)?(bv-rv_Phi(-y_val)*(bv-av)):(rv_Phi(y_val)*(bv-av)+av);
}

const tdouble RBRV_entry_RV_uniform::transform_x2y(const tdouble& x_val)
{
  if (x_val>bv || x_val<av) {
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_uniform::transform_x2y", ssV.str() );
  }
  return rv_InvPhi((x_val-av)/(bv-av));
}

const tdouble RBRV_entry_RV_uniform::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bv || x_val<av) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_uniform::calc_pdf_x", ssV.str() );
  }
  return ONE/(bv-av);
}

const tdouble RBRV_entry_RV_uniform::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bv || x_val<av) {
    if (safeCalc) {
      if (x_val<av) return ZERO;
      else return ONE;
    }
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_uniform::calc_cdf_x", ssV.str() );
  }
  return (x_val-av)/(bv-av);
}

const tdouble RBRV_entry_RV_uniform::calc_icdf_x(const tdouble p_val)
{
  return p_val*(bv-av)+av;
}

const tdouble RBRV_entry_RV_uniform::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bv || x_val<av) {
    if (safeCalc) {
      if (x_val<av) return ONE;
      else return ZERO;
    }
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_uniform::calc_sf_x", ssV.str() );
  }
  return (bv-x_val)/(bv-av);
}

const tdouble RBRV_entry_RV_uniform::calc_entropy()
{
  return log(bv-av);
}

const tdouble RBRV_entry_RV_uniform::get_mean_current_config()
{
  return (av+bv)/2;
}

const tdouble RBRV_entry_RV_uniform::get_sd_current_config()
{
  return (bv+av)/sqrt(12*ONE);
}

const bool RBRV_entry_RV_uniform::check_x(const tdouble xV)
{
  return (xV<=bv&&xV>=av);
}

const bool RBRV_entry_RV_uniform::search_circref(FlxFunction* fcr)
{
  const bool bb = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return bb || (a?(a->search_circref(fcr)):false) || (b?(b->search_circref(fcr)):false);
}

const tdouble RBRV_entry_RV_uniform::get_HPD(const tdouble p)
{
  return (ONE-p)/2;
}

const tdouble RBRV_entry_RV_uniform::Inv_cdf_x(const tdouble p)
{
  #if FLX_DEBUG
    if ( p < ZERO || p > ONE) {
      std::ostringstream ssV;
      ssV << "p (" << p << ") out of range.";
      throw FlxException("RBRV_entry_RV_uniform::Inv_cdf_x", ssV.str() );
    }
  #endif
  return av+p*(bv-av);
}

RBRV_entry_RV_Gumbel::RBRV_entry_RV_Gumbel(const std::string& name, const tuint iID, const int methID, FlxFunction* p1, FlxFunction* p2, FlxFunction* p3, FlxFunction* p4, const bool eval_once)
: RBRV_entry_RV_base(name,iID), methID(methID), p1(p1), p2(p2), p3(p3), p4(p4), eval_once(eval_once), u(ZERO), alpha(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_Gumbel::RBRV_entry_RV_Gumbel",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (p3) delete p3;
    if (p4) delete p4;
    throw;
  }
}

RBRV_entry_RV_Gumbel::RBRV_entry_RV_Gumbel(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), methID(0), p1(nullptr), p2(nullptr), p3(nullptr), p4(nullptr), eval_once(false), u(ZERO), alpha(ZERO)
{

  try {
    if (config.contains("u")) {          // 0:u,alpha
      methID = 0;
      p1 = parse_py_para("u", config);
      p2 = parse_py_para("alpha", config);
    } else if (config.contains("mu")) {   // 1:mean,sd;
      methID = 1;
      p1 = parse_py_para("mu", config);
      p2 = parse_py_para("sd", config);
    } else if (config.contains("pr_2")) {   // 2:val_1,pr_1,val_2,pr_2
      methID = 2;
      p1 = parse_py_para("val_1", config);
      p2 = parse_py_para("pr_1", config);
      p3 = parse_py_para("val_2", config);
      p4 = parse_py_para("pr_2", config);
    } else {
      throw FlxException_NeglectInInteractive("RBRV_entry_RV_lognormal::RBRV_entry_RV_lognormal_70", "Required parameters to define distribution not found in Python <dict>.");
    }

    eval_once = parse_py_para_as_bool("eval_once", config, false, false);

    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_Gumbel::RBRV_entry_RV_Gumbel_99",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (p3) delete p3;
    if (p4) delete p4;
    throw;
  }
}

RBRV_entry_RV_Gumbel::~RBRV_entry_RV_Gumbel()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (p3) delete p3;
  if (p4) delete p4;
}

void RBRV_entry_RV_Gumbel::eval_para()
{
  if (!eval_once || (eval_once&&p1) ) {
    switch (methID) {
      case 0:
      {
        const tdouble p1d = p1->calc();
        const tdouble p2d = p2->cast2positive();
        u = p1d;
        alpha = p2d;
        break;
      }
      case 1:
      {
        const tdouble mu = p1->calc();
        const tdouble sigma = p2->cast2positive();
        alpha = PI/(std::sqrt(6.)*sigma);
        u = mu - GAMMA/alpha;
        break;
      }
      case 2:
      {
        // evaluate the parameters
          const tdouble x1 = p1->calc();
          const tdouble Pr1 = p2->cast2positive();
          if (Pr1>=ONE) {
            std::ostringstream ssV;
            ssV << "Expected a probability, which must not be larger than one (" << GlobalVar.Double2String(Pr1) << ").";
            throw FlxException("RBRV_entry_read_Gumbel::get_pars_01", ssV.str() );
          }
          const tdouble x2 = p3->calc();
          const tdouble Pr2 = p4->cast2positive();
          if (Pr2>=ONE) {
            std::ostringstream ssV;
            ssV << "Expected a probability, which must not be larger than one (" << GlobalVar.Double2String(Pr2) << ").";
            throw FlxException("RBRV_entry_read_Gumbel::get_pars_02", ssV.str() );
          }
        // solve the system of linear equations
          const tdouble a1 = log(-log(Pr2));
          const tdouble a2 = -log(-log(Pr1));
          const tdouble det = a1+a2;
          alpha = det/(x1 - x2);
          u = (a1*x1 + a2*x2)/det;
          if (alpha<=ZERO) {
            std::ostringstream ssV;
            ssV << "Scale parameter 'alpha' of Gumbel distribution must not become negative or zero (" << GlobalVar.Double2String(alpha) << ").";
            throw FlxException("RBRV_entry_read_Gumbel::get_pars_03", ssV.str() );            
          }
        break;
      }
      default:
        throw FlxException_Crude("RBRV_entry_RV_Gumbel::get_pars_04");
    }
    if (eval_once) {
      delete p1; p1=NULL;
      delete p2; p2=NULL;
      if (p3) { delete p3; p3=NULL; }
      if (p4) { delete p4; p4=NULL; }
    }
  }
}

const tdouble RBRV_entry_RV_Gumbel::transform_y2x(const tdouble y_val)
{
  return u-(std::log((rv_Phi(y_val)==ONE)?(rv_Phi(-y_val)):(-std::log(rv_Phi(y_val)))))/alpha;
}

const tdouble RBRV_entry_RV_Gumbel::transform_x2y(const tdouble& x_val)
{
  const tdouble ep = -exp(-alpha*(x_val-u));
  const tdouble p = exp(ep);
  if (p<=0.5) {
    return rv_InvPhi_noAlert(p);
  } else {
    const tdouble q = -expm1(ep);
    if (q==ZERO) {
      return Y_INFTY_1;
    } else {
      return -rv_InvPhi_noAlert(q);
    }
  }
}

const tdouble RBRV_entry_RV_Gumbel::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  return alpha*exp(-alpha*(x_val-u)-exp(-alpha*(x_val-u)));
}

const tdouble RBRV_entry_RV_Gumbel::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  const tdouble res = log(alpha)+(-alpha*(x_val-u)-exp(-alpha*(x_val-u)));
  return res;
}

const tdouble RBRV_entry_RV_Gumbel::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return exp(-exp(-alpha*(x_val-u)));
}

const tdouble RBRV_entry_RV_Gumbel::calc_entropy()
{
  return ONE+GAMMA-log(alpha);
}

const tdouble RBRV_entry_RV_Gumbel::get_mean_current_config()
{
  return u + GAMMA/alpha;
}

const tdouble RBRV_entry_RV_Gumbel::get_sd_current_config()
{
  return PI/(sqrt(6*ONE)*alpha);
}

const tdouble RBRV_entry_RV_Gumbel::get_median_current_config()
{
  return u - log(log(2*ONE))/alpha;
}

const tdouble RBRV_entry_RV_Gumbel::get_mode_current_config()
{
  return u;
}

const bool RBRV_entry_RV_Gumbel::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (p3?(p3->search_circref(fcr)):false) || (p4?(p4->search_circref(fcr)):false);
}

py::dict RBRV_entry_RV_Gumbel::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["u"] = u;
  res["alpha"] = alpha;
  res["mean"] = get_mean_current_config();
  res["sd"] = get_sd_current_config();
  res["entropy"] = calc_entropy();
  return res;
}

RBRV_entry_RV_normal_trunc::RBRV_entry_RV_normal_trunc(const std::string& name, const tuint iID, FlxFunction* m, FlxFunction* s, FlxFunction* a, FlxFunction* b, const bool eval_once)
:RBRV_entry_RV_base(name,iID), m(m), s(s), a(a), b(b), eval_once(eval_once), mV(ZERO), sV(ZERO), aV(ZERO), bV(ZERO), alpha(ZERO), beta(ZERO), q(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_normal_trunc::RBRV_entry_RV_normal_trunc",1);
    if (m) delete m;
    if (s) delete s;
    if (a) delete a;
    if (b) delete b;
    throw;
  }
}

RBRV_entry_RV_normal_trunc::~RBRV_entry_RV_normal_trunc()
{
  if (m) delete m;
  if (s) delete s;
  if (a) delete a;
  if (b) delete b;
}

void RBRV_entry_RV_normal_trunc::eval_para()
{
  if (!eval_once || (eval_once&&m) ) {
    mV = m->calc();
    sV = s->calc();
    aV = a?(a->calc()):(mV-1e5*sV);
    bV = b?(b->calc()):(mV+1e5*sV);
    alpha = (aV-mV)/sV;
    beta = (bV-mV)/sV;
    q = rv_Phi_diff(alpha,beta);
    if (q<1e-100) {
      throw FlxException("RBRV_entry_RV_normal_trunc::get_pars","Parametrization numerically instable.");
    }
    if (eval_once) {
      delete m; m=NULL;
      delete s; s=NULL;
      if (a) { delete a; a=NULL; }
      if (b) { delete b; b=NULL; }
    }
  }
}

const tdouble RBRV_entry_RV_normal_trunc::transform_y2x(const tdouble y_val)
{
  tdouble res;
  if (y_val<=ZERO) {
    if (alpha<=ZERO) {
      const tdouble p_low = rv_Phi(y_val)*q+rv_Phi(alpha);
      res = (sV*rv_InvPhi_noAlert(p_low) + mV);
    } else {
      const tdouble p_up = rv_Phi(-alpha) - rv_Phi(y_val)*q;
      res =  (mV - sV*rv_InvPhi_noAlert(p_up));
    }
  } else {
    if (beta<=ZERO) {
      const tdouble p_low = rv_Phi(beta) - rv_Phi(-y_val)*q;
      res =  (sV*rv_InvPhi_noAlert(p_low) + mV);
    } else {
      const tdouble p_up = rv_Phi(-y_val)*q+rv_Phi(-beta);
      res =  (mV - sV*rv_InvPhi_noAlert(p_up));
    }
  }
  if (res<aV) {
    if (aV-res>GlobalVar.TOL()*sV) {
      throw FlxException_Crude("RBRV_entry_RV_normal_trunc::transform_y2x_01");
    }
    res = aV + GlobalVar.TOL()*sV;
  } else if (res>bV) {
    if (res-bV>GlobalVar.TOL()*sV) {
      throw FlxException_Crude("RBRV_entry_RV_normal_trunc::transform_y2x_02");
    }
    res = bV - GlobalVar.TOL()*sV;
  }
  return res;
}

const tdouble RBRV_entry_RV_normal_trunc::transform_x2y(const tdouble& x_val)
{
  if (x_val>bV || x_val<aV) {
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(aV) << ";" << GlobalVar.Double2String(bV) << "].";
    throw FlxException("RBRV_entry_RV_normal_trunc::transform_x2y", ssV.str() );
  }
  return rv_InvPhi_noAlert((rv_Phi((x_val-mV)/sV)-rv_Phi(alpha))/q);
}

const tdouble RBRV_entry_RV_normal_trunc::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bV || x_val<aV) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(aV) << ";" << GlobalVar.Double2String(bV) << "].";
    throw FlxException("RBRV_entry_RV_normal_trunc::calc_pdf_x", ssV.str() );
  }
  if (q==ZERO) return ZERO;
  return rv_phi((x_val-mV)/sV)/(sV*q);
}

const tdouble RBRV_entry_RV_normal_trunc::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bV || x_val<aV) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(aV) << ";" << GlobalVar.Double2String(bV) << "].";
    throw FlxException("RBRV_entry_RV_normal_trunc::calc_pdf_x", ssV.str() );
  }
  if (q==ZERO) return(log(ZERO));
  return rv_phi_log((x_val-mV)/sV)-log(sV*q);
}

const tdouble RBRV_entry_RV_normal_trunc::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bV || x_val<aV) {
    if (safeCalc) {
      if (x_val<aV) return ZERO;
      else return ONE;
    }
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(aV) << ";" << GlobalVar.Double2String(bV) << "].";
    throw FlxException("RBRV_entry_RV_normal_trunc::calc_cdf_x", ssV.str() );
  }
  return (rv_Phi((x_val-mV)/sV)-rv_Phi(alpha))/q;
}

const tdouble RBRV_entry_RV_normal_trunc::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bV || x_val<aV) {
    if (safeCalc) {
      if (x_val<aV) return ONE;
      else return ZERO;
    }
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(aV) << ";" << GlobalVar.Double2String(bV) << "].";
    throw FlxException("RBRV_entry_RV_normal_trunc::calc_sf_x", ssV.str() );
  }
  return (rv_Phi((mV-x_val)/sV)-rv_Phi(beta))/q;
}

const tdouble RBRV_entry_RV_normal_trunc::get_mean_current_config()
{
  return mV+(rv_phi(alpha)-rv_phi(beta))/q*sV;
}

const tdouble RBRV_entry_RV_normal_trunc::get_sd_current_config()
{
  return sV*sqrt(
    ONE
    + (alpha*rv_phi(alpha)-beta*rv_phi(beta))/q
    - pow2((rv_phi(alpha)-rv_phi(beta))/q)
  );
}

const tdouble RBRV_entry_RV_normal_trunc::get_median_current_config()
{
  return mV+rv_InvPhi((rv_Phi(alpha)+rv_Phi(beta))/2)*sV;
}

const tdouble RBRV_entry_RV_normal_trunc::get_mode_current_config()
{
  if (mV<alpha) return alpha;
  if (mV>beta) return beta;
  return mV;
}

const tdouble RBRV_entry_RV_normal_trunc::calc_entropy()
{
  return log(sqrt(2*PI*exp(ONE))*sV*q)+(alpha*rv_phi(alpha)-beta*rv_phi(beta))/(2*q);
}

const bool RBRV_entry_RV_normal_trunc::check_x(const tdouble xV)
{
  return (xV<=bV&&xV>=aV);
}

const bool RBRV_entry_RV_normal_trunc::search_circref(FlxFunction* fcr)
{
  const bool bb = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return bb || (m?(m->search_circref(fcr)):false) || (s?(s->search_circref(fcr)):false) || (a?(a->search_circref(fcr)):false) || (b?(b->search_circref(fcr)):false);
}

py::dict RBRV_entry_RV_normal_trunc::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["m"] = mV;
  res["s"] = sV;
  res["a"] = aV;
  res["b"] = bV;
  res["alpha"] = alpha;
  res["beta"] = beta;
  res["q"] = q;
  return res;
}

RBRV_entry_RV_beta::RBRV_entry_RV_beta(const std::string& name, const tuint iID, const bool is_mean, FlxFunction* p1, FlxFunction* p2, FlxFunction* a, FlxFunction* b, const bool eval_once)
: RBRV_entry_RV_base(name,iID), is_mean(is_mean), p1(p1), p2(p2), a(a), b(b), eval_once(eval_once), alpha(ZERO), beta(ZERO), av(ZERO), bv(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_beta::RBRV_entry_RV_beta",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (a) delete a;
    if (b) delete b;
    throw;
  }
}

RBRV_entry_RV_beta::RBRV_entry_RV_beta(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), is_mean(false), p1(nullptr), p2(nullptr), a(nullptr), b(nullptr), eval_once(false), alpha(ZERO), beta(ZERO), av(ZERO), bv(ZERO)
{
  try {
    if (config.contains("mu")) {          // 0:lambda,zeta;
      is_mean = true;
      p1 = parse_py_para("mu", config);
      p2 = parse_py_para("sd", config);
    } else if (config.contains("alpha")) {   // 1:mean,sd;
      is_mean = false;
      p1 = parse_py_para("alpha", config);
      p2 = parse_py_para("beta", config);
    } else {
      throw FlxException_NeglectInInteractive("RBRV_entry_RV_beta::RBRV_entry_RV_beta_70", "Required parameters to define distribution not found in Python <dict>.");
    }

    a =  parse_py_para("a", config, false);
    b =  parse_py_para("b", config, false);

    eval_once = parse_py_para_as_bool("eval_once", config, false, false);

    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_beta::RBRV_entry_RV_beta_99",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (a) delete a;
    if (b) delete b;
    throw;
  }
}

RBRV_entry_RV_beta::~RBRV_entry_RV_beta()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (a) delete a;
  if (b) delete b;
}

void RBRV_entry_RV_beta::eval_para()
{
  if (!eval_once || (eval_once&&p1) ) {
    const tdouble p1d = (is_mean?p1->calc():p1->cast2positive());
    const tdouble p2d = p2->cast2positive();
    av = (a?(a->calc()):ZERO);
    bv = (b?(b->calc()):ONE);
    if (bv<=av) {
      std::ostringstream ssV;
      ssV << "'b' (" << GlobalVar.Double2String(bv) << ") must be larger than 'a' (" << GlobalVar.Double2String(av) << ").";
      throw FlxException("RBRV_entry_RV_beta::get_pars_1", ssV.str() );
    }
    if (is_mean) {
      const tdouble& mu = p1d;
      const tdouble& sigma = p2d;
      if (mu < av || mu > bv) {
        std::ostringstream ssV;
        ssV << "'mu' (" << GlobalVar.Double2String(mu) << ") must be within the bounds of 'a' (" << GlobalVar.Double2String(av) << ") and 'b' (" << GlobalVar.Double2String(bv) << ").";
        throw FlxException("RBRV_entry_RV_beta::get_pars_2", ssV.str() );
      }
      if (pow2(sigma)>=(mu-av)*(bv-mu)) {
        std::ostringstream ssV;
        ssV << this->name << ": 'sigma^2' (" << GlobalVar.Double2String(sigma) << "²=" << GlobalVar.Double2String(pow2(sigma)) << ") must be smaller than '(mu-a)*(b-mu)' (" << GlobalVar.Double2String((bv-mu)*(mu-av)) << ") (mu=" << GlobalVar.Double2String(mu) << ").";
        throw FlxException("RBRV_entry_RV_beta::get_pars_3", ssV.str() );
      }
      alpha = (mu-av)/(bv-av)*((mu-av)*(bv-mu)/pow2(sigma)-ONE);
      beta =  (bv-mu)/(bv-av)*((mu-av)*(bv-mu)/pow2(sigma)-ONE);
    } else {
      alpha = p1d;
      beta = p2d;
    }
    if (eval_once) {
      delete p1; p1 = NULL;
      delete p2; p2 = NULL;
      if (a) { delete a; a = NULL; }
      if (b) { delete b; b = NULL; }
    }
  }
}

const tdouble RBRV_entry_RV_beta::transform_y2x(const tdouble y_val)
{
  if (y_val<=-Y_INFTY_2) {
    return av;
  }
  if (y_val>=Y_INFTY_2) {
    return bv;
  }
  const tdouble res = iBeta_reg_inv(alpha,beta,rv_Phi(y_val));
  return (res*(bv-av)+av);
}

const tdouble RBRV_entry_RV_beta::transform_x2y(const tdouble& x_val)
{
  if (x_val>bv || x_val<av) {
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_beta::transform_x2y", ssV.str() );
  }
  const tdouble xx = (x_val-av)/(bv-av);        // scale x to [0;1]
  return rv_InvPhi(iBeta_reg(alpha,beta,xx));
}

const tdouble RBRV_entry_RV_beta::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bv || x_val<av) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_beta::calc_pdf_x", ssV.str() );
  }
  const tdouble xx = (x_val-av)/(bv-av);        // scale x to [0;1]
  const tdouble lnbetafun = BetaFunLn(alpha,beta);
  return exp(
    (alpha-ONE)*log(xx)+(beta-ONE)*log(ONE-xx)
    -lnbetafun
  )/(bv-av);
}

const tdouble RBRV_entry_RV_beta::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bv || x_val<av) {
    if (safeCalc) return log(ZERO);
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_beta::calc_pdf_x", ssV.str() );
  }
  const tdouble xx = (x_val-av)/(bv-av);        // scale x to [0;1]
  const tdouble lnbetafun = BetaFunLn(alpha,beta);
  return (
    (alpha-ONE)*log(xx)+(beta-ONE)*log(ONE-xx)
    -lnbetafun
  )-log(bv-av);
}

const tdouble RBRV_entry_RV_beta::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bv || x_val<av) {
    if (safeCalc) {
      if (x_val<av) return ZERO;
      else return ONE;
    }
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_beta::calc_cdf_x", ssV.str() );
  }
  const tdouble xx = (x_val-av)/(bv-av);        // scale x to [0;1]
  return iBeta_reg(alpha,beta,xx);
}

const tdouble RBRV_entry_RV_beta::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bv || x_val<av) {
    if (safeCalc) {
      if (x_val<av) return ONE;
      else return ZERO;
    }
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_beta::calc_sf_x", ssV.str() );
  }
  const tdouble xx = (bv-x_val)/(bv-av);        // scale x to [0;1]
  return iBeta_reg(beta,alpha,xx);
}

const tdouble RBRV_entry_RV_beta::calc_entropy()
{
  const tdouble lnbetafun = BetaFunLn(alpha,beta);
  const tdouble entropy01 = lnbetafun 
          - (alpha-ONE)*flxdigamma(alpha)
          - (beta-ONE)*flxdigamma(beta)
          + (alpha+beta-2*ONE)*flxdigamma(alpha+beta);
  return entropy01/(bv-av) + log(bv-av);
}

const tdouble RBRV_entry_RV_beta::get_mean_current_config()
{
  return alpha/(alpha+beta)*(bv-av)+av;
}

const tdouble RBRV_entry_RV_beta::get_sd_current_config()
{
  return sqrt(alpha*beta/(pow2(alpha+beta)*(alpha+beta+ONE))) * (bv-av);
}

const tdouble RBRV_entry_RV_beta::get_median_current_config()
{
  return transform_y2x(ZERO);
}

const tdouble RBRV_entry_RV_beta::get_mode_current_config()
{
  if (alpha>ONE && beta>ONE) {
    return (alpha-ONE)/(alpha+beta-2*ONE)*(bv-av)+av;
  }
  if (beta>ONE && alpha<=ONE) return ZERO;
  if (alpha>ONE && beta<=ONE) return ONE;
  throw FlxException_NotImplemented("RBRV_entry_RV_beta::get_mode_current_config");
}

const bool RBRV_entry_RV_beta::check_x(const tdouble xV)
{
  return (xV<=bv&&xV>=av);
}

const bool RBRV_entry_RV_beta::search_circref(FlxFunction* fcr)
{
  const bool bb = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return bb || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (a?(a->search_circref(fcr)):false) || (b?(b->search_circref(fcr)):false);
}

py::dict RBRV_entry_RV_beta::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["alpha"] = alpha;
  res["beta"] = beta;
  res["a"] = av;
  res["b"] = bv;
  res["mean"] = get_mean_current_config();
  res["sd"] = get_sd_current_config();
  res["entropy"] = calc_entropy();
  return res;
}

const tdouble RBRV_entry_RV_beta::Inv_cdf_x(const tdouble p)
{
  #if FLX_DEBUG
    if ( p < ZERO || p > ONE) {
      std::ostringstream ssV;
      ssV << "p (" << p << ") out of range.";
      throw FlxException("RBRV_entry_RV_beta::Inv_cdf_x", ssV.str() );
    }
  #endif
  const tdouble x = iBeta_reg_inv(alpha,beta,p);
  return x*(bv-av)+av;        // scale x to [a;b]
}

RBRV_entry_RV_exponential::RBRV_entry_RV_exponential(const std::string& name, const tuint iID, FlxFunction* lambda, FlxFunction* epsilon)
: RBRV_entry_RV_base(name,iID), lambda(lambda), epsilon(epsilon), lambdaV(ZERO), eps(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_exponential::RBRV_entry_RV_exponential",1);
    if (lambda) delete lambda;
    if (epsilon) delete epsilon;
    throw;
  }
}

RBRV_entry_RV_exponential::RBRV_entry_RV_exponential(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), lambda(nullptr), epsilon(nullptr), lambdaV(ZERO), eps(ZERO)
{
  try {
    lambda = parse_py_para("lambda", config);
    epsilon = parse_py_para("epsilon", config, false);
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_exponential::RBRV_entry_RV_exponential_99",1);
    if (lambda) delete lambda;
    if (epsilon) delete epsilon;
    throw;
  }
}

RBRV_entry_RV_exponential::~RBRV_entry_RV_exponential()
{
  if (lambda) delete lambda;
  if (epsilon) delete epsilon;
}

void RBRV_entry_RV_exponential::eval_para()
{
  lambdaV = lambda->cast2positive();
  eps = (epsilon?(epsilon->calc()):ZERO);
}

const tdouble RBRV_entry_RV_exponential::transform_y2x(const tdouble y_val)
{
  tdouble res = -ONE*std::log((y_val>ZERO)?(rv_Phi(-y_val)):(ONE-rv_Phi(y_val)))/lambdaV;
  res += eps;
  return res;
}

const tdouble RBRV_entry_RV_exponential::transform_x2y(const tdouble& x_val)
{
  if (x_val<eps) {
    std::ostringstream ssV;
    ssV << "A negative value (" << GlobalVar.Double2String(x_val) << ") is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_exponential::transform_x2y", ssV.str() );
  }
  return rv_InvPhi_noAlert(-expm1(-lambdaV*(x_val-eps)));
}

const tdouble RBRV_entry_RV_exponential::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<eps) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A negative value (" << GlobalVar.Double2String(x_val) << ") is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_exponential::calc_pdf_x", ssV.str() );
  }
  return lambdaV*exp(-lambdaV*(x_val-eps));
}

const tdouble RBRV_entry_RV_exponential::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<eps) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A negative value (" << GlobalVar.Double2String(x_val) << ") is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_exponential::calc_pdf_x", ssV.str() );
  }
  return log(lambdaV)+(-lambdaV*(x_val-eps));
}

const tdouble RBRV_entry_RV_exponential::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return ONE-this->calc_sf_x(x_val,safeCalc);
}

const tdouble RBRV_entry_RV_exponential::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<eps) {
    if (safeCalc) return ONE;
    std::ostringstream ssV;
    ssV << "A negative value (" << GlobalVar.Double2String(x_val) << ") is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_exponential::calc_sf_x", ssV.str() );
  }
  return exp(-lambdaV*(x_val-eps));
}

const tdouble RBRV_entry_RV_exponential::calc_entropy()
{
  return ONE - log(lambdaV);
}

const tdouble RBRV_entry_RV_exponential::get_mean_current_config()
{
  return eps+ONE/lambdaV;
}

const tdouble RBRV_entry_RV_exponential::get_sd_current_config()
{
  return ONE/lambdaV;
}

const tdouble RBRV_entry_RV_exponential::get_median_current_config()
{
  return eps+log(2*ONE)/lambdaV;
}

const tdouble RBRV_entry_RV_exponential::get_mode_current_config()
{
  return eps;
}

const bool RBRV_entry_RV_exponential::check_x(const tdouble xV)
{
  return (xV>=eps);
}

const bool RBRV_entry_RV_exponential::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (lambda?(lambda->search_circref(fcr)):false) || (epsilon?(epsilon->search_circref(fcr)):false);
}

py::dict RBRV_entry_RV_exponential::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["lambda"] = lambdaV;
  res["epsilon"] = eps;
  res["mean"] = get_mean_current_config();
  res["sd"] = get_sd_current_config();
  res["entropy"] = calc_entropy();
  return res;
}


RBRV_entry_RV_gamma::RBRV_entry_RV_gamma(const std::string& name, const tuint iID, const bool is_mean, FlxFunction* p1, FlxFunction* p2, FlxFunction* epsilon, const bool eval_once)
: RBRV_entry_RV_base(name,iID), is_mean(is_mean), p1(p1), p2(p2), epsilon(epsilon), eval_once(eval_once), k(ZERO), lambda(ZERO), eps(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_gamma::RBRV_entry_RV_gamma",1);
    if (p2) delete p2;
    if (p1) delete p1;
    if (epsilon) delete epsilon;
    throw;
  }
}

RBRV_entry_RV_gamma::~RBRV_entry_RV_gamma()
{
  if (p2) delete p2;
  if (p1) delete p1;
  if (epsilon) delete epsilon;
}

void RBRV_entry_RV_gamma::eval_para()
{
  if (!eval_once || (eval_once&&p1) ) {
    // compute parameters
      if (epsilon) {
        eps = epsilon->calc();
      } else {
        eps = ZERO;
      }
      if (is_mean) {
        const tdouble m = p1->cast2positive();
        const tdouble s2 = pow2(p2->cast2positive());
        if (m<=eps) {
          std::ostringstream ssV;
          ssV << "The mean (" << GlobalVar.Double2String(m) << ") must not be smaller than epsilon (" << GlobalVar.Double2String(eps) << ").";
          throw FlxException("RBRV_entry_RV_gamma::get_paras_1", ssV.str() );
        }
        lambda = (m-eps)/s2;
        k = pow2(m-eps)/s2;
      } else {
        k = p1->cast2positive();
        lambda = p2->cast2positive();
      }
    if (eval_once) {
      delete p1; p1 = NULL;
      delete p2; p2 = NULL;
      if (epsilon) { delete epsilon; epsilon=NULL; }
    }
  }
}

const tdouble RBRV_entry_RV_gamma::transform_y2x(const tdouble y_val)
{
  tdouble res;
  if (y_val<=ZERO) {
    res = flxgamma_rl_inv(k,rv_Phi(y_val));
  } else {
    res = flxgamma_ru_inv(k,rv_Phi(-y_val));
  }
  res /= lambda;
  res += eps;
  return res;
}

const tdouble RBRV_entry_RV_gamma::transform_x2y(const tdouble& x_val)
{
  if (x_val<=eps) {
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(eps) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_gamma::transform_x2y", ssV.str() );
  }
  if (x_val<=eps+k/lambda) {
    return rv_InvPhi_noAlert(flxgamma_rl(k,lambda*(x_val-eps)));
  } else {
    return -rv_InvPhi_noAlert(flxgamma_ru(k,lambda*(x_val-eps)));
  }
}

const tdouble RBRV_entry_RV_gamma::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than " << GlobalVar.Double2String(eps) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_gamma::calc_cdf_x_01", ssV.str() );
  }
  const tdouble res = flxgamma_rl(k,lambda*(x_val-eps));
  if (std::isnan(res)) {
    std::ostringstream ssV;
    ssV << "CDF of gamma for k=" << GlobalVar.Double2String(k) << " and lambda=" << GlobalVar.Double2String(lambda) << " (with eps=" << GlobalVar.Double2String(eps) << ") returned 'nan'.";
    throw FlxException("RBRV_entry_RV_gamma::calc_cdf_x_02", ssV.str() );
  }
  return res;
}

const tdouble RBRV_entry_RV_gamma::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return ONE;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than " << GlobalVar.Double2String(eps) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_gamma::calc_sf_x_01", ssV.str() );
  }
  const tdouble res = flxgamma_ru(k,lambda*(x_val-eps));
  if (std::isnan(res)) {
    std::ostringstream ssV;
    ssV << "Survival function of gamma for k=" << GlobalVar.Double2String(k) << " and lambda=" << GlobalVar.Double2String(lambda) << " (with eps=" << GlobalVar.Double2String(eps) << ") returned 'nan'.";
    throw FlxException("RBRV_entry_RV_gamma::calc_sf_x_02", ssV.str() );
  }
  return res;
}

const tdouble RBRV_entry_RV_gamma::calc_entropy()
{
  return k-log(lambda)+GammaLn(k)+(ONE-k)*flxdigamma(k);
}

const tdouble RBRV_entry_RV_gamma::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(eps) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_gamma::calc_pdf_x", ssV.str() );
  }
  return (pow(x_val-eps,k-ONE)*exp(-(x_val-eps)*lambda)*pow(lambda,k))/flxgamma(k);
}

const tdouble RBRV_entry_RV_gamma::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return log(ZERO);
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(eps) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_gamma::calc_pdf_x", ssV.str() );
  }
  return ((k-ONE)*log(x_val-eps)+(-(x_val-eps)*lambda)+k*log(lambda))-GammaLn(k);
}

const tdouble RBRV_entry_RV_gamma::get_mean_current_config()
{
  return eps+k/lambda;
}

const tdouble RBRV_entry_RV_gamma::get_sd_current_config()
{
  return sqrt(k)/lambda;
}

const tdouble RBRV_entry_RV_gamma::get_median_current_config()
{
  return transform_y2x(ZERO);
}

const tdouble RBRV_entry_RV_gamma::get_mode_current_config()
{
  if (k>=ONE) return eps+(k-ONE)/lambda;
  throw FlxException_NotImplemented("RBRV_entry_RV_gamma::get_mode_current_config");
}

const bool RBRV_entry_RV_gamma::check_x(const tdouble xV)
{
  return (xV>eps);
}

const bool RBRV_entry_RV_gamma::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (epsilon?(epsilon->search_circref(fcr)):false);
}

py::dict RBRV_entry_RV_gamma::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["k"] = k;
  res["lambda"] = lambda;
  res["epsilon"] = eps;
  res["mean"] = get_mean_current_config();
  res["sd"] = get_sd_current_config();
  res["entropy"] = calc_entropy();
  return res;
}

RBRV_entry_RV_Poisson::RBRV_entry_RV_Poisson(const std::string& name, const tuint iID, FlxFunction* mean)
: RBRV_entry_RV_base(name,iID), mean(mean), meanV(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_Poisson::RBRV_entry_RV_Poisson",1);
    if (mean) delete mean;
    throw;
  }
}

RBRV_entry_RV_Poisson::~RBRV_entry_RV_Poisson()
{
  delete mean;
}

void RBRV_entry_RV_Poisson::eval_para()
{
  meanV = mean->cast2positive();
}

const tdouble RBRV_entry_RV_Poisson::transform_y2x(const tdouble y_val)
{
  const tdouble p = rv_Phi(y_val);
  const tuint ns = meanV*2;        // number of steps
  // determin the relevant interval
    tuint lc = 0;                // lower bound (in the set)
    while (true) {
      tdouble t = flxgamma_ru(lc+ns,meanV);
      if (t>=p) {
        break;
      } else {
        lc += ns;
      }
    }
    tuint uc = lc + ns;                // upper bound (not in the set!)
    tuint dc = ns;                // distance of bounds
  // search for p
    while (dc > 1) {
      tuint dh = dc / 2;
      tdouble t = flxgamma_ru(lc+dh,meanV);
      if (t>=p) {        // decrease upper bound
        uc = lc + dh;
        dc = dh;
      } else {                // increase lower bound
        lc = lc + dh;
        dc = uc - lc;
      }
    }
  return lc;
}

const tdouble RBRV_entry_RV_Poisson::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return flxgamma_ru(std::floor(x_val)+1,meanV);
}

const tdouble RBRV_entry_RV_Poisson::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  return flxgamma_rl(std::floor(x_val)+1,meanV);
}

const tdouble RBRV_entry_RV_Poisson::get_mean_current_config()
{
  return meanV;
}

const tdouble RBRV_entry_RV_Poisson::get_sd_current_config()
{
  return sqrt(meanV);
}

const bool RBRV_entry_RV_Poisson::check_x(const tdouble xV)
{
  throw FlxException_NotImplemented("RBRV_entry_RV_Poisson::check_x");
}

const bool RBRV_entry_RV_Poisson::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (mean?(mean->search_circref(fcr)):false);
}


RBRV_entry_RV_Binomial::RBRV_entry_RV_Binomial(const std::string& name, const tuint iID, FlxFunction* p, FlxFunction* N, const bool eval_once)
: RBRV_entry_RV_base(name,iID), p(p), N(N), eval_once(eval_once), _p(ZERO), _N(0)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_Binomial::RBRV_entry_RV_Binomial",1);
    if (p) delete p;
    if (N) delete N;
    throw;
  }
}

RBRV_entry_RV_Binomial::~RBRV_entry_RV_Binomial()
{
  if (p) delete p;
  if (N) delete N;
}

void RBRV_entry_RV_Binomial::eval_para()
{
  if (!eval_once || (eval_once&&p) ) {
    _p = p->cast2positive_or0();
    if (_p>ONE) {
      std::ostringstream ssV;
      ssV << "'p' (" << GlobalVar.Double2String(_p) << ") denotes a probability and must be smaller or equal than 1.";
      throw FlxException("RBRV_entry_RV_Binomial::get_pars", ssV.str() );
    }
    _N = N->cast2tuintW0();
    if (eval_once) {
      delete p; p=NULL;
      delete N; N=NULL;
    }
  }
}

const tdouble RBRV_entry_RV_Binomial::transform_y2x(const tdouble y_val)
{
  const tdouble cdf_p = rv_Phi(y_val);
  if (_N==0) {
    return ZERO;
  }
  tuint uk = _N;                // upper k
  tuint lk = 0;                // lower k
    tdouble lp = rv_cdf_Binomial(_p,_N,lk);
  if (cdf_p<=lp) {
    return ZERO;
  }
  while (uk-lk>1) {
    tuint dh = lk+ ((uk-lk)/2);
    tdouble dhp = rv_cdf_Binomial(_p,_N,dh);
    if (cdf_p<=dhp) {        // take lower part
      uk = dh;
    } else {                // take upper part
      lk = dh;
      lp = dhp;
    }
  }
  return uk;
}

const tdouble RBRV_entry_RV_Binomial::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<ZERO) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than 0 is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Binomial::calc_cdf_x_1", ssV.str() );
  } else if (x_val>=tdouble(_N)) {
    if (safeCalc) return ONE;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") larger than " << _N << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Binomial::calc_cdf_x_2", ssV.str() );
  }
  const tdouble x_floor = floor(x_val);
  return iBeta_reg(tdouble(_N)-x_floor,x_floor+ONE,ONE-_p);
}

const tdouble RBRV_entry_RV_Binomial::get_mean_current_config()
{
  return _p*_N;
}

const tdouble RBRV_entry_RV_Binomial::get_sd_current_config()
{
  return sqrt(_p*_N*(ONE-_p));
}

const bool RBRV_entry_RV_Binomial::check_x(const tdouble xV)
{
  throw FlxException_NotImplemented("RBRV_entry_RV_Binomial::check_x");
}

const bool RBRV_entry_RV_Binomial::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p?(p->search_circref(fcr)):false) || (N?(N->search_circref(fcr)):false);
}


RBRV_entry_RV_Cauchy::RBRV_entry_RV_Cauchy(const std::string& name, const tuint iID, FlxFunction* loc, FlxFunction* scale)
: RBRV_entry_RV_base(name,iID), loc(loc), scale(scale), locv(ZERO), scalev(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_Cauchy::RBRV_entry_RV_Cauchy",1);
    delete loc;
    delete scale;
    throw;
  }
}

RBRV_entry_RV_Cauchy::~RBRV_entry_RV_Cauchy()
{
  delete loc;
  delete scale;
}

void RBRV_entry_RV_Cauchy::eval_para()
{
  locv = loc->calc();
  scalev = scale->cast2positive();
}

const tdouble RBRV_entry_RV_Cauchy::get_mean_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_Cauchy::get_mean_current_config");
}

const tdouble RBRV_entry_RV_Cauchy::get_sd_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_Cauchy::get_sd_current_config");
}

const tdouble RBRV_entry_RV_Cauchy::get_median_current_config()
{
  return locv;
}

const tdouble RBRV_entry_RV_Cauchy::get_mode_current_config()
{
  return locv;
}

const bool RBRV_entry_RV_Cauchy::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (loc?(loc->search_circref(fcr)):false) || (scale?(scale->search_circref(fcr)):false);
}


const tdouble RBRV_entry_RV_Cauchy::transform_y2x(const tdouble y_val)
{
  const tdouble p = (y_val>ZERO)?(ONE/2-rv_Phi(-y_val)):(rv_Phi(y_val)-ONE/2);
  return locv+scalev*tan(PI*(p));
}

const tdouble RBRV_entry_RV_Cauchy::transform_x2y(const tdouble& x_val)
{
  const tdouble p = std::atan((x_val-locv)/scalev)/PI;
  return (p<=ZERO)?(rv_InvPhi_noAlert(p+ONE/2)):(-rv_InvPhi_noAlert(ONE/2-p));
}

const tdouble RBRV_entry_RV_Cauchy::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  return (scalev/(pow2(x_val-locv)+pow2(scalev)))/PI;
}

const tdouble RBRV_entry_RV_Cauchy::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return std::atan((x_val-locv)/scalev)/PI+ONE/2;
}


RBRV_entry_RV_Weibull::RBRV_entry_RV_Weibull(const std::string& name, const tuint iID, const bool is_mean, FlxFunction* p1, FlxFunction* p2, FlxFunction* epsilon, const bool eval_once)
: RBRV_entry_RV_base(name,iID), is_mean(is_mean), p1(p1), p2(p2), epsilon(epsilon), eval_once(eval_once), k(ZERO), lambda(ZERO), eps(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_Weibull::RBRV_entry_RV_Weibull",1);
    if (p2) delete p2;
    if (p1) delete p1;
    if (epsilon) delete epsilon;
    throw;
  }
}

RBRV_entry_RV_Weibull::~RBRV_entry_RV_Weibull()
{
  if (p2) delete p2;
  if (p1) delete p1;
  if (epsilon) delete epsilon;
}

void RBRV_entry_RV_Weibull::eval_para()
{
  if (!eval_once || (eval_once&&p1) ) {
    // compute parameters
      if (epsilon) {
        eps = epsilon->calc();
      } else {
        eps = ZERO;
      }
      if (is_mean) {
        const tdouble m = p1->cast2positive();
        const tdouble s = p2->cast2positive();
        if (m<=eps) {
          std::ostringstream ssV;
          ssV << "The mean (" << GlobalVar.Double2String(m) << ") must not be smaller than epsilon (" << GlobalVar.Double2String(eps) << ").";
          throw FlxException("RBRV_entry_RV_Weibull::get_paras_01", ssV.str() );
        }
        const tdouble cov = s/m;
        // root search
          lambda = ONE;        // dummy-assignment (to eliminate lambda in equations)
          // initial value for k
            k = ONE;        
            tdouble cv0 = get_cov_help()-cov;
            tdouble k0 = log(k);
          // second value for k
            tdouble cv1 = cv0;
            tdouble k1 = k0;
            tuint N = 0;
            do {
              if (cv0>=ZERO) {
                k *= 2;
              } else {
                k /= 2;
              }
              cv1 = get_cov_help()-cov;
              k1 = log(k);
              ++N;
            } while (cv0*cv1>ZERO && N<100);
            #if FLX_DEBUG
              if (N>=100 && cv0*cv1>ZERO) {
                throw FlxException_Crude("RBRV_entry_RV_Weibull::get_pars_02");
              }
            #endif
          // start the iteration (regula falsi)
            N = 0;
            tdouble k2, cv2;
            do {
              k2 = (k0*cv1-k1*cv0)/(cv1-cv0);
              k = exp(k2);
              const tdouble cov_of_k = get_cov_help();
              if (!std::isfinite(cov_of_k)) {
                throw FlxException("RBRV_entry_RV_Weibull::get_pars_03","The value is not finite.");
              }
              cv2 = cov_of_k-cov;
              // Anderson/Björck algorithm
                tdouble m = ONE - cv2/cv1;
                if (m<=ZERO) m = ONE/2;
              //GlobalVar.slogcout(1) << "     regula-falsi " << N << "   " << k0 << "   " << k1 << "   " << k2 << "   " << cv0 << "   " << cv1 << "   " << cv2 << "   " << m << std::endl;
              if (cv2*cv1<=0) {
                k0 = k1;
                cv0 = cv1;
                k1 = k2;
                cv1 = cv2;
              } else {
                cv0 *= m;
                k1 = k2;
                cv1 = cv2;
              }
              ++N;
            } while (fabs(cv2)>=1e-8 && N<100 );
            #if FLX_DEBUG
              if (N>=100) {
                  throw FlxException_Crude("RBRV_entry_RV_Weibull::get_pars_03");
              }
            #endif
          // detect lambda
            const tdouble l1 = m/get_mean_help();
            const tdouble l2 = s/get_sd_help();
            lambda = (l1+l2)/2;
      } else {
        k = p1->cast2positive();
        lambda = p2->cast2positive();
      }
    if (eval_once) {
      delete p1; p1 = NULL;
      delete p2; p2 = NULL;
      if (epsilon) { delete epsilon; epsilon=NULL; }
    }
  }
}

const tdouble RBRV_entry_RV_Weibull::get_mean_current_config()
{
  return get_mean_help();
}

const tdouble RBRV_entry_RV_Weibull::get_sd_current_config()
{
  return get_sd_help();
}

const tdouble RBRV_entry_RV_Weibull::get_median_current_config()
{
  return lambda*pow(log(2*ONE),ONE/k);
}

const tdouble RBRV_entry_RV_Weibull::get_mode_current_config()
{
  if (k>ONE) return lambda*pow((k-ONE)/k,ONE/k);
  return ZERO;
}

const tdouble RBRV_entry_RV_Weibull::get_mean_help()
{
  return lambda*flxgamma(ONE+ONE/k)+eps;
}

const tdouble RBRV_entry_RV_Weibull::get_sd_help()
{
  return lambda*sqrt(flxgamma(ONE+2/k)-pow2(flxgamma(ONE+ONE/k)));
}

const tdouble RBRV_entry_RV_Weibull::get_cov_help()
{
  if (fabs(eps)>GlobalVar.TOL()) {
    return get_sd_help()/get_mean_help();
  }
  const tdouble log_m = 2*GammaLn(ONE+ONE/k);
  const tdouble log_s1 = GammaLn(ONE+2/k);
  const tdouble log_s2 = 2*GammaLn(ONE+ONE/k);
  const tdouble log_f = max(log_s1, log_s2);
  const tdouble r1 = exp((log_f - log_m)/2);
  const tdouble r2 = sqrt(exp(log_s1-log_f) - exp(log_s2-log_f));
  return r1*r2;
}

const bool RBRV_entry_RV_Weibull::check_x(const tdouble xV)
{
  return (xV>eps);
}

const bool RBRV_entry_RV_Weibull::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (epsilon?(epsilon->search_circref(fcr)):false);
}

py::dict RBRV_entry_RV_Weibull::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["k"] = k;
  res["lambda"] = lambda;
  res["epsilon"] = eps;
  res["mean"] = get_mean_current_config();
  res["sd"] = get_sd_current_config();
  res["entropy"] = calc_entropy();
  return res;
}

const tdouble RBRV_entry_RV_Weibull::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(eps) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Weibull::calc_pdf_x", ssV.str() );
  }
  const tdouble d = (x_val-eps)/lambda;
  return (k/lambda)*pow(d,k-ONE)*exp(-pow(d,k));
}

const tdouble RBRV_entry_RV_Weibull::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return log(ZERO);
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(eps) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Weibull::calc_pdf_x_log", ssV.str() );
  }
  const tdouble d = (x_val-eps)/lambda;
  return log(k)-log(lambda)+(k-ONE)*log(d)-pow(d,k);
}

const tdouble RBRV_entry_RV_Weibull::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(eps) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Weibull::calc_cdf_x", ssV.str() );
  }
  const tdouble x = pow((x_val-eps)/lambda,k);
  if (x<=1e-7) {
    return x;
  } else {
    return ONE-exp(-x);
  };
}

const tdouble RBRV_entry_RV_Weibull::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=eps) {
    if (safeCalc) return ONE;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(eps) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Weibull::calc_sf_x", ssV.str() );
  }
  const tdouble x = pow((x_val-eps)/lambda,k);
  if (x<=1e-7) {
    return ONE-x;
  } else {
    return exp(-x);
  };
}

const tdouble RBRV_entry_RV_Weibull::calc_entropy()
{
  return GAMMA*(ONE-ONE/k)+log(lambda/k)+ONE;
}

const tdouble RBRV_entry_RV_Weibull::transform_y2x(const tdouble y_val)
{
  return pow(-log(rv_Phi(-y_val)),ONE/k)*lambda+eps;
}

const tdouble RBRV_entry_RV_Weibull::transform_x2y(const tdouble& x_val)
{
  return -rv_InvPhi_noAlert(exp(-pow((x_val-eps)/lambda,k)));
}

RBRV_entry_RV_ChiSquared::RBRV_entry_RV_ChiSquared(const std::string& name, const tuint iID, FlxFunction* p1, const bool eval_once)
: RBRV_entry_RV_base(name,iID), p1(p1), eval_once(eval_once)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_ChiSquared::RBRV_entry_RV_ChiSquared",1);
    if (p1) delete p1;
    throw;
  }
}

RBRV_entry_RV_ChiSquared::~RBRV_entry_RV_ChiSquared()
{
  if (p1) delete p1;
}

void RBRV_entry_RV_ChiSquared::eval_para()
{
  if (!eval_once || (eval_once&&p1) ) {
    // compute parameters
      dof = p1->cast2positive();
    if (eval_once) {
      delete p1; p1 = NULL;
    }
  }
}

const bool RBRV_entry_RV_ChiSquared::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false);
}

const bool RBRV_entry_RV_ChiSquared::check_x(const tdouble xV)
{
  return (xV>=ZERO);
}

const tdouble RBRV_entry_RV_ChiSquared::get_mean_current_config()
{
  return dof;
}

const tdouble RBRV_entry_RV_ChiSquared::get_sd_current_config()
{
  return sqrt(2*dof);
}

const tdouble RBRV_entry_RV_ChiSquared::get_median_current_config()
{
  return transform_y2x(ZERO);
}

const tdouble RBRV_entry_RV_ChiSquared::get_mode_current_config()
{
  if (dof-2>ZERO) return dof-2;
  return ZERO;
}

const tdouble RBRV_entry_RV_ChiSquared::calc_entropy()
{
  const tdouble dofh = dof/2;
  return dofh 
        + std::log(2*flxgamma(dofh)) 
        + (ONE-dofh)*flxdigamma(dofh);
}

const tdouble RBRV_entry_RV_ChiSquared::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=ZERO) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(ZERO) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_ChiSquared::calc_pdf_x", ssV.str() );
  }
  const tdouble dofh = dof/2;
  return pow(x_val,dofh-ONE)*exp(-x_val/2)/(pow(tdouble(2),dofh)*flxgamma(dofh));
}

const tdouble RBRV_entry_RV_ChiSquared::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=ZERO) {
    if (safeCalc) return log(ZERO);
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(ZERO) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_ChiSquared::calc_pdf_x_log", ssV.str() );
  }
  const tdouble dofh = dof/2;
  return (dofh-ONE)*log(x_val)-x_val/2-dofh*log(2)-GammaLn(dofh);
}

const tdouble RBRV_entry_RV_ChiSquared::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<ZERO) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(ZERO) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_ChiSquared::calc_cdf_x", ssV.str() );
  }
  return flxgamma_rl(dof/2,x_val/2);
}

const tdouble RBRV_entry_RV_ChiSquared::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<ZERO) {
    if (safeCalc) return ONE;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(ZERO) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_ChiSquared::calc_cdf_x", ssV.str() );
  }
  return flxgamma_ru(dof/2,x_val/2);
}

const tdouble RBRV_entry_RV_ChiSquared::transform_y2x(const tdouble y_val)
{
  const tdouble dofh = dof/2;
  if (y_val<=ZERO) {
    return flxgamma_rl_inv(dofh,rv_Phi(y_val))*2;
  } else {
    return flxgamma_ru_inv(dofh,rv_Phi(-y_val))*2;
  }
}

const tdouble RBRV_entry_RV_ChiSquared::transform_x2y(const tdouble& x_val)
{
  if (x_val <= dof) {
    return rv_InvPhi_noAlert(flxgamma_rl(dof/2,x_val/2));
  } else {
    return -rv_InvPhi_noAlert(flxgamma_ru(dof/2,x_val/2));
  }
}

RBRV_entry_RV_Chi::RBRV_entry_RV_Chi(const std::string& name, const tuint iID, FlxFunction* p1, const bool eval_once)
: RBRV_entry_RV_base(name,iID), p1(p1), eval_once(eval_once)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_Chi::RBRV_entry_RV_Chi",1);
    if (p1) delete p1;
    throw;
  }
}

RBRV_entry_RV_Chi::~RBRV_entry_RV_Chi()
{
  if (p1) delete p1;
}

void RBRV_entry_RV_Chi::eval_para()
{
  if (!eval_once || (eval_once&&p1) ) {
    // compute parameters
      dof = p1->cast2positive();
    if (eval_once) {
      delete p1; p1 = NULL;
    }
  }
}

const tdouble RBRV_entry_RV_Chi::transform_y2x(const tdouble y_val)
{
  const tdouble dofh = dof/2;
  if (y_val<=ZERO) {
    return sqrt(flxgamma_rl_inv(dofh,rv_Phi(y_val))*2);
  } else {
    return sqrt(flxgamma_ru_inv(dofh,rv_Phi(-y_val))*2);
  }
}

const tdouble RBRV_entry_RV_Chi::transform_x2y(const tdouble& x_val)
{
  if (x_val <= dof) {
    return rv_InvPhi_noAlert(flxgamma_rl(dof/2,pow2(x_val)/2));
  } else {
    return -rv_InvPhi_noAlert(flxgamma_ru(dof/2,pow2(x_val)/2));
  }
}

const tdouble RBRV_entry_RV_Chi::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=ZERO) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(ZERO) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Chi::calc_pdf_x", ssV.str() );
  }
  const tdouble dofh = dof/2;
  return pow(2*ONE,ONE-dofh)*pow(x_val,dof-ONE)*exp(-pow2(x_val)/2)/flxgamma(dofh);
}

const tdouble RBRV_entry_RV_Chi::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=ZERO) {
    if (safeCalc) return log(ZERO);
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(ZERO) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Chi::calc_pdf_x_log", ssV.str() );
  }
  const tdouble dofh = dof/2;
  return (ONE-dofh)*log(2)+(dof-ONE)*log(x_val)-pow2(x_val)/2-GammaLn(dofh);
}

const tdouble RBRV_entry_RV_Chi::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<ZERO) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(ZERO) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Chi::calc_cdf_x", ssV.str() );
  }
  return flxgamma_rl(dof/2,pow2(x_val)/2);
}

const tdouble RBRV_entry_RV_Chi::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<ZERO) {
    if (safeCalc) return ONE;
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller than (" << GlobalVar.Double2String(ZERO) << " is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_Chi::calc_sf_x", ssV.str() );
  }
  return flxgamma_ru(dof/2,pow2(x_val)/2);
}

const tdouble RBRV_entry_RV_Chi::calc_entropy()
{
  const tdouble dofh = dof/2;
  return GammaLn(dofh)
         + (dof-log(2)-(dof-ONE)*flxdigamma(dofh))/2;
}

const tdouble RBRV_entry_RV_Chi::get_mean_current_config()
{
  return sqrt(2*ONE)*flxgamma((dof+ONE)/2)/flxgamma(dof/2);
}

const tdouble RBRV_entry_RV_Chi::get_sd_current_config()
{
  const tdouble mu = get_mean_current_config();
  return sqrt(dof-pow2(mu));
}

const tdouble RBRV_entry_RV_Chi::get_median_current_config()
{
  return transform_y2x(ZERO);
}

const tdouble RBRV_entry_RV_Chi::get_mode_current_config()
{
  if (dof>=ONE) return sqrt(dof-ONE);
  throw FlxException_NotImplemented("RBRV_entry_RV_Chi::get_mode_current_config");
}

const bool RBRV_entry_RV_Chi::check_x(const tdouble xV)
{
  return (xV>=ZERO);
}

const bool RBRV_entry_RV_Chi::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false);
}


RBRV_entry_RV_StudentsT::RBRV_entry_RV_StudentsT(const std::string& name, const tuint iID, FlxFunction* p1, const bool eval_once)
: RBRV_entry_RV_base(name,iID), p1(p1), eval_once(eval_once)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_StudentsT::RBRV_entry_RV_StudentsT_100",1);
    if (p1) delete p1;
    throw;
  }
}

RBRV_entry_RV_StudentsT::RBRV_entry_RV_StudentsT(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), p1(nullptr), eval_once(false)
{
  try {
    p1 = parse_py_para("dof", config);

    eval_once = parse_py_para_as_bool("eval_once", config, false, false);
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_StudentsT::RBRV_entry_RV_StudentsT_99",1);
    if (p1) delete p1;
    throw;
  }
}

RBRV_entry_RV_StudentsT::~RBRV_entry_RV_StudentsT()
{
  if (p1) delete p1;
}

void RBRV_entry_RV_StudentsT::eval_para()
{
  if (!eval_once || (eval_once&&p1) ) {
    // compute parameters
      dof = p1->cast2positive();
    if (eval_once) {
      delete p1; p1 = NULL;
    }
  }
}

const tdouble RBRV_entry_RV_StudentsT::transform_y2x(const tdouble y_val)
{
  if (y_val<=ZERO) {
    return rv_InvCDF_Studentst(dof,rv_Phi(y_val));
  } else {
    return -rv_InvCDF_Studentst(dof,rv_Phi(-y_val));
  }
}

const tdouble RBRV_entry_RV_StudentsT::transform_x2y(const tdouble& x_val)
{
  if (x_val<=ZERO) {
    return rv_InvPhi_noAlert(rv_cdf_Studentst(dof,x_val));
  } else {
    return -rv_InvPhi_noAlert(rv_cdf_Studentst(dof,-x_val));
  }
}

const tdouble RBRV_entry_RV_StudentsT::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  return rv_pdf_Studentst(dof,x_val);
}

const tdouble RBRV_entry_RV_StudentsT::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return rv_cdf_Studentst(dof,x_val);
}

const tdouble RBRV_entry_RV_StudentsT::calc_entropy()
{
  const tdouble h1 = (dof+1)/2;
  return h1*(flxdigamma(h1)-flxdigamma(dof/2)) + log(sqrt(dof)*BetaFun(dof/2,ONE/2));
}

const tdouble RBRV_entry_RV_StudentsT::get_sd_current_config()
{
  if (dof<2) return std::numeric_limits<tdouble>::infinity();
  else return sqrt(dof/(dof-2));
}

const bool RBRV_entry_RV_StudentsT::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false);
}

const tdouble RBRV_entry_RV_StudentsT::get_HPD(const tdouble p)
{
  return (ONE-p)/2;
}

py::dict RBRV_entry_RV_StudentsT::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["dof"] = dof;
  res["mean"] = get_mean_current_config();
  res["sd"] = get_sd_current_config();
  return res;
}


RBRV_entry_RV_StudentsT_generalized::RBRV_entry_RV_StudentsT_generalized(const std::string& name, const tuint iID, FlxFunction* nu, FlxFunction* locf, FlxFunction* scalef)
:RBRV_entry_RV_base(name,iID), pid(0), p1(nu),p2(locf),p3(scalef), dof(ZERO), loc(ZERO), scale(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_StudentsT_generalized::RBRV_entry_RV_StudentsT_generalized_100",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (p3) delete p3;
    if (p4) delete p4;
    throw;
  }
}

RBRV_entry_RV_StudentsT_generalized::RBRV_entry_RV_StudentsT_generalized(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), pid(0), p1(nullptr), p2(nullptr), p3(nullptr), p4(nullptr), dof(ZERO), loc(ZERO), scale(ZERO)
{
  try {
    if (config.contains("scale")) {          // dof, loc, scale
      pid = 0;
      p1 = parse_py_para("dof", config);
      p2 = parse_py_para("loc", config);
      p3 = parse_py_para("scale", config);
    } else if (config.contains("pr_1")) {   // dof, loc, val_1, pr_1
      pid = 1;
      p1 = parse_py_para("dof", config);
      p2 = parse_py_para("loc", config);
      p3 = parse_py_para("val_1", config);
      p4 = parse_py_para("pr_1", config);
    } else {
      throw FlxException_NeglectInInteractive("RBRV_entry_RV_StudentsT_generalized::RBRV_entry_RV_StudentsT_generalized_01", "Required parameters to define distribution not found in Python <dict>.");
    }
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_StudentsT_generalized::RBRV_entry_RV_StudentsT_generalized_99",1);
    if (p1) delete p1;
    if (p2) delete p2;
    if (p3) delete p3;
    if (p4) delete p4;
    throw;
  }
}

RBRV_entry_RV_StudentsT_generalized::~RBRV_entry_RV_StudentsT_generalized()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (p3) delete p3;
  if (p4) delete p4;
}

tdouble RV_StudentsT_generalized_pid1_root_search_fun(const tdouble scale, void* dp) {
  // convert data-pointer
    tdouble *dp_ = (tdouble *)dp;
    const tdouble dof = dp_[0];
    const tdouble loc = dp_[1];
    const tdouble val_1 = dp_[2];
    const tdouble pr_1 = dp_[3];
  // evalute function expression for root
  const tdouble x_ = (val_1-loc)/scale;
  return (rv_cdf_Studentst(dof,x_) - pr_1)/pr_1;
}

void RBRV_entry_RV_StudentsT_generalized::eval_para()
{
  switch (pid)
  {
    case 0:
      dof = p1->cast2positive();
      loc = p2->calc();
      scale = p3->cast2positive();
      return;
    case 1:
    {
      // evaluate and prepare parameter values
        dof = p1->cast2positive();
        loc = p2->calc();
        const tdouble val_1 = p3->calc();
        const tdouble pr_1 = p4->cast2positive();
        tdouble dp[] = {dof, loc, val_1, pr_1};
        if (pr_1>=ONE) {
          throw FlxException("RBRV_entry_RV_StudentsT_generalized::get_pars_10", "A value larger or equal than one is not allowed.");
        }
      // find initial guess (for root search)
        tdouble mu_ = ZERO;
        tdouble sd_ = ZERO;
        RBRV_entry_RV_normal::get_para_from_quantile(mu_,sd_,val_1,pr_1,loc,0.5);
        const tdouble start_ub = sd_;
        const tdouble start_lb = 0.1*sd_;
      // perform root search
        std::ostringstream ssV;
        scale = flx_RootSearch_RegulaFalsi(&RV_StudentsT_generalized_pid1_root_search_fun,dp,start_lb,start_ub);  // ,1e-6,1e-8,&ssV);  // TODO remove cout after debug
        // GlobalVar.slogcout(1) << ssV.str() << std::endl << std::flush;
      return;
    }
    default:
      throw FlxException_Crude("RBRV_entry_RV_StudentsT_generalized::get_pars_99");
  }
}

const tdouble RBRV_entry_RV_StudentsT_generalized::transform_y2x(const tdouble y_val)
{
  if (y_val<=ZERO) {
    return loc + rv_InvCDF_Studentst(dof,rv_Phi(y_val))*scale;
  } else {
    return loc - rv_InvCDF_Studentst(dof,rv_Phi(-y_val))*scale;
  }
}

const tdouble RBRV_entry_RV_StudentsT_generalized::transform_x2y(const tdouble& x_val)
{
  const tdouble x_ = (x_val-loc)/scale;
  if (x_<=ZERO) {
    return rv_InvPhi_noAlert(rv_cdf_Studentst(dof,x_));
  } else {
    return -rv_InvPhi_noAlert(rv_cdf_Studentst(dof,-x_));
  }
}

const tdouble RBRV_entry_RV_StudentsT_generalized::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  const tdouble x_ = (x_val-loc)/scale;
  return rv_pdf_Studentst(dof,x_)/scale;
}

const tdouble RBRV_entry_RV_StudentsT_generalized::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  const tdouble x_ = (x_val-loc)/scale;
  return rv_cdf_Studentst(dof,x_);
}

const tdouble RBRV_entry_RV_StudentsT_generalized::calc_entropy()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_StudentsT_generalized::calc_entropy");
  // in principle, this should be not too difficult to implement, if integration by substitution is used.
}

const tdouble RBRV_entry_RV_StudentsT_generalized::get_mean_current_config()
{
  return loc;
}

const tdouble RBRV_entry_RV_StudentsT_generalized::get_sd_current_config()
{
  if (dof<2) return std::numeric_limits<tdouble>::infinity();
  else return sqrt(dof/(dof-2))*scale;
}

const tdouble RBRV_entry_RV_StudentsT_generalized::get_median_current_config()
{
  return loc;
}

const tdouble RBRV_entry_RV_StudentsT_generalized::get_mode_current_config()
{
  return loc;
}

const bool RBRV_entry_RV_StudentsT_generalized::search_circref(FlxFunction* fcr)
{
  bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  b = b || p1->search_circref(fcr) || p2->search_circref(fcr) || p3->search_circref(fcr);
  if (p4) {
    b = b || p4->search_circref(fcr);
  }
  return b;
}

const tdouble RBRV_entry_RV_StudentsT_generalized::get_HPD(const tdouble p)
{
  return (ONE-p)/2;
}

py::dict RBRV_entry_RV_StudentsT_generalized::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["dof"] = dof;
  res["loc"] = loc;
  res["scale"] = scale;
  res["sd"] = get_sd_current_config();
  return res;
}



RBRV_entry_RV_logt::RBRV_entry_RV_logt(const std::string& name, const tuint iID, FlxFunction* nu, FlxFunction* locf, FlxFunction* scalef)
:RBRV_entry_RV_StudentsT_generalized(name, iID, nu, locf, scalef)
{
  value = exp(value);
}

RBRV_entry_RV_logt::RBRV_entry_RV_logt(const std::string& name, const tuint iID, py::dict config)
:RBRV_entry_RV_StudentsT_generalized(name, iID, config)
{
  value = exp(value);
}

RBRV_entry_RV_logt::~RBRV_entry_RV_logt()
{

}

const tdouble RBRV_entry_RV_logt::transform_y2x(const tdouble y_val)
{
  const tdouble z = RBRV_entry_RV_StudentsT_generalized::transform_y2x(y_val);
  return exp(z);
}

const tdouble RBRV_entry_RV_logt::transform_x2y(const tdouble& x_val)
{
  const tdouble z = log(x_val);
  return RBRV_entry_RV_StudentsT_generalized::transform_x2y(z);
}

const tdouble RBRV_entry_RV_logt::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  return RBRV_entry_RV_StudentsT_generalized::calc_pdf_x(log(x_val),safeCalc)/x_val;
}

const tdouble RBRV_entry_RV_logt::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return RBRV_entry_RV_StudentsT_generalized::calc_cdf_x(log(x_val),safeCalc);
}

const tdouble RBRV_entry_RV_logt::calc_entropy()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_logt::calc_entropy");
  // in principle, this should be not too difficult to implement, if integration by substitution is used.
}

const tdouble RBRV_entry_RV_logt::get_mean_current_config()
{
  return std::numeric_limits<tdouble>::infinity();
}

const tdouble RBRV_entry_RV_logt::get_sd_current_config()
{
  return std::numeric_limits<tdouble>::infinity();
}

const tdouble RBRV_entry_RV_logt::get_median_current_config()
{
  return exp(RBRV_entry_RV_StudentsT_generalized::get_median_current_config());
}

const tdouble RBRV_entry_RV_logt::get_mode_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_logt::get_mode_current_config");
}

const tdouble RBRV_entry_RV_logt::get_HPD(const tdouble p)
{
  throw FlxException_NotImplemented("RBRV_entry_RV_logt::get_HPD");
}

py::dict RBRV_entry_RV_logt::info()
{
  return RBRV_entry_RV_StudentsT_generalized::info();
}


RBRV_entry_RV_Laplace::RBRV_entry_RV_Laplace(const std::string& name, const tuint iID, FlxFunction* locf, FlxFunction* scalef)
:RBRV_entry_RV_base(name,iID), locf(locf),scalef(scalef), loc(ZERO), scale(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_Laplace::RBRV_entry_RV_Laplace",1);
    if (locf) delete locf;
    if (scalef) delete scalef;
    throw;
  }
}

RBRV_entry_RV_Laplace::~RBRV_entry_RV_Laplace()
{
  if (locf) delete locf;
  if (scalef) delete scalef;
}

void RBRV_entry_RV_Laplace::eval_para()
{
  loc = locf->calc();
  scale = scalef->cast2positive();
}

const tdouble RBRV_entry_RV_Laplace::transform_y2x(const tdouble y_val)
{
  if (y_val<=ZERO) {
    return loc + log(2*rv_Phi(y_val))*scale;
  } else {
    return loc - log(2*rv_Phi(-y_val))*scale;
  }
}

const tdouble RBRV_entry_RV_Laplace::transform_x2y(const tdouble& x_val)
{
  const tdouble x_ = (x_val-loc)/scale;
  if (x_<=ZERO) {
    const tdouble z = exp(x_)/2;
    return rv_InvPhi_noAlert(z);
  } else {
    return -rv_InvPhi_noAlert(exp(-x_)/2);
  }
}

const tdouble RBRV_entry_RV_Laplace::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  const tdouble x_ = fabs(x_val-loc)/scale;
  return exp(-x_)/(2*scale);
}

const tdouble RBRV_entry_RV_Laplace::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  const tdouble x_ = fabs(x_val-loc)/scale;
  return -x_ - log(2*scale);
}

const tdouble RBRV_entry_RV_Laplace::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=loc) {
    return exp((x_val-loc)/scale)/2;
  } else {
    return ONE - exp((loc-x_val)/scale)/2;
  }
}

const tdouble RBRV_entry_RV_Laplace::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<=loc) {
    return ONE - exp((x_val-loc)/scale)/2;
  } else {
    return exp((loc-x_val)/scale)/2;
  }
}

const tdouble RBRV_entry_RV_Laplace::calc_entropy()
{
  return ONE+log(2*scale);
}

const tdouble RBRV_entry_RV_Laplace::get_mean_current_config()
{
  return loc;
}

const tdouble RBRV_entry_RV_Laplace::get_sd_current_config()
{
  return sqrt(2.)*scale;
}

const tdouble RBRV_entry_RV_Laplace::get_median_current_config()
{
  return loc;
}

const tdouble RBRV_entry_RV_Laplace::get_mode_current_config()
{
  return loc;
}

const bool RBRV_entry_RV_Laplace::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || locf->search_circref(fcr) || scalef->search_circref(fcr);
}

const tdouble RBRV_entry_RV_Laplace::get_HPD(const tdouble p)
{
  return (ONE-p)/2;
}


RBRV_entry_RV_genpareto::RBRV_entry_RV_genpareto(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), xif(nullptr), locf(nullptr), scalef(nullptr), xi(ZERO), loc(ZERO), scale(ZERO), eval_once(false)
{
  try {
    xif = parse_py_para("xi",config,true);
    locf = parse_py_para("loc",config,false);
    scalef = parse_py_para("scale",config,false);

    eval_once = parse_py_para_as_bool("eval_once", config, false, false);
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_genpareto::RBRV_entry_RV_genpareto_99",1);
    free_mem();
    throw;
  }
}

RBRV_entry_RV_genpareto::~RBRV_entry_RV_genpareto()
{
  free_mem();
}

void RBRV_entry_RV_genpareto::free_mem()
{
  if (xif) delete xif;
  if (locf) delete locf;
  if (scalef) delete scalef;
}

void RBRV_entry_RV_genpareto::eval_para()
{
  if (!eval_once || (eval_once&&xif) ) {
    // compute parameters
      xi = xif->calc();
      if (scalef) {
        scale = scalef->cast2positive();
      } else {
        scale = ONE;
      }
      if (locf) {
        loc = locf->cast2positive();
      } else {
        loc = ZERO;
      }
    if (eval_once) {
      delete xif; xif = nullptr;
      if (scalef) {
        delete scalef;
        scalef = nullptr;
      }
      if (locf) {
        delete locf;
        locf = nullptr;
      }
    }
  }
}

const tdouble RBRV_entry_RV_genpareto::transform_y2x(const tdouble y_val)
{
  const tdouble p_ = rv_Phi(-y_val); // NOTE this is p_=1-p !!!
  if (xi==ZERO) {
    return loc-scale*log(p_);
  } else {
    return loc+scale/xi*(pow(p_,-xi)-ONE);
  }
}

const tdouble RBRV_entry_RV_genpareto::transform_x2y(const tdouble& x_val)
{
  if (x_val<=get_median_current_config()) {
    return rv_InvPhi_noAlert(ONE-eval_cdf_help(x_val));
  } else {
    return -rv_InvPhi_noAlert(eval_cdf_help(x_val));
  }
}

const tdouble RBRV_entry_RV_genpareto::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  const tdouble x_ = (x_val-loc)/scale;
  return pow(ONE+xi*x_,-(ONE/xi+ONE))/scale;
}

const tdouble RBRV_entry_RV_genpareto::eval_cdf_help(const tdouble x_val)
{
  const tdouble x_ = (x_val-loc)/scale;
  if (xi==ZERO) {
    return exp(-x_);
  } else {
    return pow(ONE+xi*x_,-ONE/xi);
  }
}

const tdouble RBRV_entry_RV_genpareto::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return ONE-eval_cdf_help(x_val);
}

const tdouble RBRV_entry_RV_genpareto::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  return eval_cdf_help(x_val);
}

const tdouble RBRV_entry_RV_genpareto::calc_entropy()
{
  return log(scale)+xi+ONE;
}

const tdouble RBRV_entry_RV_genpareto::get_mean_current_config()
{
  if (xi>=ONE) return std::numeric_limits<tdouble>::infinity();
  return loc+scale/(ONE-xi);
}

const tdouble RBRV_entry_RV_genpareto::get_sd_current_config()
{
  if (xi>=ONE/2) return std::numeric_limits<tdouble>::infinity();
  return scale/((ONE-xi)*sqrt(ONE-2*xi));
}

const tdouble RBRV_entry_RV_genpareto::get_median_current_config()
{
  return loc + scale*(pow(2*ONE,xi)-ONE)/xi;
}

const tdouble RBRV_entry_RV_genpareto::get_mode_current_config()
{
  return loc;
}

const bool RBRV_entry_RV_genpareto::check_x(const tdouble xV)
{
  if (xi>=ZERO) {
    return xV>=loc;
  } else {
    return (xV>=loc) && (loc-scale/xi);
  }
}

const bool RBRV_entry_RV_genpareto::search_circref(FlxFunction* fcr)
{
  bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  b = b || xif->search_circref(fcr);
  if (locf) {
    b = b || locf->search_circref(fcr);
  }
  if (scalef) {
    b = b || scalef->search_circref(fcr);
  }
  return b;
}

py::dict RBRV_entry_RV_genpareto::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["xi"] = xi;
  res["loc"] = loc;
  res["scale"] = scale;
  return res;
}


RBRV_entry_RV_UserTransform::RBRV_entry_RV_UserTransform(const std::string& name, const tuint iID, const bool is_z2x, FlxFunction* t1, FlxFunction* t2, FlxFunction* dh, FlxFunction* checkXf, RBRV_entry_RV_base* rv_z, const bool manage_z)
: RBRV_entry_RV_base(name,iID), is_z2x(is_z2x), t1(t1), t2(t2), dh(dh), checkXf(checkXf), rv_z(rv_z), manage_z(manage_z),
  tPL(1), tPLp(&tPL[0])
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_UserTransform::RBRV_entry_RV_UserTransform",1);
    if (t1) delete t1;
    if (t2) delete t2;
    if (dh) delete dh;
    if (checkXf) delete checkXf;
    if (manage_z) delete rv_z;
    throw;
  }
}

RBRV_entry_RV_UserTransform::~RBRV_entry_RV_UserTransform()
{
  if (t1) delete t1;
  if (t2) delete t2;
  if (dh) delete dh;
  if (checkXf) delete checkXf;
  if (manage_z) delete rv_z;
}

const tdouble RBRV_entry_RV_UserTransform::eval_para_fun(FlxFunction* fp, const tdouble value)
{
  const tdouble* const tT = FunPara::ParaList;
  const tuint tTS = FunPara::ParaListSize;
  FunPara::ParaList = tPLp;
  FunPara::ParaListSize = 1;
  *tPLp = value;
  const tdouble res = fp->calc();
  FunPara::ParaList = tT;
  FunPara::ParaListSize = tTS;
  return res;
}

const tdouble RBRV_entry_RV_UserTransform::transform_x2y(const tdouble& x_val)
{
  if (t2==NULL) {
    std::ostringstream ssV;
    ssV << "Parameter '" << (is_z2x?"x2z":"z2y") << "' not set.";
    throw FlxException("RBRV_entry_RV_UserTransform::transform_x2y", ssV.str() );
  }
  if (is_z2x) {
    const tdouble z = eval_para_fun(t2,x_val);
    return rv_z->transform_x2y(z);
  } else {
    const tdouble z = rv_z->transform_x2y(x_val);
    return eval_para_fun(t2,z);
  }
}

const tdouble RBRV_entry_RV_UserTransform::transform_y2x(const tdouble y_val)
{
  if (is_z2x) {
    const tdouble z = rv_z->transform_y2x(y_val);
    return eval_para_fun(t1,z);
  } else {
    const tdouble z = eval_para_fun(t1,y_val);
    return rv_z->transform_y2x(z);
  }
}

const tdouble RBRV_entry_RV_UserTransform::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  // perform some checks
    if (dh==NULL) {
      std::ostringstream ssV;
      ssV << "Parameter '" << (is_z2x?"dx2z":"dz2y") << "' not set.";
      throw FlxException("RBRV_entry_RV_UserTransform::calc_pdf_x_01", ssV.str() );
    }
    if (t2==NULL) {
      std::ostringstream ssV;
      ssV << "Parameter '" << (is_z2x?"x2z":"z2y") << "' not set.";
      throw FlxException("RBRV_entry_RV_UserTransform::calc_pdf_x_02", ssV.str() );
    }
    if (check_x(x_val)==false) {
      if (safeCalc) return ZERO;        // TODO: critical, as it could also be ONE, depending on the value of 'x_val'
      std::ostringstream ssV;
      ssV << "The value (" << GlobalVar.Double2String(x_val) << ") is outside of the support of this random variable.";
      throw FlxException("RBRV_entry_RV_UserTransform::calc_pdf_x_03", ssV.str() );
    }
  if (is_z2x) {
    const tdouble z = eval_para_fun(t2,x_val);
    const tdouble dzdx = eval_para_fun(dh,x_val);
    return rv_z->calc_pdf_x(z) * fabs(dzdx);
  } else {
    const tdouble z = rv_z->transform_x2y(x_val);
    const tdouble y = eval_para_fun(t2,z);
    const tdouble dydz = eval_para_fun(dh,z);
    return rv_phi(y) * fabs(dydz) * rv_z->calc_pdf_x(x_val)/rv_phi(z);
  }
}

const tdouble RBRV_entry_RV_UserTransform::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  // perform some checks
    if (dh==NULL) {
      std::ostringstream ssV;
      ssV << "Parameter '" << (is_z2x?"dx2z":"dz2y") << "' not set.";
      throw FlxException("RBRV_entry_RV_UserTransform::calc_pdf_x_log_01", ssV.str() );
    }
    if (t2==NULL) {
      std::ostringstream ssV;
      ssV << "Parameter '" << (is_z2x?"x2z":"z2y") << "' not set.";
      throw FlxException("RBRV_entry_RV_UserTransform::calc_pdf_x_log_02", ssV.str() );
    }
    if (check_x(x_val)==false) {
      if (safeCalc) return log(ZERO);        // TODO: critical, as it could also be ONE, depending on the value of 'x_val'
      std::ostringstream ssV;
      ssV << "The value (" << GlobalVar.Double2String(x_val) << ") is outside of the support of this random variable.";
      throw FlxException("RBRV_entry_RV_UserTransform::calc_pdf_x_log_03", ssV.str() );
    }
  if (is_z2x) {
    const tdouble z = eval_para_fun(t2,x_val);
    const tdouble dzdx = eval_para_fun(dh,x_val);
    return rv_z->calc_pdf_x_log(z) + log(fabs(dzdx));
  } else {
    const tdouble z = rv_z->transform_x2y(x_val);
    const tdouble y = eval_para_fun(t2,z);
    const tdouble dydz = eval_para_fun(dh,z);
    return rv_phi_log(y) + log(fabs(dydz)) + rv_z->calc_pdf_x_log(x_val) - rv_phi_log(z);
  }
}

const tdouble RBRV_entry_RV_UserTransform::eval_cdf_sf(const bool is_cdf, const tdouble& x_val, const bool safeCalc)
{
  // perform some checks
    if (t2==NULL) {
      std::ostringstream ssV;
      ssV << "Parameter '" << (is_z2x?"x2z":"z2y") << "' not set.";
      throw FlxException("RBRV_entry_RV_UserTransform::calc_cdf_x_01", ssV.str() );
    }
    if (check_x(x_val)==false) {
      if (safeCalc) return ZERO;
      std::ostringstream ssV;
      ssV << "The value (" << GlobalVar.Double2String(x_val) << ") is outside of the support of this random variable.";
      throw FlxException("RBRV_entry_RV_UserTransform::calc_cdf_x_01", ssV.str() );
    }
  if (is_z2x) {
    const tdouble z = eval_para_fun(t2,x_val);
    if (is_cdf) {
      return rv_z->calc_cdf_x(z,safeCalc);
    } else {
      return rv_z->calc_sf_x(z,safeCalc);
    }
  } else {
    const tdouble y = transform_x2y(x_val);
    if (is_cdf) {
      return rv_Phi(y);
    } else {
      return rv_Phi(-y);
    }
  }
}

const tdouble RBRV_entry_RV_UserTransform::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return eval_cdf_sf(true,x_val,safeCalc);
}

const tdouble RBRV_entry_RV_UserTransform::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  return eval_cdf_sf(false,x_val,safeCalc);
}

const tdouble RBRV_entry_RV_UserTransform::get_mean_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_UserTransform::get_mean_current_config");
}

const tdouble RBRV_entry_RV_UserTransform::get_sd_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_UserTransform::get_sd_current_config");
}

const bool RBRV_entry_RV_UserTransform::check_x(const tdouble xV)
{
  if (is_z2x) {
    return checkXf?(eval_para_fun(checkXf,xV)>ZERO):true;
  } else {
    return rv_z->check_x(xV);
  }  
}

const bool RBRV_entry_RV_UserTransform::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (t1?(t1->search_circref(fcr)):false) || (t2?(t2->search_circref(fcr)):false) || rv_z->search_circref(fcr);
}

py::dict RBRV_entry_RV_UserTransform::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  if (is_z2x) {
    res["z2x"] = t1->write();
    res["x2z"] = t2->write();
    res["dx2z"] = dh->write();
    res["sd"] = get_sd_current_config();
    res["rv_z"] = rv_z->info();
  }
  return res;
}

void RBRV_entry_RV_UserTransform::replace_rv_z(RBRV_entry_RV_base* rv_z_)
{
  if (manage_z) throw FlxException_Crude("RBRV_entry_RV_UserTransform::replace_rv_z");
  rv_z = rv_z_;
}

RBRV_entry_RV_Truncated::RBRV_entry_RV_Truncated(const std::string& name, const tuint iID, FlxFunction* a, FlxFunction* b, RBRV_entry_RV_base* rv_z, const bool manage_z)
: RBRV_entry_RV_base(name,iID), a(a), b(b), rv_z(rv_z), manage_z(manage_z), aV(ZERO), bV(ZERO), q(ONE), aV_cdf(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_Truncated::RBRV_entry_RV_Truncated",1);
    if (a) delete a;
    if (b) delete b;
    if (manage_z) delete rv_z;
    throw;
  }
}

RBRV_entry_RV_Truncated::~RBRV_entry_RV_Truncated()
{
  if (a) delete a;
  if (b) delete b;
  if (manage_z) delete rv_z;
}

void RBRV_entry_RV_Truncated::eval_para()
{
  const tdouble ybound = 1e5;
  aV = a?(a->calc()):rv_z->transform_y2x(-ybound);
  bV = b?(b->calc()):rv_z->transform_y2x(ybound);
  aV_cdf = a?(rv_z->calc_cdf_x(aV)):ZERO;
  q = (b?(rv_z->calc_cdf_x(bV)):ONE)-aV_cdf;
  if (aV_cdf>ONE/2 || (q<GlobalVar.TOL()&&aV_cdf>GlobalVar.TOL())) {
    const tdouble y_a = a?(rv_z->transform_x2y(aV)):-ybound;
    const tdouble y_b = b?(rv_z->transform_x2y(bV)):ybound;
    q = rv_Phi_diff(y_a,y_b);
  }
}


const tdouble RBRV_entry_RV_Truncated::transform_y2x(const tdouble y_val)
{
  if (y_val<=ZERO || aV_cdf<0.5) {
    const tdouble p_orig = rv_Phi(y_val);
    const tdouble p_transf = p_orig*q+aV_cdf;
    if (p_transf<0.95) {
      const tdouble y_transf = rv_InvPhi_noAlert(p_transf);
      return rv_z->transform_y2x(y_transf);
    }
  }
  const tdouble p_orig = rv_Phi(-y_val);
  tdouble p_transf = p_orig*q;
  if (b) {
    const tdouble b2y = rv_z->transform_x2y(bV);
    p_transf += rv_Phi(-b2y);
  }
  const tdouble y_transf = -rv_InvPhi_noAlert(p_transf);
  return rv_z->transform_y2x(y_transf);

}

const tdouble RBRV_entry_RV_Truncated::transform_x2y(const tdouble& x_val)
{
  return rv_InvPhi_noAlert(calc_cdf_x(x_val,false));
}

const tdouble RBRV_entry_RV_Truncated::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bV || x_val<aV) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(aV) << ";" << GlobalVar.Double2String(bV) << "].";
    throw FlxException("RBRV_entry_RV_Truncated::calc_pdf_x", ssV.str() );
  }
  if (q==ZERO) return ZERO;
  return rv_z->calc_pdf_x(x_val,safeCalc)/q;
}

const tdouble RBRV_entry_RV_Truncated::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bV || x_val<aV) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(aV) << ";" << GlobalVar.Double2String(bV) << "].";
    throw FlxException("RBRV_entry_RV_Truncated::calc_pdf_x", ssV.str() );
  }
  if (q==ZERO) return log(ZERO);
  return rv_z->calc_pdf_x_log(x_val,safeCalc)-log(q);
}

const tdouble RBRV_entry_RV_Truncated::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val>bV || x_val<aV) {
    if (safeCalc) {
      if (x_val<aV) return ZERO;
      else return ONE;
    }
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(aV) << ";" << GlobalVar.Double2String(bV) << "].";
    throw FlxException("RBRV_entry_RV_normal_trunc::calc_cdf_x", ssV.str() );
  }
  return (rv_z->calc_cdf_x(x_val)-aV_cdf)/q;
}

const tdouble RBRV_entry_RV_Truncated::get_mean_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_Truncated::get_mean_current_config");
}

const tdouble RBRV_entry_RV_Truncated::get_sd_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_Truncated::get_sd_current_config");
}

const bool RBRV_entry_RV_Truncated::check_x(const tdouble xV)
{
  if (xV<=bV&&xV>=aV) {
    return rv_z->check_x(xV);
  } else {
    return false;
  }
}

const bool RBRV_entry_RV_Truncated::search_circref(FlxFunction* fcr)
{
  const bool bb = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return bb || (a?(a->search_circref(fcr)):false) || (b?(b->search_circref(fcr)):false) || rv_z->search_circref(fcr);
}

py::dict RBRV_entry_RV_Truncated::info()
{
  py::dict res = RBRV_entry_RV_base::info();
  res["lower"] = aV;
  res["upper"] = bV;
  res["q"] = q;
  res["aV_cdf"] = aV_cdf;
  res["rv_z"] = rv_z->info();
  return res;
}

RBRV_entry_RV_maxminTransform::RBRV_entry_RV_maxminTransform(const std::string& name, const tuint iID, const bool is_max, FlxFunction* n, RBRV_entry_RV_base* rv_z)
: RBRV_entry_RV_base(name,iID), is_max(is_max), n(n), rv_z(rv_z), nV(ZERO)
{
  try {
    this->init();
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_maxminTransform::RBRV_entry_RV_maxminTransform",1);
    if (n) delete n;
    delete rv_z;
    throw;
  }
}

RBRV_entry_RV_maxminTransform::~RBRV_entry_RV_maxminTransform()
{
  if (n) delete n;
  delete rv_z;
}

void RBRV_entry_RV_maxminTransform::eval_para()
{
  nV = n->cast2positive();
}

const tdouble RBRV_entry_RV_maxminTransform::transform_y2x(const tdouble y_val)
{
  if (is_max) {
    const tdouble px = pow(rv_Phi(y_val),ONE/nV);
    const tdouble yx = rv_InvPhi_noAlert(px);
    return rv_z->transform_y2x(yx);
  } else {
    const tdouble px = pow(rv_Phi(-y_val),ONE/nV);
    const tdouble yx = -rv_InvPhi_noAlert(px);
    return rv_z->transform_y2x(yx);
  }
}

const tdouble RBRV_entry_RV_maxminTransform::transform_x2y(const tdouble& x_val)
{
  const tdouble px = rv_z->calc_cdf_x(x_val);
  if (is_max) {
    const tdouble py = pow(px,nV);
    return rv_InvPhi_noAlert(py);
  } else {
    const tdouble py = pow(ONE-px,nV);
    return -rv_InvPhi_noAlert(py);
  }
}

const tdouble RBRV_entry_RV_maxminTransform::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  const tdouble px = rv_z->calc_cdf_x(x_val);
  const tdouble qx = rv_z->calc_pdf_x(x_val);
  if (is_max) {
    return nV*qx*pow(px,nV-ONE);
  } else {
    return nV*qx*pow(ONE-px,nV-ONE);
  }
}

const tdouble RBRV_entry_RV_maxminTransform::eval_cdf_sf(const bool is_cdf, const tdouble& x_val, const bool safeCalc)
{
  if (is_max) {
    const tdouble py = pow(rv_z->calc_cdf_x(x_val),nV);
    if (is_cdf) {
      return py;
    } else {
      return ONE-py;
    }
  } else {
    const tdouble py = pow(rv_z->calc_sf_x(x_val),nV);
    if (is_cdf) {
      return ONE-py;
    } else {
      return py;
    }
  }
}

const tdouble RBRV_entry_RV_maxminTransform::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  return eval_cdf_sf(true, x_val, safeCalc);
}

const tdouble RBRV_entry_RV_maxminTransform::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  return eval_cdf_sf(false, x_val, safeCalc);
}

const tdouble RBRV_entry_RV_maxminTransform::get_mean_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_maxminTransform::get_mean_current_config");
}

const tdouble RBRV_entry_RV_maxminTransform::get_sd_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_maxminTransform::get_sd_current_config");
}

const bool RBRV_entry_RV_maxminTransform::check_x(const tdouble xV)
{
    return rv_z->check_x(xV);
}

const bool RBRV_entry_RV_maxminTransform::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (n?(n->search_circref(fcr)):false) || rv_z->search_circref(fcr);
}



RBRV_entry_RV_base* extract_tail_rv(py::dict& dict_tail)
{
  const std::string tail_model = parse_py_para_as_string("use_model", dict_tail, true);
  if (tail_model=="None") {
    return nullptr;
  }
  if (dict_tail.contains("models")==false) {
    throw FlxException("extract_tail_rv_01", "Key 'models' not found in tail configuration.");
  }
  py::dict dict_models = parse_py_obj_as_dict(dict_tail["models"], "Key 'models' in 'config'");
  if (dict_models.contains(tail_model.c_str())==false) {
    throw FlxException("extract_tail_rv_02", "Model '" + tail_model + "' not found in tail configuration.");
  }
  py::dict config = parse_py_obj_as_dict(dict_models[tail_model.c_str()], "Model '" + tail_model + "' in 'config'");
  return parse_py_obj_as_rv(config,false,0,"","internal");
}

RBRV_entry_RV_quantiles::RBRV_entry_RV_quantiles(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), N_bins(0), p_vec(nullptr), q_vec(nullptr), bin_rv_params(nullptr), N_vec(nullptr), tail_up(nullptr), tail_low(nullptr),
  interpol_type(interpol_type_t::uniform)
{
  try {
    N_bins = parse_py_para_as_tuintNo0("N_bins",config,true);
    if (N_bins<5) {
      throw FlxException("RBRV_entry_RV_quantiles::RBRV_entry_RV_quantiles_01", "Not enough data.");
    }
    // import p_vec
      pv.resize(N_bins+1);
      p_vec = &pv[0];
      flxVec pv_tmp(p_vec,N_bins+1);
      pv_tmp.assign_save(parse_py_para_as_flxVec("p_vec",config,true));
    // import q_vec
      qv.resize(N_bins+1);
      q_vec = &qv[0];
      flxVec qv_tmp(q_vec,N_bins+1);
      qv_tmp.assign_save(parse_py_para_as_flxVec("q_vec",config,true));
    // tail fit
      const bool use_tail_fit = parse_py_para_as_bool("use_tail_fit",config,false,false);
      if (use_tail_fit) {
        // lower tail
          if (config.contains("tail_lower")) {
            py::dict dict_tail = parse_py_obj_as_dict(config["tail_lower"], "'tail_lower' in 'config'");
            tail_low = extract_tail_rv(dict_tail);
          }
        // upper tail
          if (config.contains("tail_upper")) {
            py::dict dict_tail = parse_py_obj_as_dict(config["tail_upper"], "'tail_upper' in 'config'");
            tail_up = extract_tail_rv(dict_tail);
          }
      }
    // type for interpolation
      if (config.contains("interpol")) {
        const std::string str_interpol = parse_py_para_as_word("interpol",config,true,true);
        if (str_interpol=="uniform") interpol_type = interpol_type_t::uniform;
        else if (str_interpol=="pchip") interpol_type = interpol_type_t::pchip;
        else if (str_interpol=="bin_beta") interpol_type = interpol_type_t::bin_beta;
        else if (str_interpol=="bin_linear") interpol_type = interpol_type_t::bin_linear;
        else if (str_interpol=="pdf_linear") interpol_type = interpol_type_t::pdf_linear;
        else {
          throw FlxException("RBRV_entry_RV_quantiles::RBRV_entry_RV_quantiles_02","Unknown identifier for type of interpolation: " + str_interpol);
        }
      }

    // beta-fit of bins
      if (interpol_type == interpol_type_t::bin_beta) {
        if (config.contains("bin_rvbeta_params")) {
          bin_rv_params = new tdouble[N_bins*2];
          flxVec brvb_tmp(bin_rv_params,N_bins*2);
          brvb_tmp.assign_save(parse_py_para_as_flxVec("bin_rvbeta_params",config,true));
        } else {
          throw FlxException("RBRV_entry_RV_quantiles::RBRV_entry_RV_quantiles_03", "'bin_rvbeta_params' required in config.");
        }
      }
    // linear-fit of bins
      if (interpol_type == interpol_type_t::bin_linear) {
        if (config.contains("bin_rvlinear_params")) {
          bin_rv_params = new tdouble[N_bins];
          flxVec brvb_tmp(bin_rv_params,N_bins);
          brvb_tmp.assign_save(parse_py_para_as_flxVec("bin_rvlinear_params",config,true));
        } else {
          throw FlxException("RBRV_entry_RV_quantiles::RBRV_entry_RV_quantiles_04", "'bin_rvlinear_params' required in config.");
        }
      }
    // linear fit of entire pdf
      if (interpol_type == interpol_type_t::pdf_linear) {
        if (config.contains("pdf_vec")) {
          bin_rv_params = new tdouble[N_bins+1];
          flxVec brvb_tmp(bin_rv_params,N_bins+1);
          brvb_tmp.assign_save(parse_py_para_as_flxVec("pdf_vec",config,true));
        } else {
          throw FlxException("RBRV_entry_RV_quantiles::RBRV_entry_RV_quantiles_05", "'pdf_vec' required in config.");
        }
        // correct q_vec
          // remember tail probabilities
            tdouble Pr_tail_low = ZERO;
            tuint i_start = 0;
            if (tail_low) {
              Pr_tail_low = p_vec[1];
              i_start = 1;
            }
            tdouble Pr_tail_up = ZERO;
            tuint i_end = N_bins;
            if (tail_up) {
              Pr_tail_up = ONE - p_vec[N_bins-1];
              --i_end;
            }
          // reset p_vec
            pv_tmp = ZERO;
          // restore tail probabilities
            p_vec[1] = Pr_tail_low;
          //  iterate over bins
            pdouble Pr_sum = Pr_tail_low;
            for (tuint i=i_start;i<i_end;++i) {
              Pr_sum += (bin_rv_params[i]+bin_rv_params[i+1])/2*(q_vec[i+1]-q_vec[i]);
              p_vec[i+1] = Pr_sum.cast2double();
            }
          // restore tail properties
            if (tail_up) {
              Pr_sum += Pr_tail_up;
              p_vec[N_bins] = Pr_sum.cast2double();
            }
          // try to 'normalize out' round-off errors
            if (fabs(ONE-Pr_sum.cast2double())>1e-10) {
              throw FlxException("RBRV_entry_RV_quantiles::RBRV_entry_RV_quantiles_06", GlobalVar.Double2String(ONE-Pr_sum.cast2double()));
            }
            for (tuint i=1;i<N_bins;++i) {
              p_vec[i] /= Pr_sum.cast2double();
            }
            p_vec[N_bins] = ONE;
      }

    // N_vec = new tuint[N_total]; TODO
    // import pchip
      if (interpol_type == interpol_type_t::pchip) {
        // NOTE not very memory efficient. One could use » std::make_shared<std::vector<tdouble>>(x);
        pchip_cdf.emplace( std::vector<tdouble>(qv), std::vector<tdouble>(pv) );
        pchip_icdf.emplace( std::vector<tdouble>(pv), std::vector<tdouble>(qv) );
      }
  } catch (FlxException& e) {
    free_mem();
    throw;
  }
}


RBRV_entry_RV_quantiles::~RBRV_entry_RV_quantiles()
{
    free_mem();
}

void RBRV_entry_RV_quantiles::free_mem()
{
    if (bin_rv_params) delete [] bin_rv_params;
    if (N_vec) delete [] N_vec;
    if (tail_up) delete tail_up;
    if (tail_low) delete tail_low;
}

void RBRV_entry_RV_quantiles::eval_para()
{

}

const tdouble RBRV_entry_RV_quantiles::transform_y2x(const tdouble y_val)
{
  const tdouble p = rv_Phi(y_val);
  if (p<p_vec[1]) {
    if (tail_low) {
      const tdouble p_transformed = ONE - p/p_vec[1];
      const tdouble x_transformed = tail_low->calc_icdf_x(p_transformed);
      return q_vec[1] - x_transformed;
    }
  } else if (p>p_vec[N_bins-1]) {
    if (tail_up) {
      const tdouble p_transformed = (p-p_vec[N_bins-1])/(ONE-p_vec[N_bins-1]);
      const tdouble x_transformed = tail_up->calc_icdf_x(p_transformed);
      return x_transformed + q_vec[N_bins-1];
    }
  }
  switch (interpol_type) {
    case interpol_type_t::uniform:
      return flx_interpolate_linear(p, p_vec, q_vec, N_bins+1);
    case interpol_type_t::pchip:
      return (*pchip_icdf)(p);
    case interpol_type_t::bin_beta:
    {
      const size_t firstPleq = flx_interpolate_find_larger_eq(p,p_vec,N_bins+1);
      if (p==p_vec[firstPleq]) return q_vec[firstPleq];
      const tdouble alpha = bin_rv_params[(firstPleq-1)*2];
      const tdouble beta = bin_rv_params[(firstPleq-1)*2+1];
      const tdouble p_ = (p-p_vec[firstPleq-1])/(p_vec[firstPleq]-p_vec[firstPleq-1]);
      const tdouble q_bin = iBeta_reg_inv(alpha,beta,p_);
      return q_bin * (q_vec[firstPleq]-q_vec[firstPleq-1]) + q_vec[firstPleq-1];
    }
    case interpol_type_t::bin_linear:
    {
      const size_t firstPleq = flx_interpolate_find_larger_eq(p,p_vec,N_bins+1);
      if (p==p_vec[firstPleq]) return q_vec[firstPleq];
      const tdouble p_ = (p-p_vec[firstPleq-1])/(p_vec[firstPleq]-p_vec[firstPleq-1]);
      const tdouble m = bin_rv_params[firstPleq-1];
      const tdouble xx = (fabs(p_)<1e-12)?p_:((-(ONE-m)+sqrt(pow2(ONE-m)+4*m*p_))/(2*m));
      return xx * (q_vec[firstPleq]-q_vec[firstPleq-1]) + q_vec[firstPleq-1];
    }
    case interpol_type_t::pdf_linear:
    {
      const size_t firstPleq = flx_interpolate_find_larger_eq(p,p_vec,N_bins+1);
      const tdouble delta_pr = p-p_vec[firstPleq-1];
      const tdouble dx = q_vec[firstPleq]-q_vec[firstPleq-1];
      const tdouble a = bin_rv_params[firstPleq-1];
      const tdouble b = bin_rv_params[firstPleq];
      const tdouble delta_x = dx/(b-a)*(sqrt(pow2(a)+2*(b-a)/dx*delta_pr)-a);
      return delta_x + q_vec[firstPleq-1];
    }
    default:
      throw FlxException("RBRV_entry_RV_quantiles::transform_y2x");
  }
}

const tdouble RBRV_entry_RV_quantiles::transform_x2y(const tdouble& x_val)
{
  return rv_InvPhi_noAlert(this->calc_cdf_x(x_val));
}

const tdouble RBRV_entry_RV_quantiles::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  const tdouble x_ref = q_vec[N_bins/2];
  const tdouble stepf = 1e-6;
  const tdouble dx = stepf*x_ref;
  const tdouble pr_low = this->calc_cdf_x(x_val-dx);
  const tdouble pr_up  = this->calc_cdf_x(x_val+dx);
  return (pr_up-pr_low)/(2*dx);
}

const tdouble RBRV_entry_RV_quantiles::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  if (x_val<q_vec[1]) {
    if (tail_low) {
      const tdouble x_transformed = q_vec[1] - x_val;
      const tdouble p_transformed = tail_low->calc_cdf_x(x_transformed, safeCalc);
      return (ONE-p_transformed)*p_vec[1];
    }
  } else if (x_val>q_vec[N_bins-1]) {
    if (tail_up) {
      const tdouble x_transformed = x_val-q_vec[N_bins-1];
      const tdouble p_transformed = tail_up->calc_cdf_x(x_transformed, safeCalc);
      return p_transformed*(ONE-p_vec[N_bins-1])+p_vec[N_bins-1];
    }
  }
  switch (interpol_type) {
    case interpol_type_t::uniform:
      return flx_interpolate_linear(x_val, q_vec, p_vec, N_bins+1);
    case interpol_type_t::pchip:
      return (*pchip_cdf)(x_val);
    case interpol_type_t::bin_beta:
    {
      const size_t firstXleq = flx_interpolate_find_larger_eq(x_val,q_vec,N_bins+1);
      if (x_val==q_vec[firstXleq]) return p_vec[firstXleq];
      const tdouble alpha = bin_rv_params[(firstXleq-1)*2];
      const tdouble beta = bin_rv_params[(firstXleq-1)*2+1];
      const tdouble xx = (x_val-q_vec[firstXleq-1])/(q_vec[firstXleq]-q_vec[firstXleq-1]);
      const tdouble p_bin = iBeta_reg(alpha,beta,xx);
      return p_bin * (p_vec[firstXleq]-p_vec[firstXleq-1]) + p_vec[firstXleq-1];
    }
    case interpol_type_t::bin_linear:
    {
      const size_t firstXleq = flx_interpolate_find_larger_eq(x_val,q_vec,N_bins+1);
      if (x_val==q_vec[firstXleq]) return p_vec[firstXleq];
      const tdouble m = bin_rv_params[firstXleq-1];
      const tdouble delta_x = x_val-q_vec[firstXleq-1];
      const tdouble xx = delta_x/(q_vec[firstXleq]-q_vec[firstXleq-1]);   // normalized to [0,1]
      const tdouble p1 = ONE-m;
      const tdouble p_x = ONE+m*(2*xx-ONE);
      const tdouble delta_p = (p1+p_x)/2*xx*(p_vec[firstXleq]-p_vec[firstXleq-1]);
      return p_vec[firstXleq-1] + delta_p;
    }
    case interpol_type_t::pdf_linear:
    {
      const size_t firstXleq = flx_interpolate_find_larger_eq(x_val,q_vec,N_bins+1);
      const tdouble dx = (q_vec[firstXleq]-q_vec[firstXleq-1]);
      const tdouble xx = (x_val-q_vec[firstXleq-1])/dx;
      const tdouble pdfx = bin_rv_params[firstXleq-1] + xx*(bin_rv_params[firstXleq]-bin_rv_params[firstXleq-1]);
      return (bin_rv_params[firstXleq-1]+pdfx)/2*(xx*dx) + p_vec[firstXleq-1];
    }
    default:
      throw FlxException_Crude("RBRV_entry_RV_quantiles::calc_cdf_x");
  }
}

const tdouble RBRV_entry_RV_quantiles::get_mean_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_quantiles::get_mean_current_config");
}

const tdouble RBRV_entry_RV_quantiles::get_median_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_quantiles::get_median_current_config");
}

const tdouble RBRV_entry_RV_quantiles::get_mode_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_quantiles::get_mode_current_config");
}

const tdouble RBRV_entry_RV_quantiles::get_sd_current_config()
{
  throw FlxException_NotImplemented("RBRV_entry_RV_quantiles::get_sd_current_config");
}

const bool RBRV_entry_RV_quantiles::check_x(const tdouble xV)
{
  if (xV<q_vec[0]) {
    if (tail_low) {
      const tdouble x_transformed = q_vec[1] - xV;
      return tail_low->check_x(x_transformed);
    } else {
      return false;
    }
  } else if (xV>q_vec[N_bins]) {
    if (tail_up) {
      const tdouble x_transformed = xV - q_vec[N_bins-1];
      return tail_up->check_x(x_transformed);
    } else {
      return false;
    }
  } else {
    return true;
  }
}

const bool RBRV_entry_RV_quantiles::search_circref(FlxFunction* fcr)
{
  bool b = false;
  if (tail_low) {
    b = b || tail_low->search_circref(fcr);
  }
  if (tail_up) {
    b = b || tail_up->search_circref(fcr);
  }
  return b;
}



tdouble calc_expectation_numerical_1D::calc_expectation(const tulong N_Interv, const tulong N, const tdouble pRedistribute, RBRV_entry_RV_base& rnd, const tdouble LB, const tdouble UB)
{
  // remember original setting
  tdouble orig_x = ZERO;
  bool orig_x_valid = false;
  try {
    orig_x = rnd.get_value();
    orig_x_valid = true;
  } catch (FlxException& e) {
    orig_x_valid = false;
  }
  
  rnd.eval_para();
  tdouble res = ZERO;
  if ( N<=1 || N_Interv<=1 ) {
    res = calc_Interval(UB,LB,(N<=1)?N_Interv:N,rnd).cast2double();
  } else {
    const tdouble step=(UB-LB)/N_Interv;
    flxpVec I_vec(N_Interv);
    std::valarray<tulong> N_Vec(1,N_Interv);
    pdouble Integ=ZERO, t;
    tdouble x_u = LB, x_o = x_u+step;
    
    // first crude sampling
      for (tulong i=0;i<N_Interv;++i) {
        t = calc_Interval(x_o,x_u,1,rnd);
        I_vec[i] = t;
        Integ += sqrt(fabs(t.cast2double()));
        x_u=x_o;
        x_o+=step;
      } 
      tdouble tt = Integ.cast2double();
    if (std::isfinite(tt)==false) {
      res = tt;
    } else if (tt==ZERO) {
      res = ZERO;
    } else {
      //redistribution
        // adaptive
          tulong max = 0;
          tulong maxI = 0;
          tulong count = 0;
          for (tulong i=0;i<N_Interv;++i) {
            N_Vec[i] = tulong(std::floor((sqrt(fabs((I_vec[i]).cast2double()))/tt)*N));
            if (N_Vec[i]>max) {
              maxI = i;
              max = N_Vec[i];
            }
            count+=N_Vec[i];
          }
          N_Vec[maxI] += N - count;
        // redistribute (diffusion)
          count = 0;
          tulong rd = tulong(N_Vec[0]*pRedistribute/2);
          N_Vec[0]-=rd;
          count +=N_Vec[0];
          for (tulong i=1;i<N_Interv-1;++i) {
            N_Vec[i]+=rd;
            rd = tulong((N_Vec[i]-rd)*pRedistribute/2);
            N_Vec[i]-=2*rd;
            N_Vec[i-1]+=rd;
            count += N_Vec[i]+rd;
          }
          N_Vec[N_Interv-1]+=rd;
          rd = tulong((N_Vec[N_Interv-1]-rd)*pRedistribute/2);
          N_Vec[N_Interv-1]-=rd;
          N_Vec[N_Interv-2]+=rd;
          count += N_Vec[N_Interv-1]+rd;
          N_Vec[maxI] += N - count;
      // resampling
        Integ=ZERO;
        x_u = LB;
        x_o = x_u+step;
        for (tulong i=0;i<N_Interv;++i) {
          if (N_Vec[i]<=1) {
            Integ+=I_vec[i];
          } else {
            Integ += calc_Interval(x_o,x_u,N_Vec[i],rnd);
          }
          x_u=x_o;
          x_o+=step;
        }
      res = Integ.cast2double();
    }
    
  }
  // restore original setting
  if (orig_x_valid) {
    #if FLX_DEBUG
      rnd.set_is_valid(false);
    #endif
    rnd.set_x(orig_x);    
  }
  return res;
}

pdouble calc_expectation_numerical_1D::calc_Interval(const tdouble& x_o, const tdouble& x_u, const tulong& N, RBRV_entry_RV_base& rnd)
{
  const tdouble step=(x_o-x_u)/N;
  
  tulong NpI;
  if (N>=1e4) {
    NpI = tulong(floor(sqrt(tdouble(N))));
  } else {
    NpI = N;
  }
  
  tdouble y = x_u;
  tdouble CDF_new, CDF_old;
  tdouble f_new, f_old;
  pdouble Integ=ZERO, IntegI;
  // assign values
    #if FLX_DEBUG
      rnd.set_is_valid(false);
    #endif
    rnd.set_x(rnd.transform_y2x(y));
    f_old = fun->calc();
    bool b_comp;
    if (y>ZERO) {
      b_comp = true;
      CDF_old = rv_Phi(-y);
    } else {
      b_comp = false;
      CDF_old = rv_Phi(y);
    }
    y+=step;
  tulong j_start=1;
  for (int k=(NpI==N)?2:1; k<=2; ++k) {
    for (tulong j = (NpI==N)?NpI:j_start; j <= NpI; ++j) {
      IntegI=ZERO;
      for (tulong i = 0; i < NpI; ++i) {
        CDF_new = (b_comp)?rv_Phi(-y):rv_Phi(y);
        #if FLX_DEBUG
          rnd.set_is_valid(false);
        #endif
        rnd.set_x(rnd.transform_y2x(y));
        f_new = fun->calc();
        IntegI+= ((b_comp)?(CDF_old-CDF_new):(CDF_new-CDF_old))/2*(f_new + f_old);
        y += step;
        f_old = f_new;
        CDF_old = CDF_new;
      }
      Integ+=IntegI;
    }
    if (N!=NpI) {
      NpI = N-pow2(NpI);
      j_start = NpI;
    }
  }
  return Integ;
}



tdouble calc_expectation_numerical_MCI::calc_expectation(FunBase* fun, RBRV_constructor& RndBoxN, const tulong N_Interv, const tulong N_mciI_1, const tulong N_mciI_2, const tdouble pRedistribute, const tdouble Bound)
{
  funR = fun;
  // remember x and y
    flxVec orig_X(RndBoxN.get_NOX());
    RndBoxN.get_x_Vec(orig_X.get_tmp_vptr());
    flxVec orig_Y(RndBoxN.get_NRV());
    RndBoxN.get_y_Vec(orig_Y.get_tmp_vptr());
    
  pdouble Integ(0), IntegI, IntegI2;
  flxpVec Ivec(N_Interv);
  
  // get dimension of the problem
    DIM = RndBoxN.get_NRV();
    flxVec tempVec(DIM);

  const tdouble step=Bound/N_Interv;
  tdouble h_inv = Bound;
  
  tdouble x_u, x_o;
  tulong N_mciI_T, N_mciI_T_I;
  std::valarray<tulong> N_Vec(N_mciI_1,N_Interv);
  FunNumber* lb = new FunNumber(-ONE);
  FunNumber* ub = new FunNumber(ONE);
  uD = new RBRV_entry_RV_uniform("p",0,new FlxFunction(lb),new FlxFunction(ub),false);
  for (int k=1;k<=2;++k) {
    x_u=ZERO; x_o = step;
    for (tulong i=0;i<N_Interv;++i) {
      // change properties of uniform distribution
        lb->set_thenumber(x_u);
        ub->set_thenumber(x_o);
      IntegI = ZERO;
      N_mciI_T = N_Vec[i];
      if (N_mciI_T>=1e4) {
        N_mciI_T_I = tulong(floor(sqrt(tdouble(N_mciI_T))));
        for (int l=1;l<=2;++l) {
          for (tulong j=(l==2)?N_mciI_T_I-1:1;j<N_mciI_T_I;++j) {
            IntegI2 = ZERO;
            for (tulong m=0;m<N_mciI_T_I;++m) {
              IntegI2+=Integrationstep(tempVec,RndBoxN);
            }
            IntegI+=IntegI2;
          }
          N_mciI_T_I = N_mciI_T - pow2(N_mciI_T_I) + N_mciI_T_I;
        }
      } else {
        for (tulong j=0;j<N_mciI_T;++j) {
          IntegI+=Integrationstep(tempVec,RndBoxN);
        }
      }
      Ivec[i]+=IntegI;
      Integ+=sqrt(fabs(IntegI.cast2double()));        // only relevant for k=1
      x_u = x_o;
      x_o += step;
    }
    if (k==1 && N_Interv > 1) {        //redistribution
      // adaptive 
        tulong max = 0;
        tulong maxI = 0;
        tulong count = 0;
        tdouble tt = Integ.cast2double();
        for (tulong i=0;i<N_Interv;++i) {
          N_Vec[i] = tulong(std::floor((sqrt(fabs((Ivec[i]).cast2double()))/tt)*N_mciI_2));
          if (N_Vec[i]>max) {
            maxI = i;
            max = N_Vec[i];
          }
          count+=N_Vec[i];
        }
        N_Vec[maxI] += N_mciI_2 - count;
      // redistribute (diffusion)
        count = 0;
        tulong rd = tulong(N_Vec[0]*pRedistribute/2);
        N_Vec[0]-=rd;
        count +=N_Vec[0];
        for (tulong i=1;i<N_Interv-1;++i) {
          N_Vec[i]+=rd;
          rd = tulong((N_Vec[i]-rd)*pRedistribute/2);
          N_Vec[i]-=2*rd;
          N_Vec[i-1]+=rd;
          count += N_Vec[i]+rd;
        }
        N_Vec[N_Interv-1]+=rd;
        rd = tulong((N_Vec[N_Interv-1]-rd)*pRedistribute/2);
        N_Vec[N_Interv-1]-=rd;
        N_Vec[N_Interv-2]+=rd;
        count += N_Vec[N_Interv-1]+rd;
        N_Vec[maxI] += N_mciI_2 - count;
    } else if (k==1 && N_Interv == 1) {
      N_Vec[0] = N_mciI_2;
    }
  }
  delete uD;
  
  // sum of everything
    Integ = ZERO;
    tulong j = N_Interv*N_mciI_1+N_mciI_2;
    tdouble jm = tdouble(j)/tdouble(N_Interv);
    for (tulong i=0;i<N_Interv;++i) {
      Integ += Ivec[i] * (jm/(N_mciI_1+N_Vec[i]));
    }
    tdouble result = Integ.cast2double() * h_inv / j;
  
  // restore x and y
    RndBoxN.set_smp_y(orig_Y);
    RndBoxN.set_smp_x(orig_X);
  // free memory
    funR = NULL;
  return result;
}

const tdouble calc_expectation_numerical_MCI::Integrationstep(flxVec& tempVec, RBRV_constructor& RndBoxN)
{
  // generate y-vector
    flxVec &y = tempVec;
    RndBoxN.propose_y(y);

  // transform to z-space
    const tdouble ly2 = y.get_Norm2_NOroot();
    const tdouble lz2 = uD->Inv_cdf_x(rv_cdf_ChiSquare(DIM,ly2));
    flxVec &z = tempVec;
    z*= sqrt(lz2/ly2);
  
  // calculate f
    const tdouble f = rv_pdf_ChiSquare(DIM,lz2);
    
  // update standard normal space
    RndBoxN.set_smp(z);

  // return function * foverh;
    return funR->calc()*f;        // h has to be multiplied afterwards
}





