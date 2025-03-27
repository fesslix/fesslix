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

#include "flxrbrv_rvs.h"

#include "flxdata.h"

FlxFunction * parse_py_para(const std::string& para_name, py::dict config)
{
  if (config.contains(para_name.c_str()) == false) {
    std::ostringstream ssV;
    ssV << "Key '" << para_name << "' not found in Python <dict>.";
    throw FlxException_NeglectInInteractive("parse_py_para_01", ssV.str());
  }
  return parse_function(config[para_name.c_str()], "key '"+para_name+"' in Python <dict>");
}

const bool parse_py_para_as_bool(const std::string& para_name, py::dict config, const bool required, const bool def_val)
{
  if (config.contains(para_name.c_str())) {
    try {
      return py::cast<bool>(config[para_name.c_str()]);
    } catch (const py::cast_error &e) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_bool_01", "Key '"+para_name+"' in Python <dict> cannot be cast into type 'bool'.");
    }
  } else {
    if (required) {
      throw FlxException_NeglectInInteractive("parse_py_para_as_bool_02", "Key '" + para_name + "' not found in Python <dict>.");
    } else {
      return def_val;
    }
  }
}


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

void RBRV_entry_RV_stdN::info(std::ostream& os)
{
  os << "standard Normal distribution" << std::endl;
  os << "  mean    = " << GlobalVar.Double2String(get_mean_current_config()) << std::endl;
  os << "  std.dev = " << GlobalVar.Double2String(get_sd_current_config()) << std::endl;
  os << "  entropy = " << GlobalVar.Double2String(calc_entropy()) << std::endl;
}

const tdouble RBRV_entry_RV_stdN::get_HPD(const tdouble p)
{
  return (ONE-p)/2;
}


RBRV_entry_RV_normal::RBRV_entry_RV_normal(const std::string& name, const tuint iID, const int pid, FlxFunction* p1v, FlxFunction* p2v, FlxFunction* p3v, FlxFunction* p4v, const bool eval_once)
: RBRV_entry_RV_base(name,iID), pid(pid), p1(p1v), p2(p2v), p3(p3v), p4(p4v), eval_once(eval_once), mu(ZERO), sigma(ZERO)
{

}

RBRV_entry_RV_normal::RBRV_entry_RV_normal(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), pid(0), p1(nullptr), p2(nullptr), p3(nullptr), p4(nullptr), eval_once(false), mu(ZERO), sigma(ZERO)
{
  try {
    if (config.contains("mu")) {          // mean, standard deviation
      pid = 0;
      p1 = parse_py_para("mu", config);
      p2 = parse_py_para("sd", config);
    } else if (config.contains("cov")) {   // // C.o.V. and quantile value
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

void RBRV_entry_RV_normal::get_paras()
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
  get_paras();
  return mu+sigma*y_val;
}

const tdouble RBRV_entry_RV_normal::transform_x2y(const tdouble& x_val)
{
  get_paras();
  return (x_val-mu)/sigma;
}

const tdouble RBRV_entry_RV_normal::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_paras();
  return rv_phi((x_val-mu)/sigma)/sigma;
}

const tdouble RBRV_entry_RV_normal::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  get_paras();
  return rv_phi_log((x_val-mu)/sigma)-log(sigma);
}

const tdouble RBRV_entry_RV_normal::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_paras();
  return rv_Phi((x_val-mu)/sigma);
}

const tdouble RBRV_entry_RV_normal::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  get_paras();
  return rv_Phi((mu-x_val)/sigma);
}

const tdouble RBRV_entry_RV_normal::calc_entropy()
{
  get_paras();
  const tdouble s = sigma;
  return log(2*PI*exp(ONE)*pow2(s))/2;
}

const bool RBRV_entry_RV_normal::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || p1->search_circref(fcr) || p2->search_circref(fcr) || (p3?(p3->search_circref(fcr)):false) || (p4?(p4->search_circref(fcr)):false);
}

void RBRV_entry_RV_normal::info(std::ostream& os)
{
  get_paras();
  os << "Normal distribution" << std::endl;
  os << "  mean    = " << GlobalVar.Double2String(get_mean_current_config()) << std::endl;
  os << "  std.dev = " << GlobalVar.Double2String(get_sd_current_config()) << std::endl;
  os << "  entropy = " << GlobalVar.Double2String(calc_entropy()) << std::endl;
}

const tdouble RBRV_entry_RV_normal::get_HPD(const tdouble p)
{
  get_paras();
  return (ONE-p)/2;
}

RBRV_entry_RV_lognormal::~RBRV_entry_RV_lognormal()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (p3) delete p3;
  if (p4) delete p4;
  if (epsilon) delete epsilon;
}

void RBRV_entry_RV_lognormal::get_paras()
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
  get_paras();
  return exp(y_val*zeta+lambda)+eps;
}

const tdouble RBRV_entry_RV_lognormal::transform_x2y(const tdouble& x_val)
{
  get_paras();
  if (x_val<=eps) {
    std::ostringstream ssV;
    ssV << "A value (" << GlobalVar.Double2String(x_val) << ") smaller or equal than '" << GlobalVar.Double2String(eps) << "' is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_lognormal::transform_x2y", ssV.str() );
  }
  return (log(x_val-eps)-lambda)/zeta;
}

const tdouble RBRV_entry_RV_lognormal::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_paras();
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
  get_paras();
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
  get_paras();
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
  get_paras();
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
  get_paras();
  return (ONE+log(2*PI*pow2(zeta)))/2+lambda;
}

const tdouble RBRV_entry_RV_lognormal::get_mean_current_config()
{
  get_paras();
  return exp(lambda+pow2(zeta)/(ONE*2))+eps;
}

const tdouble RBRV_entry_RV_lognormal::get_sd_current_config()
{
  get_paras();
  return sqrt(exp(pow2(zeta))-ONE)*exp(lambda+pow2(zeta)/(ONE*2));
}

const tdouble RBRV_entry_RV_lognormal::get_median_current_config()
{
  get_paras();
  return exp(lambda)+eps;
}

const tdouble RBRV_entry_RV_lognormal::get_mode_current_config()
{
  get_paras();
  return exp(lambda-pow2(zeta))+eps;
}

const bool RBRV_entry_RV_lognormal::check_x(const tdouble xV)
{
  get_paras();
  return (xV>eps);
}

const bool RBRV_entry_RV_lognormal::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (p3?(p3->search_circref(fcr)):false) || (p4?(p4->search_circref(fcr)):false) || (epsilon?(epsilon->search_circref(fcr)):false);
}

void RBRV_entry_RV_lognormal::info(std::ostream& os)
{
  get_paras();
  os << "log-Normal distribution" << std::endl;
  os << "  lambda  = " << GlobalVar.Double2String(lambda) << std::endl;
  os << "  zeta    = " << GlobalVar.Double2String(zeta) << std::endl;
  os << "  epsilon = " << GlobalVar.Double2String(eps) << std::endl;
  os << "  mean    = " << GlobalVar.Double2String(get_mean_current_config()) << std::endl;
  os << "  std.dev = " << GlobalVar.Double2String(get_sd_current_config()) << std::endl;
  os << "  entropy = " << GlobalVar.Double2String(calc_entropy()) << std::endl;
}

const tdouble RBRV_entry_RV_lognormal::get_CoeffOfVar_withoutEpsilon()
{
  get_paras();
  return get_sd_current_config()/(get_mean_current_config()-eps);
}


RBRV_entry_RV_uniform::~RBRV_entry_RV_uniform()
{
  if (a) delete a;
  if (b) delete b;
}

void RBRV_entry_RV_uniform::get_paras()
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
  get_paras();
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
  get_paras();
  if (x_val>bv || x_val<av) {
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_uniform::transform_x2y", ssV.str() );
  }
  return rv_InvPhi((x_val-av)/(bv-av));
}

const tdouble RBRV_entry_RV_uniform::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_paras();
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
  get_paras();
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

const tdouble RBRV_entry_RV_uniform::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  get_paras();
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
  get_paras();
  return log(bv-av);
}

const tdouble RBRV_entry_RV_uniform::get_mean_current_config()
{
  get_paras();
  return (av+bv)/2;
}

const tdouble RBRV_entry_RV_uniform::get_sd_current_config()
{
  get_paras();
  return (bv+av)/sqrt(12*ONE);
}

const bool RBRV_entry_RV_uniform::check_x(const tdouble xV)
{
  get_paras();
  return (xV<=bv&&xV>=av);
}

const bool RBRV_entry_RV_uniform::search_circref(FlxFunction* fcr)
{
  const bool bb = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return bb || (a?(a->search_circref(fcr)):false) || (b?(b->search_circref(fcr)):false);
}

const tdouble RBRV_entry_RV_uniform::get_HPD(const tdouble p)
{
  get_paras();
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
  get_paras();
  return av+p*(bv-av);
}


RBRV_entry_RV_Gumbel::~RBRV_entry_RV_Gumbel()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (p3) delete p3;
  if (p4) delete p4;
}

void RBRV_entry_RV_Gumbel::get_pars()
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
  get_pars();
  return u-(std::log((rv_Phi(y_val)==ONE)?(rv_Phi(-y_val)):(-std::log(rv_Phi(y_val)))))/alpha;
}

const tdouble RBRV_entry_RV_Gumbel::transform_x2y(const tdouble& x_val)
{
  get_pars();
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
  get_pars();
  return alpha*exp(-alpha*(x_val-u)-exp(-alpha*(x_val-u)));
}

const tdouble RBRV_entry_RV_Gumbel::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
  const tdouble res = log(alpha)+(-alpha*(x_val-u)-exp(-alpha*(x_val-u)));
  return res;
}

const tdouble RBRV_entry_RV_Gumbel::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
  return exp(-exp(-alpha*(x_val-u)));
}

const tdouble RBRV_entry_RV_Gumbel::calc_entropy()
{
  get_pars();
  return ONE+GAMMA-log(alpha);
}

const tdouble RBRV_entry_RV_Gumbel::get_mean_current_config()
{
  get_pars();
  return u + GAMMA/alpha;
}

const tdouble RBRV_entry_RV_Gumbel::get_sd_current_config()
{
  get_pars();
  return PI/(sqrt(6*ONE)*alpha);
}

const tdouble RBRV_entry_RV_Gumbel::get_median_current_config()
{
  get_pars();
  return u - log(log(2*ONE))/alpha;
}

const tdouble RBRV_entry_RV_Gumbel::get_mode_current_config()
{
  get_pars();
  return u;
}

const bool RBRV_entry_RV_Gumbel::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (p3?(p3->search_circref(fcr)):false) || (p4?(p4->search_circref(fcr)):false);
}

void RBRV_entry_RV_Gumbel::info(std::ostream& os)
{
  get_pars();
  os << "Gumbel distribution" << std::endl;
  os << "  u       = " << GlobalVar.Double2String(u) << std::endl;
  os << "  alpha   = " << GlobalVar.Double2String(alpha) << std::endl;
  os << "  mean    = " << GlobalVar.Double2String(get_mean_current_config()) << std::endl;
  os << "  std.dev = " << GlobalVar.Double2String(get_sd_current_config()) << std::endl;
  os << "  entropy = " << GlobalVar.Double2String(calc_entropy()) << std::endl;
}

RBRV_entry_RV_normal_trunc::~RBRV_entry_RV_normal_trunc()
{
  if (m) delete m;
  if (s) delete s;
  if (a) delete a;
  if (b) delete b;
}

void RBRV_entry_RV_normal_trunc::get_pars()
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
  get_pars();
  
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
  get_pars();
  if (x_val>bV || x_val<aV) {
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(aV) << ";" << GlobalVar.Double2String(bV) << "].";
    throw FlxException("RBRV_entry_RV_normal_trunc::transform_x2y", ssV.str() );
  }
  return rv_InvPhi_noAlert((rv_Phi((x_val-mV)/sV)-rv_Phi(alpha))/q);
}

const tdouble RBRV_entry_RV_normal_trunc::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
  return mV+(rv_phi(alpha)-rv_phi(beta))/q*sV;
}

const tdouble RBRV_entry_RV_normal_trunc::get_sd_current_config()
{
  get_pars();
  return sV*sqrt(
    ONE
    + (alpha*rv_phi(alpha)-beta*rv_phi(beta))/q
    - pow2((rv_phi(alpha)-rv_phi(beta))/q)
  );
}

const tdouble RBRV_entry_RV_normal_trunc::get_median_current_config()
{
  get_pars();
  return mV+rv_InvPhi((rv_Phi(alpha)+rv_Phi(beta))/2)*sV;
}

const tdouble RBRV_entry_RV_normal_trunc::get_mode_current_config()
{
  get_pars();
  if (mV<alpha) return alpha;
  if (mV>beta) return beta;
  return mV;
}

const tdouble RBRV_entry_RV_normal_trunc::calc_entropy()
{
  get_pars();
  return log(sqrt(2*PI*exp(ONE))*sV*q)+(alpha*rv_phi(alpha)-beta*rv_phi(beta))/(2*q);
}

const bool RBRV_entry_RV_normal_trunc::check_x(const tdouble xV)
{
  get_pars();
  return (xV<=bV&&xV>=aV);
}

const bool RBRV_entry_RV_normal_trunc::search_circref(FlxFunction* fcr)
{
  const bool bb = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return bb || (m?(m->search_circref(fcr)):false) || (s?(s->search_circref(fcr)):false) || (a?(a->search_circref(fcr)):false) || (b?(b->search_circref(fcr)):false);
}

void RBRV_entry_RV_normal_trunc::info(std::ostream& os)
{
  get_pars();
  os << "truncated-Normal distribution" << std::endl;
  os << "  m       = " << GlobalVar.Double2String(mV) << std::endl;
  os << "  s       = " << GlobalVar.Double2String(sV) << std::endl;
  os << "  a       = " << GlobalVar.Double2String(aV) << std::endl;
  os << "  b       = " << GlobalVar.Double2String(bV) << std::endl;
  os << "  alpha   = " << GlobalVar.Double2String(alpha) << std::endl;
  os << "  beta    = " << GlobalVar.Double2String(beta) << std::endl;
  os << "  q       = " << GlobalVar.Double2String(q) << std::endl;
}

RBRV_entry_RV_beta::~RBRV_entry_RV_beta()
{
  if (p1) delete p1;
  if (p2) delete p2;
  if (a) delete a;
  if (b) delete b;
}

void RBRV_entry_RV_beta::get_pars()
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
      if (pow2(sigma)>=(ONE-mu)*mu) {
        std::ostringstream ssV;
        ssV << this->name << ": 'sigma^2' (" << GlobalVar.Double2String(sigma) << "²=" << GlobalVar.Double2String(pow2(sigma)) << ") must be smaller than 'mu*(1.-mu)' (" << GlobalVar.Double2String((ONE-mu)*mu) << ") (mu=" << GlobalVar.Double2String(mu) << ").";
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
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
  if (x_val>bv || x_val<av) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_beta::calc_cdf_x", ssV.str() );
  }
  const tdouble xx = (x_val-av)/(bv-av);        // scale x to [0;1]
  return iBeta_reg(alpha,beta,xx);
}

const tdouble RBRV_entry_RV_beta::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
  if (x_val>bv || x_val<av) {
    if (safeCalc) return ZERO;
    std::ostringstream ssV;
    ssV << "Value (" << GlobalVar.Double2String(x_val) << ") is not within the valid bounds [" << GlobalVar.Double2String(av) << ";" << GlobalVar.Double2String(bv) << "].";
    throw FlxException("RBRV_entry_RV_beta::calc_sf_x", ssV.str() );
  }
  const tdouble xx = (bv-x_val)/(bv-av);        // scale x to [0;1]
  return iBeta_reg(beta,alpha,xx);
}

const tdouble RBRV_entry_RV_beta::calc_entropy()
{
  get_pars();
  const tdouble lnbetafun = BetaFunLn(alpha,beta);
  const tdouble entropy01 = lnbetafun 
          - (alpha-ONE)*flxdigamma(alpha)
          - (beta-ONE)*flxdigamma(beta)
          + (alpha+beta-2*ONE)*flxdigamma(alpha+beta);
  return entropy01/(bv-av) + log(bv-av);
}

const tdouble RBRV_entry_RV_beta::get_mean_current_config()
{
  get_pars();
  return alpha/(alpha+beta)*(bv-av)+av;
}

const tdouble RBRV_entry_RV_beta::get_sd_current_config()
{
  get_pars();
  return sqrt(alpha*beta/(pow2(alpha+beta)*(alpha+beta+ONE))) * (bv-av);
}

const tdouble RBRV_entry_RV_beta::get_median_current_config()
{
  return transform_y2x(ZERO);
}

const tdouble RBRV_entry_RV_beta::get_mode_current_config()
{
  get_pars();
  if (alpha>ONE && beta>ONE) {
    return (alpha-ONE)/(alpha+beta-2*ONE)*(bv-av)+av;
  }
  if (beta>ONE && alpha<=ONE) return ZERO;
  if (alpha>ONE && beta<=ONE) return ONE;
  throw FlxException_NotImplemented("RBRV_entry_RV_beta::get_mode_current_config");
}

const bool RBRV_entry_RV_beta::check_x(const tdouble xV)
{
  get_pars();
  return (xV<=bv&&xV>=av);
}

const bool RBRV_entry_RV_beta::search_circref(FlxFunction* fcr)
{
  const bool bb = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return bb || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (a?(a->search_circref(fcr)):false) || (b?(b->search_circref(fcr)):false);
}

void RBRV_entry_RV_beta::info(std::ostream& os)
{
  get_pars();
  os << "beta distribution" << std::endl;
  os << "  alpha   = " << GlobalVar.Double2String(alpha) << std::endl;
  os << "  beta    = " << GlobalVar.Double2String(beta) << std::endl;
  os << "  a       = " << GlobalVar.Double2String(av) << std::endl;
  os << "  b       = " << GlobalVar.Double2String(bv) << std::endl;
  os << "  mean    = " << GlobalVar.Double2String(get_mean_current_config()) << std::endl;
  os << "  std.dev = " << GlobalVar.Double2String(get_sd_current_config()) << std::endl;
  os << "  entropy = " << GlobalVar.Double2String(calc_entropy()) << std::endl;
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
  get_pars();
  const tdouble x = iBeta_reg_inv(alpha,beta,p);
  return x*(bv-av)+av;        // scale x to [a;b]
}

RBRV_entry_RV_exponential::~RBRV_entry_RV_exponential()
{
  delete lambda;
  if (epsilon) delete epsilon;
}

const tdouble RBRV_entry_RV_exponential::transform_y2x(const tdouble y_val)
{
  const tdouble lambdaV = lambda->cast2positive();
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
  tdouble res = -ONE*std::log((y_val>ZERO)?(rv_Phi(-y_val)):(ONE-rv_Phi(y_val)))/lambdaV;
  res += eps;
  return res;
}

const tdouble RBRV_entry_RV_exponential::transform_x2y(const tdouble& x_val)
{
  const tdouble lambdaV = lambda->cast2positive();
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
  if (x_val<eps) {
    std::ostringstream ssV;
    ssV << "A negative value (" << GlobalVar.Double2String(x_val) << ") is not allowed at this point.";
    throw FlxException("RBRV_entry_RV_exponential::transform_x2y", ssV.str() );
  }
  return rv_InvPhi_noAlert(-expm1(-lambdaV*(x_val-eps)));
}

const tdouble RBRV_entry_RV_exponential::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  const tdouble lambdaV = lambda->cast2positive();
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
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
  const tdouble lambdaV = lambda->cast2positive();
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
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
  const tdouble lambdaV = lambda->cast2positive();
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
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
  return ONE - log(lambda->cast2positive());
}

const tdouble RBRV_entry_RV_exponential::get_mean_current_config()
{
  const tdouble lambdaV = lambda->cast2positive();
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
  return eps+ONE/lambdaV;
}

const tdouble RBRV_entry_RV_exponential::get_sd_current_config()
{
  const tdouble lambdaV = lambda->cast2positive();
  return ONE/lambdaV;
}

const tdouble RBRV_entry_RV_exponential::get_median_current_config()
{
  const tdouble lambdaV = lambda->cast2positive();
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
  return eps+log(2*ONE)/lambdaV;
}

const tdouble RBRV_entry_RV_exponential::get_mode_current_config()
{
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
  return eps;
}

const bool RBRV_entry_RV_exponential::check_x(const tdouble xV)
{
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
  return (xV>=eps);
}

const bool RBRV_entry_RV_exponential::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (lambda?(lambda->search_circref(fcr)):false) || (epsilon?(epsilon->search_circref(fcr)):false);
}

void RBRV_entry_RV_exponential::info(std::ostream& os)
{
  const tdouble eps = (epsilon?(epsilon->calc()):ZERO);
  os << "exponential distribution" << std::endl;
  os << "  lambda  = " << GlobalVar.Double2String(lambda->calc()) << std::endl;
  os << "  epsilon = " << GlobalVar.Double2String(eps) << std::endl;
  os << "  mean    = " << GlobalVar.Double2String(get_mean_current_config()) << std::endl;
  os << "  std.dev = " << GlobalVar.Double2String(get_sd_current_config()) << std::endl;
  os << "  entropy = " << GlobalVar.Double2String(calc_entropy()) << std::endl;
}



RBRV_entry_RV_gamma::~RBRV_entry_RV_gamma()
{
  if (p2) delete p2;
  if (p1) delete p1;
  if (epsilon) delete epsilon;
}

void RBRV_entry_RV_gamma::get_pars()
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
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
  return k-log(lambda)+GammaLn(k)+(ONE-k)*flxdigamma(k);
}

const tdouble RBRV_entry_RV_gamma::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
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
  get_pars();
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
  get_pars();
  return eps+k/lambda;
}

const tdouble RBRV_entry_RV_gamma::get_sd_current_config()
{
  get_pars();
  return sqrt(k)/lambda;
}

const tdouble RBRV_entry_RV_gamma::get_median_current_config()
{
  return transform_y2x(ZERO);
}

const tdouble RBRV_entry_RV_gamma::get_mode_current_config()
{
  get_pars();
  if (k>=ONE) return eps+(k-ONE)/lambda;
  throw FlxException_NotImplemented("RBRV_entry_RV_gamma::get_mode_current_config");
}

const bool RBRV_entry_RV_gamma::check_x(const tdouble xV)
{
  get_pars();
  return (xV>eps);
}

const bool RBRV_entry_RV_gamma::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (epsilon?(epsilon->search_circref(fcr)):false);
}

void RBRV_entry_RV_gamma::info(std::ostream& os)
{
  get_pars();
  os << "Gamma distribution" << std::endl;
  os << "  k       = " << GlobalVar.Double2String(k) << std::endl;
  os << "  lambda  = " << GlobalVar.Double2String(lambda) << std::endl;
  os << "  epsilon = " << GlobalVar.Double2String(eps) << std::endl;
  os << "  mean    = " << GlobalVar.Double2String(get_mean_current_config()) << std::endl;
  os << "  std.dev = " << GlobalVar.Double2String(get_sd_current_config()) << std::endl;
  os << "  entropy = " << GlobalVar.Double2String(calc_entropy()) << std::endl;
}


RBRV_entry_RV_Poisson::~RBRV_entry_RV_Poisson()
{
  delete mean;
}

const tdouble RBRV_entry_RV_Poisson::transform_y2x(const tdouble y_val)
{
  const tdouble p = rv_Phi(y_val);
  const tdouble meanV = mean->cast2positive();
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
  const tdouble lambdat = mean->cast2positive();
  return flxgamma_ru(std::floor(x_val)+1,lambdat);
}

const tdouble RBRV_entry_RV_Poisson::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  const tdouble lambdat = mean->cast2positive();
  return flxgamma_rl(std::floor(x_val)+1,lambdat);
}

const tdouble RBRV_entry_RV_Poisson::get_mean_current_config()
{
  return mean->cast2positive();
}

const tdouble RBRV_entry_RV_Poisson::get_sd_current_config()
{
  return sqrt(mean->cast2positive());
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


RBRV_entry_RV_Binomial::~RBRV_entry_RV_Binomial()
{
  if (p) delete p;
  if (N) delete N;
}

void RBRV_entry_RV_Binomial::get_pars()
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
  // get the parameters of the distribution
    get_pars();
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
  get_pars();
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
  get_pars();
  return _p*_N;
}

const tdouble RBRV_entry_RV_Binomial::get_sd_current_config()
{
  get_pars();
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


RBRV_entry_RV_Cauchy::~RBRV_entry_RV_Cauchy()
{
  delete loc;
  delete scale;
}

void RBRV_entry_RV_Cauchy::get_paras(tdouble& locv, tdouble& scalev) const
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
  tdouble locv, scalev;
  get_paras(locv,scalev);
  return locv;
}

const tdouble RBRV_entry_RV_Cauchy::get_mode_current_config()
{
  tdouble locv, scalev;
  get_paras(locv,scalev);
  return locv;
}

const bool RBRV_entry_RV_Cauchy::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (loc?(loc->search_circref(fcr)):false) || (scale?(scale->search_circref(fcr)):false);
}


const tdouble RBRV_entry_RV_Cauchy::transform_y2x(const tdouble y_val)
{
  tdouble locv, scalev;
  get_paras(locv,scalev);
  const tdouble p = (y_val>ZERO)?(ONE/2-rv_Phi(-y_val)):(rv_Phi(y_val)-ONE/2);
  return locv+scalev*tan(PI*(p));
}

const tdouble RBRV_entry_RV_Cauchy::transform_x2y(const tdouble& x_val)
{
  tdouble locv, scalev;
  get_paras(locv,scalev);
  const tdouble p = std::atan((x_val-locv)/scalev)/PI;
  return (p<=ZERO)?(rv_InvPhi_noAlert(p+ONE/2)):(-rv_InvPhi_noAlert(ONE/2-p));
}

const tdouble RBRV_entry_RV_Cauchy::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  tdouble locv, scalev;
  get_paras(locv,scalev);
  return (scalev/(pow2(x_val-locv)+pow2(scalev)))/PI;
}

const tdouble RBRV_entry_RV_Cauchy::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  tdouble locv, scalev;
  get_paras(locv,scalev);
  return std::atan((x_val-locv)/scalev)/PI+ONE/2;
}


RBRV_entry_RV_Weibull::~RBRV_entry_RV_Weibull()
{
  if (p2) delete p2;
  if (p1) delete p1;
  if (epsilon) delete epsilon;
}

void RBRV_entry_RV_Weibull::get_pars()
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
  get_pars();
  return get_mean_help();
}

const tdouble RBRV_entry_RV_Weibull::get_sd_current_config()
{
  get_pars();
  return get_sd_help();
}

const tdouble RBRV_entry_RV_Weibull::get_median_current_config()
{
  get_pars();
  return lambda*pow(log(2*ONE),ONE/k);
}

const tdouble RBRV_entry_RV_Weibull::get_mode_current_config()
{
  get_pars();
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
  get_pars();
  return (xV>eps);
}

const bool RBRV_entry_RV_Weibull::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false) || (p2?(p2->search_circref(fcr)):false) || (epsilon?(epsilon->search_circref(fcr)):false);
}

void RBRV_entry_RV_Weibull::info(std::ostream& os)
{
  get_pars();
  os << "Weibull distribution" << std::endl;
  os << "  k       = " << GlobalVar.Double2String(k) << std::endl;
  os << "  lambda  = " << GlobalVar.Double2String(lambda) << std::endl;
  os << "  epsilon = " << GlobalVar.Double2String(eps) << std::endl;
  os << "  mean    = " << GlobalVar.Double2String(get_mean_current_config()) << std::endl;
  os << "  std.dev = " << GlobalVar.Double2String(get_sd_current_config()) << std::endl;
  os << "  entropy = " << GlobalVar.Double2String(calc_entropy()) << std::endl;
}

const tdouble RBRV_entry_RV_Weibull::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
  return GAMMA*(ONE-ONE/k)+log(lambda/k)+ONE;
}

const tdouble RBRV_entry_RV_Weibull::transform_y2x(const tdouble y_val)
{
  get_pars();
  return pow(-log(rv_Phi(-y_val)),ONE/k)*lambda+eps;
}

const tdouble RBRV_entry_RV_Weibull::transform_x2y(const tdouble& x_val)
{
  get_pars();
  return -rv_InvPhi_noAlert(exp(-pow((x_val-eps)/lambda,k)));
}

RBRV_entry_RV_ChiSquared::~RBRV_entry_RV_ChiSquared()
{
  if (p1) delete p1;
}

void RBRV_entry_RV_ChiSquared::get_pars()
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
  get_pars();
  return (xV>=ZERO);
}

const tdouble RBRV_entry_RV_ChiSquared::get_mean_current_config()
{
  get_pars();
  return dof;
}

const tdouble RBRV_entry_RV_ChiSquared::get_sd_current_config()
{
  get_pars();
  return sqrt(2*dof);
}

const tdouble RBRV_entry_RV_ChiSquared::get_median_current_config()
{
  return transform_y2x(ZERO);
}

const tdouble RBRV_entry_RV_ChiSquared::get_mode_current_config()
{
  get_pars();
  if (dof-2>ZERO) return dof-2;
  return ZERO;
}

const tdouble RBRV_entry_RV_ChiSquared::calc_entropy()
{
  get_pars();
  const tdouble dofh = dof/2;
  return dofh 
        + std::log(2*flxgamma(dofh)) 
        + (ONE-dofh)*flxdigamma(dofh);
}

const tdouble RBRV_entry_RV_ChiSquared::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
  const tdouble dofh = dof/2;
  if (y_val<=ZERO) {
    return flxgamma_rl_inv(dofh,rv_Phi(y_val))*2;
  } else {
    return flxgamma_ru_inv(dofh,rv_Phi(-y_val))*2;
  }
}

const tdouble RBRV_entry_RV_ChiSquared::transform_x2y(const tdouble& x_val)
{
  get_pars();
  if (x_val <= dof) {
    return rv_InvPhi_noAlert(flxgamma_rl(dof/2,x_val/2));
  } else {
    return -rv_InvPhi_noAlert(flxgamma_ru(dof/2,x_val/2));
  }
}

RBRV_entry_RV_Chi::~RBRV_entry_RV_Chi()
{
  if (p1) delete p1;
}

void RBRV_entry_RV_Chi::get_pars()
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
  get_pars();
  const tdouble dofh = dof/2;
  if (y_val<=ZERO) {
    return sqrt(flxgamma_rl_inv(dofh,rv_Phi(y_val))*2);
  } else {
    return sqrt(flxgamma_ru_inv(dofh,rv_Phi(-y_val))*2);
  }
}

const tdouble RBRV_entry_RV_Chi::transform_x2y(const tdouble& x_val)
{
  get_pars();
  if (x_val <= dof) {
    return rv_InvPhi_noAlert(flxgamma_rl(dof/2,pow2(x_val)/2));
  } else {
    return -rv_InvPhi_noAlert(flxgamma_ru(dof/2,pow2(x_val)/2));
  }
}

const tdouble RBRV_entry_RV_Chi::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
  const tdouble dofh = dof/2;
  return GammaLn(dofh)
         + (dof-log(2)-(dof-ONE)*flxdigamma(dofh))/2;
}

const tdouble RBRV_entry_RV_Chi::get_mean_current_config()
{
  get_pars();
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
  get_pars();
  if (dof>=ONE) return sqrt(dof-ONE);
  throw FlxException_NotImplemented("RBRV_entry_RV_Chi::get_mode_current_config");
}

const bool RBRV_entry_RV_Chi::check_x(const tdouble xV)
{
  get_pars();
  return (xV>=ZERO);
}

const bool RBRV_entry_RV_Chi::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || (p1?(p1->search_circref(fcr)):false);
}



RBRV_entry_RV_StudentsT::RBRV_entry_RV_StudentsT(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), p1(nullptr), eval_once(false)
{
  try {
    p1 = parse_py_para("dof", config);

    eval_once = parse_py_para_as_bool("eval_once", config, false, false);
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

void RBRV_entry_RV_StudentsT::get_pars()
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
  get_pars();
  if (y_val<=ZERO) {
    return rv_InvCDF_Studentst(dof,rv_Phi(y_val));
  } else {
    return -rv_InvCDF_Studentst(dof,rv_Phi(-y_val));
  }
}

const tdouble RBRV_entry_RV_StudentsT::transform_x2y(const tdouble& x_val)
{
  get_pars();
  if (x_val<=ZERO) {
    return rv_InvPhi_noAlert(rv_cdf_Studentst(dof,x_val));
  } else {
    return -rv_InvPhi_noAlert(rv_cdf_Studentst(dof,-x_val));
  }
}

const tdouble RBRV_entry_RV_StudentsT::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
  return rv_pdf_Studentst(dof,x_val);
}

const tdouble RBRV_entry_RV_StudentsT::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
  return rv_cdf_Studentst(dof,x_val);
}

const tdouble RBRV_entry_RV_StudentsT::calc_entropy()
{
  get_pars();
  const tdouble h1 = (dof+1)/2;
  return h1*(flxdigamma(h1)-flxdigamma(dof/2)) + log(sqrt(dof)*BetaFun(dof/2,ONE/2));
}

const tdouble RBRV_entry_RV_StudentsT::get_sd_current_config()
{
  get_pars();
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
  get_pars();
  return (ONE-p)/2;
}


RBRV_entry_RV_StudentsT_generalized::RBRV_entry_RV_StudentsT_generalized(const std::string& name, const tuint iID, FlxFunction* nu, FlxFunction* locf, FlxFunction* scalef)
:RBRV_entry_RV_base(name,iID), nu(nu),locf(locf),scalef(scalef), dof(ZERO), loc(ZERO), scale(ZERO)
{

}

RBRV_entry_RV_StudentsT_generalized::RBRV_entry_RV_StudentsT_generalized(const std::string& name, const tuint iID, py::dict config)
: RBRV_entry_RV_base(name,iID), nu(nullptr), locf(nullptr), scalef(nullptr), dof(ZERO), loc(ZERO), scale(ZERO)
{
  try {
    if (config.contains("dof")) {          // dof, loc, scale
      nu = parse_py_para("dof", config);
      locf = parse_py_para("loc", config);
      scalef = parse_py_para("scale", config);
    } else {
      throw FlxException_NeglectInInteractive("RBRV_entry_RV_StudentsT_generalized::RBRV_entry_RV_StudentsT_generalized_01", "Required parameters to define distribution not found in Python <dict>.");
    }
  } catch (FlxException& e) {
    FLXMSG("RBRV_entry_RV_StudentsT_generalized::RBRV_entry_RV_StudentsT_generalized_99",1);
    if (nu) delete nu;
    if (locf) delete locf;
    if (scalef) delete scalef;
    throw;
  }
}

RBRV_entry_RV_StudentsT_generalized::~RBRV_entry_RV_StudentsT_generalized()
{
  if (nu) delete nu;
  if (locf) delete locf;
  if (scalef) delete scalef;
}

void RBRV_entry_RV_StudentsT_generalized::get_pars()
{
  dof = nu->cast2positive();
  loc = locf->calc();
  scale = scalef->cast2positive();
}

const tdouble RBRV_entry_RV_StudentsT_generalized::transform_y2x(const tdouble y_val)
{
  get_pars();
  if (y_val<=ZERO) {
    return loc + rv_InvCDF_Studentst(dof,rv_Phi(y_val))*scale;
  } else {
    return loc - rv_InvCDF_Studentst(dof,rv_Phi(-y_val))*scale;
  }
}

const tdouble RBRV_entry_RV_StudentsT_generalized::transform_x2y(const tdouble& x_val)
{
  get_pars();
  const tdouble x_ = (x_val-loc)/scale;
  if (x_<=ZERO) {
    return rv_InvPhi_noAlert(rv_cdf_Studentst(dof,x_));
  } else {
    return -rv_InvPhi_noAlert(rv_cdf_Studentst(dof,-x_));
  }
}

const tdouble RBRV_entry_RV_StudentsT_generalized::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
  const tdouble x_ = (x_val-loc)/scale;
  return rv_pdf_Studentst(dof,x_)/scale;
}

const tdouble RBRV_entry_RV_StudentsT_generalized::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
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
  get_pars();
  return loc;
}

const tdouble RBRV_entry_RV_StudentsT_generalized::get_sd_current_config()
{
  get_pars();
  if (dof<2) return std::numeric_limits<tdouble>::infinity();
  else return sqrt(dof/(dof-2))*scale;
}

const tdouble RBRV_entry_RV_StudentsT_generalized::get_median_current_config()
{
  get_pars();
  return loc;
}

const tdouble RBRV_entry_RV_StudentsT_generalized::get_mode_current_config()
{
  get_pars();
  return loc;
}

const bool RBRV_entry_RV_StudentsT_generalized::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || nu->search_circref(fcr) || locf->search_circref(fcr) || scalef->search_circref(fcr);
}

const tdouble RBRV_entry_RV_StudentsT_generalized::get_HPD(const tdouble p)
{
  get_pars();
  return (ONE-p)/2;
}


RBRV_entry_RV_Laplace::RBRV_entry_RV_Laplace(const std::string& name, const tuint iID, FlxFunction* locf, FlxFunction* scalef)
:RBRV_entry_RV_base(name,iID), locf(locf),scalef(scalef), loc(ZERO), scale(ZERO)
{

}

RBRV_entry_RV_Laplace::~RBRV_entry_RV_Laplace()
{
  if (locf) delete locf;
  if (scalef) delete scalef;
}

void RBRV_entry_RV_Laplace::get_pars()
{
  loc = locf->calc();
  scale = scalef->cast2positive();
}

const tdouble RBRV_entry_RV_Laplace::transform_y2x(const tdouble y_val)
{
  get_pars();
  if (y_val<=ZERO) {
    return loc + log(2*rv_Phi(y_val))*scale;
  } else {
    return loc - log(2*rv_Phi(-y_val))*scale;
  }
}

const tdouble RBRV_entry_RV_Laplace::transform_x2y(const tdouble& x_val)
{
  get_pars();
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
  get_pars();
  const tdouble x_ = fabs(x_val-loc)/scale;
  return exp(-x_)/(2*scale);
}

const tdouble RBRV_entry_RV_Laplace::calc_pdf_x_log(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
  const tdouble x_ = fabs(x_val-loc)/scale;
  return -x_ - log(2*scale);
}

const tdouble RBRV_entry_RV_Laplace::calc_cdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
  if (x_val<=loc) {
    return exp((x_val-loc)/scale)/2;
  } else {
    return ONE - exp((loc-x_val)/scale)/2;
  }
}

const tdouble RBRV_entry_RV_Laplace::calc_sf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
  if (x_val<=loc) {
    return ONE - exp((x_val-loc)/scale)/2;
  } else {
    return exp((loc-x_val)/scale)/2;
  }
}

const tdouble RBRV_entry_RV_Laplace::calc_entropy()
{
  get_pars();
  return ONE+log(2*scale);
}

const tdouble RBRV_entry_RV_Laplace::get_mean_current_config()
{
  get_pars();
  return loc;
}

const tdouble RBRV_entry_RV_Laplace::get_sd_current_config()
{
  get_pars();
  return sqrt(2.)*scale;
}

const tdouble RBRV_entry_RV_Laplace::get_median_current_config()
{
  get_pars();
  return loc;
}

const tdouble RBRV_entry_RV_Laplace::get_mode_current_config()
{
  get_pars();
  return loc;
}

const bool RBRV_entry_RV_Laplace::search_circref(FlxFunction* fcr)
{
  const bool b = (corr_valF==NULL)?false:(corr_valF->search_circref(fcr));
  return b || locf->search_circref(fcr) || scalef->search_circref(fcr);
}

const tdouble RBRV_entry_RV_Laplace::get_HPD(const tdouble p)
{
  get_pars();
  return (ONE-p)/2;
}

RBRV_entry_RV_UserTransform::RBRV_entry_RV_UserTransform(const std::string& name, const tuint iID, const bool is_z2x, FlxFunction* t1, FlxFunction* t2, FlxFunction* dh, FlxFunction* checkXf, RBRV_entry_RV_base* rv_z, const bool manage_z)
: RBRV_entry_RV_base(name,iID), is_z2x(is_z2x), t1(t1), t2(t2), dh(dh), checkXf(checkXf), rv_z(rv_z), manage_z(manage_z),
  tPL(1), tPLp(&tPL[0])
{
  
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

void RBRV_entry_RV_UserTransform::info(std::ostream& os)
{
  os << "user-transform distribution" << std::endl;
  if (is_z2x) {
    if (t1)      { os << "  z2x     = " << t1->write() << std::endl; }
    if (t2)      { os << "  x2z     = " << t2->write() << std::endl; }
    if (dh)      { os << "  dx2z    = " << dh->write() << std::endl; }
    //if (checkXf) { os << "  checkx  = " << checkXf->write() << std::endl; }
    os << "  Distribution of Z (" << rv_z->get_type() << "):" << std::endl;
    rv_z->info(os);
  } else {
    throw FlxException_NotImplemented("RBRV_entry_RV_UserTransform::info");
  }

}

void RBRV_entry_RV_UserTransform::replace_rv_z(RBRV_entry_RV_base* rv_z_)
{
  if (manage_z) throw FlxException_Crude("RBRV_entry_RV_UserTransform::replace_rv_z");
  rv_z = rv_z_;
}

RBRV_entry_RV_Truncated::RBRV_entry_RV_Truncated(const std::string& name, const tuint iID, FlxFunction* a, FlxFunction* b, RBRV_entry_RV_base* rv_z, const bool manage_z)
: RBRV_entry_RV_base(name,iID), a(a), b(b), rv_z(rv_z), manage_z(manage_z), aV(ZERO), bV(ZERO), q(ONE), aV_cdf(ZERO)
{

}

RBRV_entry_RV_Truncated::~RBRV_entry_RV_Truncated()
{
  if (a) delete a;
  if (b) delete b;
  if (manage_z) delete rv_z;
}

void RBRV_entry_RV_Truncated::get_pars()
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
  get_pars();
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
  get_pars();
  return rv_InvPhi_noAlert(calc_cdf_x(x_val,false));
}

const tdouble RBRV_entry_RV_Truncated::calc_pdf_x(const tdouble& x_val, const bool safeCalc)
{
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
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

void RBRV_entry_RV_Truncated::info(std::ostream& os)
{
  get_pars();
  os << "truncated distribution" << std::endl;
  os << "  lower   = " << GlobalVar.Double2String(aV) << std::endl;
  os << "  upper   = " << GlobalVar.Double2String(bV) << std::endl;
  os << "  q       = " << GlobalVar.Double2String(q) << std::endl;
  os << "  aV_cdf  = " << GlobalVar.Double2String(aV_cdf) << std::endl;
  os << "  Distribution of Z (" << rv_z->get_type() << "):" << std::endl;
  rv_z->info(os);
}

RBRV_entry_RV_maxminTransform::RBRV_entry_RV_maxminTransform(const std::string& name, const tuint iID, const bool is_max, FlxFunction* n, RBRV_entry_RV_base* rv_z)
: RBRV_entry_RV_base(name,iID), is_max(is_max), n(n), rv_z(rv_z), nV(ZERO)
{

}

RBRV_entry_RV_maxminTransform::~RBRV_entry_RV_maxminTransform()
{
  if (n) delete n;
  delete rv_z;
}

void RBRV_entry_RV_maxminTransform::get_pars()
{
  nV = n->cast2positive();
}

const tdouble RBRV_entry_RV_maxminTransform::transform_y2x(const tdouble y_val)
{
  get_pars();
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
  get_pars();
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
  get_pars();
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
  get_pars();
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





