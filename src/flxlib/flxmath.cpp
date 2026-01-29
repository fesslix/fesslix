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
#include "flxmath_rnd.h"

#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>


/**
* @brief Array Size of FactorialLn-Function (stores calculated values)
*/
const int FactorialLnArraySize = 100;
/**
* @brief Array Size of Factorial-Function (stores calculated values up to !(FactorialArraySize+1) )
*/
const int FactorialArraySize = 33;
/**
* @brief Array Size of DoubleFactorial-Function (stores calculated values up to !(DoubleFactorialArraySize+1) )
*/
const int DoubleFactorialArraySize = 33;


const tdouble GammaLn(const tdouble &xx)
{
  static const tdouble cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  tdouble y=xx;
  const tdouble &x=xx;
  tdouble tmp=x+5.5;
  tmp-=(x+0.5)*std::log(tmp);
  tdouble ser=1.000000000190015;
  for (int j=0;j<6;j++) ser+=cof[j]/++y;
  return -tmp+std::log(2.5066282746310005*ser/x);
}
// {
//   return boost::math::lgamma(xx);
// }

const tdouble GammaLn_diff(const tdouble alpha, const tdouble beta)
{
  static const tdouble cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  // for beta
  tdouble tmp_beta = beta+5.5;
  tmp_beta -= (beta+0.5)*std::log(tmp_beta);
  tdouble beta_ser=1.000000000190015;
  tdouble y_beta = beta;
  for (int j=0;j<6;j++) beta_ser+=cof[j]/++y_beta;
  // for alpha + beta
  const tdouble log_ab = std::log(5.5+alpha+beta);
  const tdouble tmp_1 = 5.5+alpha-(0.5+alpha)*log_ab;
  const tdouble tmp_2 = beta*(ONE-log_ab);
  tdouble ab_ser=1.000000000190015;
  tdouble y_ab = alpha+beta;
  for (int j=0;j<6;j++) ab_ser+=cof[j]/++y_ab;
  return ((std::log(2.5066282746310005*ab_ser/(alpha+beta)) - std::log(2.5066282746310005*beta_ser/beta)) - tmp_1) + (tmp_beta-tmp_2);
}

const tdouble BetaFunLn(const tdouble &z, const tdouble &w) {
  // VARIANT: to avoid numerical problems for either very large z or w
  const tdouble v1 = max(z,w);
  const tdouble v2 = min(z,w);
  return GammaLn(v2)+(GammaLn(v1)-GammaLn(v1+v2));
  // ORIGINAL variant
  // return GammaLn(z)+GammaLn(w)-GammaLn(z+w);
}

const tdouble flxgamma(const tdouble x)
{
  return tgamma(x);
//   return boost::math::tgamma(x);
}

const tdouble flxgamma(const tdouble a, const tdouble z)
{
  return boost::math::tgamma(a,z);
}

const tdouble flxgamma_l(const tdouble a, const tdouble z)
{
  return boost::math::tgamma_lower(a,z);
}

const tdouble flxgamma_ru(const tdouble a, const tdouble z)
{
  return boost::math::gamma_q(a,z);
}

const tdouble flxgamma_rl(const tdouble a, const tdouble z)
{
  return boost::math::gamma_p(a,z);
}

const tdouble flxgamma_ru_inv(const tdouble a, const tdouble q)
{
  return boost::math::gamma_q_inv(a,q);
}

const tdouble flxgamma_rl_inv(const tdouble a, const tdouble p)
{
  return boost::math::gamma_p_inv(a,p);
}

const tdouble flxerf_inv(const tdouble p)
{
  return boost::math::erf_inv(p);
}

const tdouble iBeta_reg(const tdouble alpha, const tdouble beta, const tdouble x)
{
  return boost::math::ibeta(alpha,beta,x);
}

const tdouble iBeta_regc(const tdouble alpha, const tdouble beta, const tdouble x)
{
  return boost::math::ibetac(alpha,beta,x);
}

const tdouble iBeta_reg_inv(const tdouble alpha, const tdouble beta, const tdouble p)
{
  // check for approximate solution to avoid iteration errors
  const tdouble mu = alpha/(alpha+beta);
  const tdouble sd = sqrt(alpha*beta/(pow2(alpha+beta)*(alpha+beta+ONE)));
  if (alpha+beta>1e6 && sd/mu<1e-4) {
    return rv_InvPhi(p)*sd+mu;
  }
  try {
    const tdouble res = boost::math::ibeta_inv(alpha,beta,p);
    if (std::isnan(res)) {
      // shomehow this error occurs for small rv_Phi in combination with Valgrind
      using boost::math::policies::make_policy;
      using boost::math::policies::pole_error;
      using boost::math::policies::domain_error;
      using boost::math::policies::overflow_error;
      using boost::math::policies::evaluation_error;
      using boost::math::policies::throw_on_error;
      using boost::math::policies::digits10;
      const tdouble res_2 = boost::math::ibeta_inv(alpha,beta,p,make_policy(
            domain_error<throw_on_error>(),
            pole_error<throw_on_error>(),
            overflow_error<throw_on_error>(),
            evaluation_error<throw_on_error>(),
            digits10<20>()
         ));
      std::ostringstream ssV;
      ssV << "alpha=" << GlobalVar.Double2String(alpha) << std::endl
          << "beta="  << GlobalVar.Double2String(beta) << std::endl
          << "p=" << GlobalVar.Double2String(p) << std::endl
          << "res_2=" << GlobalVar.Double2String(res_2) << std::endl ;
      throw FlxException("iBeta_reg_inv_01",ssV.str());
    }
    return res;
  } catch (const boost::math::evaluation_error& e) {
      *(GlobalVar.get_true_cerr()) << std::endl << "WARNING: iBeta_reg_inv" << std::endl;
      throw FlxException("iBeta_reg_inv_02");
  }
}

const tdouble iBetac_reg_inv(const tdouble alpha, const tdouble beta, const tdouble p)
{
  return boost::math::ibetac_inv(alpha,beta,p);
}

const tdouble flxdigamma(const tdouble x)
{
  return boost::math::digamma(x);
}

const tdouble Factorial(const int n) {
  static int ntop=4;
  static tdouble a[FactorialArraySize]={2.0,6.0,24.0};
  if (n<=1) return ONE;
  else if (n>FactorialArraySize+1) return exp(GammaLn(tdouble(n)+ONE));
  while (ntop<n) {
    const int j=ntop++;
    a[ntop-2]=a[j-2]*ntop;
  }
  return a[n-2];
}

const tdouble DoubleFactorial(const int n)
{
  static int ntop=5;
  static tdouble a[FactorialArraySize]={2.0,3.0,8.0,15.0};
  if (n<=1) return ONE;
  else if (n>FactorialArraySize+1) {
    if (n%2 == 0) {
      const int k = n/2;
      return pow(2.,k)*Factorial(k);                        // TODO efficiency of pow?
    } else {
      const int k = (n+1)/2;
      return Factorial(2*k)/(pow(2.,k)*Factorial(k));        // TODO efficiency of pow?
    }
  }
  else {
    while (ntop<n) {
      int j=ntop++;
      a[ntop-2]=a[j-3]*ntop;
    }
    return a[n-2];
  }
}

const tdouble FactorialLn(const int n) {
  static tdouble a[FactorialLnArraySize];
  if (n<=1) return ZERO; 
  if (n<=FactorialLnArraySize+1) return (a[n-2]!=ZERO?a[n-2]:(a[n-2]=GammaLn(tdouble(n)+ONE)));
  else return GammaLn(tdouble(n)+ONE);
}


// ------------------------------------ optimization routines ----------------------------------


void flx_optim(tdouble x_lb, tdouble x_ub, tdouble& x_min, flx_math_fun_ptr fun, void* dp, const bool use_brent, const bool use_initial, const tuint NITER, const tuint NEX, const tdouble eps1, const tdouble eps2, ostreamp osp)
{
  try {
  // check the intervals
    if (x_ub<=x_lb) {
      std::ostringstream ssV;
      ssV << "Error with optimization boundaries.";
      throw FlxException("flx_optim_01", ssV.str() );
    }
    if (eps1 < GlobalVar.TOL()) {
      GlobalVar.alert.alert("flx_optim_02a","The 'eps1' threshold might be too restrictive.");
    }
    if (eps2 < GlobalVar.TOL()) {
      GlobalVar.alert.alert("flx_optim_02b","The 'eps2' threshold might be too restrictive.");
    }

  // initialize variables
    if (osp) {
      *osp << "  =====================================" << std::endl;
      *osp << "  1D-optimization (minimization)" << std::endl;
      *osp << "  =====================================" << std::endl;
    }
    tdouble fa=fun(x_lb,dp);
    if (fa!=fa) throw FlxException("flx_optim_nan_01a", "NaN not allowed as function result.");
    tdouble fb = -1e9;                // dummy initialization
    tdouble fc=fun(x_ub,dp);
    if (fc!=fc) throw FlxException("flx_optim_nan_01b", "NaN not allowed as function result.");
    bool do_bracket = true;
    if (use_initial) {
      if (x_min<x_lb || x_min>x_ub) {
        throw FlxException("flx_optim_02", "The 'x_min' must be within ['x_lb','x_ub']");
      }
      fb = fun(x_min,dp);
      if (fb!=fb) throw FlxException("flx_optim_nan_03a", "NaN not allowed as function result.");
      if (fa==fb && fb==fc) {  // TODO: maybe include the case when values are near-constant
        throw FlxException("flx_optim_nan_03b", "Function appears to be constant.");
      }
      do_bracket = !((fa-fb)>ZERO&&(fc-fb)>ZERO);
    }
    if (osp) {
      *osp << "  fa=f(" << x_lb << ")=" << fa << std::endl;
      if (use_initial) *osp << "  fb=f(" << x_min << ")=" << fb << std::endl;
      *osp << "  fc=f(" << x_ub << ")=" << fc << std::endl;
    }
  // try to bracket the minimum: preliminary search (Brent's exponential search)
    if (do_bracket) {
      if (osp) *osp << "  bracketing the mimimum" << std::endl;
      const tdouble phi = (ONE + sqrt(tdouble(5.0))) / (2*ONE);
      const tdouble resphi = 2*ONE - phi;
      const tdouble phih = (2*ONE)/(sqrt(tdouble(5.0))-ONE);
      if (!use_initial) {
        x_min = x_lb + resphi*(x_ub-x_lb);
        fb = fun(x_min,dp);
        if (fb!=fb) throw FlxException("flx_optim_nan_04", "NaN not allowed as function result.");
        if (osp) *osp << "  fb=f(" << x_min << ")=" << fb << std::endl;
      }
      if (fa > fb && fb > fc) {
        for (tuint i=0;(i<NEX)&&(fa > fb && fb > fc);++i) {
          const tdouble x = x_ub + phih*(x_ub-x_min);
          x_lb = x_min;
          x_min = x_ub;
          x_ub = x;
          fa = fb;
          fb = fc;
          fc = fun(x_ub,dp);
          if (fc!=fc) throw FlxException("flx_optim_nan_05", "NaN not allowed as function result.");
          if (osp) *osp << "    " << i << ": fc=f(" << x_ub << ")=" << fc << std::endl;
        }
      } else if (fa < fb && fb < fc) {
        for (tuint i=0;(i<NEX)&&(fa < fb && fb < fc);++i) {
          const tdouble x = x_lb - phih*(x_min-x_lb);
          x_ub = x_min;
          x_min = x_lb;
          x_lb = x;
          fc = fb;
          fb = fa;
          fa = fun(x_lb,dp);
          if (fa!=fa) throw FlxException("flx_optim_nan_06", "NaN not allowed as function result.");
          if (osp) *osp << "    " << i << ": fa=f(" << x_lb << ")=" << fa << std::endl;
        }
      } else if (fa < fb && fb > fc) {
        std::ostringstream ssV;
        ssV << "There appears to be a maximum and not a minimum within the chosen interval.";
        throw FlxException("flx_optim_04", ssV.str() );
      } else if (fa == fb && fb == fc) {
        std::ostringstream ssV;
        ssV << "The function seems to be constant on the specified interval.";
        throw FlxException("flx_optim_05", ssV.str() );
      }
      if (!((fa-fb)>ZERO&&(fc-fb)>ZERO)) {
        std::ostringstream ssV;
        ssV << "Maximum number of iterations (" << NEX << ") reached when bracketing the minimum." << std::endl;
        ssV << "        fa(" << x_lb << ")=" << fa << ", fb(" << x_min << ")=" << fb << ", fc(" << x_ub << ")=" << fc;
        throw FlxException("flx_optim_06", ssV.str() );
      }
    }
  // employ the GSL-library
    // initialize the function to optimize
      gsl_function F;
      F.function = fun;
      F.params = dp;
    const gsl_min_fminimizer_type *T = use_brent?gsl_min_fminimizer_brent:gsl_min_fminimizer_goldensection;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (T);
    int status = gsl_min_fminimizer_set_with_values (s, &F, x_min,fb, x_lb,fa, x_ub,fc);
    if (status!=GSL_SUCCESS) throw FlxException_Crude("flx_optim_07");
    int iter = 0, max_iter = 100;
    // initial output
      if (osp) {
        *osp << "  using " << gsl_min_fminimizer_name (s) << " method" << std::endl;
        *osp << "  iter [lower, upper] xmin fmin err(est)" << std::endl;
        *osp << "  " << iter << "\t" << x_lb << "\t" << x_ub << "\t" << x_min << "\t" << x_ub-x_lb << std::endl;
      }
    do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);
      x_min = gsl_min_fminimizer_x_minimum (s);
      x_lb = gsl_min_fminimizer_x_lower (s);
      x_ub = gsl_min_fminimizer_x_upper (s);
      status = gsl_min_test_interval (x_lb, x_ub, eps1, eps2);
      if (osp) {
        *osp << "  " << iter << "\t" << x_lb << "\t" << x_ub << "\t" << x_min << "\t" << x_ub-x_lb << std::endl;
      }
      if (status == GSL_SUCCESS) {
        if (osp) {
          *osp << "  Optimization converged" << std::endl;
        }
        gsl_min_fminimizer_free (s);
        return;
      }
      if (iter>=max_iter) {
        gsl_min_fminimizer_free (s);
        throw FlxException("flx_optim_08", "Optimization did not converge in specified number of steps.");
      }
    } while (status == GSL_CONTINUE);
    gsl_min_fminimizer_free (s);
    throw FlxException("flx_optim_09", "A problem occurred in the optimization algorithm.");
  } catch (FlxException &e) {
    if (osp) *osp << std::endl << e.what();
    throw;
  }
}


// ------------------------------------ root search routines ----------------------------------

const tdouble flx_RootSearch_Bisec(flx_math_fun_ptr fun, void* dp, tdouble start, tdouble end, const tdouble dx, const tdouble dy, ostreamp streamp)
{
  const int fixW = 15;                // output-width for logging
  if (streamp) {
    *streamp << std::endl;
    *streamp << "*************************************" << std::endl;
    *streamp << "* root-search: bisection method     *" << std::endl;
    *streamp << "*************************************" << std::endl;
    *streamp << " dx = " << GlobalVar.Double2String(dx) << std::endl;
    *streamp << " dy = " << GlobalVar.Double2String(dy) << std::endl;    
  }
  // check the intervals
    if (end<=start) {
      std::ostringstream ssV;
      ssV << "Error with root-search boundaries.";
      throw FlxException("FlxFun_RootSearch_Bisec_1", ssV.str() );
    }
    if (dx<ZERO) {
      std::ostringstream ssV;
      ssV << "'dx' threshold must not be negativ";
      throw FlxException("FlxFun_RootSearch_Bisec_2", ssV.str() );
    }
    if (dy<ZERO) {
      std::ostringstream ssV;
      ssV << "'dy' threshold must not be negativ";
      throw FlxException("FlxFun_RootSearch_Bisec_3", ssV.str() );
    }
    if (dx < GlobalVar.TOL() && dy < GlobalVar.TOL()) {
      GlobalVar.alert.alert("FlxFun_RootSearch_Bisec_4","The 'dx' and 'dy' thresholds might be too restrictive.");
    }
  bool found = false;
  // start - configuration
    tdouble res = ZERO;  // dummy initialization to prevent compiler warning
    tdouble f1 = fun(start,dp);
    tdouble f2 = fun(end,dp);
    if (fabs(f1)<=GlobalVar.TOL()) {
      res = start;
      found = true;
    } else if (fabs(f2)<=GlobalVar.TOL()) {
      res = end;
      found = true;
    } else if (f2*f1>=ZERO) {
      std::ostringstream ssV;
      ssV << "Function has the same sign at both boundaries.";
      throw FlxException("FlxFun_RootSearch_Bisec_5", ssV.str() );
    }
  if (streamp) {
    *streamp << " f(" << GlobalVar.Double2String(start,false,-1,fixW) << ") = " << GlobalVar.Double2String(f1,false,-1,fixW) << std::endl;
    *streamp << " f(" << GlobalVar.Double2String(end,false,-1,fixW) << ") = "   << GlobalVar.Double2String(f2,false,-1,fixW) << std::endl;
  }
  if (!found) {
    while (end-start>dx) {
      res = (start+end)/2;
      const tdouble f3 = fun(res,dp);
      if (streamp) {
        *streamp << " f(" << GlobalVar.Double2String(res,false,-1,fixW) << ") = " << GlobalVar.Double2String(f3,false,-1,fixW) << std::endl;
      }
      if (fabs(f3)<=dy) {
        found = true;
        break;
      }
      if ( f1*f3>ZERO ) {
        start = res;
        f1 = f3;
      } else {
        end = res;
        f2 = f3;
      }
    }
  }
  if (!found) {
    res = (start+end)/2;
  }
  if (streamp) {
    if (!found) {
      *streamp << " ALERT: root could not be localized for the specified dy-error-bound." << std::endl;
    }
    *streamp << " result: " << GlobalVar.Double2String(res,false,-1,fixW) << std::endl;
    *streamp << "*************************************" << std::endl;
  }
  return res;
}

const tdouble flx_RootSearch_RegulaFalsi(flx_math_fun_ptr fun, void* dp, tdouble start, tdouble end, const tdouble dx, const tdouble dy, ostreamp streamp)
{
  const int fixW = 15;                // output-width for logging
  if (streamp) {
    *streamp << std::endl;
    *streamp << "*************************************" << std::endl;
    *streamp << "* root-search: regula falsi         *" << std::endl;
    *streamp << "*************************************" << std::endl;
    *streamp << " dx = " << GlobalVar.Double2String(dx) << std::endl;
    *streamp << " dy = " << GlobalVar.Double2String(dy) << std::endl;   
  }
  const tuint NMAX = tuint(1e4);
  // check the intervals
    if (end<=start) {
      std::ostringstream ssV;
      ssV << "Error with root-search boundaries.";
      throw FlxException("FlxFun_RootSearch_RegulaFalsi_1", ssV.str() );
    }
    if (dx<ZERO) {
      std::ostringstream ssV;
      ssV << "'dx' threshold must not be negativ";
      throw FlxException("FlxFun_RootSearch_RegulaFalsi_2", ssV.str() );
    }
    if (dy<ZERO) {
      std::ostringstream ssV;
      ssV << "'dy' threshold must not be negativ";
      throw FlxException("FlxFun_RootSearch_RegulaFalsi_3", ssV.str() );
    }
    if (dx < GlobalVar.TOL() && dy < GlobalVar.TOL()) {
      GlobalVar.alert.alert("FlxFun_RootSearch_RegulaFalsi_4","The 'dx' and 'dy' thresholds might be too restrictive.");
    }
  bool found = false;
  // start - configuration
    tdouble res = ZERO;
    tdouble f1 = fun(start,dp);
    tdouble f2 = fun(end,dp);
    if (streamp) {
      *streamp << " f(" << GlobalVar.Double2String(start,false,-1,fixW) << ") = " << GlobalVar.Double2String(f1,false,-1,fixW) << std::endl;
      *streamp << " f(" << GlobalVar.Double2String(end,false,-1,fixW) << ") = "   << GlobalVar.Double2String(f2,false,-1,fixW) << std::endl;
    }
    if (fabs(f1)<=GlobalVar.TOL()) {
      res = start;
      found = true;
    } else if (fabs(f2)<=GlobalVar.TOL()) {
      res = end;
      found = true;
    }
  if (!found) {
    tuint c = 0;
    while (fabs(end-start)>dx && c<NMAX) {
      if (fabs(f1-f2)<=GlobalVar.TOL()) {
        std::ostringstream ssV;
        ssV << "Function value is the same at both root-search boundaries: " << GlobalVar.Double2String(f1);
        throw FlxException("FlxFun_RootSearch_RegulaFalsi_5", ssV.str() );
      }
      ++c;
      tdouble f3;
      if (f1*f2<=ZERO) {        // Pegasus-algorithm
        res = (start*f2-end*f1)/(f2-f1);
        f3 = fun(res,dp);
        if (streamp) {
          *streamp << " f(" << GlobalVar.Double2String(res,false,-1,fixW) << ") = " << GlobalVar.Double2String(f3,false,-1,fixW) << " (Pegasus-algorithm)" << std::endl;
        }
        if ( f2*f3<ZERO ) {
          start = end;
          f1 = f2;
          end = res;
          f2 = f3;
        } else {
          const tdouble m = f2/(f2+f3);
          f1 = m*f1;
          end = res;
          f2 = f3;
        }
      } else {                        // Sekantenverfahren
        res = end - ((end-start)/(f2-f1))*f2;
        f3 = fun(res,dp);
        if (streamp) {
          *streamp << " f(" << GlobalVar.Double2String(res,false,-1,fixW) << ") = " << GlobalVar.Double2String(f3,false,-1,fixW) << " (Secant-algorithm)" << std::endl;
        }
        f1 = f2;
        start = end;
        f2 = f3;
        end = res;
      }
      if (fabs(f3)<=dy) {
        found = true;
        break;
      }
      #if FLX_DEBUG
        if (std::isnan(start) || std::isnan(end) ) throw FlxException_Crude("FlxFun_RootSearch_RegulaFalsi_6");
      #endif
    }
    if (c>=NMAX) {
      res = flx_RootSearch_Bisec(fun,dp,start,end,dx,dy,streamp);
      found = true;
    }
  }
  if (!found) {
    res = (start+end)/2;
  }
  if (streamp) {
    if (!found) {
      *streamp << " ALERT: root could not be localized for the specified error bounds." << std::endl;
    }
    *streamp << " result: " << GlobalVar.Double2String(res,false,-1,fixW) << std::endl;
    *streamp << "*************************************" << std::endl;
  }
  #if FLX_DEBUG
    if (std::isnan(res)) {
      throw FlxException_Crude("FlxFun_RootSearch_RegulaFalsi_7");
    }
  #endif
  return res;
}





