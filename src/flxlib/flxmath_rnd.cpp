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

#include "flxmath_rnd.h"

#include <boost/math/distributions.hpp>
#include <boost/math/special_functions.hpp>

rng_type randgen;

std::normal_distribution<tdouble> pd_normal(ZERO, ONE);
std::uniform_real_distribution<tdouble> pd_uniform(ZERO, ONE);

boost::math::normal ndist(ZERO,ONE);

void rv_initialize(bool startup, bool user_seed, tuint user_seed_int, tuint user_init_calls, rng_type* rngp, const bool do_output)
{
  if (rngp==NULL) rngp = &randgen;
  if (startup) {
    user_seed = GlobalVar.MT19937_init_seed;
    user_seed_int = GlobalVar.MT19937_init_seedvalue;
    user_init_calls = GlobalVar.MT19937_init_calls;
  } 
  
  if (user_seed) {
    if (do_output) GlobalVar.slogcout(3) << "  Random Number Generator: MT19937 - initialized with user seed (" << user_seed_int << ")" << std::endl;
    rngp->seed(static_cast<unsigned int>(user_seed_int));
  } else {
    if (GlobalVar.MT19937_init_RAND) {
      tuint rV = static_cast<unsigned int>(rand());
      if (do_output) GlobalVar.slogcout(3) << "Random Number Generator: MT19937 - initialized with rand()=" << rV << ";" << std::endl;
      rngp->seed(rV);
    } else {
      if (do_output) GlobalVar.slogcout(3) << "Random Number Generator: MT19937 - initialized with time (" << std::time(0) << ")" << std::endl;
      rngp->seed(static_cast<unsigned int>(std::time(0)));
    }
  }

  if (do_output) GlobalVar.slogcout(3) << "Random Number Generator: MT19937 - initialized with " << user_init_calls << " initial calls." << std::endl;
  for (tuint i = 0; i < user_init_calls; ++i) {
    rv_normal(*rngp);
  }
}

rng_type& get_rng()
{
  return randgen;
}

const tdouble rv_uniform()
{
  return pd_uniform(randgen);
}

const tdouble rv_uniform(rng_type& rng)
{
  return pd_uniform(rng);
}

const tdouble rv_normal()
{
  return pd_normal(randgen);

}

const tdouble rv_normal(rng_type& rng)
{
    return pd_normal(rng);
}

void rv_normal(flxVec& y)
{
  rv_normal(y,randgen);
}

void rv_normal(flxVec& y,rng_type& rng)
{
  const tuint S = y.get_N();
  tdouble* yp = y.get_tmp_vptr();
  for (tuint i = 0; i<S;++i) {
    yp[i] = pd_normal(rng);
  }
}


const tdouble rv_Phi(const tdouble& y)
{
  #if FLX_DEBUG
    #ifdef __unix__
      if (std::isnan(y)) {
        throw FlxException_Crude("rv_Phi_a");
      }
    #endif
  #endif
  return boost::math::cdf(ndist,y);
}

const tdouble rv_Phi_diff(const tdouble& y_a, const tdouble& y_b)
{
  #if FLX_DEBUG
    if (y_b<y_a) {
      throw FlxException_Crude("rv_Phi_b");
    }
  #endif
  const tdouble p_a = rv_Phi(y_a);
  const tdouble p_b = rv_Phi(y_b);
  if (p_b<ONE/2) {
    return p_b-p_a;
  }
  const tdouble q_a = rv_Phi(-y_a);
  const tdouble q_b = rv_Phi(-y_b);
  if (p_a>ONE/2) {
    return q_a-q_b;
  }
  return ((p_b-(q_b+p_a))+q_a)/2;
}

const tdouble rv_InvPhi(const tdouble& p)
{
  try {
    return boost::math::quantile( ndist,p );
  } catch ( std::exception const& e ) {
    FLXMSG("rv_InvPhi_1",1);
    GlobalVar.alert.alert("rv_InvPhi_2", "  for p="+GlobalVar.Double2String(p)+"\n"+e.what() );
    if (p < ONE/2) { return -Y_INFTY_1; }
    else { return Y_INFTY_1; }
  }
}

const tdouble rv_InvPhi_noAlert(const tdouble& p)
{
  try {
    if (p<=ZERO) return -Y_INFTY_2;
    else if (p>=ONE) return Y_INFTY_2;
    return boost::math::quantile( ndist,p );
  } catch ( std::exception const& e ) {
    FLXMSG("rv_InvPhi_noAlert",5);
    if (p < ONE/2) { 
      return -Y_INFTY_2;
    }
    else { return Y_INFTY_2; }
  }
}

// -----------------------------------------------------------------------------------------------------

const tdouble rv_cdf_ChiSquare(const tuint& v, const tdouble& x)
{
  return (x<=ZERO?ZERO:boost::math::gamma_p(v/2.0, x/2.0));
}

const tdouble rv_pdf_ChiSquare(const tuint& v, const tdouble& x)
{
  return (x<=ZERO?ZERO:
    (pow(x,v/(2*ONE)-1)*exp(-x/(2*ONE)))/(pow(2*ONE,v/(2*ONE))*boost::math::tgamma(v/2,0))
  );
}

const tdouble rv_InvCDF_ChiSquare(const tuint& v, const tdouble& p)
{
  return 2*ONE*boost::math::gamma_p_inv(v/(2*ONE), p);
}

// -----------------------------------------------------------------------------------------------------

const tdouble rv_cdf_Studentst(const tdouble& v, const tdouble& x)
{
  boost::math::students_t dist(v);
  return boost::math::cdf(dist,x);
}

const tdouble rv_pdf_Studentst(const tdouble& v, const tdouble& x)
{
  boost::math::students_t dist(v);
  return boost::math::pdf(dist,x);
}

const tdouble rv_InvCDF_Studentst(const tdouble& v, const tdouble& p)
{
  try {
    boost::math::students_t dist(v);
    return boost::math::quantile(dist,p);
  } catch ( std::exception const& e ) {
    FLXMSG("rv_InvCDF_Studentst_01",1);
    GlobalVar.alert.alert("rv_InvCDF_Studentst_02", "  for p="+GlobalVar.Double2String(p)+"\n"+e.what() );
    if (p < ONE/2) { return -Y_INFTY_1; }
    else { return Y_INFTY_1; }
  }
}


// -----------------------------------------------------------------------------------------------------

const tdouble rv_cdf_Binomial(const tdouble p, const tuint N, const tuint k)
{
  return boost::math::cdf(boost::math::binomial(N,p),k);
}

// -----------------------------------------------------------------------------------------------------

const tdouble rv_pmf_beta_binomial_ln(const tuint M, const tuint N, const tdouble alpha, const tdouble beta)
{
  const tdouble t1 = BinomCoeff_ln(N,M);

  // VARIANT 1: (numerically instable)
  // const tdouble t2 = BetaFunLn(alpha+M,beta+(N-M));
  // const tdouble t3 = BetaFunLn(alpha,beta);
  // const tdouble res = t1 + (t2 - t3);

  // VARIANT 2: (numerically instable)
  // const tdouble t23_ = GammaLn(alpha+M)-GammaLn(alpha)
  //    +GammaLn(beta+N-M)-GammaLn(beta+N-M+alpha+M)
  //    +GammaLn(beta+alpha)-GammaLn(beta);
  // const tdouble res = t1 + t23_;

  // VARIANT 3: (appears to be robust)
  const tdouble t23 = GammaLn_diff(M,alpha)+GammaLn_diff(alpha,beta)-GammaLn_diff(alpha+M,beta+N-M);
  const tdouble res = t1 + t23;

  // TODO ... ideally, we should avoid res>ZERO ... however, this seems to be numercially challenging
  // if (res>ZERO) {
  //   GlobalVar.slogcout(1) << "DEBUG_rv_pmf_beta_binomial_ln " << res << "  " << t1 << "  " << t23 << "  " << M << "  " << N << "  " << alpha << "  " << beta << "  " << alpha/(alpha+beta) << std::endl;
  //   throw FlxException_Crude("rv_pmf_beta_binomial_ln");
  // }
  return res;
}

