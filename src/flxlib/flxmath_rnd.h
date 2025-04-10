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

#include "flxVec.h"

// #include <boost/random.hpp>
// typedef boost::random::mt19937 rng_type;

#include <random>
typedef std::mt19937 rng_type;


/**
* @brief initializes the random number generator
*/
FLXLIB_EXPORT void rv_initialize(bool startup, bool user_seed=false, tuint user_seed_int = 0, tuint user_init_calls = 0, rng_type* rngp=NULL, const bool do_output=true);

rng_type& get_rng();

/**
* @brief realization of a random variable with a uniform distribution [0;1]
*/
FLXLIB_EXPORT const tdouble rv_uniform();
FLXLIB_EXPORT const tdouble rv_uniform(rng_type& rng);

/**
* @brief realization of a random variable with a standard normal distribution
*/
FLXLIB_EXPORT const tdouble rv_normal();
FLXLIB_EXPORT const tdouble rv_normal(rng_type& rng);
FLXLIB_EXPORT void rv_normal(flxVec& y);
FLXLIB_EXPORT void rv_normal(flxVec& y,rng_type& rng);

/**
* @brief realization of a random variable with a lognormal distribution
*/
inline const tdouble rv_lognormal(tdouble &lambda, tdouble &xi, tdouble &epsilon) {
  return exp(rv_normal()*xi+lambda)+epsilon;
}

/**
* @brief probability density function of the standard normal distribution: phi
*/
inline const tdouble rv_phi(const tdouble &y) { return exp(-0.5*y*y) / sqrt(2*PI) ;  }
inline const tdouble rv_phi_log(const tdouble &y) { return (-y*y - log(2*PI))/2; }

/**
* @brief cumulative distribution function of the standard normal distribution: Phi
*/
FLXLIB_EXPORT const tdouble rv_Phi(const tdouble& y);
/**
* @brief returns rv_Phi(y_b)-rv_Phi(y_a)
*/
FLXLIB_EXPORT const tdouble rv_Phi_diff(const tdouble& y_a,const tdouble& y_b);

/**
* @brief inverse of the cumulative distribution function of the standard normal distribution: Phi^-1
*/
FLXLIB_EXPORT const tdouble rv_InvPhi(const tdouble& p);
/**
* @brief inverse of the cumulative distribution function of the standard normal distribution: Phi^-1
* ... no alerts
*/
FLXLIB_EXPORT const tdouble rv_InvPhi_noAlert(const tdouble& p);

// -----------------------------------------------------------------------------------------------------

/**
* @brief cdf of the Chi Square distribution
*/
const tdouble rv_cdf_ChiSquare(const tuint& v, const tdouble& x);

/**
* @brief pdf of the Chi Square distribution
*/
const tdouble rv_pdf_ChiSquare(const tuint& v, const tdouble& x);

/**
* @brief inverse of the cdf of the Chi Square distribution
*/
const tdouble rv_InvCDF_ChiSquare(const tuint& v, const tdouble& p);

// -----------------------------------------------------------------------------------------------------

/**
* @brief cdf of the Chi Square distribution
*/
const tdouble rv_cdf_Studentst(const tdouble& v, const tdouble& x);

/**
* @brief pdf of the Chi Square distribution
*/
const tdouble rv_pdf_Studentst(const tdouble& v, const tdouble& x);

/**
* @brief inverse of the cdf of the Chi Square distribution
*/
const tdouble rv_InvCDF_Studentst(const tdouble& v, const tdouble& p);

// -----------------------------------------------------------------------------------------------------

/**
* @brief cdf of the Binomial distribution with parameters p, N at k
*/
const tdouble rv_cdf_Binomial(const tdouble p, const tuint N, const tuint k);

// -----------------------------------------------------------------------------------------------------

/**
* @brief PMF of the beta-binomial distribution (log-transform)
*/
const tdouble rv_pmf_beta_binomial_ln(const tuint M, const tuint N, const tdouble alpha, const tdouble beta);


