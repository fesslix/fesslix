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

#ifndef fesslix_flxgp_relmeth_H
#define fesslix_flxgp_relmeth_H

#include "flxgp_kernel.h"





// ------------------------------------------------------------------------------------------------
class flxGP_data_ptr;

class PYBIND11_EXPORT flxGP_MCI {
  protected:
    flxGPProj_base& gp;
    const tuint Ndim;             // number of basic random variables of the problem
    const tuint user_seed_int;    // seed value for the random number generator
    const tuint user_init_calls;  // number of initial calls to the random number generator
    const tdouble tqi_val;        // quantile value for evaluating stopping criterion (0.99 by default)
    std::vector<tdouble> dmV;     // training set: X
    std::vector<tdouble> doV;     // training set: Y
    std::vector<tdouble> mcs_pi;  // vector of 0-1 probabilities in the Monte Carlo simulation
    flxVec tqi_vec;                // vector needed to evaluate the variance of the tqi
    flxVec tqi_vec_rv_u;           // vector to store pseudo-random uniform samples -> NOT modified after initialization in flxGP_MCI::flxGP_MCI
    tulong id_next_point;          // id of the next sample for model evaluation
    std::stringstream logStream;

    tdouble static_sum;
    tdouble last_m;
    tdouble last_n;

    void init_RNG();

    /**
    * @brief tqi_val-quantile of posterior beta distribution for 'm' hits in 'n' samples (with uniform prior)
    *
    * i.e., this function returns the expected sampling uncertainty WITHOUT Kriging
    */
    const tdouble tqi_eval(const tdouble m, const tdouble n) const;
    const tdouble tqi_eval_covar() const;
    const tdouble tqi_eval_pr(const tdouble p) const;
    const tdouble get_mean_tqi(const tdouble ref_m, const tulong n, const tulong* skip_id=NULL, const tdouble nrep=ONE);

    void generate_sample(flxVec &uvec_);
  public:
    /**
    * @brief MCI with a GP-surrogate model
    * @param Ndim number of uncertain model parameters
    * @param Nreserve intented (maximum) number of runs of the 'actual' model
    */
    flxGP_MCI(flxGPProj_base& gp, const tuint Nreserve, const tuint user_seed_int, const tuint user_init_calls, const tdouble tqi_val);
    virtual ~flxGP_MCI() {}

    void assemble_lh_samples(flxVec& lh_samples);
    const bool is_point_unique(const flxVec& uvec_) const;
    void get_next_point(flxVec& uvec_);
    void register_sample(const tdouble lsfval, const flxVec& uvec_);
    void condition_on_data(const bool register_pvec, const bool learn_noise);
    void optimize_gp_para(const tuint iterMax);

    const tuint get_N_model_calls() const { return doV.size(); }
    const tuint get_Ndim() const { return Ndim; }

    py::dict simulate_GP_mci(const tulong Nsmpls, tdouble& err, int& proposed_action_id);
    void output_summary();

    friend class flxGP_data_ptr;
};

// ------------------------------------------------------------------------------------------------

class FLXLIB_EXPORT flxGP_data_ptr : public flxGP_data_base {
  private:
    flxGP_MCI *akmcs;
  public:
    flxGP_data_ptr(flxGP_MCI *akmcs) : akmcs(akmcs) { }
    virtual ~flxGP_data_ptr() {}

    virtual flxGP_data_base* copy() const;
    virtual const tdouble* get_data_ptr_mtx(tuint &N_obsv, const tuint N_dim);
    virtual const tdouble* get_data_ptr_vec(const tuint N_obsv);
};









#endif // fesslix_flxgp_relmeth_H

