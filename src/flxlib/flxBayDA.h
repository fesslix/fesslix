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
#include "flxmtxfun_data.h"
#include "flxrbrv_rvs.h"
#include "flxrbrv_rvs_read.h"
#include "flxBayUp.h"


class flxBayDA_likeli {
  protected:
    /**
    * @brief type of the distribution
    */
    tuint type;
    /**
    * @brief dimension of the parameter space
    */
    tuint M;
    /**
    * @brief parameter vector
    */
    flxVec* pvec;
    /**
    * @brief distribution object
    */
    RBRV_entry_RV_base* rv;
    /**
    * @brief distribution object (to work with externally)
    *
    * considers internal data-transformations
    */
    RBRV_entry_RV_base* rv_ext;
    /**
    * @brief pointer to data vector
    */
    flxVec* data_vec;
    /**
    * @brief defines the domain of the data (see flxBayDA)
    */
    tuint id_transform;
    tuint id_transform2;
    /**
    * @brief value for correcting the transformed likelihood
    */
    tdouble tcl;
    FlxRndCreator& RndCreator;
    /**
    * @brief Number of chains
    */
    const tuint Nchain;
    /**
    * @brief total number of posterior samples to generate
    */
    const tuint Npost;
    /**
    * @brief array for random shuffling of posterior samples
    */
    tuint* ashuffle;
    /**
    * @brief parameter vector for current chain positions
    */
    flxVec* pchain;
    /**
    * @brief parameter vector for posterior samples
    */
    flxVec* ppost;
    /**
    * @brief step vector (for MCMC step)
    */
    flxVec* pstep;
    /**
    * @brief vector of likelihood computed for each chain
    */
    flxVec* lcvec;
    /**
    * @brief number of posterior samples currently stored
    */
    tuint N;
    /**
    * @brief adapt the step-size after how many chain moves?
    */
    const tuint N_adapt;
    /**
    * @brief counter - when was the step adapted the last time
    */
    tuint N_adapt_;
    /**
    * @brief counter for adapting pstep
    */
    tuint* N_adapt_count;
    /**
    * @brief sample of the dimensions - order of sampling
    */
    tuint* dim_sample;
    /**
    * @brief estimates the posterior data fit
    */
    vdouble lsum;
    vdouble* chain_stats;

    void apply_transform2(const tuint idt);
    FlxFunction* gen_para_fun(const tuint ftype, const tuint pid);
  public:
    flxBayDA_likeli(tuint type, flxVec& data_vec_, const tuint id_transform, FlxRndCreator& RndCreator, tuint Nchain, tuint Npost, tuint N_adapt);
    ~flxBayDA_likeli();

    const tuint get_M() const { return M; }
    flxVec& get_pvec() { return *pvec; }

    const tdouble calc_likeli();
    void initialize_chains();
    void move_chains(const bool is_tuning_phase);
    void fill_post_samples();
    void get_post_sample(flxVec* sample_ptr=NULL);

    const tdouble get_fit() { return lsum.get_mean(); }
    const tdouble get_fit_sd() { return sqrt(lsum.get_variance()); }
    void get_posterior_mean(flxVec& mv, flxVec& sv) const;
    const flxVec& get_pvec() const { return *pvec; }
    RBRV_entry_RV_base* get_rv() { return rv_ext; }
    const tdouble eval_Gelman_Rubin_convergence(const tuint m);
};



/**
* @brief main class for Bayesian data analysis
*/
class flxBayDA : public FlxDataBase {
  protected:
    /**
    * @brief name-ID of the Bayesian data analysis object
    */
    std::string name;
    /**
    * @brief defines the domain of the data
    *
    * 0: no-id_transform
    * 1: only positive values allowed
    */
    tuint id_transform;
    /**
    * @brief the observed data - for which a distribution is to be derived
    *
    * values are transformed according to id_transform, such that they could potentially span the entire real line
    */
    flxVec data_vec;
    FlxRndCreator& RndCreator;
    /**
    * @brief Number of chains
    */
    const tuint Nchain;
    /**
    * @brief length of burn-in period (for each chain)
    */
    const tuint Nburn;
    /**
    * @brief length of tuning period (for each chain)
    *
    * must be smaller than Nburn
    */
    const tuint Ntune;
    /**
    * @brief total number of posterior samples to generate
    */
    const tuint Npost;
    /**
    * @brief adapt the step-size after how many chain moves?
    */
    const tuint N_adapt;
    /**
    * @brief random variable associated with the problem
    */
    RBRV_entry_RV_UserTransform* rv;
    RBRV_entry_RV_stdN rv_dummy;
    flxVec* model_pr;
    flxBayDA_likeli** models;
    /**
    * @brief threshold for deciding whether or not a model is unplausible
    */
    const tdouble tplaus;
    /**
    * @brief random variable associated with the problem
    */
    std::valarray<int> dtfv;
    /**
    * @brief random variable associated with the problem
    */
    const std::string pvec_name;
    const std::string distid_name;

    const tdouble find_MLE(flxBayDA_likeli& likeli, const tdouble step_size, const bool return_ini);
    void free_models();
  public:
    flxBayDA(const std::string& name, const tuint id_transform, FlxSMtx& smtx, FlxRndCreator& RndCreator, tuint Nchain, tuint Nburn, tuint Ntune, tuint Npost, tuint N_adapt, const tdouble tplaus, std::valarray<int> dtfv, const std::string& pvec_name, const std::string& distid_name);
    virtual ~flxBayDA();

    void gen_samples();
    void sample();
    RBRV_entry_RV_UserTransform* get_rv_ptr() { return rv; }

};


class FunRVcheckX : public FunBaseFun_onePara {
  private:
    RBRV_entry* rv;
  public:
    FunRVcheckX (RBRV_entry* rv, FunBase *child_1) : FunBaseFun_onePara(child_1), rv(rv) {};
    const tdouble calc();
    const std::string write_v();
};

class FlxBayUpBox;
class RBRV_entry_read_bayDA : public RBRV_entry_read_base {
  protected:
    FlxString* name_bayDA;
  public:
    RBRV_entry_read_bayDA(const bool readName, const bool readBrakets=true);
    virtual ~RBRV_entry_read_bayDA();

    static FlxBayUpBox* BayUpBox;

    virtual RBRV_entry* generate_entry(const std::string& family, tuint& running_iID);

};



