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

#include "flxdata.h"
#include "flxrandom.h"
#include "flxrbrv_rvs.h"
#include "flxStatBox.h"
#include "flxBayDA.h"

class FlxBayUP_csm_base;
class flxBayUp;



/**
* @brief controls the adaptive process (base class)
*/
class flxBayUp_adaptive_ctrl_base {
  protected:
    FlxFunction* updatesAfterNsamples;
    const tuint smpl_order;
    
    /**
    * @brief returns the number of adaptive updates per conditioning step
    */
    const tuint get_maxUpdatesPerCStep();
  public:
    flxBayUp_adaptive_ctrl_base(FlxFunction* updatesAfterNsamples, const tuint smpl_order);
    virtual ~flxBayUp_adaptive_ctrl_base();
    /**
    * @brief returns a copy of this element
    */
    virtual flxBayUp_adaptive_ctrl_base* copy() = 0;
    
    virtual const bool is_adaptive() const = 0;
    /**
    * @brief controls the adaptive process
    */
    virtual void eval() = 0;
    /**
    * @brief checks if an adaptive adaption is required
    */
    virtual void requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm) = 0;
    /**
    * @param Nc_now number of chains/seeds in the conditioning step
    * @param returns after how many seeds an adaptive update is performed (0 if no adaptive step is desired)
    */
    const tuint get_updatesAfterNsamples();
    /**
    * @brief returns the proposed order of samples
    * 0: (small to large g) no reordering
    * 1: (large to small g)
    * 2: random
    * 3: (increasing length)
    */
    const tuint get_smpl_order() const;
    
    virtual void print_info(std::ostream& sout) const;
};


/**
* @brief controls the adaptive process (enforce bounds)
*/
class flxBayUp_adaptive_ctrl_bounds : public flxBayUp_adaptive_ctrl_base {
  private:
    tdouble factor_d;
    tdouble lower_d;
    tdouble upper_d;
    FlxFunction* factor;
    FlxFunction* lower;
    FlxFunction* upper;
  public:
    flxBayUp_adaptive_ctrl_bounds(FlxFunction* factor, FlxFunction* lower, FlxFunction* upper, FlxFunction* maxUpdatesPerCStep, const tuint smpl_order);
    ~flxBayUp_adaptive_ctrl_bounds();
    virtual flxBayUp_adaptive_ctrl_bounds* copy();
    
    virtual const bool is_adaptive() const { return factor_d>GlobalVar.TOL(); }
    virtual void eval();
    virtual void requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm);
    virtual void print_info(std::ostream& sout) const;
};

/**
* @brief controls the adaptive process (enforce bounds)
*/
class flxBayUp_adaptive_ctrl_log : public flxBayUp_adaptive_ctrl_base {
  private:
    FlxFunction* f1;
    FlxFunction* f2;
    FlxFunction* ftacr;
  public:
    flxBayUp_adaptive_ctrl_log(FlxFunction* f1, FlxFunction* f2, FlxFunction* ftacr, FlxFunction* maxUpdatesPerCStep, const tuint smpl_order);
    ~flxBayUp_adaptive_ctrl_log();
    virtual flxBayUp_adaptive_ctrl_log* copy();
    
    virtual const bool is_adaptive() const { return true; }
    virtual void eval() {}
    virtual void requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm);
    virtual void print_info(std::ostream& sout) const;
};


const tuint acdcs_dim = 10;
class FlxBayUP_csm_dcs_MCMC;
class FlxBayUP_csm_csus_MCMC;
class flxBayUp_adaptive_ctrl_dcs : public flxBayUp_adaptive_ctrl_base {
  protected:
    FlxBayUP_csm_dcs_MCMC* csm_dcs;
    FlxBayUP_csm_csus_MCMC* csm_csus;
    tdouble factor_d;
    FlxFunction* factor;
    tdouble psd_max_d;
    FlxFunction* psd_max;
    
    // history of past samples
      tuint smpl_i;        // current position in list
      tuint smpl_N;        // total number of entries in list
      tuint smpl_Nmax;        // maximum number of entries in list
      tdouble* smpl_list;        // acdcs_dim x smpl_Nmax
        // 0: yR  (u-transform)
        // 1: yW
        // 2: sdR
        // 3: sdW
        // 4: seed_radius
        // 5: seed_radius (y-transform)
        // 6: stepR
        // 7: cosw
        // 8: step_size
        // 9: sample_direction
        // a+0: measure to optimize
      bool* smpl_acc;
      tdouble* smpl_weightsP;        // smpl_Nmax-vector that contains importance weights
      tdouble* smpl_resP;                // smpl_Nmax-vector that contains quantities of interst in LSF
    // current configuration  ( transformed to standard-Normal space)
      tdouble cur_sdR_ut;
      tdouble cur_sdW_ut;
      tdouble cur_sdWS_ut;
      tdouble cur_pSD;
      tuint shift;                // the shift to use when evaluating the LSF
    // average quantities
      tdouble omega_sum;
      tuint omega_N;
      tdouble sdR_sum;
      tuint sdR_N;
      tdouble sdW_sum;
      tuint sdW_N;
      tdouble sdWS_sum;
      tuint sdWS_N;
      tdouble pSD_sum;
      tuint pSD_N;
      
      tuint callN;
      
    const tdouble smpl_mean(const tuint shiftV, const bool consider_acc, const bool dois);
    
    void do_gsl_opti(const tuint mode);  // mode: 0 full; 1 radius; 2 angle
    void plot_shift();
    void plot_smpls();
    const tdouble adopt_to_acr(const tdouble acr, const tdouble cur_sd_ut, const tdouble prev_sd);
    
    void requires_adptv_step_dcs(const tdouble acr, FlxBayUP_csm_base& csm);
    void requires_adptv_step_csus(const tdouble acr, FlxBayUP_csm_base& csm);
  public:
    flxBayUp_adaptive_ctrl_dcs(FlxFunction* maxUpdatesPerCStep, FlxFunction* factor, FlxFunction* psd_max, const tuint smpl_order);
    virtual ~flxBayUp_adaptive_ctrl_dcs();
    virtual flxBayUp_adaptive_ctrl_dcs* copy();
    
    virtual const bool is_adaptive() const { return true; }
    virtual void eval();
    virtual void requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm);
    
    // interface to communicate with FlxBayUP_csm_dcs_MCMC
      virtual void register_csm(FlxBayUP_csm_dcs_MCMC* csm_dcsV);                // register MCMC algorithm initially
      virtual void register_csm(FlxBayUP_csm_csus_MCMC* csm_csusV);                // register MCMC algorithm initially
      void append_smpl(const flxVec& lastS, const bool accpetedV);        // append a proposed sample
      void write_adaptive_info(std::ostream& sout);
      
    const tdouble LSF(const tdouble sdR, const tdouble sdW, const tuint mode);
    const tdouble get_cur_sdR_ut() const { return cur_sdR_ut; }
};


const tuint acvelo_dim = 3;
class flxBayUp_adaptive_ctrl_velo : public flxBayUp_adaptive_ctrl_base {
  private:
    tdouble vspreadh_d;
    FlxRndCreator& RndCreator;
    FlxFunction* vspread;
    // history of past samples
      tuint smpl_i;        // current position in list
      tuint smpl_N;        // total number of entries in list
      tuint smpl_Nmax;        // maximum number of entries in list
      tdouble* smpl_list;        // acvelo_dim x smpl_Nmax
        // 0: sd
        // 1: velo2
        // 2: acc
    // current configuration  ( transformed to standard-Normal space)
      tdouble cur_sd_ut;
    // average quantities
      tdouble sd_sum;
      tuint sd_N;
      int velo_seq_count;
    const tdouble get_dynamic_spread() const;
  public:
    flxBayUp_adaptive_ctrl_velo(FlxRndCreator& RndCreator, FlxFunction* vspread, FlxFunction* maxUpdatesPerCStep, const tuint smpl_order);
    ~flxBayUp_adaptive_ctrl_velo();
    virtual flxBayUp_adaptive_ctrl_velo* copy();
    
    virtual const bool is_adaptive() const { return true; }
    virtual void eval();
    virtual void requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm);
    void append_smpl(const flxVec& lastS);        // append a proposed sample
    virtual void print_info(std::ostream& sout) const;
    const tdouble get_working_sd() const;
    void write_adaptive_info(std::ostream& sout);
};


const tuint spread_hist_N_max = 10;
const tuint acopti_jump_dim = 6;
class flxBayUp_adaptive_ctrl_opti_jump : public flxBayUp_adaptive_ctrl_base {
  private:
    FlxRndCreator& RndCreator;
    FlxFunction* acr_min;        // the smallest allowed acceptance rate
    FlxFunction* esjd_scale;        // find the spread that decreases the optimal ESJD by this factor
    FlxFunction* pw_p1;                // 1. weight according to number of adaptive step; 0. equal weight for each step 
    FlxFunction* pw_p2;                // pw_p1=0 ==> weighting factor for each adaptive step
    FlxFunction* aeps;                // C.o.V. of optimum to deactivate adaptive algorithm
    FlxFunction* Nmax_fun;        // evaluates to smpl_Nmax;
    tuint M_dim;                // number of random variables in the problem
    tdouble sd_prev_min;        // minimum from last optimization
    // history of past samples
      tuint smpl_i;        // current position in list
      tuint smpl_N;        // total number of entries in list
      tuint smpl_Nmax;        // maximum number of entries in list
      tdouble* smpl_list;        // acopti_jump_dim x smpl_Nmax
        // 0: lseed
        // 1: a
        // 2: b^2
        // 3: spread*
        // 4: acr (ONE or ZERO)
        // 5: measure to optimize
      tdouble* swl;                // weight list for the seeds
      tdouble* twl;                // weight list for samples
      flxVec spread_hist_vec;
      tuint spread_hist_N;
      bool deactivated;                // true if the adaptive procedure is deactivated
    tdouble limit_acr;                // value of acr_min
    flx_interp ipds;                // used to interpolate function to optimize
    const tdouble proposal_pdf_ln(const tdouble* c, const tdouble sd) const;
    void compute_seed_weights();
    const tdouble compute_overall_acr();
    const bool skip_adaptive_step(const tdouble acr);
    const tdouble get_pweight();
  public:
    flxBayUp_adaptive_ctrl_opti_jump(FlxRndCreator& RndCreator, FlxFunction* acr_min, FlxFunction* esjd_scale, FlxFunction* pw_p1, FlxFunction* pw_p2, FlxFunction* aeps, FlxFunction* Nmax_fun, FlxFunction* maxUpdatesPerCStep, const tuint smpl_order);
    ~flxBayUp_adaptive_ctrl_opti_jump();
    virtual flxBayUp_adaptive_ctrl_opti_jump* copy();
    void initialize(const tuint M_dim_v, const tuint Nfinal);
    
    virtual const bool is_adaptive() const { return true; }
    virtual void eval();
    virtual void requires_adptv_step(const tdouble acr, FlxBayUP_csm_base& csm);
    virtual void print_info(std::ostream& sout) const;
    const tdouble perfFun(const tdouble sd);
    void write_adaptive_info(std::ostream& sout);
    
    // interface to communicate with FlxBayUP_csm_dcs_MCMC
      void append_smpl(const flxVec& lastS);        // append a proposed sample
};


class Flx_SuS_CLevelStat;
class Flx_SuS_Control;
/**
* @brief manages the seed list
*/
class FlxBayUp_Update_List : public FlxBoxBaseR {
  public:
    flxBayUp& parent;
  private:
    RBRV_constructor& RndBoxRB;
    const tuint N_RV;        // number of random variables
    const tuint N_OX;        // dimension of the x-vector
    const tuint Nc;        // number of chains
    const tuint Ncl;        // number of samples per chain
    /**
    * @brief the number of seeds to generate in the final updating stage
    */
    const tuint Ns_final;
    const tuint Nburn;
    /**
    * @brief stores some extented output .. for logging purposes
    */
    std::ostringstream ext_out;
  public:
    // methods for Bayesian updating supported by this class
    enum MethType { BUS, UBUS, BUST, ABCSUBSIM, RASUBSIM, TMCMC, MHRS, RS, ABCRS, RAMCI, LS };
    // technique used to randomize seed values
    enum randomizeTec { NONE, INIT, RNDPICK };
    /**
    * @brief ID of the method to perform the Bayesian updating
    */
    const MethType meth_id;
  
  private:
    tdouble& c;                // the scaling constant we are working with (ideally max_L) - log-transform
    tdouble max_L;        // the largest observed value of the Likelihood function - log-transform
    tdouble s_thr;        // the threshold value of the LSF
    bool fullList;        // true: full list is used (up to Ns_final)
    tdouble* y_list;        // contains Ns_final realizations of N_RVs
    flxVec yh_vec;                 // for temporary storing a y-vector
    tdouble* x_list;        // contains Ns_final realizations of N_OXs
    flxVec xh_vec;                // for temporary storing a x-vector
    tdouble* p_list;        // contains Ns_final last entries of y_list        (not 'uniform' but standard Normal)
    tdouble* s_list;        // the 'limit-state' function -> needs to be updated if c is changed
    tdouble* sh_list;                // for temporary sorting of s_list
    tdouble* L_list;        // contains Ns_final Likelihood values
    bool* multpleVec_help;// temporary vector for finding multiple samples 
    /**
    * @brief 
    * before updating:
    *   -2: invalid -> outsid of 'failure' domain
    *   -1: end of list
    *   0: free entry
    *   1: occupied
    *   2: seed value
    * after updating:
    *   -2: invalid -> outside of 'failure' domain
    *   0: seed that has not been used before
    *   1: seed that was used before -> draw next value from chain
    */
    int* i_list;
    /**
    * @brief contains the indices of seeds (no ordered!) 
    * @note the list must end with an entry Ns_final+1
    */
    tuint* seed_idx;
    /**
    * @brief length of the individual chains
    */
    tuint* chain_length;
    const tuint max_runs;
    const bool log_LSF;                // extensive logging of results
    // remember previous c-values
      tuint oldC_N;
      tdouble* oldC_list;
      tdouble* oldS_list;
    flxBayUp_adaptive_ctrl_base* adpt_ctrl;
    const randomizeTec randomizeID;
    
    bool finalized;
    /**
    * for updating: the current ID in the list
    * for post-processing: the current ID in the list
    */
    tuint curID;
    /**
    * for updating: the ID of the current seed 
    * for post-processing: the current ID if randomizeID==INIT
    */
    tuint curSID;
    
    // ==================================================================================================
    // Methods used (exclusively) in the updating step
    // ==================================================================================================
    /**
    * @brief returns the percentile of LSF-values smaller or equal than zero for given 's_thr'
    */
    const tdouble get_perc_BUST();
    /**
    * @brief stores the original seed ids to g_seed_ID!
    */
    void prepare_for_gamma_comp_1(Flx_SuS_CLevelStat& curLevelStat_next) const;
    void prepare_for_gamma_comp_2(Flx_SuS_CLevelStat& curLevelStat_next) const;
    
    // ==================================================================================================
    // Methods used (exclusively) in the posterior step
    // ==================================================================================================
    
    /**
    * @brief the sample at curID introduces a new maxL -> fix the statistics
    * @param smpl_part_of_posterior: if false, we just know the the likelihood used is too small
    *       a comment on smpl_part_of_posterior=false:
    *                 (based on a sample NOT directly generated from the prior & MCMC)
    *           -> we realized this just by checking whether an (arbitrary) sample is in the failure domain
    * @param lMay_: needs to be set if smpl_part_of_posterior=false
    * @returns true if sample at curID was replaced (only relevant for smpl_part_of_posterior=true)
    */
    const bool update_c_posterior(const bool smpl_part_of_posterior=true, const tdouble lMax_ = ZERO);
    
    // ==================================================================================================
    // Methods used in the updating and in the posterior step
    // ==================================================================================================
    
    /**
    * @brief sets curID to the next free entry (or stops at the end)
    * @param init true: only if to be initialized for the start of the conditioning level
    */
    void set_next(const bool init=false);
    
    const tdouble propose_qlnL(const tdouble pa_maxL);
    void swap_smpls(const tuint id1, const tuint id2);
    void update_LSF_vals(const tuint N_in_list, const bool do_check=true);
    void assign_new_p_vals();
    const tdouble MHRS_uBUS_help(const tdouble maxLc, const flxVec& ut_vec, flxVec& Lt_vec);
    void MHRS_uBUS(tdouble& p_mod);
    /**
    * @param L the value of the likelihood function (in log-transform)
    * @param ci the value of the scaling constant (<=maxL) (in log-transform)
    * @param hi the threshold value of the current level
    * @returns log-transform of distribution
    */
    const tdouble eval_L_thresh(const tdouble L, const tdouble ci, const tdouble hi) const { const tdouble Lt=L+hi; return (Lt>ci)?ci:Lt; }
    
  public:
    FlxBayUp_Update_List(flxBayUp& parent, const tuint Nc, const tuint Ncl, const tuint Ns_finalV, const tuint Nburn, const randomizeTec randomizeID, flxBayUp_adaptive_ctrl_base* adpt_ctrl, const tuint max_runs, const bool use_cStart, const MethType meth_id, const bool log_LSF, const bool find_multiple);
    ~FlxBayUp_Update_List();
    
    // ==================================================================================================
    // General methods
    // ==================================================================================================
    
    const tuint get_Nrv() const { return N_RV; }
    const tuint get_NOX() const { return N_OX; }
    const tuint get_Ns_final() const { return Ns_final; }
    /**
    * @brief returns the proposed number of chains
    */
    const tuint get_Nc() const { return Nc; }
    const tuint get_Nburn() const { return Nburn; }
    
    RBRV_constructor& get_RndBox() { return RndBoxRB; }
    
    /**
    * @brief add all seeds to sbox (used by a MCMC method)
    */
    void fill_sbox(FlxStatBox& sbox);
    void fill_slist(std::vector<tdouble*>& slist);
    
    flxBayUp_adaptive_ctrl_base& get_adpt_ctrl() { return *adpt_ctrl; }
    
    void print_ext_out(std::ostream& os);
    
    // ==================================================================================================
    // Methods used (exclusively) in the updating step
    // ==================================================================================================
    
    const tuint get_Ncl() const { return Ncl; }
    const tuint get_max_runs() const { return max_runs; }
    /**
    * @brief computes the 'velocity' of samples in the set
    * @param NcNow Number of seed values (if=0: no seeds in list)
    *   primarily to be used for the updating step
    */
    const tdouble get_velo(const tuint NcNow = 0) const;
    /**
    * @returns the current threshold (for subset simulation)
    */
    const tdouble get_thr() const { return s_thr; } 
    tuint& get_curSID_ref();
    tdouble* get_seed_y_list();        // returns the y-vector of the current seed
    #if FLX_DEBUG
      const tdouble get_seed_L_list() const;
      const tdouble get_seed_p_list() const;
      const tdouble get_seed_s_list() const;
    #endif
    
    /**
    * @brief updates s_thr and marks seeds 
    * @note for BUST s_thr needs to be set beforehand!!!
    * @returns returns the number of chains (seeds) in the next run
    */
    const tuint update_thr(Flx_SuS_CLevelStat& curLevelStat_next, Flx_SuS_CLevelStat& curLevelStat);
    /**
    * @brief for BUST: needs to be executed instead of update_thr
    *   (calls update_thr internally)
    */
    const tuint update_thr_BUST(tdouble& pnow, Flx_SuS_CLevelStat& curLevelStat_next, Flx_SuS_CLevelStat& curLevelStat);
    /**
    * @brief sets curSID to the next entry (for the updating step)
    */
    void set_next_seed() { ++curSID; }
    /**
    * @brief returns true if all Ns_final samples are used deliberately
    */
    const bool is_fullList() const { return fullList; }
    /**
    * @brief returns true if the threshold for the limit-state reached zero (important during subset simulation)
    */
    const bool is_gt_zero() const;
    const tdouble get_maxL() const { return max_L; }
    /**
    * @brief updates the LSF - if max_L has changed
    * @note do not use this function if you are drawing samples from the posterior!
    * @param is_first true: we are in the zeroth level of Subset Simulation
    * @param zeroLSF true: we are in the last level of Subset Simulation (actual failure domain)
    * @returns true, sampling distribution is not affected by change in c
    */
    const bool update_c(tdouble& p_mod, const bool is_first=false);
    /**
    * @brief check if some intermediate c-LSF are nested and can be removed
    */
    void nested_c();
    /**
    * @brief returns the length of the chain associated with the current seed
    */
    const tuint get_cur_chain_length() const { return chain_length[curSID]; }
    /**
    * @brief computes the factor for chain correlation
    */
    void compute_gamma(Flx_SuS_CLevelStat& curLevelStat, const Flx_SuS_Control& sctrl) const;
    /**
    * returns the mean of the log-transformed likelihood values
    */
    const tdouble expectation_likelihood() const;
    
    /**
    * @brief prepares the list for drawing posterior realizations
    * @note indicates the end of the updating step
    */
    const tuint finalize();
    const bool is_finalized() { return finalized; }
    void reset_finalized();
    
    tdouble* TMCMC_get_y_list() { return y_list; }
    tdouble* TMCMC_get_x_list() { return x_list; }
    tdouble* TMCMC_get_L_list() { return L_list; }
    
    // ==================================================================================================
    // Methods used (exclusively) in the posterior step
    // ==================================================================================================
    
    const tdouble* get_cur_y_list();        // returns the y-vector of the current sample 
    const tdouble* get_cur_x_list();        // returns the x-vector of the current samples
    const tdouble get_cur_L();                // returns the current likelihood value
    
    /**
    * @brief chooses the next curID (for drawing from the posterior)
    * @note before changing curID, the i_list-entry of the sample is marked with 1
    */
    void set_next_draw();    
    
    // ==================================================================================================
    // Methods used in the updating and in the posterior step
    // ==================================================================================================
    
    const int get_cur_i_list() const;        // returns the identifier of the current entry
    
    /**
    * @brief inserts a new entry at the current position (curID)
    * @param is_first always accept the sample
    * @param rejectIt do not accept the sample
    * @param is_posterior updating step is complete; this sample is drawn as a posterior sample (after finalize!)
    * @param just_accpt_chk just check whether the sample can be accepted (do not store it) ... is_posterior must be true
    * @param tLik value of the likelihood (only relevant in case of RS)
    * @return true, if proposal has been accepted
    * @note if is_posterior==false&&meth==SubsetBased 
    *           if sample accepted: seed_idx[curSID] is set to the current sample!
    *           if sample rejected: then seed_idx[curSID] has to temporarily point to the previously used sample 
    */
    const bool insert_entry(const bool is_first, bool rejectIt, const bool is_posterior, const bool just_accpt_chk, std::ofstream* os_smpl, const tdouble tLik=ZERO, tdouble* acr=NULL);
    /**
    * @brief evaluates the limit-state function of the BUS problem
    * @param p standard normal realization
    * @param L the Likelihood function - log-transform
    */
    const tdouble eval_LSF(const tdouble p, const tdouble L) const;
    
    void write_smpl(const tuint ID2write, std::ofstream& os_smpl);
    
};


class FlxBayUP_csm_base : public FlxDataBase {
  protected:
    FlxRndCreator& RndCreator;
    FlxFunction* h_fun;
    // for adaptively optimizing the ESJD
      flxBayUp_adaptive_ctrl_velo* adpt_velo;
      flxVec lastVeloInfo;
    
  public:
    FlxBayUP_csm_base(FlxRndCreator& RndCreator, FlxFunction* h_fun);
    virtual ~FlxBayUP_csm_base() { if (h_fun) delete h_fun; }
    
    virtual void register_adpt_ctrl(flxBayUp_adaptive_ctrl_base* adpt_ctrl);
    /**
    * @brief prepare the chain at the beginning of a conditioned sampling step
    */
    virtual void prepare() {}
    /**
    * @brief prepare for a new Markov-chain seed (where y is the seed)
    */
    virtual void prepareS(const flxVec& y) {};
    /**
    * @brief propose a new sample based on the previous sample
    * @returns false if the y_prop = y_prev
    */
    virtual const bool propose(flxVec& y_prop, const flxVec& y_prev) = 0;
    /**
    * @param was_accepted = true, if proposed sample was accpeted; false, otherwise
    */
    virtual void acceptance_feedback(const bool was_accepted);
    virtual void adptv_spread_multiply(const tdouble f) = 0;
    virtual const std::string print_info() = 0;
    virtual void write_adaptive_info(std::ostream& sout, const bool is_adaptive) = 0;
    /**
    * @brief acceptance rate of a 'single-component' in the pre-canditate state
    */
    virtual const tdouble get_ac1d() const = 0;
};

class FlxBayUP_csm_cwmh_MCMC : public FlxBayUP_csm_base {
  private:
    FlxRndKernel_base* kernel;
    size_t N1Dacc;        // for get_ac1d()
    size_t N1Dtotal;        // for get_ac1d()
  public:
    FlxBayUP_csm_cwmh_MCMC(FlxRndCreator& RndCreator, const std::string& kernelName, const tdouble h_value, FlxFunction* h_fun);
    ~FlxBayUP_csm_cwmh_MCMC() { delete kernel; }
    
    virtual void prepare();
    const bool propose(flxVec& y_prop, const flxVec& y_prev);
    virtual void adptv_spread_multiply(const tdouble f);
    const std::string print_info();
    void write_adaptive_info(std::ostream& sout, const bool is_adaptive);
    virtual const tdouble get_ac1d() const { return tdouble(N1Dacc)/tdouble(N1Dtotal); }
};

class FlxBayUP_csm_cov_MCMC : public FlxBayUP_csm_base {
  private:
    const tuint M;        // dimension of the problem
    tdouble h;                // h-value 
    tdouble p;                // fraction of samples to consider for sample-covariance
    tuint Nmax;                // maximum number of samples to consider for sample-covariance
    tdouble p_single;        // fraction of samples to consider for sample-covariance
    tuint Nmax_single;        // maximum number of samples to consider for sample-covariance
    flxVec sde;                // standard deviations in eigen-directions
    flxVec hvN1;                // this is a vector that acts as a temporary storage
    flxVec hvN2;                
    flxVec hvM;
    iVec ivN;
    FlxMtxSym covmtx;        // sample covariance matrix
    FlxMtx Tinv;        // matrix that transforms from eigen-space to standard normal space
    FlxStatBox sbox;        
    FlxBayUp_Update_List& list;
    std::vector<flxVec> Eigenvectors;
    FlxRndKernel_base* kernel;
    size_t N1Dacc;        // for get_ac1d()
    size_t N1Dtotal;        // for get_ac1d()
  public:
    FlxBayUP_csm_cov_MCMC(FlxRndCreator& RndCreator, const tuint M, const std::string& kernelName, const tdouble h_value, FlxFunction* h_fun, const tdouble p, const tuint Nmax, const tdouble p_single, const tuint Nmax_single, FlxBayUp_Update_List& listV);
    ~FlxBayUP_csm_cov_MCMC() { delete kernel; }
    
    void prepare();
    virtual void prepareS(const flxVec& y);
    const bool propose(flxVec& y_prop, const flxVec& y_prev);
    virtual void adptv_spread_multiply(const tdouble f);
    const std::string print_info();
    void write_adaptive_info(std::ostream& sout, const bool is_adaptive);
    virtual const tdouble get_ac1d() const { return tdouble(N1Dacc)/tdouble(N1Dtotal); }
};

class FlxBayUP_csm_csus_MCMC : public FlxBayUP_csm_base {
  private:
    tdouble rho;
    tdouble sD;
    flxBayUp_adaptive_ctrl_opti_jump* adpt_ctrl;
    // descriptor of last proposed state
      flxVec lastS;
        // 0: length^2 of standard Normal vector
        // 1: N_RV
        // 2: sd
        // 4: seed radius
        // 7: cos(omega)
        // 8: squared step_size
  public:
    FlxBayUP_csm_csus_MCMC(FlxRndCreator& RndCreator, const tdouble sD, FlxFunction* h_fun);
    
    virtual void register_adpt_ctrl(flxBayUp_adaptive_ctrl_base* adpt_ctrlV);
    void prepare();
    const bool propose(flxVec& y_prop, const flxVec& y_prev);
    void acceptance_feedback(const bool was_accepted);
    virtual void adptv_spread_multiply(const tdouble f);
    const std::string print_info();
    void write_adaptive_info(std::ostream& sout, const bool is_adaptive);
    virtual const tdouble get_ac1d() const { return ONE; }
    // interface to communicate with flxBayUp_adaptive_ctrl_dcs
      void get_cur_spread(tdouble& sdV) const;
      void set_cur_spread(const tdouble& sdV);
};

class FlxBayUP_csm_dcs_MCMC : public FlxBayUP_csm_base {
  private:
    tdouble rhoR;
    tdouble sdR;
    tdouble sdW;        // std.dev. of angle-proposal in case of isotropic direction
    tdouble sdWS;        // std.dev. of angle-proposal in case of seed-direction 
    tdouble pSD;        // probability of generating sample based on seeds
    flxBayUp_adaptive_ctrl_dcs* adpt_ctrl;
    // descriptor of last proposed state
      flxVec lastS;
        // 0: yR (u-transform)
        // 1: yW
        // 2: sdR
        // 3: sdW
        // 4: seed_radius
        // 5: seed_radius (u-transform)
        // 6: stepR
        // 7: cosw
        // 8: squared step_size
        // 9: sample_direction
    FlxBayUp_Update_List& list;
    std::vector<tdouble*> seedVec;
  public:
    FlxBayUP_csm_dcs_MCMC(FlxRndCreator& RndCreator, const tdouble sdV, const tdouble pSD, FlxFunction* h_fun, FlxBayUp_Update_List& list);
    
    virtual void register_adpt_ctrl(flxBayUp_adaptive_ctrl_base* adpt_ctrlV);
    void prepare();
    const bool propose(flxVec& y_prop, const flxVec& y_prev);
    void acceptance_feedback(const bool was_accepted);
    virtual void adptv_spread_multiply(const tdouble f);
    const std::string print_info();
    void write_adaptive_info(std::ostream& sout, const bool is_adaptive);
    virtual const tdouble get_ac1d() const { return ONE; }
    // interface to communicate with flxBayUp_adaptive_ctrl_dcs
      void get_cur_spread(tdouble& sdRV, tdouble& sdWV, tdouble& sdWSV, tdouble& pSDV) const;
      void set_cur_spread(const tdouble& sdRV, const tdouble& sdWV, const tdouble& sdWSV, const tdouble& pSDV);
};

class FlxBayUP_csm_TMCMC : public FlxBayUP_csm_base {
  private:
    tdouble beta;
    FlxMtxLTri Acov;
    bool is_set;
  public:
    FlxBayUP_csm_TMCMC(FlxRndCreator& RndCreator, const tuint M, const tdouble beta, FlxFunction* h_fun);
    
    void prepare(const bool is_posterior=false);
    void prepare(const flxVec& covMtx);
    const bool propose(flxVec& y_prop, const flxVec& y_prev);
    virtual void adptv_spread_multiply(const tdouble f);
    const std::string print_info();
    void write_adaptive_info(std::ostream& sout, const bool is_adaptive);
    virtual const tdouble get_ac1d() const { return ONE; }
};

class Flx_SuS_CLevelStat : public FlxDataBase {
  private:
    const tuint get_MaxLevelDepth() const;        // maximum number of levels for common seeds
    // correlation of seeds
      const tuint get_seed_N_groups() const;        // total number of groups
      const tuint get_seed_group_size(tuint depth) const;
      const tuint get_seed_group_depth(const tuint seed_group) const;
      // get group -> from delta_level & delta_pos
      const tuint get_seed_group(const tuint delta_level, const tuint delta_pos) const;
      //const tuint get_seed_group_l1(const tuint depth) const { return (depth+1)/2; }
      //const tuint get_seed_group_p1(const tuint depth, const tuint l1p) const { return l1p*(1+depth-l1p); }
      //const tuint get_seed_group_l2(tuint depth, const tuint l1) const { return depth-l1; }
      // finds the first ID of cID in seed_chainID
      const tuint find_start_in_seed_chainID(const tuint cID) const;
      void add2seedCorr_2group(const tuint cID1, const tuint cID2, const tuint group, const bool isSeedCorr);
      /**
      * @brief allocates memory for computing seed correlation (or level correlation)
      * @param isSeedCorr true: seed correlation; false: level correlation
      */
      void allocate_g_seed_corrE(const bool isSeedCorr);
      void deallocate_g_seed_corrE(const bool isSeedCorr);
    // correlation of pis
      const tuint get_pic_N_groups() const;        // total number of groups
      const tuint get_pic_group_depth(const tuint pic_group) const;
      void add2piCorr_2group0(const tuint cID1, const tuint cID2, const tuint delta_pos);
    /**
    *   ordering of seedID is acording to acutal ordering - not original ordering
    * @param delta_level must initially be set to zero
    */
    const bool find_common_seed( const tuint chainID_this, const tuint posInChain_this, const tuint level_rhs, const tuint chainID_rhs, const tuint posInChain_rhs, tuint& delta_level, tuint& delta_pos, tuint mLdp=0 );
    /**
    * @brief empirical seed-Gamma (contribution of seed correlation to gamma)
    * @param lid - distance to the current level
    * @param f_l - vector with fractions of the different levels
    */
    void empirical_seedGamma(const tuint lid, tdouble* e_l, const tuint mld, const tdouble NsmpX, const tdouble* para);
    void empirical_corrLevel(const tuint lid, tdouble* e_l, const tuint mld, const tdouble NsmpX, const tdouble* para);
    
    tuint level;
    tuint maxLevelDepth;
  public:
    Flx_SuS_CLevelStat* prev;
    tdouble pi;                // condt. probability
    tuint Nsamples;        // total number of samples
    tuint Nchains;        // number of initial chains
    tuint Nfailures;        // number of chains retrieved!!!
    tdouble g_t;        // threshold value picked based on samples of this level
    // all related to gamma
      tdouble gamma;                // factor of chain correlation
      tdouble gamma_chain;        
      tdouble gamma_from_seed;        // factor of chain correlation (that comes from the seeds)
      tuint g_Ncl_max;                // the maximum chain length in this level
      tuint* g_seed_ID_original;        // seeds how they appeared in previous list of samples
      tuint* g_seed_ID;                // up to rho is computed: IDs of the initial Nchain chains
                                // after rho is computed: mapping of chain to parent seed in previous level
      tuint* g_chain_length;        // length of the initial Nchain chains
      tuint* g_N;                // g_N_entries*[total_N,N_first=1;N_1st&&N_2nd==1]
      tuint g_N_entries;        // sets in g_N
      tdouble* pi_chain;        // the estim. cond. prob. based on the position of the chain [g_Ncl_max]
      tdouble p_00, p_11;        
      tdouble lag1_corr;        // lag-1 autocorrelation
      private:
      tuint** g_seed_corrE;        // similar to g_N for seed-correlation
      public:
      tdouble eff;                // efficiency of the chain (in percent)
      tdouble eff_Gelman;        // efficiency-estimate based on Gelman
    // history of seeds
      tdouble corr_pi_prev;        // estimated correlation with previous conditional probabilities
      tdouble corr_pi_prev_negFrac;        // negative fraction when evaluating corr_pi_prev ....
      tdouble frac_rel_chains;        // fraction of chains from different levels that are correlated
      tuint* seed_chainID;        // from which chain does the seed come from (seed is for the next level!) [dimension: Nfailures]
      tuint* seed_chainPos;        // at which possition in the chain did the seed occurr? (seed is for the next level!) [dimension: Nfailures]
      tuint* find_multiples;
      
    Flx_SuS_CLevelStat(const tuint level, Flx_SuS_CLevelStat* prev);
    Flx_SuS_CLevelStat(const Flx_SuS_CLevelStat &rhs);
    ~Flx_SuS_CLevelStat();
    
    Flx_SuS_CLevelStat& operator=( const Flx_SuS_CLevelStat& rhs);
    
    const tuint get_level() const { return level; }
    
    void empirical_Corr();
    void add2seedCorr();
    void add2piCorr();
};

class Flx_SuS_Control {
  public:
    /**
    * @brief how to estimate the credible intervals of the failure probability
    *   none:        credible intervals are not estimated
    *   simpleBayes: use a Bayesian estimator that neglects chain correlation
    *   ccorrBayes:  use a Bayesian estimator that considers chain correlation (but regards seeds as uncorrelated)
    */
    enum credibleEstim { none, simpleBayes, ccorrBayes, fcorrBayes, icorrBayes };
    
    bool prt_alert; 
    bool verbose;
    credibleEstim credEst;
    FlxMtxConstFun* pc;                // contains the probabilities for the credible intervals
    tuint N_cred_smpl;                // samples used to evaluate credible intervals
    bool comp_gamma;                // compute the chain correlation
    bool consider_seed_corr;        // consider correlation of seeds when computing gamma
    bool consider_pi_corr;        // consider correlation of levels
    bool empirical_corr;        // determine seed_corr & pi_corr empirically
    bool find_multiples;        // find repeated samples
    FlxString* os_samples;        // file to write the samples to
    FlxFunction* TMCMC_target_COV;        // target coefficient of variation in TMCMC
    tuint TMCMC_update_weights;        // original TMCMC (false); improved version (true)
    FlxFunction* TMCMC_alpha;        // factor to approximate the intermediate distributions
    FlxMtxConstFun* LS_SPNT;        // start-point of the Line sampling
    FlxFunction* LS_tol;        // tolerance parameter of line search
    FlxFunction* LS_max_iter;        // maximum number of iterations in line search
    FlxFunction* pa_maxL;                // work with a reduced maxL in BUS
  
    Flx_SuS_Control() : prt_alert(true), verbose(false), credEst(none), pc(NULL), N_cred_smpl(0), comp_gamma(false), consider_seed_corr(false), consider_pi_corr(false), empirical_corr(false), find_multiples(false), os_samples(NULL), TMCMC_target_COV(NULL), TMCMC_update_weights(1), TMCMC_alpha(NULL), LS_SPNT(NULL), LS_tol(NULL), LS_max_iter(NULL), pa_maxL(NULL) {}
    Flx_SuS_Control(const Flx_SuS_Control &rhs);
    ~Flx_SuS_Control();
    
    Flx_SuS_Control& operator=( const Flx_SuS_Control& rhs);
    
    static credibleEstim parse_credibleEstim(const std::string& strCredEst);
    static const std::string get_credibleStr(const credibleEstim cID);
};


class flx_LS_line_prop {
  private:
    bool topo_set;
    bool only_in;
    tdouble l_out;
    tdouble l_in;
    tdouble u_out;
    tdouble u_in;
    std::stack<tdouble> *ostack;
    
    void set_topo();
  public:
    flx_LS_line_prop() : topo_set(false), only_in(true), l_out(-1e5),l_in(1e5),u_out(1e5),u_in(-1e5), ostack(NULL) {}
    flx_LS_line_prop(const flx_LS_line_prop& rhs);
    flx_LS_line_prop& operator=( const flx_LS_line_prop& rhs);
    ~flx_LS_line_prop();
    
    void register_in(const tdouble c);
    void register_out(const tdouble c);
    void register_c(const tdouble c, const tdouble g);
    const tdouble get_upper_Pr(const tdouble betaNorm) const;
    const tdouble get_mean_Pr(const tdouble betaNorm) const;
    const tdouble propose(const tdouble pr, const tdouble betaNorm) const;
    const bool is_topo_set() const { return topo_set; }
    const bool is_only_in() const { return only_in; }
    void force_topo(const tdouble c1, const tdouble g1, const tdouble c2, const tdouble g2);
    const tdouble get_lower_scale() const;
    void print_topo(std::ostream& os) const ;
};


class flxBayUp_RBRV_set;
class FlxBayUp_Update : public FlxDataBase {
  private:
    FlxBayUp_Update_List* list;
    /**
    * @brief the class managing the MCMC method
    */
    FlxBayUP_csm_base* csm;
    FlxRndCreator& RndCreator;
    flxBayUp_RBRV_set* burbrvs;
    tdouble& iadpt;
    // variables of the updating step
      tdouble PrMod;                 // stores the 'probability' of the model (after the updating was performed) (log-transform)
      tdouble DataFit;                // data fit term (after updating was performed)
      tdouble relEntr;                // relative Entropy (after updating was performed)
      tulong N_tot_LSFc_sim;        // sotres number of model calls (right after the updating was performed)
      Flx_SuS_Control susControl;
    // variables of the posterior step
      tuint post_adpt_calls;        // how often was the adaptive scheme already called?
      tuint post_adptcount_N;        // number of samples per adaptive step (0 if no adpt. step is desired)
      tuint post_adptcount;        // number of total samples since the last adaptive step
      tuint post_adpt_accsmpl;        // number of samples accepted since the last adaptive step
    
    std::vector<Flx_SuS_CLevelStat*> CLevelStat;
    
    void output_forwardEstimators();
    void output_credibleIntervals();
    std::ofstream* open_smpl_file4write();
    
    void TMCMC_weight_vec(const tdouble q_prev, const tdouble q_now, const flxVec& Lvec, flxVec& Wvec);
    const tdouble TMCMC_COV(const tdouble q_prev, const tdouble q_now, const flxVec& Lvec, flxVec& Wvec);
    const tdouble TMCMC_new_q(const tdouble q_prev, const tdouble target_cov, const flxVec& Lvec, flxVec& Wvec);
    void TMCMC_assemble_smplCOV(flxVec& y_covV, const flxVec& W_vec, const tdouble* y_list_prev, const tuint Nc, const tuint NRV, flxpVec& y_meanP, flxVec& y_mean, flxVec& y_tmp );
    
    const tdouble line_search_LSF_call(const tdouble c, const flxVec &rv_base, flxVec &rv_prop, const flxVec &betaVec, tulong& N_LSF_calls, flx_LS_line_prop& lsp, tdouble& L) const;
    const flx_LS_line_prop perform_line_search(const flxVec &rv_base, flxVec &rv_prop, const flxVec &betaVec, const tdouble tol, const tuint iter_max,tulong& N_LSF_calls) const;
  public:
    FlxBayUp_Update(FlxRndCreator& RndCreator);
    ~FlxBayUp_Update();
    
    FlxRndCreator& get_RndCreator() { return RndCreator; }
    /**
    * @brief performs the Bayesian updating 
    *         -> generate sample from the posterior starting with samples from the prior
    * @param listV (is deallocated even if an error is thrown)
    * @param csmV the method to perform the conditioned sampling
    * @param use_cStart only relevant if update is executed more than once ... reset cStart to initial value?!
    */
    void update(FlxBayUp_Update_List* listV, FlxBayUP_csm_base* csmV, const bool use_cStart, const Flx_SuS_Control& susControlV);
    /**
    * @brief returns a new realization - based on existing seeds (posterior sampling!)
    * @param smpl is the sample drawn from the posterior - the size of the vector must match the number of random variables in the updating problem
    */
    void draw_realization(flxVec& smpl);
    /**
    * @brief in the posterior stage: check if the current sample can be accepted
    * @note the sample to check MUST be in the Rnd-environment (as current sample in y-space)
    * @returns true if sample can be accepted
    */
    const bool chk_accept_cur_smpl();
    /**
    * @brief returns the probability of the model
    */
    const tdouble get_PrMod() const { return exp(PrMod); }
    const tdouble get_LogPrMod() const { return PrMod; }
    const tdouble get_LogDataFit() const { return DataFit; }
    const tdouble get_LogRelEntr() const { return relEntr; }
    /**
    * @brief total number of model evaluations !during! the updating step
    */
    const tulong get_total_LSFc_sim() const { return N_tot_LSFc_sim; }
    
    const FlxBayUp_Update_List& get_list() const;
    void reset_finalized_smpls();
    void get_sus_level_info(const std::string vecs, const tuint pid, const tuint pid2);
    
    static void define_constants();
};


/**
* @brief Manages the posterior distribution
*/
class flxBayUp_RBRV_set : public RBRV_set_base {
  private:
    flxBayUp& parent;
    const std::vector<RBRV_set_base*>& setvec;                // a list with all the relevant sets (the order does matter)
    const tuint NRV;
    const tuint NOX;
    const tuint Nsets;
    bool proposed;        // true, if x corresponds to y (because y was actually proposed)
    
  public:
    flxBayUp_RBRV_set(const bool internal, flxBayUp& parentV);
    virtual ~flxBayUp_RBRV_set() {}
    
    virtual const tuint get_NRV() const { return NRV; }
    virtual const tuint get_NOX() const { return NOX; }
    virtual const tuint get_NRV_only_this() const { return 0; }
    virtual const tuint get_NOX_only_this() const { return 0; }
    virtual void set_is_valid(const bool set_is_valid);
    virtual const flxVec& propose_y();
    virtual void transform_y2x();
    virtual const bool allow_x2y() const { return false; }
    /**
    * @brief sets y without setting x
    */
    virtual void set_y(const tdouble* const y_vec);
    virtual void set_y_only_this(const tdouble* const y_vec);
    virtual void get_y(tdouble* const y_vec);
    virtual void get_y_only_this(tdouble* const y_vec);
    virtual flxVec& get_y();
    /**
    * @brief sets x without setting y first
    */
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec);
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec);
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec);
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec);
    virtual void find_dependent_sets(std::vector<RBRV_set_base*>& setvecV );
    virtual const tuint group_dependent_sets(std::vector<RBRV_set_base*>& setvecV, const tuint pos_this );
    void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    
    // needed for weighting model probabilities 
    const flxBayUp& get_parent() const { return parent; }
};


/**
* @brief Uncertain observations
*/
class flxBayUp_uncertobsv_set : public RBRV_set_base {
  private:
    RBRV_set* singleObsv;        // a single observation
    FlxFunction* flxLike;        // the Likelihood function
    const tuint Nobserv;        // total number of observations
    const tuint paraN;                // number of parameters per observation
    tdouble* dataPtr;
    tdouble Likelihood;                // the Likelihood value (assume independency)
    const bool is_log;
    
  public:
    flxBayUp_uncertobsv_set(const std::string& name, RBRV_set* singleObsvV, FlxFunction* flxLike, const tuint NobservV, const tuint paraN, FlxIstream_vector *isv, const bool is_log);
    virtual ~flxBayUp_uncertobsv_set();
    
    virtual const tuint get_NRV() const { return singleObsv->get_NRV()*Nobserv; }
    virtual const tuint get_NOX() const { return 0; }
    virtual const tuint get_NRV_only_this() const { return Nobserv*singleObsv->get_NRV_only_this(); }
    virtual const tuint get_NOX_only_this() const { return 0; }
    virtual void set_is_valid(const bool set_is_valid);
    virtual void transform_y2x();
    virtual const bool allow_x2y() const { return false; }
    virtual void set_y_only_this(const tdouble* const y_vec) { set_y(y_vec); }
    virtual void get_y_only_this(tdouble* const y_vec) { get_y(y_vec); }
    /**
    * @brief sets x without setting y first
    */
    virtual void set_x(const tdouble* const x_vec) {}
    virtual void set_x_only_this(const tdouble* const x_vec) {}
    virtual void get_x(tdouble* const x_vec) {}
    virtual void get_x_only_this(tdouble* const x_vec) {}
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec);
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec);
    virtual void find_dependent_sets(std::vector<RBRV_set_base*>& setvec );
    void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
    
    const tdouble& get_logLike() const { return Likelihood; }
};


class flxBayUp : public FlxDataBase {
  public:
    enum MethCategory { BUS, TMCMC, ABC, RA, UNDEFINED };
    const bool is_subsetRel;        // true if object is not used for updating, but for subset simulation
  private:
    const tdouble scaleconst;                // the likelihood function is scaled by this constant (default: 1.0) - log-transform
    /**
    * @brief the maximum observed value of the total Likelihood function - log-Transform
    */
    tdouble cStart;
    tdouble cStart_init;        // the inital value
    tdouble pa_maxL;                // cStart can be reduced
    /**
    * @brief contains the RBRV-sets needed in the analysis
    */
    std::vector<RBRV_set_base*> setvec;
    /**
    * @brief individual Likelihood functions
    * the values can be accessed by nameID::[index] where [index] is the index in the list
    */
    std::vector<RBRV_entry*> LklVec;
    tuint N_LklVec;
    /**
    * @brief if !NULL: use this; otherwise multiply all local likelihood functions
    */
    FlxFunction* global_Likelihood;
    bool global_Likelihood_is_log;
    MethCategory methCat;
    /**
    * @brief contains the entries of LklVec -> if !NULL, box is locked
    */
    RBRV_set* LklSet;                // the last RV is 'p'
      RBRV_entry_RV_stdN* p_rv; // direct reference to 'p'
    RBRV_constructor* RndBox;
    const std::string nameID;
    
  public:
    /**
    * @brief NOTE: cStartV must be log-transformed
    */
    flxBayUp(const std::string& nameID, const tdouble& scaleconst, const tdouble& cStartV, const std::string& parentsetstr);
    flxBayUp(const std::string& parentsetstr);
    ~flxBayUp();
    
    FlxBayUp_Update updater;
    
    const std::string& get_name() const { return nameID; }
    const MethCategory get_methCat() const { return methCat; }
    /**
    * @brief adds a local Likelihood function
    * make sure to delete lklEntry if an exception is thrown
    */
    void add_localLkl(RBRV_entry* lklEntry);
    void add_localLkl(flxBayUp_uncertobsv_set* ts);
    
    const tuint get_N_localLkl() const { return N_LklVec; }
    /**
    * @brief sets the global Likelihood function
    * a copy of this function is generated
    */
    void set_globalLkl(FlxFunction& global_LikelihoodV, const bool is_log, const MethCategory methCatV);
    
    /**
    * @brief defines LklSet -> no more local Likelihood functions can be added afterwards
    */
    void freeze();
    
    /**
    * @brief evaluates the current value of the likelihood function
    * @return log(Likelihood)
    */
    const tdouble eval_Likelihood();
    /**
    * @brief evaluates the current value of the limit-state function (for reliability analysis)
    * needed for RAmci and RAsubSim
    * @return LSF -> basically the 'likelihood' without the log
    */
    const tdouble eval_RAlsf();
    tdouble& get_cStart() { return cStart; }
    tdouble& get_pa_maxL() { return pa_maxL; }
    const tdouble get_cStart_init() const { return cStart_init; }
    const bool cStart_has_changed() { return  cStart!=cStart_init; /*fabs(cStart-cStart_init)>GlobalVar.TOL();*/ }
    const tdouble get_p();
    void set_TMCMC();
    
    RBRV_constructor& get_RndBox();
    const std::vector<RBRV_set_base*>& get_setvec() const { return setvec; }
};



class flxBayDA;
/**
* @brief A class for storing Bayesian updating objects
*/
class FlxBayUpBox {
  private:
    #if FLX_DEBUG
      /**
      * @brief number of instances of this class - only one class may exist
      */
      static int Cinst;
    #endif
    std::map<std::string, flxBayUp*> box;

    std::map<std::string, flxBayDA*> box_DA;
  public:
    FlxBayUpBox();
    ~FlxBayUpBox();
    /**
    * @brief Insert a Bayesian updating object.
    * @param nameID an unique name (lowercase!!!)
    * @param value the object to store
    * @note an existing object can be redefined. However, the old reference is not vailid anymore.
    */
    void insert ( const std::string& nameID, flxBayUp* obj, const bool errSerious=true);
    /**
    * @brief Get an Bayesian updating object by name.
    * @param name the nameID of the object (lowercase!!!).
    * @return flxBayUp_base&
    */
    flxBayUp& get ( const std::string& nameID );

    /**
    * @brief Insert a Bayesian updating object.
    * @param nameID an unique name (lowercase!!!)
    * @param value the object to store
    * @note an existing object cannot be redefined.
    */
    void insert_DA ( const std::string& nameID, flxBayDA* obj, const bool errSerious=true);
    /**
    * @brief Get an Bayesian updating object by name.
    * @param name the nameID of the object (lowercase!!!).
    * @return flxBayDA&
    */
    flxBayDA& get_DA ( const std::string& nameID );
};



//--------------------------------------------------------------
// Combine models - corresponding to their posterior probability
//--------------------------------------------------------------

class RBRV_entry_value : public RBRV_entry {
  protected:
    const std::string get_type() const { return "value"; }
  public:
    RBRV_entry_value(const std::string& name) : RBRV_entry(name) {}
    virtual ~RBRV_entry_value() {}
    void transform_y2x(const tdouble* const y_vec);
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr);
};

/**
* @brief draw samples - considering the posterior probability of models
*/
class flxBayUp_mProb_set : public RBRV_set_base, public FlxDataBase {
  private:
    const tuint Nmodels;
    flxBayUp_RBRV_set** modelVec;        // the different models (posterior)
    flxVec postPrVec;                        // posterior probability of each model
    tdouble sumPrVec;
    RBRV_entry_RV_uniform p;
    
    std::vector <RBRV_set_base* > setvec;        // contains all sets that underly this object
    tuint NRV_total;                        // total number of random variables (not counting p)
    tuint NOX_total;                        
    flxVec* y_total;
    
    const tuint N_model_res;
    RBRV_entry_value** model_res_list;
    FlxFunction** model_res_map;
    
    void free_mem();                        // deallocates the dynamic memory of this object
    /**
    * @brief returns the ID of the model that belongs to the current realization of p
    */
    const tuint get_model_ID() const;
    void update_model_res(const tuint id);
  public:
    flxBayUp_mProb_set(const bool internal, const std::string& name, const tuint Nmodels, flxBayUp_RBRV_set** modelVec, flxVec priorPrVec, const tuint N_model_res, std::vector<std::string>& model_res_list_Str, FlxFunction** model_res_map);
    virtual ~flxBayUp_mProb_set();
    
    virtual const tuint get_NRV() const { return NRV_total+1; }
    virtual const tuint get_NOX() const { return NOX_total+1+N_model_res; }
    virtual const tuint get_NRV_only_this() const { return 1; }
    virtual const tuint get_NOX_only_this() const { return N_model_res+1; }
    virtual void set_is_valid(const bool is_valid);
    virtual const flxVec& propose_y();
    virtual void transform_y2x();
    virtual const bool allow_x2y() const { return false; }
    /**
    * @brief sets y without setting x
    */
    virtual void set_y(const tdouble* const y_vec);
    virtual void set_y_only_this(const tdouble* const y_vec);
    virtual void get_y(tdouble* const y_vec);
    virtual void get_y_only_this(tdouble* const y_vec);
    virtual flxVec& get_y();
    /**
    * @brief sets x without setting y first
    */
    virtual void set_x(const tdouble* const x_vec);
    virtual void set_x_only_this(const tdouble* const x_vec);
    virtual void get_x(tdouble* const x_vec);
    virtual void get_x_only_this(tdouble* const x_vec);
    virtual const bool check_xVec(const tdouble* xp);
    virtual void get_mean(tdouble* const m_vec);
    virtual void get_mean_only_this(tdouble* const m_vec);
    virtual void get_sd(tdouble* const s_vec);
    virtual void get_sd_only_this(tdouble* const s_vec);
    virtual void find_dependent_sets(std::vector<RBRV_set_base*>& setvecV );
    virtual const tuint group_dependent_sets(std::vector<RBRV_set_base*>& setvecV, const tuint pos_this );
    virtual void print(std::ostream& sout, const std::string prelim, tuint& counter, const bool printID);
};




//------------------------------------------------------------


class RBRV_entry_fun_data : public RBRV_entry_fun {
  protected:
    const tuint paraN;
    const tuint dataN;
    tdouble* dataPtr;
    const tdouble is_log;
    
    const std::string get_type() const { return "fun_data"; }
  public:
    RBRV_entry_fun_data(const std::string& name, FlxFunction* fun, const tuint paraN, FlxIstream_vector *isv, const tdouble is_log);
    virtual ~RBRV_entry_fun_data();
    
    virtual void eval_para() {}
    virtual void transform_y2x(const tdouble* const y_vec);
    virtual const tdouble get_value_log() const { return value; }
};

class RBRV_entry_ref_log : public RBRV_entry {
  protected:
    const tdouble& reflog;
    
    const std::string get_type() const { return "ref_log"; }
  public:
    RBRV_entry_ref_log(const std::string& name, const tdouble& reflog) : RBRV_entry(name), reflog(reflog) {}
    
    virtual void transform_y2x(const tdouble* const y_vec);
    virtual const tdouble get_value_log() const { return reflog; }
    virtual const tdouble get_mean_current_config();
    virtual const tdouble get_sd_current_config();
    virtual const bool check_x(const tdouble xV);
    virtual const bool search_circref(FlxFunction* fcr) { return false; }
};

class RBRV_entry_fun_log : public RBRV_entry_fun {
  protected:

    const std::string get_type() const { return "fun_log"; }
  public:
    RBRV_entry_fun_log(const std::string& name, FlxFunction* fun) : RBRV_entry_fun(name,fun) {}
    virtual ~RBRV_entry_fun_log() { }
    
    virtual const tdouble get_value_log() const { return value; }
};


