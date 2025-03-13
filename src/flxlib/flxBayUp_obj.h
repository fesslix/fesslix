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

#include "flxBayUp.h"

#include "flxobjects.h"
#include "flxrbrv_rvs_read.h"
#include "flxobjrandom.h"


class FLXLIB_EXPORT FlxCreateObjReaders_BU : public FlxCreateObjReaders {
  public:
    void createObjReaders (FlxObjectReadBox* objReadBox );
    void createFunReaders (FlxData* dataBox );
};


#ifdef fesslix_flxBayUp_obj_CPP
  FlxBayUpBox* BayUpBox = NULL;
#else
  extern FlxBayUpBox* BayUpBox;
#endif
  


//----------------------------------------------------------------------------------------------

/**
* @brief object class: creates a new Bayesian-updating object
*/
class FlxObjBayUp_new : public FlxObjOutputBase {
  private:
    FlxString* nameID;
    FlxString* rbrvsets;
    FlxFunction* cStart;
    FlxFunction* scaleconst;
    const bool cStart_isLog;
    void task();
  public:
    FlxObjBayUp_new ( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxString* rbrvsets, FlxFunction* cStart, FlxFunction* scaleconst, const bool cStart_isLog );
    virtual ~FlxObjBayUp_new();
};

/**
* @brief object read class: for FlxObjBayUp_new
*/
class FlxObjReadBayUp_new : public FlxObjReadOutputBase {
  public:
    FlxObjReadBayUp_new();
    FlxObjBase* read ();
};


/**
* @brief object class: explicit formulation of the likelihood function
*/
class FlxObjBayUp_likelihood : public FlxObjOutputBase {
  private:
    FlxString* nameID;
    FlxFunction* lfun;
    const bool is_log;
    void task();
  public:
    FlxObjBayUp_likelihood ( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* lfun, const bool is_log );
    virtual ~FlxObjBayUp_likelihood();
};

/**
* @brief object class: explicit formulation of the likelihood function - data form an ivstream
*/
class FlxObjBayUp_likelihood_data : public FlxObjOutputBase {
  private:
    FlxString* nameID;
    const tuint paraN;
    /**
    * @brief the name of the input stream
    */
    FlxString* isname;
    FlxFunction* lfun;
    const bool is_log;
    void task();
  public:
    FlxObjBayUp_likelihood_data ( const bool dolog, const std::string& ostreamV, FlxString* nameID, const tuint paraN, FlxString* isname, FlxFunction* lfun, const bool is_log );
    virtual ~FlxObjBayUp_likelihood_data();
};

/**
* @brief object read class: for FlxObjBayUp_likelihood
*/
class FlxObjReadBayUp_likelihood : public FlxObjReadOutputBase {
  public:
    FlxObjReadBayUp_likelihood();
    FlxObjBase* read ();
};


/**
* @brief object class: explicit formulation of the likelihood function - uncertain observations form an ivstream
*/
class FlxObjBayUp_uncertobsv : public FlxObjOutputBase {
  private:
    FlxString* nameID_BU;
    const tuint paraN;
    /**
    * @brief the name of the input stream
    */
    FlxString* isname;
    FlxFunction* lfun;
    FlxString* set_name;
    std::vector<RBRV_entry_read_base*> set_entries;
    const bool is_log;
    
    void task();
  public:
    FlxObjBayUp_uncertobsv ( const bool dolog, const std::string& ostreamV, FlxString* nameID_BU, const tuint paraN, FlxString* isname, FlxFunction* lfun, FlxString* set_name, std::vector<RBRV_entry_read_base*>& set_entries, const bool is_log );
    virtual ~FlxObjBayUp_uncertobsv();
};

/**
* @brief object read class: for FlxObjBayUp_uncertobsv
*/
class FlxObjReadBayUp_uncertobsv : public FlxObjReadOutputBase {
  public:
    FlxObjReadBayUp_uncertobsv();
    FlxObjBase* read ();
};


/**
* @brief object class: explicit formulation of the likelihood function
*/
class FlxObjBayUp_glbllikelihood : public FlxObjOutputBase {
  private:
    FlxString* nameID;
    FlxFunction* lfun;
    const bool is_log;
    const flxBayUp::MethCategory methCat;
    void task();
  public:
    FlxObjBayUp_glbllikelihood ( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* lfun, const bool is_log, const flxBayUp::MethCategory methCat );
    virtual ~FlxObjBayUp_glbllikelihood();
};

/**
* @brief object read class: for FlxObjBayUp_glbllikelihood
*/
class FlxObjReadBayUp_glbllikelihood : public FlxObjReadOutputBase {
  private:
    /**
    * @brief if true; the log-option is not available
    */
    const flxBayUp::MethCategory methCat;
  public:
    FlxObjReadBayUp_glbllikelihood(const flxBayUp::MethCategory methCat);
    FlxObjBase* read ();
};


/**
* @brief object class: performs the Bayesian updating
*/
class FlxObjBayUp_update : public FlxObjSuS {
  private:
    FlxString* nameID;
    FlxFunction* Ns_final;        // final number of samples (seeds!)
    FlxFunction* Nburn;
    const bool use_cStart;
    const FlxBayUp_Update_List::MethType meth_id;
    const bool log_LSF;
    
    void task();
  public:
    FlxObjBayUp_update ( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* Nc, FlxFunction* Ncl, FlxFunction* Nburn, FlxFunction* Ns_final, FlxFunction* max_runs, const FlxBayUp_Update_List::randomizeTec randomize, flxBayUp_adaptive_ctrl_base* adpt_ctrl, const bool use_cStart, const FlxBayUp_Update_List::MethType meth_id, const bool log_LSF, const Flx_SuS_Control& susControl, SuS_csm_evalStorage* csm_eval );
    virtual ~FlxObjBayUp_update();
};

/**
* @brief object read class: for FlxObjBayUp_update
*/
class FlxObjReadBayUp_update : public FlxObjReadSuS {
  public:
    FlxObjReadBayUp_update();
    FlxObjBase* read ();
};


/**
* @brief object class: changes the method that draws samples from the posterior
*/
class FlxObjBayUp_Set : public FlxObjBase {
  private:
    FlxString* setname;
    std::vector<FlxString*> BUname_vec;
    std::vector<FlxFunction*> priorWeight_vec;
    const tuint Nmodels;
    
    std::vector<FlxString*> model_res_list_Str; 
    std::vector<FlxFunction**> model_res_map;    
    const tuint Nvalues;
    
    void task();
  public:
    FlxObjBayUp_Set ( const bool dolog, FlxString* setname, const std::vector<FlxString*>& BUname_vec, const std::vector<FlxFunction*>& priorWeight_vec, std::vector<FlxString*> model_res_list_Str, std::vector<FlxFunction**> model_res_map );
    virtual ~FlxObjBayUp_Set();
};

/**
* @brief object read class: for FlxObjBayUp_Set
*/
class FlxObjReadBayUp_Set : public FlxObjReadBase {
  public:
    FlxObjBase* read();
};


/**
* @brief object class: changes the method that draws samples from the posterior
*/
class FlxObjBayUp_Reset_Smpls : public FlxObjBase {
  private:
    FlxString* nameID;
    
    void task();
  public:
    FlxObjBayUp_Reset_Smpls ( const bool dolog, FlxString* nameID ) : FlxObjBase(dolog), nameID(nameID) {}
    virtual ~FlxObjBayUp_Reset_Smpls() { delete nameID; }
};

/**
* @brief object read class: for FlxObjBayUp_Reset_Smpls
*/
class FlxObjReadBayUp_Reset_Smpls : public FlxObjReadBase {
  public:
    FlxObjBase* read();
};


/**
* @brief arithmetic class: numerical integration of a function
*/
class FunBayUp_Prop : public FunBase {
  private:
    flxBayUp& buo;
    FunBase* pidf;
    
    const tdouble calc_help(const tuint PID);
  public:
    FunBayUp_Prop(flxBayUp& buo, FunBase* pidf) : buo(buo), pidf(pidf)  {}
    ~FunBayUp_Prop();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return pidf->search_circref(fcr); }
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return pidf->dependOn_Const(thenumber); }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
    virtual const bool evalw() { return pidf->evalw(); }
};

class FunReadFunBayUp_Prop : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};


/**
* @brief arithmetic class: returns evaulates lsf of problem (if defined)
*/
class FunBayUp_lsf : public FunBase {
  private:
    flxBayUp& buo;
  public:
    FunBayUp_lsf(flxBayUp& buo) : buo(buo)  {}
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return false; }
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadFunBayUp_lsf : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};



/**
* @brief arithmetic class: convolution integral - solution with IPS
*/
class FunConvExp : public FunBase, public FlxDataBase {
  private:
    RBRV_set_base* dist_set;        // the modelling error
    RBRV_set_base* pulse_set;        // the measurement error
    FlxString* GM_name;                // the model output
    FlxString* DO_name;                // the observed output
    FlxString* name_dist_set;
    FlxString* name_pulse_set;
    const tuint seed;                // the seed to initialize the rng before each execution
    const tuint Ninteg;                // number of samples used for integration
    const tdouble eps;                // stopping criteria (size of the interval)
    const tuint Nrgir;                // number of initial runs of the random number generator
    
    tuint N;                        // dimension of the problem
    // NOTE: orignial space is dist_set!
    flxVec* cv;                        // the center of the IS density (original space)
    flxVec* y;                        // the current sample (in standard normal space)
    flxVec* eps_Mo;                // the current sample in original space
    flxVec* eps_me;                // = Odif-eps_Mo: the current measurment error
    flxVec* Odif;                // = DOp-GMp: difference between observation and model output (in original space)
    tdouble* DOp;                // temporary pointer -> Observation
    tdouble* GMp;                // temporary pointer -> Model output
    
    const tdouble get_pulse_log();        // returns the pulse of the current eps_Mo
    const tdouble compute_cv();
    
  public:
    FunConvExp ( FlxString* GM_name, FlxString* DO_name, FlxString* name_dist_set, FlxString* name_pulse_set, const tuint seed, const tuint Ninteg, const tdouble eps, const tuint Nrgir);
    virtual ~FunConvExp();
    
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr);
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber);
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

/**
* @brief object read class: for FunConvExp
*/
class FunReadFunConvExp : public FunReadFunBase {
  public:
    FunBase* read( bool errSerious );
};



//----------------------------------------------------------------------------------------------

class FlxObjBayDA_new : public FlxObjOutputBase {
  protected:
    FlxString* nameID;                       // name associated with Bayesian data analysis
    FlxMtxConstFun* mcf_data;                // where to store
    FlxFunction* id_transform;
    FlxFunction* Nchain;
    FlxFunction* Nburn;
    FlxFunction* Ntune;
    FlxFunction* Npost;
    FlxFunction* N_adapt;
    FlxFunction* t_plaus;
    FlxMtxConstFun* dtf;
    FlxString* pvec_str;
    FlxString* distid_str;

    void task();
  public:
    FlxObjBayDA_new ( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxMtxConstFun* mcf_data, FlxFunction* id_transform,
                        FlxFunction* Nchain, FlxFunction* Nburn, FlxFunction* Ntune, FlxFunction* Npost, FlxFunction* N_adapt, FlxFunction* t_plaus, FlxMtxConstFun* dtf, FlxString* pvec_str, FlxString* distid_str);
    virtual ~FlxObjBayDA_new();
};

class FlxObjReadBayDA_new : public FlxObjReadOutputBase {
  public:
    FlxObjReadBayDA_new();
    FlxObjBase* read ();
};


class FlxObjBayDA_sample : public FlxObjOutputBase {
  protected:
    FlxString* nameID;                       // name associated with Bayesian data analysis

    void task();
  public:
    FlxObjBayDA_sample ( const bool dolog, const std::string& ostreamV, FlxString* nameID);
    virtual ~FlxObjBayDA_sample() { delete nameID; }
};

class FlxObjReadBayDA_sample : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};



