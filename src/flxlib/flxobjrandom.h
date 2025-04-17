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

#include "flxobjects.h"
#include "flxBayUp.h"

class FLXLIB_EXPORT FlxCreateObjReaders_RND : public FlxCreateObjReaders {
  public:
    void createObjReaders (FlxObjectReadBox* objReadBox );
    void createFunReaders (FlxData* dataBox );
};

// #################################################################################
// post-processors
// #################################################################################

class PYBIND11_EXPORT post_proc_base {
  public:
    virtual ~post_proc_base() {}
    post_proc_base& operator=(const post_proc_base& rhs) = delete;

    virtual void append_data(const flxVec& vec_full) = 0;
    virtual py::dict eval() = 0;
};

class PYBIND11_EXPORT post_proc_mean_double : public post_proc_base {
  private:
    tdouble sum;
    tulong N;
    const tuint colID;
  public:
    post_proc_mean_double(const tuint colID);

    virtual void append_data(const flxVec& vec_full);
    virtual py::dict eval();
};


// #################################################################################
// dataBox
// #################################################################################

class PYBIND11_EXPORT flxDataBox {
  private:
    const tuint M_in;         // dimension of model input
    const tuint M_out;        // dimension of model output
    const tuint M;            // total dimension (M_in+M_out)
  public:
    flxVec vec_full;
    flxVec vec_out;
    flxVec vec_in;
  private:
    // for writing into a file
      std::ofstream *fstream;
      tuint fs_N_col;
      tuint* fs_cols;
      bool fs_binary;
    // for storing data in memory
      tuint mem_N_reserved;   // total number of data that can be stored
      tuint mem_N;            // total number of data in memory
      tfloat* mem_ptr;
      tuint mem_N_col;
      tuint* mem_cols;
    // for managing post-processors
      std::vector<post_proc_base*> pp_vec;

    tuint* process_col_input(tuint& N_col, py::dict config);
    const tuint extract_colID(py::object col);
    const tuint extract_colID_(py::dict config);
  public:
    flxDataBox(const tuint M_in, const tuint M_out);
    flxDataBox() = delete;
    flxDataBox(flxDataBox& rhs) = delete;
    ~flxDataBox();

    flxDataBox& operator=(const flxDataBox& rhs) = delete;

    // for storing data in memory
      void write2mem(py::dict config);
      py::array_t<tfloat> extract_col_from_mem(py::object col);
      void free_mem();

    // for writing into a file
      void write2file(py::dict config);
      void close_file();

    // for managing post-processors
      post_proc_base& register_post_processor(py::dict config);

    // for consistency checks
      const tuint get_M_in() const { return M_in; }
      const tuint get_M_out() const { return M_out; }
      /**
      * @brief ensure M == M_
      */
      void ensure_M(const tuint M_) const;
      /**
      * @brief ensure M_in == M_
      */
      void ensure_M_in(const tuint M_) const;
      /**
      * @brief ensure M_in == M_
      */
      void ensure_M_out(const tuint M_) const;

    /**
    * @brief registers a data point in the box
    *
    * @note before calling this function, the input/output needs to be assigned to vec_full or (vec_out & vec_in) !!!
    */
    void append_data();


};




//======================== Monte Carlo Integration ===========================


/**
* @brief perform a Monte Carlo Simulation
* @param N number of samples
* @param model the model to sample
* @param sampler for generating samples
* @return Python dictionary with results of the Monte Carlo simulation
*/
FLXLIB_EXPORT void flx_perform_MCS(const tulong N, FlxMtxFun_base& model, RBRV_constructor& sampler, flxDataBox& dbox);


/**
* @brief object class: Monte Carlo Integration
*
* MCI (const-name, Np, g) { ... };
*/
class FlxObjMCI : public FlxObjBase {
  protected:
    /**
    * @brief number of sampling points
    */
    tulong Np;
    /**
    * @brief pointer to the variable storing the result of the calculation (the integral)
    */
    tdouble& theResult;
    /**
    * @brief the function used to calculate Np
    */
    FlxFunction* funNp;
    /**
    * @brief the function to integrate
    */
    FlxFunction* fung;
    /**
    * @brief true: split the sampling into intervals
    */
    bool interv;
    /**
    * @brief log verbose output
    */
    bool verboseLog;
    /**
    * @brief Integration in order to obtain reliability 
    */
    bool reliability;
    /**
    * @brief probability of exceedance (used for confidence studies)
    */
    FlxMtxConstFun* pc;
    FlxString* rbrvsets;
    /**
    * @brief do one Integration-Step (needed for deriving Importance Sampling)
    */
    #if FLX_KAHAN_MCI
      virtual void Integrationstep(pdouble& Integ, tdouble &hits, RBRV_constructor& RndBox) { 
    #else
      virtual void Integrationstep(tdouble& Integ, tdouble &hits, RBRV_constructor& RndBox) { 
    #endif
        RndBox.gen_smp(); const tdouble t=fung->calc(); Integ+=t; if (t>0) ++hits; }
    /**
    * @brief needed to calculate FlxMtxFun (faster) (needed for deriving Importance Sampling)
    */
    virtual void FirstThingsFirst(RBRV_constructor& RndBox);
    /**
    * @brief needed for deriving Importance Sampling
    */
    virtual void LastThingsLast() {}
    /**
    * @brief log additional information after integration
    */
    virtual void log_AddResInfo(std::ostream &sout, const tdouble hits, const tdouble Npd);
    void output_Bayesian_credible_reliablity(std::ostream& sout, const tdouble hits, const tdouble Npd);
    
    void task();
  public:
    FlxObjMCI ( bool dolog, tdouble& theResult, FlxFunction* funNp, FlxFunction* fung, bool interv, bool verboseLog, bool reliability, FlxMtxConstFun* pc, FlxString* rbrvsets ) 
      : FlxObjBase(dolog),  Np(0),theResult(theResult), funNp(funNp), fung(fung),interv(interv), verboseLog(verboseLog), reliability(reliability), pc(pc), rbrvsets(rbrvsets) {}
    virtual ~FlxObjMCI();
};

/**
* @brief object read class: for FlxObjMCI
*/
class FlxObjReadMCI : public FlxObjReadLogBase {
  protected:
    /**
    * @brief read the settings of the Monte Carlo integration
    */
    void read_MCIblock(tdoublePtr& theResult, FlxFunctionPtr& funNp, FlxFunctionPtr& fung, bool errSerious=true);
  public:
    FlxObjReadMCI();
    FlxObjBase* read ();
};

//======================== Importance Sampling ===========================

/**
* @brief object class: Importance Sampling
*
* IPS (const-name, Np, g) SMPL_SPACE { ... };
*/
class FlxObjIpS : public FlxObjMCI {
  protected:
    FlxRndSamplingSpace_Generator_base* sspaceG;
    FlxRndSamplingSpace_base* sspace;
    vdouble sum_weights_hit;
    vdouble sum_weights_nohit;
    vdouble sum_weights_total;
    #if FLX_KAHAN_MCI
      virtual void Integrationstep(pdouble& Integ, tdouble &hits, RBRV_constructor& RndBox);
    #else
      virtual void Integrationstep(tdouble& Integ, tdouble &hits, RBRV_constructor& RndBox);
    #endif
    void FirstThingsFirst(RBRV_constructor& RndBox);
    void LastThingsLast() {delete sspace; }
    void log_AddResInfo(std::ostream &sout, const tdouble hits, const tdouble Npd);
  public:
    FlxObjIpS ( bool dolog, tdouble& theResult, FlxFunction* funNp, FlxFunction* fung, bool interv, bool verboseLog, bool reliability, FlxMtxConstFun* pc, FlxRndSamplingSpace_Generator_base* sspaceG, FlxString* rbrvsets) 
      : FlxObjMCI(dolog, theResult,funNp,fung,interv,verboseLog,reliability,pc,rbrvsets), sspaceG(sspaceG), sum_weights_hit(), sum_weights_nohit(), sum_weights_total() {}
    ~FlxObjIpS() { delete sspaceG; }
};

/**
* @brief object read class: for FlxObjIpS
*/
class FlxObjReadIpS : public FlxObjReadMCI {
  public:
    FlxObjBase* read ();
};


//======================== Line Sampling ===========================

/**
* @brief object class: Monte Carlo Integration
*
* MCI (const-name, Np, g) { ... };
*/
class FlxObjLineSmpl : public FlxObjBase {
  protected:
    /**
    * @brief pointer to the variable storing the result of the calculation (the integral)
    */
    tdouble& theResult;
    /**
    * @brief the function used to calculate Np
    */
    FlxFunction* funNp;
    /**
    * @brief the function to integrate
    */
    FlxFunction* fung;
    FlxMtxConstFun* LS_SPNT;        // start-point of the Line sampling
    FlxFunction* LS_tolF;        // tolerance parameter of line search
    FlxFunction* LS_max_iterF;        // maximum number of iterations in line search
    const bool extended_ls;
    /**
    * @brief log verbose output
    */
    const bool verboseLog;
    FlxString* rbrvsets;
    const bool use_bisec;
    std::vector< std::pair<tdouble,tdouble> > line_hist;
    
    RBRV_constructor* RndBoxp;
    
    void hist_push(const tdouble c, const tdouble g);
    const tdouble hist_eval(const tdouble betaNorm);
    
    const tdouble LSF_call(const tdouble c, const flxVec &rv_base, flxVec &rv_prop, const flxVec &betaVec, tulong& N_LSF_calls);
    /**
    * @brief performs the actual line-serach
    * @returns 'c' of the root-serach & fdright
    * @param fdright true: failure domain right of 'c'; false: failure domain left of 'c'
    */
    const tdouble perform_line_search(const flxVec &rv_base, flxVec &rv_prop, const flxVec &betaVec, const tdouble tol, const tuint iter_max,tulong& N_LSF_calls, bool& fdright, bool& found, const tdouble startV, const tdouble endV);
    const tdouble perform_line_search_rgfsi(const flxVec &rv_base, flxVec &rv_prop, const flxVec &betaVec, const tdouble tol, const tuint iter_max,tulong& N_LSF_calls, bool& fdright, bool& found, const tdouble startV, const tdouble endV);
    const tdouble perform_line_search_bisec(const flxVec &rv_base, flxVec &rv_prop, const flxVec &betaVec, const tdouble tol, const tuint iter_max,tulong& N_LSF_calls, bool& fdright, bool& found, const tdouble startV, const tdouble endV);
    void task();
  public:
    FlxObjLineSmpl ( bool dolog, tdouble& theResult, FlxFunction* funNp, FlxFunction* fung,FlxMtxConstFun* LS_SPNT,FlxFunction* LS_tolF,FlxFunction* LS_max_iterF, const bool extended_ls, const bool verboseLog, FlxString* rbrvsets, const bool use_bisec ) 
      : FlxObjBase(dolog),  theResult(theResult), funNp(funNp), fung(fung), LS_SPNT(LS_SPNT), LS_tolF(LS_tolF), LS_max_iterF(LS_max_iterF), extended_ls(extended_ls), verboseLog(verboseLog), rbrvsets(rbrvsets), use_bisec(use_bisec), RndBoxp(NULL) {}
    virtual ~FlxObjLineSmpl() { delete funNp; delete fung; delete LS_SPNT; delete LS_tolF; delete LS_max_iterF; delete rbrvsets; }
};

/**
* @brief object read class: for FlxObjMCI
*/
class FlxObjReadLineSmpl : public FlxObjReadLogBase {
  public:
    FlxObjReadLineSmpl();
    FlxObjBase* read ();
};



//======================== Subset Simulation ===========================



class SuS_csm_evalStorage : public FlxDataBase {
  private:
    FlxFunction* h_value;        // ... for the kernel
    FlxString* kernel_name;
    FlxString* MCMCmethS;
    FlxFunction* csm_p;
    FlxFunction* csm_nmax;
    FlxFunction* csm_p_single;
    FlxFunction* csm_nmax_single;
    FlxFunction* dcs_pSD;
  public:
    SuS_csm_evalStorage ( FlxFunction* h_value, FlxString* kernel_name, FlxString* MCMCmethS, FlxFunction* csm_p, FlxFunction* csm_nmax, FlxFunction* csm_p_single, FlxFunction* csm_nmax_single, FlxFunction* dcs_pSD );
    ~SuS_csm_evalStorage();
    
    /**
    * @brief evaluates and returns the csm-method
    * @note deallocates list if an error is thrown
    */
    FlxBayUP_csm_base* eval (FlxBayUp_Update_List* list);
};


/**
* @brief object class: performs Subset simulation
*/
class FlxObjSuS : public FlxObjOutputBase {
  protected:
    FlxFunction* Nc;                // number of chains
    FlxFunction* Ncl;                // length of a chain
    FlxFunction* max_runs;        // maximum number of SSS-iterations
    const FlxBayUp_Update_List::randomizeTec randomize;
    flxBayUp_adaptive_ctrl_base* adpt_ctrl;
    const Flx_SuS_Control susControl;
    SuS_csm_evalStorage* csm_eval;
    
    FlxString* rbrvsets;
    FlxFunction* lsf;
    
    void task();
  public:
    FlxObjSuS ( const bool dolog, const std::string& ostreamV, FlxFunction* Nc, FlxFunction* Ncl, FlxFunction* max_runs, const FlxBayUp_Update_List::randomizeTec randomize, flxBayUp_adaptive_ctrl_base* adpt_ctrl, const Flx_SuS_Control& susControl, SuS_csm_evalStorage* csm_eval, FlxString* rbrvsets, FlxFunction* lsf);
    virtual ~FlxObjSuS();
};

/**
* @brief object read class: for FlxObjSus (Subset Simulation)
*/
class FlxObjReadSuS : public FlxObjReadOutputBase {
  protected:
    FlxBayUp_Update_List::randomizeTec get_randomize_id();
    flxBayUp_adaptive_ctrl_base* get_adpt_ctrl();
    SuS_csm_evalStorage* get_csm_eval();
    const Flx_SuS_Control get_susControl();
  public:
    FlxObjReadSuS();
    ~FlxObjReadSuS();
    FlxObjBase* read ();
    
    static flxBayUp* lastSuS;
};

class FlxObjSus_level_info : public FlxObjBase {
  private:
    FlxMtxConstFun* VecStr;
    FlxString* nameID;
    FlxFunction* pidf;
    FlxFunction* pidf2;
    
    void task();
  public:
    FlxObjSus_level_info ( const bool dolog, FlxMtxConstFun* VecStr, FlxString* nameID, FlxFunction* pidf, FlxFunction* pidf2 ) : FlxObjBase(dolog), VecStr(VecStr), nameID(nameID), pidf(pidf), pidf2(pidf2) {}
    virtual ~FlxObjSus_level_info() { delete VecStr; if (nameID) delete nameID; delete pidf; if (pidf2) delete pidf2;  }
};

class FlxObjReadSus_level_info : public FlxObjReadBase {
  public:
    FlxObjBase* read();
};

//======================== FORM ===========================


class FlxObjFORM_base : public FlxObjOutputBase {
  protected:
    tuint DIM;  
    FlxFunction* LSF;
    FlxFunction* fdstep;
    FlxFunction* epsdyfF;
    FlxFunction* eps1;
    FlxFunction* eps2;
    FlxFunction* iHLRF_lambda_start;
    FlxFunction* iHLRF_epsilon;
    FlxFunction* iHLRF_reduce;
    FlxMtxConstFun* xstart;
    FlxMtxConstFun* dx_min;
    unsigned int maxIter;
    /**
    * @brief id of the finite difference method
    * 1: forward
    * 2: central
    * 3: backward
    */
    int fd_method;
    /**
    * @brief id of the optimization method
    * 1: HL-RF
    * 2: iHL-RF
    * 3: iGP
    */
    int opt_method;
    bool dxdyAnalytical;
    const bool verboseLog;
    FlxString* rbrvsets;
        
    RBRV_constructor* RndBox;
    
    void update_Start();
    void eval_xStart(flxVec &x);
    flxVec do_FORM(flxVec &x, flxVec &y, tdouble &beta_new, tuint& LSFcalls,const bool only_partial=false);
  public:
    FlxObjFORM_base ( bool dolog, FlxFunction* LSFv, FlxFunction* fdstepV, FlxFunction* epsdyfF, FlxFunction* eps1, FlxFunction* eps2, FlxFunction* iHLRF_lambda_start, FlxFunction* iHLRF_epsilon, FlxFunction* iHLRF_reduce, unsigned int maxIter, bool verbose, std::string ostreamV, bool dxdyAnalytical, FlxMtxConstFun* xstart, FlxMtxConstFun* dx_min, int fd_method, int opt_method, FlxString* rbrvsets ) 
     : FlxObjOutputBase(dolog,ostreamV), LSF(LSFv), fdstep(fdstepV),epsdyfF(epsdyfF), eps1(eps1), eps2(eps2), iHLRF_lambda_start(iHLRF_lambda_start),iHLRF_epsilon(iHLRF_epsilon), iHLRF_reduce(iHLRF_reduce), xstart(xstart), dx_min(dx_min), maxIter(maxIter), fd_method(fd_method),opt_method(opt_method), dxdyAnalytical(dxdyAnalytical), verboseLog(verbose), rbrvsets(rbrvsets), RndBox(NULL) {}
    ~FlxObjFORM_base();
};


/**
* @brief object class
*
* form ( Y, X, LSF );
*/
class FlxObjFORM : public FlxObjFORM_base {
  private:
    const std::string cnamey;
    const std::string cnamex;
    FlxString* betaDP;
    const bool only_partial;

    void task();
  public:
    FlxObjFORM ( const bool dolog, const std::string cnamey, const std::string cnamex, FlxMtxConstFun* xstart, FlxFunction* LSFv, FlxFunction* fdstepV, FlxFunction* epsdyfF, FlxFunction* eps1, FlxFunction* eps2, FlxFunction* iHLRF_lambda_start, FlxFunction* iHLRF_epsilon, FlxFunction* iHLRF_reduce, unsigned int maxIter, FlxString* betaDP, bool verbose, bool dxdyAnalytical, FlxMtxConstFun* dx_min, int fd_method, int opt_method, FlxString* rbrvsets, const bool only_partial ) 
     : FlxObjFORM_base(dolog,LSFv,fdstepV,epsdyfF,eps1,eps2,iHLRF_lambda_start,iHLRF_epsilon,iHLRF_reduce,maxIter,verbose,"cout",dxdyAnalytical,xstart,dx_min,fd_method,opt_method,rbrvsets), cnamey(cnamey), cnamex(cnamex), betaDP(betaDP), only_partial(only_partial) {}
    ~FlxObjFORM();
    
    static void sensitivities(const flxVec &y, RBRV_constructor& RndBox, std::ostream& slog, flxVec* svp = NULL);
};

class FlxObjReadFORM_base : public FlxObjReadOutputBase {
  public:
    FlxObjReadFORM_base();
};

/**
* @brief object read class: for FlxObjFORM
*/
class FlxObjReadFORM : public FlxObjReadFORM_base {
  private:
    const bool only_partial;
  public:
    FlxObjReadFORM(const bool only_partial=false);
    FlxObjBase* read ();
};

/**
* @brief object class
*
* form_pdf (  );
*/
class FlxObjFORM_pdf : public FlxObjFORM_base {
  private:
    FlxFunction* rvfun;
    FlxFunction* lboundF;
    FlxFunction* uboundF;
    FlxFunction* nintervalF;
    void task();
  public:
    FlxObjFORM_pdf ( bool dolog, FlxFunction* rvfun, FlxFunction* lboundF, FlxFunction* uboundF, FlxFunction* nintervalF, FlxMtxConstFun* xstart, FlxFunction* fdstepV, FlxFunction* epsdyfF, FlxFunction* eps1, FlxFunction* eps2, FlxFunction* iHLRF_lambda_start, FlxFunction* iHLRF_epsilon, FlxFunction* iHLRF_reduce, unsigned int maxIter, bool verbose, std::string ostreamV, bool dxdyAnalytical, FlxMtxConstFun* dx_min, int fd_method, int opt_method, FlxString* rbrvsets );
    ~FlxObjFORM_pdf();
};

/**
* @brief object read class: for FlxObjFORM_pdf
*/
class FlxObjReadFORM_pdf : public FlxObjReadFORM_base {
  public:
    FlxObjReadFORM_pdf();
    FlxObjBase* read ();
};


class FlxObjFORM_betaSensitivities : public FlxObjOutputBase {
  private:
    FlxMtxConstFun* sv;
    FlxMtxConstFun* rvy;
    FlxString* rvsets;
    void task();
  public:
    FlxObjFORM_betaSensitivities ( const bool dolog, const std::string& ostreamV,  FlxMtxConstFun* sv, FlxMtxConstFun* rvy, FlxString* rvsets ) : FlxObjOutputBase(dolog,ostreamV), sv(sv), rvy(rvy), rvsets(rvsets) {}
    ~FlxObjFORM_betaSensitivities() { delete sv; delete rvy; delete rvsets; }
};

/**
* @brief object read class: for FlxObjFORM
*/
class FlxObjReadFORMbetaSensitivities : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};


// -------------------------- Kernel Density Estimation ----------------------------------------------------------------------

/**
* @brief object read class: for FunPlot-Objects
*/
class FlxObjReadKDE : public FlxObjReadOutputBase {
  public:
    FlxObjReadKDE();
    FlxObjBase* read();
};

/**
* @brief object class: plot some functions (for GNUplot)
*
* KDE ( N, KERNEL, BANDWIDTH ) { ... };
*/
class FlxObjKDE : public FlxObjOutputBase {
  protected:
    FlxFunction* funR;
    FlxString* rbrvsets;
    FlxFunction* N;
    FlxRndKernel_base* kernel;
    FlxFunction* h;
    FlxFunction* lbound;
    FlxFunction* ubound;
    FlxFunction* Ninterval;
    bool do_cdf;
    virtual void task();
  public:
    FlxObjKDE( bool dolog, FlxFunction* funR, FlxString* rbrvsets, FlxFunction* N, FlxRndKernel_base* kernel, FlxFunction* h, FlxFunction* lbound, FlxFunction* ubound, FlxFunction* Ninterval, bool do_cdf, std::string ostreamV );
    ~FlxObjKDE();
};


// -------------------------- Statistical Commands ----------------------------------------------------------------------


/**
* @brief object read class: for statsmp-Objects
*/
class FlxObjReadMCSsensi : public FlxObjReadOutputBase {
  public:
    FlxObjReadMCSsensi();
    FlxObjBase* read();
};

struct kvEl
{
    double key;
    double val;
    
    const bool operator<(const kvEl& ref) const {
        return key<ref.key;
    }
};

/**
* @brief object class: descriptive statistics of a sample set
*/
class FlxObjMCSsensi : public FlxObjOutputBase {
  protected:
    FlxMtxConstFun* mcn;
    FlxString* istrm;
    FlxFunction* rvN;
    FlxFunction* nb;
    
    void task();
  public:
    FlxObjMCSsensi( bool dolog, std::string ostreamV, FlxMtxConstFun* mcn, FlxString* istrm, FlxFunction* rvN, FlxFunction* nb );
    ~FlxObjMCSsensi();
};

/**
* @brief object read class: for statsmp-Objects
*/
class FlxObjReadStatSmp : public FlxObjReadOutputBase {
  public:
    FlxObjReadStatSmp();
    FlxObjBase* read();
};

/**
* @brief object class: descriptive statistics of a sample set
*/
class FlxObjStatSmp : public FlxObjOutputBase {
  protected:
    /**
    * @brief the name of the input stream
    */
    FlxString* isname;
    FlxIstream_vector* is_vec;
    /**
    * @brief the string to add to the name of the const-variables (to store calculated values)
    */
    FlxString* add_name;
    const std::string add_name_str;
    /**
    * @brief required for qdouble-summation
    */
    FlxFunction* Np;
    const tuint Np_val;
    /**
    * @brief the option (1: single sample set; 2: 2 correlation of 2 sample sets)
    */
    int optionP;
    /**
    * @brief output statistics about significant figures
    */
    const  bool sigfig;
    
    void sigfig_mean(const tdouble estim_mu, const tdouble estim_sd, const tulong N);
    
    void task();
  public:
    FlxObjStatSmp( bool dolog, std::string ostreamV, FlxString* isname, FlxString* add_name, FlxFunction* Np, int optionP, const bool sigfig );
    FlxObjStatSmp( bool dolog, std::string ostreamV, FlxIstream_vector* is_vec, const std::string& add_name_str, const tuint Np_val, int optionP, const bool sigfig );
    ~FlxObjStatSmp();
};


/**
* @brief object read class: for sortsmp-Objects
*/
class FlxObjReadSortSmp : public FlxObjReadOutputBase {
  public:
    FlxObjReadSortSmp();
    FlxObjBase* read();
};

/**
* @brief object class: sort numbers in a file
*/
class FlxObjSortSmp : public FlxObjOutputBase {
  protected:
    /**
    * @brief the name of the input stream
    */
    FlxString* isname;
    /**
    * @brief maxiumum block-size for sorting
    */
    FlxFunction* Np;
    
    void task();
  public:
    FlxObjSortSmp( bool dolog, std::string ostreamV, FlxString* isname, FlxFunction* Np );
    ~FlxObjSortSmp();
};

/**
* @brief object read class: for smpplot-Objects
*/
class FlxObjReadSmpPlot : public FlxObjReadOutputBase {
  public:
    FlxObjReadSmpPlot();
    FlxObjBase* read();
};

/**
* @brief object class: plotting histograms, ...
*/
class FlxObjSmpPlot : public FlxObjOutputBase {
  protected:
    /**
    * @brief the name of the input stream
    */
    FlxString* isname;
    /**
    * @brief type of the plot
    */
    FlxFunction* typeFun;
    
    bool autoBound;
    FlxFunction* xminU;
    FlxFunction* xmaxU;
    int binEstimator;
    FlxFunction* NbinsU;
    
    void task();
  public:
    FlxObjSmpPlot( const bool dolog, const std::string ostreamV, FlxString* isname, FlxFunction* typeFun, const bool autoBound, FlxFunction* xminU, FlxFunction* xmaxU, const int binEstimator, FlxFunction* NbinsU, const int prec, const int fixW );
    ~FlxObjSmpPlot();
};



class FlxObjQQplot : public FlxObjOutputBase {
  protected:
    /**
    * @brief the name of the input stream
    */
    FlxString* isname;
    /**
    * @brief distribution type for the qq-plot
    */
    RBRV_entry_RV_base* rep;
    
    void task();
  public:
    FlxObjQQplot( const bool dolog, const std::string ostreamV, FlxString* isname, RBRV_entry_RV_base* rep );
    ~FlxObjQQplot();
};

class FlxObjReadQQplot : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read();
};





// -------------------------- FlxFunction ----------------------------------------------------------------------


class FunReadFunSmpCDF : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};

/**
* @brief arithmetic class: CDF of a set of samples
*/
class FunSmpCDF : public FunBase, public FlxDataBase {
  private:
    FlxString* isname;
    FunBase* val;
    const bool inverse;
  public:
    FunSmpCDF(FlxString* isname, FunBase* val, const bool inverse) : isname(isname), val(val), inverse(inverse)  {}
    ~FunSmpCDF();
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return val->search_circref(fcr); }
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return val->dependOn_Const(thenumber); }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
    virtual const bool evalw() { return val->evalw(); }
    
    static tdouble inv_cdf(const tdouble thr, const tdouble* const tp, const tuint N);
};



// -------------------------- Sensitivity analysis ----------------------------------------------------------------------


class FlxObjSensi_s1o_new : public FlxObjOutputBase {
  protected:
    FlxString* nameID;                       // name associated with Bayesian data analysis
    FlxFunction* Nlearn;
    FlxFunction* x_dim;

    void task();
  public:
    FlxObjSensi_s1o_new ( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* Nlearn, FlxFunction* x_dim);
    virtual ~FlxObjSensi_s1o_new() { delete nameID; delete Nlearn; delete x_dim; }
};

class FlxObjReadSensi_s1o_new : public FlxObjReadOutputBase {
  public:
    FlxObjReadSensi_s1o_new();
    FlxObjBase* read ();
};




class FlxObjSensi_s1o_add : public FlxObjOutputBase {
  protected:
    FlxString* nameID;                       // name associated with Bayesian data analysis
    FlxFunction* value_x;                     // name associated with Bayesian data analysis
    FlxFunction* value_y;
    FlxMtxConstFun* xvec;

    void task();
  public:
    FlxObjSensi_s1o_add ( const bool dolog, const std::string& ostreamV, FlxString* nameID, FlxFunction* value_x, FlxMtxConstFun* xvec, FlxFunction* value_y);
    virtual ~FlxObjSensi_s1o_add();
};

class FlxObjReadSensi_s1o_add : public FlxObjReadOutputBase {
  public:
    FlxObjBase* read ();
};


class FunSensi_s1o_eval : public FunBase, public FlxDataBase {
  private:
    FlxString* nameID;
  public:
    FunSensi_s1o_eval(FlxString* nameID) : nameID(nameID)  {}
    ~FunSensi_s1o_eval() { delete nameID; };
    const tdouble calc();
    const bool search_circref(FlxFunction* fcr) { return false; }
    const std::string write();
    const tuint precedence() { return 0; }
    const bool dependOn_Const(const tdouble* const thenumber) { return false; }
    const bool optimize(FunBasePtr& optf, const Fun_OptimizeInfo &foi) { return false; }
};

class FunReadFunSensi_s1o_eval : public FunReadFunBase {
  public:
    FunBase* read ( bool errSerious );
};


/**
* @brief object class:
*/
class FlxObjSensi_s1o_dist : public FlxObjBase {
  private:
    FlxString* nameID;
    FlxMtxConstFun* MtxConstStr;
    FlxFunction* Nfun;

    void task();
  public:
    FlxObjSensi_s1o_dist ( const bool dolog, FlxString* nameID, FlxMtxConstFun* MtxConstStr, FlxFunction* Nfun ) : FlxObjBase(dolog), nameID(nameID), MtxConstStr(MtxConstStr), Nfun(Nfun) {}
    ~FlxObjSensi_s1o_dist();
};

/**
* @brief object read class: for FlxObjRBRV_vec_get
*/
class FlxObjReadSensi_s1o_dist : public FlxObjReadBase {
  public:
    FlxObjReadSensi_s1o_dist();
    FlxObjBase* read ();
};

